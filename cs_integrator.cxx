#include "FourVec.hxx"
#include "PolarFourVec.hxx"
#include "processes.hxx"
#include "argparse.hpp"
#include "Bound.hxx"

#include <iostream>
#include <thread> 
#include <string>
#include <mutex> 
#include <fstream>
#include <stdexcept> 

#include <TString.h> 
#include <TError.h> 

//#define EDCS_DEBUG
#include "DifferentialCSIntegrator.hxx"

namespace {
    constexpr int D = 4; 

    constexpr double deg = 3.1415926536 / 180.; 

    constexpr double me = 0.501;

    enum EProcessName {
        kTrident_electron =1,
        kTrident_positron, 
        kBH_photoproduction  
    };
}

void output_c_array(
    Bound<double> energy, Bound<double> cos_theta, 
    int npts_energy, int npts_cos_theta, 
    std::vector<std::vector<double>>& data, std::string path
); 

int main(int argc, char* argv[])
{
    //parse arguments
    argparse::ArgumentParser program("cs_integrator"); 

    //positional argument is the output file 
    std::string path_output;
        
    try {
        program.add_argument("process")
            .required()
            .help("specific cross section to integrate")
            .choices("trident_electron", "trident_positron", "bh_photoproduction");
            
        program.add_argument("path_output")
            .required()
            .help("Path to the output file destination")
            .store_into(path_output);
            
        //beam energy 
        program.add_argument("--beam-energy")
            .help("beam energy of incident electron")
            .scan<'g', double>()
            .default_value(2200.)
            .nargs(1);

        //energy range
        program.add_argument("-e", "--energy-range")
            .help("energy range to scan for outgoing particle")
            .scan<'g', double>()
            .default_value(std::vector<double>{1000., 2200.})
            .nargs(2);

        //beam energy 
        program.add_argument("-t", "--rel-error")
            .help("relative error of calculated cross section")
            .scan<'g', double>()
            .default_value(1e-3)
            .nargs(1);
            
        //number of scan points in energy
        program.add_argument("--n-pts-energy")
            .help("Number of scan-points in the energy range")
            .scan<'i', int>() 
            .default_value(25)
            .nargs(1);
        
        //maximum number of iterations to attempt 
        program.add_argument("--max-iterations")
            .help("Maximum number of integration steps to attempt")
            .scan<'i', unsigned int>()
            .default_value(30000000)
            .nargs(1);
            
        //range to scan in cos(theta)
        program.add_argument("--cos-theta-range")
            .help("Range to scan in cos(theta)")
            .nargs(2)
            .scan<'g', double>() 
            .default_value(std::vector<double>{1.-0.14, 1.-0.08});
            
        //number of points to scan in cos(theta)
        program.add_argument("--n-pts-cos-theta")
            .help("Number of scan points in the cosine-scan range")
            .nargs(1)
            .scan<'i', int>() 
            .default_value(25);
    
    } catch (const std::exception& e) {

        Error(__func__, "Something went wrong trying to create args. what(): %s", e.what()); 
        return 1; 
    }
    

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1; 
    }

    const std::map<std::string, EProcessName> process_map{
        {"bh_photoproduction",  kBH_photoproduction},
        {"trident_electron",    kTrident_electron},
        {"trident_positron",    kTrident_positron}
    };

    auto process_it = process_map.find(program.get<std::string>("process")); 
    if (process_it == process_map.end()) {
        Error("cs_integrator", "Process %s is not a valid name", program.get<std::string>("process").c_str());
        return 1; 
    }
    const auto process = process_it->second; 

    //relative error to shoot for 
    const double rel_error = program.get<double>("--rel-error"); 

    const long int max_iterations = program.get<unsigned int>("--max-iterations"); 

    const double beam_E = program.get<double>("--beam-energy"); 

    auto e_rng_arg = program.get<std::vector<double>>("--energy-range"); 
    const Bound<double> energy_range{ e_rng_arg[0], e_rng_arg[1] };

    const int npts_energy = program.get<int>("--n-pts-energy"); 

    auto ct_rng_arg = program.get<std::vector<double>>("--cos-theta-range");
    const Bound<double> cos_theta_range{ ct_rng_arg[0], ct_rng_arg[1] };

    const int npts_cos_theta = program.get<int>("--n-pts-cos-theta"); 

    //scan rates for energy & cos theta 
    const double dE   = (energy_range.max - energy_range.min)/((double)npts_energy-1);
    const double dCos = (cos_theta_range.max - cos_theta_range.min)/((double)npts_cos_theta-1);

    std::cout << "attempting integration..." << std::endl; 

    std::vector<std::thread> threads; 
    threads.reserve(npts_cos_theta*npts_energy);

    PolarFourVec P0; 
    P0.energy = beam_E; 
    P0.cos_theta = 1.; 
    P0.phi = 0.; 
    P0.mass2 = me*me; 

    std::vector<std::vector<double>> results(npts_cos_theta, std::vector<double>(npts_energy, 0.));
    std::vector<std::vector<double>> errors (npts_cos_theta, std::vector<double>(npts_energy, 0.));

    std::mutex save_mutex; 

    const size_t n_integration_steps = 3e8; 

    //scan over energy & cos(theta)
    double energy = energy_range.min; 
    for (int ie=0; ie<npts_energy; ie++) {

        double cos_theta = cos_theta_range.min; 
        for (int ic=0; ic<npts_cos_theta; ic++) {

            threads.emplace_back([P0, energy, cos_theta, n_integration_steps, ie,ic, &save_mutex, &results, &errors, process, rel_error]{

                PolarFourVec P1; 
                P1.energy = energy; 
                P1.cos_theta = cos_theta;
                P1.phi = 0.;  
                P1.mass2 = me*me; 

                std::vector<double> spectator_masses{}; 

                double amp = 1.; 

                int P1_ind; 

                DifferentialCSIntegrator* integrator; 
                
                switch (process) {

                    //____________________________________________________________________________
                    case kTrident_positron : {
                        integrator = new DifferentialCSIntegrator(processes::Factory::trident()); 
                        P1_ind = 3; 
                        spectator_masses.push_back(me);
                        spectator_masses.push_back(me);
                        break; 
                    }
                    
                    //____________________________________________________________________________
                    case kTrident_electron : {
                        integrator = new DifferentialCSIntegrator(processes::Factory::trident()); 
                        P1_ind = 2; 
                        spectator_masses.push_back(me);
                        spectator_masses.push_back(me);
                        break; 
                    }
                    
                    //____________________________________________________________________________
                    case kBH_photoproduction : {
                        integrator = new DifferentialCSIntegrator(processes::Factory::bh_photoproduction());
                        P1_ind = 1; 
                        spectator_masses.push_back(0.);
                        break; 
                    }

                    default : {
                        Error(__func__, "Process is unsupported.\n");
                        std::exit(1); 
                    }
                }
                auto result = integrator->Integrate(P0, P1, P1_ind, spectator_masses);

                save_mutex.lock(); 
                results[ie][ic] = result.val; 
                save_mutex.unlock(); 

                delete integrator; 
            });

            cos_theta += dCos; 
        }
        energy += dE; 
    }

    for (auto& thread : threads) thread.join(); 


    output_c_array(energy_range, cos_theta_range, npts_energy, npts_cos_theta, results, path_output);

    /*PolarFourVec P0, P1;
    
    P0.cos_theta = 1.0;
    P0.phi       = 0.; 
    P0.energy    = beam_E; 
    P0.mass2     = me*me; 
    
    //scan over theta & phi 


    P1.cos_theta = 0.95;
    P1.phi       = 0.; 
    P1.energy    = beam_E*0.95; 
    P1.mass2     = me*me; 
    
    int n_samples = 20;

    double amp{0.}, amp2{0.};
    for (int i=0; i<n_samples; i++) {
        
        double amp_i = EstimateDifferentialCS<3>(
            processes::bh_photoproduction,
            P0, 
            P1, 1, { 0. },
            3.e7, 
            Setting::kVerbose, Setting::kMISER
        );

        amp  += amp_i;
        amp2 += amp_i*amp_i; 
    }   
    
    double variance = amp2/((double)n_samples) - std::pow( amp/((double)n_samples), 2 );*/ 

    //std::printf("done. final amplitude: %.5e  +/-  %.4e\n", amp/((double)n_samples), std::sqrt(variance)); 
    
    return 0; 
}

void output_c_array(
    Bound<double> energy, Bound<double> cos_theta, 
    int npts_energy, int npts_cos_theta, 
    std::vector<std::vector<double>>& data, std::string path
)
{
    std::ofstream outfile(path, std::ios::out | std::ios::trunc);

    if (!outfile.is_open()) {
        std::cerr << "<output_c_array>: unable to open output file: '" << path << "'.\n"; 
        return; 
    }

    outfile << Form("const double energy_min = %+.8e;\n\n",energy.min); 
    outfile << Form("const double energy_max = %+.8e;\n\n",energy.max); 

    outfile << Form("const double cos_theta_min = %+.8e;\n\n",cos_theta.min); 
    outfile << Form("const double cos_theta_max = %+.8e;\n\n",cos_theta.max); 

    outfile << Form("const int npts_energy = %i;\n\n",npts_energy); 
    outfile << Form("const int npts_cos_theta = %i;\n\n",npts_cos_theta); 

    outfile << Form("const std::vector<std::vector<double>> cs_array = {\n"); 

    for (size_t i=0; i<data.size(); i++) {

        outfile << "    { ";
        
        const auto& row = data[i];
        for (size_t j=0; j<row.size(); j++) {

            outfile << Form("%+12.8e, ", row[j]); 
        }
        outfile << (i < data.size()-1 ? "},\n" : "}\n"); 

    }
    outfile << "};\n\n"; 
    
    outfile.close(); 
    std::cout << "saved output file '" << path << "'\n"; 
}

