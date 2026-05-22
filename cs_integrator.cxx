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

#include <TString.h> 

//#define EDCS_DEBUG
#include "EstimateDifferentialCS.hxx"

namespace {
    constexpr int D = 4; 

    constexpr double deg = 3.1415926536 / 180.; 

    constexpr double me = 0.501;
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
    
    std::string process; 
    program.add_argument("process")
        .required()
        .help("specific cross section to integrate")
        .choices("trident", "bh_photoproduction")
        .store_into(process); 

    std::string path_output;
    program.add_argument("path_output")
        .required()
        .help("Path to the output file destination")
        .store_into(path_output);

    //index of fixed output momentum 
    int fixed_momentum_index; 
    program.add_argument("-i", "--fixed-momentum-index")
        .default_value(1)
        .help("The index of the fixed output momentum to measure")
        .scan<'i', int>()
        .nargs(1);
        
    //beam energy 
    program.add_argument("--beam-energy")
        .help("beam energy of incident electron")
        .default_value(2200.)
        .scan<'g', double>()
        .nargs(1);

    //energy range
    program.add_argument("-e", "--energy-range")
        .help("energy range to scan for outgoing particle")
        .nargs(2)
        .scan<'g', double>()
        .default_value(std::vector<double>{1000., 2200.});
        
    //number of scan points in energy
    program.add_argument("--n-pts-energy")
        .help("Number of scan-points in the energy range")
        .nargs(1)
        .scan<'i', int>() 
        .default_value(25);

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
    

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1; 
    }

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

    const size_t n_integration_steps = 3e7; 

    //I found miser to be the most effective 
    Setting::MCStrategy integration_strategy = Setting::kMISER; 

    //scan over energy & cos(theta)
    double energy = energy_range.min; 
    for (int ie=0; ie<npts_energy; ie++) {

        double cos_theta = cos_theta_range.min; 
        for (int ic=0; ic<npts_cos_theta; ic++) {

            threads.emplace_back([P0, energy, cos_theta, integration_strategy, n_integration_steps, ie,ic, &save_mutex, &results, &errors]{

                PolarFourVec P1; 
                P1.energy = energy; 
                P1.cos_theta = cos_theta;
                P1.phi = 0.;  
                P1.mass2 = me*me; 

                double amp = EstimateDifferentialCS<3>(
                    processes::bh_photoproduction,
                    P0, 
                    P1, 1, { 0. },
                    n_integration_steps, 
                    Setting::kVerbose, Setting::kMISER
                );

                save_mutex.lock(); 
                results[ie][ic] = amp; 
                save_mutex.unlock(); 
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

