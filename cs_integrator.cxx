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
#include <cstdlib> 

#include <TString.h> 
#include <TError.h> 
#include <TStopwatch.h> 

//#define EDCS_DEBUG
#include "DifferentialCSIntegrator.hxx"

namespace {
    constexpr int D = 4; 

    constexpr double deg = 3.1415926536 / 180.; 

    constexpr double me = 0.501;

    enum EProcessName {
        kTrident_electron =1,
        kTrident_positron, 
        kBH_photoproduction, 
        kElastic
    };
}

void output_c_array(
    Bound<double> energy, Bound<double> cos_theta, 
    int npts_energy, int npts_cos_theta, 
    const std::vector<double>& data, std::string path
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
            .choices("trident_electron", "trident_positron", "bh_photoproduction", "elastic");
            
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
        {"trident_positron",    kTrident_positron}, 
        {"elastic",             kElastic}
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

    PolarFourVec P0; 
    P0.energy = beam_E; 
    P0.cos_theta = 1.; 
    P0.phi = 0.; 
    P0.mass2 = me*me; 

    std::vector<double> results(npts_cos_theta*npts_energy, 0.);
    std::vector<double> errors (npts_cos_theta*npts_energy, 0.);

    std::mutex save_mutex; 

    const size_t n_integration_steps = 3e8; 

    struct phase_space_point_t { double energy,cos_theta; int ie,ic; };

    //__________________________________________________________________________________________________________
    auto perform_integration = [P0, n_integration_steps, &save_mutex, &results, &errors, process, rel_error, npts_energy]
        (const phase_space_point_t pt)
    {
        PolarFourVec P1; 
        P1.energy    = pt.energy; 
        P1.cos_theta = pt.cos_theta;
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

            //____________________________________________________________________________
            case kElastic : {
                integrator = new DifferentialCSIntegrator(processes::Factory::elastic());
                P1_ind = 1; 
                break; 
            }

            default : {
                Error(__func__, "Process is unsupported.\n");
                std::exit(1); 
            }
        }
        //inform that this integration has started 
        save_mutex.lock(); 
        printf("starting: energy: %6.1f cos(theta): %7.5f\n",
            pt.energy, pt.cos_theta
        ); std::cout << std::flush; 
        save_mutex.unlock(); 

        integrator->SetOptions(Setting::kNone);

        TStopwatch timer; 
        auto result = integrator->Integrate(P0, P1, P1_ind, spectator_masses);
        auto elapsed = timer.RealTime(); 

        //inform that this integration has ended 
        save_mutex.lock(); 
        results[pt.ie*npts_energy + pt.ic] = result.val; 
        printf("> time: %8.3fs energy: %6.1f cos(theta): %7.5f  result: %.6e +/- %.3e\n",
            elapsed, pt.energy, pt.cos_theta, result.val, result.error
        ); std::cout << std::flush; 
        save_mutex.unlock(); 

        delete integrator; 
    };
    //__________________________________________________________________________________________________________

    std::vector<phase_space_point_t> points; points.reserve(npts_energy*npts_cos_theta);

    //scan over energy & cos(theta)
    double energy = energy_range.min; 
    for (int ie=0; ie<npts_energy; ie++) {

        double cos_theta = cos_theta_range.min; 
        for (int ic=0; ic<npts_cos_theta; ic++) {

            points.emplace_back(energy, cos_theta, ie, ic);

            cos_theta += dCos; 
        }
        energy += dE; 
    }

    //now launch our threads. 

    size_t end=0; 
    const size_t n_tasks = points.size(); 

    const char* SLURM_CPUS_PER_TASK = std::getenv("SLURM_CPUS_PER_TASK");
    size_t max_n_threads; 
    if (!SLURM_CPUS_PER_TASK) {
        max_n_threads = std::thread::hardware_concurrency(); 
        printf("hardware concurrency: %zi\n", max_n_threads);
    } else {
        max_n_threads = (size_t)std::stoi( std::string{SLURM_CPUS_PER_TASK} );
        printf("slurm cpus per task: %zi\n", max_n_threads);
    }
     
    const size_t n_threads = std::min( max_n_threads, n_tasks ); 
    std::vector<std::thread> threads; threads.reserve(n_threads);

    //divide up the integrations evenly among each thread
    for (size_t t=0; t<n_threads; t++) {
        
        size_t start = end; 
        end += (n_tasks / n_threads) + (t < n_tasks % n_threads ? 1 : 0);

        threads.emplace_back([t,n_threads, start, end, &perform_integration, &points, &save_mutex]
        {
            for (size_t ii=start; ii<end; ii++) {
                save_mutex.lock();
                std::printf("worker thread %2zi/%zi processing task %zi/%zi\n", t+1, n_threads, ii-start+1, end-start);
                save_mutex.unlock(); 

                perform_integration(points[ii]); 
            }
        });
    }

    for (auto& thread : threads) thread.join(); 

    output_c_array(energy_range, cos_theta_range, npts_energy, npts_cos_theta, results, path_output);

    return 0; 
}

void output_c_array(
    Bound<double> energy, Bound<double> cos_theta, 
    int npts_energy, int npts_cos_theta, 
    const std::vector<double>& data, std::string path
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

    outfile << "// row is energy, col is cos (theta). ends are inclusive of bounds.\n";
    outfile << Form("const std::vector<double> cs_array{\n"); 
    
    for (size_t i=0; i<npts_energy; i++) {
        
        outfile << "    ";
        for (size_t j=0; j<npts_cos_theta; j++) {

            outfile << Form("%+12.8e", data[i*npts_energy + j]);
            if (j < npts_cos_theta-1) {
                outfile << ", ";
            } else {
                if (i < npts_energy-1) outfile << ", ";  
            } 
        }
        outfile << "\n";
    }
    outfile << "};\n\n"; 
    
    outfile.close(); 
    std::cout << "saved output file '" << path << "'\n"; 
}

