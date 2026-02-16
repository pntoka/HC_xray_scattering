#include <iostream>
#include <vector>
#include <iomanip>
#include <chrono>
#include "iobs_calculator.h"
#include "iobs_parameters.h"

void benchmark_calculator(int num_points, const IOBSParameters& params) {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Benchmarking IOBSCalculator with " << num_points << " points" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Generate s_values
    std::vector<double> s_values;
    double s_min = 0.1;
    double s_max = 1.0;
    for (int i = 0; i < num_points; ++i) {
        double s = s_min + (s_max - s_min) * i / (num_points - 1);
        s_values.push_back(s);
    }
    
    std::cout << "  s range: [" << s_values.front() << ", " << s_values.back() << "]" << std::endl;
    std::cout << "  Number of points: " << s_values.size() << std::endl;
    
    // Prepare calculator
    IOBSCalculator calculator;
    std::vector<csp> csp_data = params.getCSPData();
    
    // Start timing
    auto start = std::chrono::high_resolution_clock::now();
    
    // Calculate intensities
    std::vector<double> intensities = calculator.calculate(
        s_values,
        params.useA, params.density, 
        params.absorptionCorrection,
        params.sampleThickness, params.transmission,
        params.useGradient, params.g,
        params.useCorrAutoColl, params.par_r, 
        params.par_delta, params.par_l,
        params.const1, params.const2, 
        params.useQ, params.b, params.k,
        params.cno, &csp_data,
        params.cN, params.cO, 
        params.cS, params.cH,
        params.dan, 
        static_cast<Enumerations::radiationType>(params.radiation),
        params.wavelength, params.useP, 
        params.polarizedBeam, params.polarizationDegree,
        params.coh, params.inc
    );
    
    // End timing
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    
    // Display results
    std::cout << "\nResults:" << std::endl;
    std::cout << "  Execution time: " << std::fixed << std::setprecision(6) 
              << elapsed.count() << " seconds" << std::endl;
    std::cout << "  Time per point: " << std::scientific << std::setprecision(3)
              << (elapsed.count() / num_points) << " seconds" << std::endl;
    std::cout << "  Points per second: " << std::fixed << std::setprecision(2)
              << (num_points / elapsed.count()) << std::endl;
    
    std::cout << std::string(70, '=') << std::endl;
}

int main() {
    std::cout << std::string(70, '=') << std::endl;
    std::cout << "IOBSCalculator Performance Benchmark" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // Set up parameters (same as test_fitter.cpp ground truth)
    IOBSParameters params;
    params.mu = 5.91571;
    params.beta = 2.22345;
    params.a3 = 3.2758;
    params.da3 = 0.12257;
    params.sig3 = 1.26459;
    params.u3 = 0.13175;
    params.eta = 0.91002;
    params.nu = 4;
    params.alpha = 0.21636;
    params.lcc = 1.37635;
    params.sig1 = 0.20102;
    params.q = 0.14557;
    params.g = 0.61442;
    params.const1 = 131.28153832427495;
    params.const2 = 1.1564648406287676;
    params.k = 554.665935633803;
    params.cno = 0.0;
    params.cN = 0.0;
    params.cO = 0.0;
    params.cS = 0.0;
    params.cH = 0.0;
    params.dan = 0;
    params.useA = false;
    params.density = 2.2;
    params.absorptionCorrection = 1;
    params.sampleThickness = 0.1;
    params.transmission = false;
    params.useGradient = true;
    params.useCorrAutoColl = true;
    params.par_r = 24;
    params.par_delta = 0.5;
    params.par_l = 0.7;
    params.useQ = false;
    params.b = 0.002;
    params.radiation = 0;  // X_ray
    params.wavelength = 1.5418;
    params.useP = true;
    params.polarizedBeam = false;
    params.polarizationDegree = 0;
    params.coh = true;
    params.inc = true;
    
    std::cout << "\nParameters loaded (ground truth from test_fitter.cpp)" << std::endl;
    std::cout << "  mu: " << params.mu << std::endl;
    std::cout << "  beta: " << params.beta << std::endl;
    std::cout << "  a3: " << params.a3 << std::endl;
    
    // Run benchmarks
    benchmark_calculator(10, params);
    benchmark_calculator(1000, params);
    
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Benchmark completed successfully!" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    return 0;
}
