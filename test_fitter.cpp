#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include "iobs_calculator.h"
#include "iobs_fitter.h"
#include "iobs_parameters.h"

int main() {
    std::cout << std::string(70, '=') << std::endl;
    std::cout << "IOBSFitter C++ Test - Base Case 4 with Random Start Params" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    // ========================================================================
    // STEP 1: Generate synthetic "experimental" data using ground truth params
    // ========================================================================
    std::cout << "\n[1/4] Generating synthetic experimental data..." << std::endl;
    
    // Ground truth parameters from test_case.txt
    IOBSParameters ground_truth_params;
    ground_truth_params.mu = 5.91571;
    ground_truth_params.beta = 2.22345;
    ground_truth_params.a3 = 3.2758;
    ground_truth_params.da3 = 0.12257;
    ground_truth_params.sig3 = 1.26459;
    ground_truth_params.u3 = 0.13175;
    ground_truth_params.eta = 0.91002;
    ground_truth_params.nu = 4;
    ground_truth_params.alpha = 0.21636;
    ground_truth_params.lcc = 1.37635;
    ground_truth_params.sig1 = 0.20102;
    ground_truth_params.q = 0.14557;
    ground_truth_params.g = 0.61442;
    ground_truth_params.const1 = 131.28153832427495;
    ground_truth_params.const2 = 1.1564648406287676;
    ground_truth_params.k = 554.665935633803;
    ground_truth_params.cno = 0.0;
    ground_truth_params.cN = 0.0;
    ground_truth_params.cO = 0.0;
    ground_truth_params.cS = 0.0;
    ground_truth_params.cH = 0.0;
    ground_truth_params.dan = 0;
    ground_truth_params.useA = false;
    ground_truth_params.density = 2.2;
    ground_truth_params.absorptionCorrection = 1;
    ground_truth_params.sampleThickness = 0.1;
    ground_truth_params.transmission = false;
    ground_truth_params.useGradient = true;
    ground_truth_params.useCorrAutoColl = true;
    ground_truth_params.par_r = 24;
    ground_truth_params.par_delta = 0.5;
    ground_truth_params.par_l = 0.7;
    ground_truth_params.useQ = false;
    ground_truth_params.b = 0.002;
    ground_truth_params.radiation = 0;  // X_ray
    ground_truth_params.wavelength = 1.5418;
    ground_truth_params.useP = true;
    ground_truth_params.polarizedBeam = false;
    ground_truth_params.polarizationDegree = 0;
    ground_truth_params.coh = true;
    ground_truth_params.inc = true;
    
    // Generate s_values: 1000 points from 0.1 to 1.0
    std::vector<double> s_values;
    int num_points = 1000;
    double s_min = 0.1;
    double s_max = 1.0;
    for (int i = 0; i < num_points; ++i) {
        double s = s_min + (s_max - s_min) * i / (num_points - 1);
        s_values.push_back(s);
    }
    
    // Calculate synthetic data
    IOBSCalculator calculator;
    std::vector<csp> csp_data = ground_truth_params.getCSPData();
    std::vector<double> experimental_intensities = calculator.calculate(
        s_values,
        ground_truth_params.useA, ground_truth_params.density, 
        ground_truth_params.absorptionCorrection,
        ground_truth_params.sampleThickness, ground_truth_params.transmission,
        ground_truth_params.useGradient, ground_truth_params.g,
        ground_truth_params.useCorrAutoColl, ground_truth_params.par_r, 
        ground_truth_params.par_delta, ground_truth_params.par_l,
        ground_truth_params.const1, ground_truth_params.const2, 
        ground_truth_params.useQ, ground_truth_params.b, ground_truth_params.k,
        ground_truth_params.cno, &csp_data,
        ground_truth_params.cN, ground_truth_params.cO, 
        ground_truth_params.cS, ground_truth_params.cH,
        ground_truth_params.dan, 
        static_cast<Enumerations::radiationType>(ground_truth_params.radiation),
        ground_truth_params.wavelength, ground_truth_params.useP, 
        ground_truth_params.polarizedBeam, ground_truth_params.polarizationDegree,
        ground_truth_params.coh, ground_truth_params.inc
    );
    
    std::cout << "  Data points: " << s_values.size() << std::endl;
    std::cout << "  s range: [" << s_values.front() << ", " << s_values.back() << "]" << std::endl;
    
    // Calculate intensity range
    double min_intensity = *std::min_element(experimental_intensities.begin(), 
                                            experimental_intensities.end());
    double max_intensity = *std::max_element(experimental_intensities.begin(), 
                                            experimental_intensities.end());
    std::cout << "  Intensity range: [" << min_intensity << ", " 
              << max_intensity << "]" << std::endl;
    
    // ========================================================================
    // STEP 2: Set up initial parameters from start_params_random.json
    // ========================================================================
    std::cout << "\n[2/4] Loading initial parameters..." << std::endl;
    
    IOBSParameters initial_params;
    initial_params.mu = 4.907644737883213;
    initial_params.beta = 4.176026288464299;
    initial_params.a3 = 3.8352624821765984;
    initial_params.da3 = 0.03869587848572819;
    initial_params.sig3 = 1.8922978564833286;
    initial_params.u3 = 0.13613571929243937;
    initial_params.eta =  0.8313810617909404;
    initial_params.nu = 4;
    initial_params.alpha = 0.01;
    initial_params.lcc = 1.3566137940924743;
    initial_params.sig1 = 0.44618484963827143;
    initial_params.q =  0.5288907041806619;
    initial_params.g = -0.5133963034490847;
    initial_params.const1 = 404.5894542413238;
    initial_params.const2 = -2;
    initial_params.k = 1;
    initial_params.cno = 0.0;
    initial_params.cN = 0.0;
    initial_params.cO = 0.0;
    initial_params.cS = 0.0;
    initial_params.cH = 0.0;
    initial_params.dan = 0;
    initial_params.useA = false;
    initial_params.density = 2.2;
    initial_params.absorptionCorrection = 1;
    initial_params.sampleThickness = 0.1;
    initial_params.transmission = false;
    initial_params.useGradient = true;
    initial_params.useCorrAutoColl = true;
    initial_params.par_r = 24;
    initial_params.par_delta = 0.5;
    initial_params.par_l = 0.7;
    initial_params.useQ = false;
    initial_params.b = 0.002;
    initial_params.radiation = 0;  // X_ray
    initial_params.wavelength = 1.5418;
    initial_params.useP = true;
    initial_params.polarizedBeam = false;
    initial_params.polarizationDegree = 0;
    initial_params.coh = true;
    initial_params.inc = true;
    
    std::cout << "  Parameters loaded from: start_params_random.json (hardcoded)" << std::endl;
    std::cout << "  Initial mu: " << initial_params.mu << std::endl;
    std::cout << "  Initial beta: " << initial_params.beta << std::endl;
    std::cout << "  Initial a3: " << initial_params.a3 << std::endl;
    
    // ========================================================================
    // STEP 3: Configure fitter
    // ========================================================================
    std::cout << "\n[3/4] Configuring fitter..." << std::endl;
    
    IOBSFitter fitter;
    fitter.setMaxIterations(100);  // 10 iterations
    fitter.setNumThreads(1);
    fitter.setVerbose(true);
    fitter.setFunctionTolerance(1e-15);
    fitter.setParameterTolerance(1e-15);
    fitter.setGradientTolerance(1e-15);
    
    std::cout << "  Max iterations: 10" << std::endl;
    std::cout << "  Threads: 12" << std::endl;
    std::cout << "  Verbose: true" << std::endl;
    
    // ========================================================================
    // STEP 4: Run fitting
    // ========================================================================
    std::cout << "\n[4/4] Running multi-stage fit..." << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    FitResult result = fitter.fitMultiStage(s_values, experimental_intensities, initial_params);
    
    std::cout << std::string(70, '-') << std::endl;
    
    // ========================================================================
    // Display results
    // ========================================================================
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "RESULTS" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    std::cout << "Success: " << (result.success ? "true" : "false") << std::endl;
    std::cout << "Convergence: " << result.convergence_message << std::endl;
    
    std::cout << "\nPerformance:" << std::endl;
    std::cout << "  Iterations: " << result.num_iterations << std::endl;
    std::cout << "  Function evaluations: " << result.num_function_evaluations << std::endl;
    std::cout << "  Execution time: " << std::fixed << std::setprecision(2) 
              << result.execution_time_seconds << " seconds" << std::endl;
    
    std::cout << "\nGoodness of fit:" << std::endl;
    std::cout << "  Final cost: " << std::scientific << std::setprecision(6) 
              << result.final_cost << std::endl;
    std::cout << "  RSS: " << std::scientific << std::setprecision(6) 
              << result.rss << std::endl;
    std::cout << "  R²: " << std::fixed << std::setprecision(6) 
              << result.r_squared << std::endl;
    
    std::cout << "\nFinal parameters (unscaled):" << std::endl;
    std::cout << "  mu:     " << std::fixed << std::setprecision(6) 
              << result.final_params.mu << std::endl;
    std::cout << "  beta:   " << result.final_params.beta << std::endl;
    std::cout << "  a3:     " << result.final_params.a3 << std::endl;
    std::cout << "  da3:    " << result.final_params.da3 << std::endl;
    std::cout << "  sig3:   " << result.final_params.sig3 << std::endl;
    std::cout << "  eta:    " << result.final_params.eta << std::endl;
    std::cout << "  alpha:  " << result.final_params.alpha << std::endl;
    std::cout << "  sig1:   " << result.final_params.sig1 << std::endl;
    std::cout << "  lcc:    " << result.final_params.lcc << std::endl;
    std::cout << "  q:      " << result.final_params.q << std::endl;
    std::cout << "  g:      " << result.final_params.g << std::endl;
    std::cout << "  u3:     " << result.final_params.u3 << std::endl;
    std::cout << "  const1: " << result.final_params.const1 << std::endl;
    std::cout << "  const2: " << result.final_params.const2 << std::endl;
    std::cout << "  k:      " << result.final_params.k << std::endl;
    
    std::cout << "\nFinal parameters (scaled [0,1]):" << std::endl;
    std::vector<std::string> param_names = {
        "mu", "beta", "a3", "da3", "sig3", "eta", "alpha",
        "sig1", "lcc", "q", "g", "u3", "const1", "const2", "k"
    };
    for (size_t i = 0; i < param_names.size(); ++i) {
        std::cout << "  " << std::setw(7) << std::left << param_names[i] << ": "
                  << std::fixed << std::setprecision(6) 
                  << result.final_params_scaled[i] << std::endl;
    }
    
    std::cout << "\nGround truth parameters (for comparison):" << std::endl;
    std::cout << "  mu:     " << ground_truth_params.mu << std::endl;
    std::cout << "  beta:   " << ground_truth_params.beta << std::endl;
    std::cout << "  a3:     " << ground_truth_params.a3 << std::endl;
    std::cout << "  da3:    " << ground_truth_params.da3 << std::endl;
    std::cout << "  sig3:   " << ground_truth_params.sig3 << std::endl;
    std::cout << "  eta:    " << ground_truth_params.eta << std::endl;
    std::cout << "  alpha:  " << ground_truth_params.alpha << std::endl;
    std::cout << "  sig1:   " << ground_truth_params.sig1 << std::endl;
    std::cout << "  lcc:    " << ground_truth_params.lcc << std::endl;
    std::cout << "  q:      " << ground_truth_params.q << std::endl;
    std::cout << "  g:      " << ground_truth_params.g << std::endl;
    std::cout << "  u3:     " << ground_truth_params.u3 << std::endl;
    std::cout << "  const1: " << ground_truth_params.const1 << std::endl;
    std::cout << "  const2: " << ground_truth_params.const2 << std::endl;
    std::cout << "  k:      " << ground_truth_params.k << std::endl;
    
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "Test completed successfully!" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    return 0;
}
