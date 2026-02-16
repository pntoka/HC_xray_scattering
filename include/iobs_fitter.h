#ifndef IOBS_FITTER_H
#define IOBS_FITTER_H

#include "iobs_parameters.h"
#include "iobs_parameter_scaler.h"
#include <string>
#include <vector>

struct FitStage {
    std::string name;
    std::vector<std::string> param_names;
};

struct FitResult {
    bool success = false;
    std::string convergence_message;
    
    IOBSParameters final_params;
    std::vector<double> final_params_scaled;
    std::vector<double> fitted_values;
    
    double final_cost = 0.0;
    double rss = 0.0;
    double r_squared = 0.0;
    int num_iterations = 0;
    int num_function_evaluations = 0;
    double execution_time_seconds = 0.0;
};

class IOBSFitter {
public:
    IOBSFitter();
    
    void setMaxIterations(int max_iter);
    void setFunctionTolerance(double ftol);
    void setParameterTolerance(double ptol);
    void setGradientTolerance(double gtol);
    void setNumThreads(int num_threads);
    void setVerbose(bool verbose);
    
    // Fit all parameters at once
    FitResult fit(const std::vector<double>& s_values,
                  const std::vector<double>& intensities,
                  const IOBSParameters& initial_params);

    // Multi-stage fitting with default stages (normalization, interlayer,
    // intralayer, lcc, all_parameters)
    FitResult fitMultiStage(const std::vector<double>& s_values,
                            const std::vector<double>& intensities,
                            const IOBSParameters& initial_params);

    // Multi-stage fitting with custom stages
    FitResult fitMultiStage(const std::vector<double>& s_values,
                            const std::vector<double>& intensities,
                            const IOBSParameters& initial_params,
                            const std::vector<FitStage>& stages);

private:
    int max_iterations_;
    double function_tolerance_;
    double parameter_tolerance_;
    double gradient_tolerance_;
    int num_threads_;
    bool verbose_;

    struct StageSummary {
        int iterations = 0;
        int function_evaluations = 0;
        double initial_cost = 0.0;
        double final_cost = 0.0;
        bool converged = false;
        std::string report;
    };

    // Run a single fitting stage, updating params_vec in place.
    // free_indices: which of the 15 parameters to optimise.
    StageSummary runStage(const std::vector<double>& s_values,
                          const std::vector<double>& intensities,
                          const IOBSParameters& base_params,
                          std::vector<double>& params_vec,
                          const std::vector<int>& free_indices);

    // Build final FitResult from a scaled parameter vector and accumulated summaries
    FitResult buildResult(const std::vector<double>& s_values,
                          const std::vector<double>& intensities,
                          const IOBSParameters& base_params,
                          const std::vector<double>& params_vec,
                          int total_iterations,
                          int total_function_evals,
                          bool all_converged,
                          const std::string& message,
                          double elapsed_seconds);
};

#endif
