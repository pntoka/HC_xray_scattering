#include "iobs_fitter.h"
#include "iobs_calculator.h"
#include "iobs_parameter_scaler.h"
#include <ceres/ceres.h>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <set>

// ============================================================================
// Helper: min-max scale y-data
// ============================================================================
struct YScaler {
    double y_min;
    double y_max;
    double range;

    YScaler(const std::vector<double>& y) {
        y_min = *std::min_element(y.begin(), y.end());
        y_max = *std::max_element(y.begin(), y.end());
        range = y_max - y_min;
        if (range < 1e-15) range = 1.0;  // degenerate guard
    }

    double scale(double val) const { return (val - y_min) / range; }
    double unscale(double val) const { return val * range + y_min; }

    std::vector<double> scaleVec(const std::vector<double>& v) const {
        std::vector<double> out(v.size());
        for (size_t i = 0; i < v.size(); ++i) out[i] = scale(v[i]);
        return out;
    }

    std::vector<double> unscaleVec(const std::vector<double>& v) const {
        std::vector<double> out(v.size());
        for (size_t i = 0; i < v.size(); ++i) out[i] = unscale(v[i]);
        return out;
    }
};

// ============================================================================
// Cost functor — 15 parameter blocks of size 1 each.
// Supports holding arbitrary subsets constant via Ceres problem API.
// ============================================================================
struct IOBSCostFunctor {
    const std::vector<double>& s_values;
    const std::vector<double>& intensities;
    IOBSParameters base_params;

    IOBSCostFunctor(const std::vector<double>& s,
                    const std::vector<double>& y,
                    const IOBSParameters& base)
        : s_values(s), intensities(y), base_params(base) {}

    bool operator()(const double* const* parameters, double* residuals) const {
        // Gather the 15 scaled values from individual parameter blocks
        double scaled[NUM_FIT_PARAMS];
        for (int i = 0; i < NUM_FIT_PARAMS; ++i) {
            scaled[i] = parameters[i][0];
        }

        // Unscale and evaluate
        IOBSParameters p = base_params;
        ParameterScaler::unscaleInto(scaled, p);

        IOBSCalculator calculator;
        std::vector<csp> csp_data = p.getCSPData();
        std::vector<double> model = calculator.calculate(
            s_values,
            p.useA, p.density, p.absorptionCorrection,
            p.sampleThickness, p.transmission,
            p.useGradient, p.g,
            p.useCorrAutoColl, p.par_r, p.par_delta, p.par_l,
            p.const1, p.const2, p.useQ, p.b, p.k,
            p.cno, &csp_data,
            p.cN, p.cO, p.cS, p.cH,
            p.dan,
            static_cast<Enumerations::radiationType>(p.radiation),
            p.wavelength, p.useP, p.polarizedBeam, p.polarizationDegree,
            p.coh, p.inc
        );

        for (size_t i = 0; i < s_values.size(); ++i) {
            double r = model[i] - intensities[i];
            residuals[i] = std::isfinite(r) ? r : 1e6;
        }
        return true;
    }
};

// ============================================================================
// Analytic-style cost functor that computes its own Jacobian via finite
// differences matching scipy's approach (absolute step = diff_step * max(1,|x|))
// ============================================================================
struct IOBSAnalyticCostFunctor : public ceres::CostFunction {
    const std::vector<double>& s_values;
    const std::vector<double>& intensities;  // original (unscaled) intensities
    IOBSParameters base_params;
    double diff_step;
    int num_free;
    std::vector<int> free_indices;
    YScaler y_scaler;

    IOBSAnalyticCostFunctor(const std::vector<double>& s,
                            const std::vector<double>& y,
                            const IOBSParameters& base,
                            const std::vector<int>& free_idx,
                            double step = 1e-4)
        : s_values(s), intensities(y), base_params(base),
          diff_step(step), free_indices(free_idx), y_scaler(y) {
        num_free = static_cast<int>(free_indices.size());
        mutable_parameter_block_sizes()->push_back(num_free);
        set_num_residuals(static_cast<int>(s.size()));
    }

    // Evaluate the model given a full scaled parameter vector, return Y-SCALED model
    std::vector<double> evaluateModelScaled(const double* all_scaled) const {
        IOBSParameters p = base_params;
        ParameterScaler::unscaleInto(all_scaled, p);

        IOBSCalculator calculator;
        std::vector<csp> csp_data = p.getCSPData();
        std::vector<double> raw = calculator.calculate(
            s_values,
            p.useA, p.density, p.absorptionCorrection,
            p.sampleThickness, p.transmission,
            p.useGradient, p.g,
            p.useCorrAutoColl, p.par_r, p.par_delta, p.par_l,
            p.const1, p.const2, p.useQ, p.b, p.k,
            p.cno, &csp_data,
            p.cN, p.cO, p.cS, p.cH,
            p.dan,
            static_cast<Enumerations::radiationType>(p.radiation),
            p.wavelength, p.useP, p.polarizedBeam, p.polarizationDegree,
            p.coh, p.inc
        );
        return y_scaler.scaleVec(raw);
    }

    bool Evaluate(const double* const* parameters,
                  double* residuals,
                  double** jacobians) const override {
        const double* free_params = parameters[0];
        const int m = static_cast<int>(s_values.size());

        // Build full scaled parameter vector
        auto base_scaled = ParameterScaler::scaleParams(base_params);
        double all_scaled[NUM_FIT_PARAMS];
        for (int i = 0; i < NUM_FIT_PARAMS; ++i) {
            all_scaled[i] = base_scaled[i];
        }
        for (int j = 0; j < num_free; ++j) {
            all_scaled[free_indices[j]] = free_params[j];
        }

        // Scale the target intensities
        std::vector<double> y_scaled = y_scaler.scaleVec(intensities);

        // Evaluate model in scaled y-space
        std::vector<double> model = evaluateModelScaled(all_scaled);
        for (int i = 0; i < m; ++i) {
            double r = model[i] - y_scaled[i];
            residuals[i] = std::isfinite(r) ? r : 1e6;
        }

        // Compute Jacobian if requested — central differences
        if (jacobians && jacobians[0]) {
            for (int j = 0; j < num_free; ++j) {
                int idx = free_indices[j];
                double x_val = all_scaled[idx];
                double h = diff_step * std::max(1.0, std::abs(x_val));

                double x_fwd = x_val + h;
                double x_bwd = x_val - h;
                if (x_fwd > 1.0) x_fwd = 1.0;
                if (x_bwd < 0.0) x_bwd = 0.0;
                double actual_2h = x_fwd - x_bwd;

                if (std::abs(actual_2h) < 1e-15) {
                    for (int i = 0; i < m; ++i) {
                        jacobians[0][i * num_free + j] = 0.0;
                    }
                    continue;
                }

                double perturbed_fwd[NUM_FIT_PARAMS];
                double perturbed_bwd[NUM_FIT_PARAMS];
                for (int k = 0; k < NUM_FIT_PARAMS; ++k) {
                    perturbed_fwd[k] = all_scaled[k];
                    perturbed_bwd[k] = all_scaled[k];
                }
                perturbed_fwd[idx] = x_fwd;
                perturbed_bwd[idx] = x_bwd;

                std::vector<double> model_fwd = evaluateModelScaled(perturbed_fwd);
                std::vector<double> model_bwd = evaluateModelScaled(perturbed_bwd);

                for (int i = 0; i < m; ++i) {
                    double deriv = (model_fwd[i] - model_bwd[i]) / actual_2h;
                    jacobians[0][i * num_free + j] = std::isfinite(deriv) ? deriv : 0.0;
                }
            }
        }
        return true;
    }
};

// ============================================================================
// Helper: resolve parameter names to indices
// ============================================================================
static std::vector<int> namesToIndices(const std::vector<std::string>& names) {
    std::vector<int> indices;
    indices.reserve(names.size());
    for (const auto& name : names) {
        auto it = std::find(PARAM_NAMES.begin(), PARAM_NAMES.end(), name);
        if (it != PARAM_NAMES.end()) {
            indices.push_back(static_cast<int>(std::distance(PARAM_NAMES.begin(), it)));
        }
    }
    return indices;
}

// ============================================================================
// Helper: compute fitted values from a scaled parameter vector
// ============================================================================
static std::vector<double> computeFitted(const std::vector<double>& s_values,
                                         const IOBSParameters& base_params,
                                         const std::vector<double>& params_vec) {
    IOBSParameters p = base_params;
    ParameterScaler::unscaleInto(params_vec.data(), p);

    IOBSCalculator calculator;
    std::vector<csp> csp_data = p.getCSPData();
    std::vector<double> result = calculator.calculate(
        s_values,
        p.useA, p.density, p.absorptionCorrection,
        p.sampleThickness, p.transmission,
        p.useGradient, p.g,
        p.useCorrAutoColl, p.par_r, p.par_delta, p.par_l,
        p.const1, p.const2, p.useQ, p.b, p.k,
        p.cno, &csp_data,
        p.cN, p.cO, p.cS, p.cH,
        p.dan,
        static_cast<Enumerations::radiationType>(p.radiation),
        p.wavelength, p.useP, p.polarizedBeam, p.polarizationDegree,
        p.coh, p.inc
    );
    
    // Debug: check for NaNs in final computation
    for (size_t i = 0; i < result.size(); ++i) {
        if (!std::isfinite(result[i])) {
            std::cerr << "WARNING: NaN/Inf in fitted value at index " << i 
                      << ", s=" << s_values[i] << std::endl;
            result[i] = 0.0;  // Or handle appropriately
        }
    }
    
    return result;
}

// ============================================================================
// IOBSFitter implementation
// ============================================================================

IOBSFitter::IOBSFitter()
    : max_iterations_(200),
      function_tolerance_(1e-15),
      parameter_tolerance_(1e-15),
      gradient_tolerance_(1e-15),
      num_threads_(8),
      verbose_(false) {}

void IOBSFitter::setMaxIterations(int max_iter)       { max_iterations_ = max_iter; }
void IOBSFitter::setFunctionTolerance(double ftol)     { function_tolerance_ = ftol; }
void IOBSFitter::setParameterTolerance(double ptol)    { parameter_tolerance_ = ptol; }
void IOBSFitter::setGradientTolerance(double gtol)     { gradient_tolerance_ = gtol; }
void IOBSFitter::setNumThreads(int num_threads)        { num_threads_ = num_threads; }
void IOBSFitter::setVerbose(bool verbose)              { verbose_ = verbose; }

// ============================================================================
// runStage — solve one stage with analytic Jacobian matching scipy TRF
// ============================================================================
IOBSFitter::StageSummary IOBSFitter::runStage(
        const std::vector<double>& s_values,
        const std::vector<double>& intensities,
        const IOBSParameters& base_params,
        std::vector<double>& params_vec,
        const std::vector<int>& free_indices) {

    // Extract only the free parameters into a compact vector
    const int num_free = static_cast<int>(free_indices.size());
    std::vector<double> free_params(num_free);
    for (int j = 0; j < num_free; ++j) {
        free_params[j] = params_vec[free_indices[j]];
    }

    // Build the Ceres problem with a SINGLE parameter block (free params only)
    // and our custom analytic cost function that does its own finite differences
    // matching scipy's step-size logic
    ceres::Problem problem;

    // Create base_params with current full scaled state baked in
    // so the cost functor knows the fixed parameter values
    IOBSParameters stage_base = base_params;
    // We need a version of base_params whose scaled representation matches params_vec
    // for the fixed parameters. We achieve this by passing params_vec state via
    // the base_params used for scaling. The analytic functor reads base_scaled
    // from scaleParams(base_params) and overrides free indices from the parameter block.
    // So we need base_params to produce the correct scaled values for fixed params.
    // The cleanest way: unscale the current params_vec into stage_base.
    ParameterScaler::unscaleInto(params_vec.data(), stage_base);

    auto* cost_function = new IOBSAnalyticCostFunctor(
        s_values, intensities, stage_base, free_indices, 1e-3);  // Use 1e-4 step

    problem.AddResidualBlock(cost_function, nullptr, free_params.data());

    // Set bounds [0, 1] for each free parameter
    for (int j = 0; j < num_free; ++j) {
        problem.SetParameterLowerBound(free_params.data(), j, 0.0);
        problem.SetParameterUpperBound(free_params.data(), j, 1.0);
    }

    // Configure solver — use LEVENBERG_MARQUARDT to match scipy TRF behavior
    ceres::Solver::Options options;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.linear_solver_type = ceres::DENSE_QR;
    options.max_num_iterations = max_iterations_;
    options.function_tolerance = function_tolerance_;
    options.parameter_tolerance = parameter_tolerance_;
    options.gradient_tolerance = gradient_tolerance_;
    options.num_threads = num_threads_;
    options.minimizer_progress_to_stdout = verbose_;

    // Smaller initial trust region — parameters are in [0,1]
    options.initial_trust_region_radius = 1.0;
    options.max_trust_region_radius = 10.0;

    // LM diagonal scaling
    options.min_lm_diagonal = 1e-6;
    options.max_lm_diagonal = 1e32;

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    // Write free params back into the full params_vec
    for (int j = 0; j < num_free; ++j) {
        params_vec[free_indices[j]] = free_params[j];
    }

    StageSummary ss;
    ss.iterations = static_cast<int>(summary.iterations.size());
    ss.function_evaluations = summary.num_residual_evaluations;
    ss.initial_cost = summary.initial_cost;
    ss.final_cost = summary.final_cost;
    ss.converged = (summary.termination_type == ceres::CONVERGENCE);
    ss.report = summary.BriefReport();
    return ss;
}

// ============================================================================
// buildResult — assemble a FitResult from final state
// ============================================================================
FitResult IOBSFitter::buildResult(const std::vector<double>& s_values,
                                  const std::vector<double>& intensities,
                                  const IOBSParameters& base_params,
                                  const std::vector<double>& params_vec,
                                  int total_iterations,
                                  int total_function_evals,
                                  bool all_converged,
                                  const std::string& message,
                                  double elapsed_seconds) {
    FitResult result;

    // Unscale final parameters
    result.final_params = base_params;
    ParameterScaler::unscaleInto(params_vec.data(), result.final_params);
    result.final_params_scaled.assign(params_vec.begin(), params_vec.end());

    // Compute fitted values (unscaled, in original intensity space)
    result.fitted_values = computeFitted(s_values, base_params, params_vec);

    // RSS and R² in ORIGINAL (unscaled) intensity space
    double ss_res = 0.0, ss_tot = 0.0, mean_y = 0.0;
    for (double y : intensities) mean_y += y;
    mean_y /= intensities.size();
    for (size_t i = 0; i < intensities.size(); ++i) {
        double d = intensities[i] - result.fitted_values[i];
        ss_res += d * d;
        double dm = intensities[i] - mean_y;
        ss_tot += dm * dm;
    }

    result.rss = ss_res;
    result.r_squared = (ss_tot > 0.0) ? 1.0 - (ss_res / ss_tot) : 0.0;
    result.final_cost = 0.5 * ss_res;
    result.num_iterations = total_iterations;
    result.num_function_evaluations = total_function_evals;
    result.success = all_converged;
    result.convergence_message = message;
    result.execution_time_seconds = elapsed_seconds;

    return result;
}

// ============================================================================
// fit — single stage, all parameters free (original API)
// ============================================================================
FitResult IOBSFitter::fit(const std::vector<double>& s_values,
                          const std::vector<double>& intensities,
                          const IOBSParameters& initial_params) {
    if (s_values.size() != intensities.size()) {
        FitResult r;
        r.convergence_message = "Error: s_values and intensities size mismatch";
        return r;
    }

    auto start = std::chrono::high_resolution_clock::now();

    auto scaled_arr = ParameterScaler::scaleParams(initial_params);
    std::vector<double> params_vec(scaled_arr.begin(), scaled_arr.end());

    // All indices free
    std::vector<int> all_free(NUM_FIT_PARAMS);
    for (int i = 0; i < NUM_FIT_PARAMS; ++i) all_free[i] = i;

    StageSummary ss = runStage(s_values, intensities, initial_params, params_vec, all_free);

    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();

    return buildResult(s_values, intensities, initial_params, params_vec,
                       ss.iterations, ss.function_evaluations,
                       ss.converged, ss.report, elapsed);
}

// ============================================================================
// fitMultiStage — default stages matching multi_fit_model_eval.py
// ============================================================================
FitResult IOBSFitter::fitMultiStage(const std::vector<double>& s_values,
                                    const std::vector<double>& intensities,
                                    const IOBSParameters& initial_params) {
    std::vector<FitStage> stages = {
        {"normalization", {"q", "k", "const1", "const2", "g"}},
        {"interlayer",    {"mu", "beta", "a3", "da3", "sig3", "eta"}},
        {"intralayer",    {"alpha", "sig1"}},
        {"lcc",           {"lcc"}},
        {"all_parameters", {"mu", "beta", "a3", "da3", "sig3", "eta", "alpha",
                            "sig1", "lcc", "q", "g", "u3", "const1", "const2", "k"}}
    };
    return fitMultiStage(s_values, intensities, initial_params, stages);
}

// ============================================================================
// fitMultiStage — custom stages
// ============================================================================
FitResult IOBSFitter::fitMultiStage(const std::vector<double>& s_values,
                                    const std::vector<double>& intensities,
                                    const IOBSParameters& initial_params,
                                    const std::vector<FitStage>& stages) {
    if (s_values.size() != intensities.size()) {
        FitResult r;
        r.convergence_message = "Error: s_values and intensities size mismatch";
        return r;
    }

    auto start = std::chrono::high_resolution_clock::now();

    // Scale initial parameters to [0, 1]
    auto scaled_arr = ParameterScaler::scaleParams(initial_params);
    std::vector<double> params_vec(scaled_arr.begin(), scaled_arr.end());

    int total_iterations = 0;
    int total_function_evals = 0;
    bool all_converged = true;
    std::string combined_message;

    for (size_t si = 0; si < stages.size(); ++si) {
        const auto& stage = stages[si];
        std::vector<int> free_indices = namesToIndices(stage.param_names);

        if (verbose_) {
            std::cout << "\n" << std::string(60, '=') << std::endl;
            std::cout << "Stage " << (si + 1) << "/" << stages.size()
                      << ": " << stage.name << std::endl;
            std::cout << "Fitting " << free_indices.size() << " parameters:";
            for (const auto& n : stage.param_names) std::cout << " " << n;
            std::cout << std::endl;
            std::cout << std::string(60, '-') << std::endl;
        }

        StageSummary ss = runStage(s_values, intensities, initial_params,
                                      params_vec, free_indices);

        total_iterations += ss.iterations;
        total_function_evals += ss.function_evaluations;
        if (!ss.converged) {
            all_converged = false;
        }

        if (verbose_) {
            std::cout << "Stage '" << stage.name << "': "
                      << ss.iterations << " iterations, "
                      << "cost " << std::scientific << ss.initial_cost
                      << " -> " << ss.final_cost << std::endl;

            // Print current parameter values for this stage
            for (int idx : free_indices) {
                double unscaled = ParameterScaler::unscale(
                    params_vec[idx], PARAM_BOUNDS[idx][0], PARAM_BOUNDS[idx][1]);
                std::cout << "  " << PARAM_NAMES[idx] << " = "
                          << std::fixed << std::setprecision(6) << unscaled
                          << " (scaled: " << params_vec[idx] << ")" << std::endl;
            }
        }

        combined_message += "Stage '" + stage.name + "': " +
                            ss.report + "\n";
    }

    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();

    return buildResult(s_values, intensities, initial_params, params_vec,
                       total_iterations, total_function_evals,
                       all_converged, combined_message, elapsed);
}
