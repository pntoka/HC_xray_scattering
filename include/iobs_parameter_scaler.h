#ifndef IOBS_PARAMETER_SCALER_H
#define IOBS_PARAMETER_SCALER_H

#include "iobs_parameters.h"
#include <array>
#include <string>
#include <vector>

// Number of fittable parameters
static constexpr int NUM_FIT_PARAMS = 15;

// Parameter names in order
static const std::array<std::string, NUM_FIT_PARAMS> PARAM_NAMES = {
    "mu", "beta", "a3", "da3", "sig3", "eta", "alpha",
    "sig1", "lcc", "q", "g", "u3", "const1", "const2", "k"
};

// Parameter bounds: {min, max} matching the Python param_names_scale
static constexpr std::array<std::array<double, 2>, NUM_FIT_PARAMS> PARAM_BOUNDS = {{
    {1.2, 8.0},       // mu
    {0.01, 8.0},      // beta
    {2.7, 4.0},       // a3
    {0.01, 1.0},      // da3
    {0.01, 2.0},      // sig3
    {0.6, 1.0},       // eta
    {0.01, 1.5},      // alpha
    {0.0, 0.6},       // sig1
    {1.2, 1.5},       // lcc
    {0.0, 0.75},      // q
    {-1.0, 1.0},      // g
    {0.0, 0.4},       // u3
    {-500.0, 750.0},  // const1
    {-2.0, 5.0},      // const2
    {1.0, 1000.0}     // k
}};

class ParameterScaler {
public:
    // Scale a value from [min, max] to [0, 1]
    static double scale(double value, double min_val, double max_val) {
        return (value - min_val) / (max_val - min_val);
    }

    // Unscale a value from [0, 1] to [min, max]
    static double unscale(double scaled, double min_val, double max_val) {
        return scaled * (max_val - min_val) + min_val;
    }

    // Scale IOBSParameters to a [0,1] array
    static std::array<double, NUM_FIT_PARAMS> scaleParams(const IOBSParameters& params) {
        std::array<double, NUM_FIT_PARAMS> scaled;
        std::array<double, NUM_FIT_PARAMS> values = getParamValues(params);
        for (int i = 0; i < NUM_FIT_PARAMS; ++i) {
            scaled[i] = scale(values[i], PARAM_BOUNDS[i][0], PARAM_BOUNDS[i][1]);
        }
        return scaled;
    }

    // Unscale a [0,1] array back into IOBSParameters (modifies only fit params)
    static void unscaleInto(const double* scaled, IOBSParameters& params) {
        params.mu     = unscale(scaled[0],  PARAM_BOUNDS[0][0],  PARAM_BOUNDS[0][1]);
        params.beta   = unscale(scaled[1],  PARAM_BOUNDS[1][0],  PARAM_BOUNDS[1][1]);
        params.a3     = unscale(scaled[2],  PARAM_BOUNDS[2][0],  PARAM_BOUNDS[2][1]);
        params.da3    = unscale(scaled[3],  PARAM_BOUNDS[3][0],  PARAM_BOUNDS[3][1]);
        params.sig3   = unscale(scaled[4],  PARAM_BOUNDS[4][0],  PARAM_BOUNDS[4][1]);
        params.eta    = unscale(scaled[5],  PARAM_BOUNDS[5][0],  PARAM_BOUNDS[5][1]);
        params.alpha  = unscale(scaled[6],  PARAM_BOUNDS[6][0],  PARAM_BOUNDS[6][1]);
        params.sig1   = unscale(scaled[7],  PARAM_BOUNDS[7][0],  PARAM_BOUNDS[7][1]);
        params.lcc    = unscale(scaled[8],  PARAM_BOUNDS[8][0],  PARAM_BOUNDS[8][1]);
        params.q      = unscale(scaled[9],  PARAM_BOUNDS[9][0],  PARAM_BOUNDS[9][1]);
        params.g      = unscale(scaled[10], PARAM_BOUNDS[10][0], PARAM_BOUNDS[10][1]);
        params.u3     = unscale(scaled[11], PARAM_BOUNDS[11][0], PARAM_BOUNDS[11][1]);
        params.const1 = unscale(scaled[12], PARAM_BOUNDS[12][0], PARAM_BOUNDS[12][1]);
        params.const2 = unscale(scaled[13], PARAM_BOUNDS[13][0], PARAM_BOUNDS[13][1]);
        params.k      = unscale(scaled[14], PARAM_BOUNDS[14][0], PARAM_BOUNDS[14][1]);
    }

private:
    static std::array<double, NUM_FIT_PARAMS> getParamValues(const IOBSParameters& p) {
        return {p.mu, p.beta, p.a3, p.da3, p.sig3, p.eta, p.alpha,
                p.sig1, p.lcc, p.q, p.g, p.u3, p.const1, p.const2, p.k};
    }
};

#endif
