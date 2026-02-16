#include "iobs_calculator.h"
#include <omp.h>  // Add OpenMP header

std::vector<double> IOBSCalculator::calculate(
    const std::vector<double>& s_values,
    bool useA, double density, double absorptionCorrection,
    double sampleThickness, double transmission,
    bool useGradient, double g,
    bool useCorrAutoColl, double par_r, double par_delta, double par_l,
    double const1, double const2, bool useQ, double b, double k,
    double cno, std::vector<csp>* cspData,
    double cN, double cO, double cS, double cH,
    double dan, Enumerations::radiationType radiationType,
    double wavelength, bool useP, bool polarizedBeam,
    double polarizationDegree, bool coh, bool inc) {
    
    std::vector<double> results(s_values.size());

    // Parallel loop - each thread gets its own Calculations object
    #pragma omp parallel
    {
        Calculations calculations;  // Thread-local instance
        
        #pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < s_values.size(); ++i) {
            results[i] = calculations.iObs(
                useA, density, absorptionCorrection, sampleThickness, transmission,
                useGradient, g, useCorrAutoColl, par_r, par_delta, par_l,
                const1, const2, useQ, b, k, cno, cspData, cN, cO, cS, cH,
                dan, s_values[i], radiationType, wavelength, useP, polarizedBeam,
                polarizationDegree, coh, inc);
        }
    }

    return results;
}