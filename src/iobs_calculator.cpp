#include "iobs_calculator.h"

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
    
    Calculations calculations;
    std::vector<double> results;
    results.reserve(s_values.size());

    for (double s : s_values) {
        // std::vector<csp> cspDataCopy = cspData;
        double result = calculations.iObs(
            useA, density, absorptionCorrection, sampleThickness, transmission,
            useGradient, g, useCorrAutoColl, par_r, par_delta, par_l,
            const1, const2, useQ, b, k, cno, cspData, cN, cO, cS, cH,
            dan, s, radiationType, wavelength, useP, polarizedBeam,
            polarizationDegree, coh, inc);
        results.push_back(result);
    }

    return results;
}