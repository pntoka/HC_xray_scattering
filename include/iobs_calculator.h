#ifndef IOBS_CALCULATOR_H
#define IOBS_CALCULATOR_H

#include "calculations.h"
#include "structs.h"
#include <vector>

class IOBSCalculator {
public:
    IOBSCalculator() = default;
    
    std::vector<double> calculate(
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
        double polarizationDegree, bool coh, bool inc);
};

#endif