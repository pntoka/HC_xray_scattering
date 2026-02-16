#include "iobs_calculator.h"
#include <omp.h>  // Add OpenMP header
#include "iobs_parameters.h"

std::vector<double> IOBSCalculator::calculate(
    const IOBSParameters& params,
    const std::vector<double>& s_values) {
    
    std::vector<csp> cspData = params.getCSPData();
    std::vector<double> results(s_values.size());

    #pragma omp parallel
    {
        Calculations calculations;
        
        #pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < s_values.size(); ++i) {
            results[i] = calculations.iObs(
                params.useA, params.density, params.absorptionCorrection,
                params.sampleThickness, params.transmission,
                params.useGradient, params.g,
                params.useCorrAutoColl, params.par_r, params.par_delta, params.par_l,
                params.const1, params.const2, params.useQ, params.b, params.k,
                params.cno, &cspData,
                params.cN, params.cO, params.cS, params.cH,
                params.dan, s_values[i], 
                static_cast<Enumerations::radiationType>(params.radiation),
                params.wavelength, params.useP, params.polarizedBeam,
                params.polarizationDegree, params.coh, params.inc);
        }
    }
    return results;
}