#ifndef IOBS_PARAMETERS_H
#define IOBS_PARAMETERS_H

#include "structs.h"
#include <vector>

class IOBSParameters {
public:
    // CSP parameters
    double mu = 0.0;
    double beta = 0.0;
    double a3 = 0.0;
    double da3 = 0.0;
    double sig3 = 0.0;
    double u3 = 0.0;
    double eta = 0.0;
    double nu = 0.0;
    double alpha = 0.0;
    double lcc = 0.0;
    double sig1 = 0.0;
    double q = 0.0;

    // Atomic parameters
    double cH = 0.0;
    double cN = 0.0;
    double cO = 0.0;
    double cS = 0.0;
    double dan = 0.0;
    double cno = 0.0;

    // Physical parameters
    double k = 0.0;
    double const1 = 0.0;
    double const2 = 0.0;
    bool useQ = false;
    double b = 0.0;
    bool useA = false;
    double density = 0.0;
    double sampleThickness = 0.0;
    double transmission = 0.0;
    double absorptionCorrection = 0.0;

    // Beam parameters
    bool useP = false;
    bool polarizedBeam = false;
    double polarizationDegree = 0.0;
    bool useGradient = false;
    double g = 0.0;

    // Collimation parameters
    bool useCorrAutoColl = false;
    double par_r = 0.0;
    double par_delta = 0.0;
    double par_l = 0.0;

    // Radiation parameters
    int radiation = 0;
    double wavelength = 0.0;

    // Calculation flags
    bool coh = true;
    bool inc = true;

    std::vector<csp> getCSPData() const {
        std::vector<csp> cspData(1);
        cspData[0].mu = mu;
        cspData[0].Nm = mu/beta;
        cspData[0].a3min = a3-da3;
        cspData[0].da3 = da3;
        cspData[0].sig3 = sig3;
        cspData[0].u3 = u3;
        cspData[0].eta = eta;
        cspData[0].nu = nu;
        cspData[0].lm = nu/alpha;
        cspData[0].sig1 = sig1;
        cspData[0].lcc = lcc;
        cspData[0].q = q;
        cspData[0].concn = 1;
        return cspData;
    }
    csp getCSP() const {
        csp data;
        data.mu = mu;
        data.Nm = mu/beta;
        data.a3min = a3-da3;
        data.da3 = da3;
        data.sig3 = sig3;
        data.u3 = u3;
        data.eta = eta;
        data.nu = nu;
        data.lm = nu/alpha;
        data.sig1 = sig1;
        data.lcc = lcc;
        data.q = q;
        data.concn = 1;
        return data;
    }
};

#endif