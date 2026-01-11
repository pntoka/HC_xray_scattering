#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include "calculations.h"
#include "iobs_parameters.h"
#include "structs.h"
#include "iobs_calculator.h"

namespace py = pybind11;

PYBIND11_MODULE(iobs_ngc, m) {
    m.doc() = "Python module for calculating observed intnesity from scattering of non-graphitic carbon";
    
    // Binding the radiation type enum
    py::enum_<Enumerations::radiationType>(m, "RadiationType")
        .value("X_ray", Enumerations::X_ray)
        .value("neutron", Enumerations::neutron)
        .export_values();

    // Binding the CSP struct (coherent scattering parameters)
    py::class_<csp>(m, "CSP")
        .def(py::init<>())
        .def(py::init([](const py::dict& dict) {
            csp params;
            if (dict.contains("mu")) params.mu = dict["mu"].cast<double>();
            if (dict.contains("Nm")) params.Nm = dict["Nm"].cast<double>();
            if (dict.contains("a3min")) params.a3min = dict["a3min"].cast<double>();
            if (dict.contains("da3")) params.da3 = dict["da3"].cast<double>();
            if (dict.contains("sig3")) params.da3 = dict["sig3"].cast<double>();
            if (dict.contains("u3")) params.u3 = dict["u3"].cast<double>();
            if (dict.contains("eta")) params.concn = dict["eta"].cast<double>();
            if (dict.contains("nu")) params.nu = dict["nu"].cast<double>();
            if (dict.contains("lm")) params.lm = dict["lm"].cast<double>();
            if (dict.contains("sig1")) params.sig1 = dict["sig1"].cast<double>();
            if (dict.contains("lcc")) params.lcc = dict["lcc"].cast<double>();
            if (dict.contains("q")) params.q = dict["q"].cast<double>();
            if (dict.contains("concn")) params.concn = dict["concn"].cast<double>();
            return params;
        }))
    .def_readwrite("mu", &csp::mu)
    .def_readwrite("Nm", &csp::Nm)
    .def_readwrite("a3min", &csp::a3min)
    .def_readwrite("da3", &csp::da3)
    .def_readwrite("sig3", &csp::sig3)
    .def_readwrite("u3", &csp::u3)
    .def_readwrite("eta", &csp::eta)
    .def_readwrite("nu", &csp::nu)
    .def_readwrite("lm", &csp::lm)
    .def_readwrite("sig1", &csp::sig1)
    .def_readwrite("lcc", &csp::lcc)
    .def_readwrite("q", &csp::q)
    .def_readwrite("concn", &csp::concn);

    // Bind the iObs parameters class
    py::class_<IOBSParameters>(m, "IOBSParameters")
        .def(py::init<>())
        .def(py::init([](const py::dict& dict) {
            IOBSParameters params;
            // CSP parameters
            if (dict.contains("mu")) params.mu = dict["mu"].cast<double>();
            if (dict.contains("beta")) params.beta = dict["beta"].cast<double>();
            if (dict.contains("a3")) params.a3 = dict["a3"].cast<double>();
            if (dict.contains("da3")) params.da3 = dict["da3"].cast<double>();
            if (dict.contains("sig3")) params.sig3 = dict["sig3"].cast<double>();
            if (dict.contains("u3")) params.u3 = dict["u3"].cast<double>();
            if (dict.contains("eta")) params.eta = dict["eta"].cast<double>();
            if (dict.contains("nu")) params.nu = dict["nu"].cast<double>();
            if (dict.contains("alpha")) params.alpha = dict["alpha"].cast<double>();
            if (dict.contains("lcc")) params.lcc = dict["lcc"].cast<double>();
            if (dict.contains("sig1")) params.sig1 = dict["sig1"].cast<double>();
            if (dict.contains("q")) params.q = dict["q"].cast<double>();

            // Other parameters
            if (dict.contains("cH")) params.cH = dict["cH"].cast<double>();
            if (dict.contains("cN")) params.cN = dict["cN"].cast<double>();
            if (dict.contains("cO")) params.cO = dict["cO"].cast<double>();
            if (dict.contains("cS")) params.cS = dict["cS"].cast<double>();
            if (dict.contains("dan")) params.dan = dict["dan"].cast<double>();
            if (dict.contains("cno")) params.cno = dict["cno"].cast<double>();

        
            if (dict.contains("k")) params.k = dict["k"].cast<double>();
            if (dict.contains("const1")) params.const1 = dict["const1"].cast<double>();
            if (dict.contains("const2")) params.const2 = dict["const2"].cast<double>();
            if (dict.contains("useQ")) params.useQ = dict["useQ"].cast<bool>();
            if (dict.contains("b")) params.b = dict["b"].cast<double>();
            if (dict.contains("useA")) params.useA = dict["useA"].cast<bool>();
            if (dict.contains("density")) params.density = dict["density"].cast<double>();
            if (dict.contains("sampleThickness")) params.sampleThickness = dict["sampleThickness"].cast<double>();
            if (dict.contains("transmission")) params.transmission = dict["transmission"].cast<double>();
            if (dict.contains("absorptionCorrection")) params.absorptionCorrection = dict["absorptionCorrection"].cast<double>();

            // Beam parameters
            if (dict.contains("useP")) params.useP = dict["useP"].cast<bool>();
            if (dict.contains("polarizedBeam")) params.polarizedBeam = dict["polarizedBeam"].cast<bool>();
            if (dict.contains("polarizationDegree")) params.polarizationDegree = dict["polarizationDegree"].cast<double>();
            if (dict.contains("useGradient")) params.useGradient = dict["useGradient"].cast<bool>();
            if (dict.contains("g")) params.g = dict["g"].cast<double>();

            // Collimation parameters
            if (dict.contains("useCorrAutoColl")) params.useCorrAutoColl = dict["useCorrAutoColl"].cast<bool>();
            if (dict.contains("par_r")) params.par_r = dict["par_r"].cast<double>();
            if (dict.contains("par_delta")) params.par_delta = dict["par_delta"].cast<double>();
            if (dict.contains("par_l")) params.par_l = dict["par_l"].cast<double>();

            // Radiation parameters
            if (dict.contains("radiation")) params.radiation = dict["radiation"].cast<int>();
            if (dict.contains("wavelength")) params.wavelength = dict["wavelength"].cast<double>();

            // Calculation flags
            if (dict.contains("coh")) params.coh = dict["coh"].cast<bool>();
            if (dict.contains("inc")) params.inc = dict["inc"].cast<bool>();

            return params;
        }))
        // CSP parameters
        .def_readwrite("mu", &IOBSParameters::mu)
        .def_readwrite("beta", &IOBSParameters::beta)
        .def_readwrite("a3", &IOBSParameters::a3)
        .def_readwrite("da3", &IOBSParameters::da3)
        .def_readwrite("sig3", &IOBSParameters::sig3)
        .def_readwrite("u3", &IOBSParameters::u3)
        .def_readwrite("eta", &IOBSParameters::eta)
        .def_readwrite("nu", &IOBSParameters::nu)
        .def_readwrite("alpha", &IOBSParameters::alpha)
        .def_readwrite("lcc", &IOBSParameters::lcc)
        .def_readwrite("sig1", &IOBSParameters::sig1)
        .def_readwrite("q", &IOBSParameters::q)
        // Atomic parameters
        .def_readwrite("cH", &IOBSParameters::cH)
        .def_readwrite("cN", &IOBSParameters::cN)
        .def_readwrite("cO", &IOBSParameters::cO)
        .def_readwrite("cS", &IOBSParameters::cS)
        .def_readwrite("dan", &IOBSParameters::dan)
        .def_readwrite("cno", &IOBSParameters::cno)
        // Physical parameters
        .def_readwrite("k", &IOBSParameters::k)
        .def_readwrite("const1", &IOBSParameters::const1)
        .def_readwrite("const2", &IOBSParameters::const2)
        .def_readwrite("useQ", &IOBSParameters::useQ)
        .def_readwrite("b", &IOBSParameters::b)
        .def_readwrite("useA", &IOBSParameters::useA)
        .def_readwrite("density", &IOBSParameters::density)
        .def_readwrite("sampleThickness", &IOBSParameters::sampleThickness)
        .def_readwrite("transmission", &IOBSParameters::transmission)
        .def_readwrite("absorptionCorrection", &IOBSParameters::absorptionCorrection)
        // Beam parameters
        .def_readwrite("useP", &IOBSParameters::useP)
        .def_readwrite("polarizedBeam", &IOBSParameters::polarizedBeam)
        .def_readwrite("polarizationDegree", &IOBSParameters::polarizationDegree)
        .def_readwrite("useGradient", &IOBSParameters::useGradient)
        .def_readwrite("g", &IOBSParameters::g)
        // Collimation parameters
        .def_readwrite("useCorrAutoColl", &IOBSParameters::useCorrAutoColl)
        .def_readwrite("par_r", &IOBSParameters::par_r)
        .def_readwrite("par_delta", &IOBSParameters::par_delta)
        .def_readwrite("par_l", &IOBSParameters::par_l)
        // Radiation parameters
        .def_readwrite("radiation", &IOBSParameters::radiation)
        .def_readwrite("wavelength", &IOBSParameters::wavelength)
        // Calculation flags
        .def_readwrite("coh", &IOBSParameters::coh)
        .def_readwrite("inc", &IOBSParameters::inc)
        // Methods
        .def("get_csp_data", &IOBSParameters::getCSPData)
        .def("get_csp", &IOBSParameters::getCSP);


    // Binding the iObs calculator with the IOBSParameters
    py::class_<IOBSCalculator>(m, "IOBSCalculator")
        .def(py::init<>())
        .def("calculate", [](IOBSCalculator& self, const IOBSParameters& params, std::vector<double>& s_values) {
            std::vector<csp> cspData = params.getCSPData();
            return self.calculate(
                s_values,
                params.useA, params.density, params.absorptionCorrection,
                params.sampleThickness, params.transmission,
                params.useGradient, params.g,
                params.useCorrAutoColl, params.par_r, params.par_delta, params.par_l,
                params.const1, params.const2, params.useQ, params.b, params.k,
                params.cno, &cspData,
                params.cN, params.cO, params.cS, params.cH,
                params.dan, static_cast<Enumerations::radiationType>(params.radiation),
                params.wavelength, params.useP, params.polarizedBeam, params.polarizationDegree,
                params.coh, params.inc
            );
        });

    py::class_<Calculations>(m, "Calculations")
        .def(py::init<>())
        .def("deg_to_rad", &Calculations::degToRad)
        .def("rad_to_deg", &Calculations::radToDeg)
        .def("s_to_theta", &Calculations::sToTheta)
        .def("theta_to_s", &Calculations::thetaToS, "Function to convert theta to s")
        .def("recoil", &Calculations::Recoil)
        .def("dq", &Calculations::Dq)
        .def("qabs", &Calculations::Qabs)
        .def("g0", &Calculations::g0)
        .def("I_1D", &Calculations::I1D)
        .def("iObs", &Calculations::iObs)
        .def("Ifit", &Calculations::Ifit)
        .def("Iinc", &Calculations::Iinc)
        .def("fCpara", &Calculations::fCpara)
        .def("fH", &Calculations::fH)
        .def("fN", &Calculations::fN)
        .def("fO", &Calculations::fO)
        .def("fS", &Calculations::fS)
        .def("IcomC", &Calculations::IcomC)
        .def("IcomH", &Calculations::IcomH)
        .def("IcomN", &Calculations::IcomN)
        .def("IcomO", &Calculations::IcomO)
        .def("hypgeo2F1", &Calculations::hypgeo2F1)
        .def("Jhk", &Calculations::Jhk)
        .def("JhkXprefactor", &Calculations::JhkXprefactor)
        .def("g", &Calculations::g)
        .def("Jhk_prefactor", &Calculations::Jhk_prefactor)
        .def("F", &Calculations::F)
        .def("Iintra", &Calculations::Iintra)
        .def("IcomS", &Calculations::IcomS);

}