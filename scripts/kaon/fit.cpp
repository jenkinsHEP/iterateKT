// Do a fit of K+ -> π+ π+ π- decay width and dalitz plot parameters
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#include "utilities.hpp"
#include "fitter.hpp"

#include "kaon_fitter.hpp"
#include "amplitudes/kaon.hpp"
#include "isobars/pseudoscalar.hpp"

void fit()
{
    using namespace iterateKT;

    // -----------------------------------------------------------------------
    
    // Set up general kinematics so everything knows masses
    // Assume masses are given in terms of pion mass
    kinematics kin = new_kinematics(M_KAON_PM, M_PION_PM);

    // Set up our amplitude 
    amplitude amp = new_amplitude<charged_kaon>(kin, "K⁺ → π⁺π⁺π⁻");
    
    // Add all the isobars, note the order they are added will be the order
    // the basis functions are generated
    std::vector<uint> empty = {}; // Pass empty to isobars with no sub polynomials
    isobar F0 = amp->add_isobar<I1_S0>(2,        id::I1_S0, "F0"); 
    isobar F1 = amp->add_isobar<I1_P1>(1,        id::I1_P1, "F1"); 
    isobar F2 = amp->add_isobar<I1_S2>(empty, 1, id::I1_S2, "F2");
    isobar H1 = amp->add_isobar<I2_P1>(1,        id::I2_P1, "H1"); 
    isobar H2 = amp->add_isobar<I2_S2>(empty, 1, id::I2_S2, "H2");

    // Path to precalculated isobar files
    std::string path   = "/scripts/kaon/basis_functions/basis_";
    // Import everything 
    for (auto iso : amp->get_isobars()) iso->import_iteration<11>(path+iso->name()+".dat");

    // -----------------------------------------------------------------------
    // Set up data

    double gam, dgam, g, dg, h, dh, k, dk;

    // K+ -> pi+ pi+ pi- width and dalitz parameters
    gam =  2.9590,  dgam = 218E-4;
    g   = -0.21134, dg   = 17E-5;
    h   =  0.0185,  dh   = 4E-4;
    k   = -0.00463, dk   = 14E-5;

    data_set PtoPPM;
    PtoPPM._id     = "K+ -> pi+ pi+ pi-";
    PtoPPM._z      = {gam,   g,  h,  k};
    PtoPPM._dz     = {dgam, dg, dh, dk};
    PtoPPM._N      = 4;
    PtoPPM._option = option::P_ppm;
    PtoPPM._type   = kaon::fit::kAll;

    // K+ -> pi0 pi0 pi+ dalitz parameters
    g   =  0.626,  dg   = 7E-3;
    h   =  0.052,  dh   = 8E-3;
    k   =  0.0054, dk   = 35E-4;

    data_set PtoZZP;
    PtoZZP._id     = "K+ -> pi0 pi0 pi+";
    PtoZZP._z      = { g,  h,  k};
    PtoZZP._dz     = {dg, dh, dk};
    PtoZZP._N      = 3;
    PtoZZP._option = option::P_zzp;
    PtoZZP._type   = kaon::fit::kDalitz;
    
    // -----------------------------------------------------------------------
    // Set up fitter

    // std::vector<iterateKT::complex> pars=  {1210.2084, -5965.0633, -1071.8057, 9163.4163};
    // auto ppars = kaon::fit::process_fitter_parameters(pars, amp);
    // amp->set_parameters(ppars);
    // amp->set_option(option::P_ppm);
    // // auto dpars = amp->get_dalitz_parameters(1E-6);
    // // print("Width =", amp->width());
    // // print("g =", dpars[0]);
    // // print("h =", dpars[1]);
    // // print("k =", dpars[3]);
    // print(ppars);
    // // amp->set_option(option::P_zzp);
    // // dpars = amp->get_dalitz_parameters(1E-3);
    // // print("Width =", amp->width());
    // // print(dpars[0], dpars[1], dpars[3]);
    // exit(1);

    fitter<kaon::fit> fitter(amp);

    fitter.set_print_level(1);
    
    // Add data from above
    fitter.add_data(PtoPPM);
    // fitter.add_data(PtoZZP);

    // Set up parameters (all real)
    std::vector<std::string> labels = {"alpha", "beta", "gamma", "zeta"};
    fitter.set_parameter_labels(labels);
    for (auto par : labels) fitter.make_real(par);

    std::vector<iterateKT::complex> initial_guess = {1537.7541, -8764.3614, 264.82228, 9930.3531};
    fitter.do_fit(initial_guess);
};