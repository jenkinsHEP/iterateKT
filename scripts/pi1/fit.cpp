// Fit KT amplitudes for π1(1600) decay with only P-wave to data from [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] - https://inspirehep.net/literature/1898933
// ------------------------------------------------------------------------------

#include <algorithm>
#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "plotter.hpp"
#include "fitter.hpp"

#include "amplitudes/pi1.hpp"
#include "isobars/pi1.hpp"
#include "COMPASS_pi1/fitter.hpp"
#include "COMPASS_pi1/data.hpp"

void fit()
{
    using namespace iterateKT;
    using iterateKT::complex;

    // -----------------------------------------------------------------------
    // Operating options

    int m3pibin    = 21;  // which m3pi bin to fit
    int Niter      = 10;  // Number of KT iterations

    // Import our data sets 
    data_set tbin0   = COMPASS::parse_JSON("tBin_0/dalitz_m3piBin_"+to_string(m3pibin)+"_tBin_0.json");
    tbin0._option    = option::tbin0;
    data_set tbin1   = COMPASS::parse_JSON("tBin_1/dalitz_m3piBin_"+to_string(m3pibin)+"_tBin_1.json");
    tbin1._option    = option::tbin1;
    data_set tbin2   = COMPASS::parse_JSON("tBin_2/dalitz_m3piBin_"+to_string(m3pibin)+"_tBin_2.json");
    tbin2._option    = option::tbin2;
    data_set tbin3   = COMPASS::parse_JSON("tBin_3/dalitz_m3piBin_"+to_string(m3pibin)+"_tBin_3.json");
    tbin3._option    = option::tbin3;
    
    // All should have the same m3pi
    double m3pi      = tbin0._extras["m3pi"];

    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution

    // Set up general kinematics so everything knows masses
    kinematics kin   = new_kinematics(m3pi, M_PION);
    
    // Set up our amplitude 
    // the pi1_tbins::initialize should handle adding the isobars
    amplitude amp    = new_amplitude<pi1_tbins>(kin);
    amp->set_name("π₁ → 3π");

    exit(1);
    // -----------------------------------------------------------------------
    // Set up fitter

    // 
    fitter<COMPASS::fit_all_tbins> fitter(amp, "Combined");
    fitter.set_tolerance(1E-5);
    fitter.set_print_level(3);

    // Add all bins
    fitter.add_data(tbin0);
    fitter.add_data(tbin1);
    fitter.add_data(tbin2);
    fitter.add_data(tbin3);

    // Add three t-slopes in addition to three subtraction coeffs
    fitter.add_extra_parameters(2);
    
    std::vector<std::string> labels = {"alpha", "delta", "b_alpha", "b_delta"};
    fitter.set_parameter_labels(labels);
    fitter.fix_argument("alpha", 0.); 
    fitter.make_real("b_alpha"); 
    fitter.make_real("b_delta"); 

    std::vector<complex> initial_guess  = {1., 1., 0., 0.};
    fitter.do_fit(initial_guess);
};