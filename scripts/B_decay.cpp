// Basis functions for B -> 3pi
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "timer.hpp"
#include "basis.hpp"
#include "plotter.hpp"
#include "solver.hpp"

// Identical to the structure of kaon 
#include "amplitudes/kaon.hpp"
#include "isobars/pseudoscalar.hpp"

void B_decay()
{
    using namespace iterateKT;

    // -----------------------------------------------------------------------
    
    // Set up general kinematics so everything knows masses
    kinematics kin = new_kinematics(5.27934, M_PION_PM);

    // Set up our amplitude 
    solver solver(kin);

    // Change default settings appropriate for the B decay
    settings sets = default_settings();
    sets._interpolation_points     = {150, 0, 50};
    sets._intermediate_energy      = 1.5;
    sets._cutoff                   = 4;          // Finite cutoff well below pth 
    sets._extra_cusp               = norm(0.98); // Add the KKbar cusp to special points to handle with care
    // sets._angular_integrator_depth = 0;

    // Variable step size settings
    // sets._iteration_rate_intercept = 0.1;
    // sets._iteration_rate_slope     = 0.05;

    // Use Madrid phases
    phase_args iso_0 = {"madrid/delta_00.dat", 1.69,  2, 0.33,  1.1}; // S0 -> 2π
    phase_args iso_1 = {"madrid/delta_11.dat", 1.69,  1,    1,  2.0}; // P1 -> π
    phase_args iso_2 = {"madrid/delta_02.dat", 9.99,  0,    1,  2.0}; // S2 -> 0
    sets._phase_shifts = { {id::I0_P1, iso_1}, 
                           {id::I1_S0, iso_0}, {id::I1_P1, iso_1}, {id::I1_S2, iso_2},
                           {id::I2_P1, iso_1}, {id::I2_S2, iso_2}};


    // Add all our isobars
    std::vector<uint> empty = {}; // Pass empty to isobars with no sub polynomials
    isobar F0 = solver.add_isobar<I1_S0>(3,        id::I1_S0, "F0", sets); 
    isobar F1 = solver.add_isobar<I1_P1>(1,        id::I1_P1, "F1", sets); 
    isobar F2 = solver.add_isobar<I1_S2>(empty, 1, id::I1_S2, "F2", sets);

    // Iterate N times
    int N = 3;
    solver.timed_iterate(N);
    // solver.timed_iterate(5);

    // -----------------------------------------------------------------------
    // Plot Results

    plotter plotter;

    double smin =  +0.0;
    // double smin = kin->sth();
    double smax =  +2.0;

    auto plot_basis = [&](isobar isobar, int i, std::string label)
    {
        plot p = plotter.new_plot();
        p.set_curve_points(200);
        p.set_legend(false);

        label += "(#it{s} + #it{i}#epsilon)";
        p.set_labels("#it{s} [GeV^{2}]", label);

        auto rF = [&](int j){return [j,isobar,i](double s){return real(isobar->basis_function(j, i, s+IEPS));}; };
        auto iF = [&](int j){return [j,isobar,i](double s){return imag(isobar->basis_function(j, i, s+IEPS));}; };

        for (int k = 0; k < N; k++)
        {
            p.add_curve ({smin, smax}, rF(k));
            p.add_dashed({smin, smax}, iF(k));
        }
        p.add_curve({smin, smax}, rF(N), solid (jpacColor::DarkGrey));
        p.add_curve({smin, smax}, iF(N), dashed(jpacColor::DarkGrey));
        return p;
    };

    std::vector<plot> F0s, F1s, F2s;
    F0s.emplace_back(plot_basis(F0, 0, "F_{0}^{a}"));
    F0s.emplace_back(plot_basis(F0, 1, "F_{0}^{b}"));
    F0s.emplace_back(plot_basis(F0, 2, "F_{0}^{c}"));
    F0s.emplace_back(plot_basis(F0, 3, "F_{0}^{d}"));
    plotter.combine({2,2}, F0s, "F0s.pdf");
    F1s.emplace_back(plot_basis(F1, 0, "F_{1}^{a}"));
    F1s.emplace_back(plot_basis(F1, 1, "F_{1}^{b}"));
    F1s.emplace_back(plot_basis(F1, 2, "F_{1}^{c}"));
    F1s.emplace_back(plot_basis(F1, 3, "F_{1}^{d}"));
    plotter.combine({2,2}, F1s, "F1s.pdf");
    F2s.emplace_back(plot_basis(F2, 0, "F_{2}^{a}"));
    F2s.emplace_back(plot_basis(F2, 1, "F_{2}^{b}"));
    F2s.emplace_back(plot_basis(F2, 2, "F_{2}^{c}"));
    F2s.emplace_back(plot_basis(F2, 3, "F_{2}^{d}"));
    plotter.combine({2,2}, F2s, "F2s.pdf");
};