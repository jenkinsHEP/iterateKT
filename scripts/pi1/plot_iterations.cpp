// Plot the iterations of the KT equations for pi_1 -> 3pi using a variety
// of production functions
// 
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#include "kinematics.hpp"
#include "solver.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "timer.hpp"
#include "plotter.hpp"
#include "settings.hpp"

#include "amplitudes/pi1.hpp"
#include "isobars/pi1.hpp"

void plot_iterations()
{
    using namespace iterateKT;
    using iterateKT::complex;
    
    // -----------------------------------------------------------------------
    // Operating options

    uint   N                =    9; // Number of iterations
    double t                = -0.1; // Production t
    double m3pi             = 1.40; // 3pi invariant mass
    uint   asymptotic_power =    3; // sigma^N behavior of the isobar
    
    auto   constant = [&](complex sigma){return 1.;};
    auto   linear   = [&](complex sigma){return sigma;};
    auto   omega    = [&](complex sigma){return sigma/(sigma-norm(M_OMEGA)+I*M_OMEGA*8.68E-3); };
    auto   deck     = [&](complex sigma){return pi1::deck(t, m3pi*m3pi, sigma);};

    // Which production function to use
    std::function<complex(complex)> driving_term = linear;

    // -----------------------------------------------------------------------
    // Set up, shouldnt need to change anything below this line

    kinematics kinematics = new_kinematics(m3pi, M_PION);
    solver solver(kinematics);
    
    settings sets = default_settings();
    sets._extra_cusp = norm(M_OMEGA);
    
    isobar pwave = solver.add_isobar<P_wave>(driving_term, asymptotic_power, id::P_wave, "P_wave", sets);
    
    // -----------------------------------------------------------------------
    // Plot the 0th, 1st, and last iterations

    timer timer;
    plotter plotter;

    timer.start();

    std::array<double,2> bounds = {0, 2.6};
    double A = kinematics->A();
    double B = kinematics->B();
    double C = kinematics->C();
    double D = kinematics->D();

    plot p1 = plotter.new_plot();
    p1.set_curve_points(1000);
    p1.set_xrange(bounds);
    p1.set_labels("#sigma   [GeV^{2}]", "#it{F}(#it{t}, #it{m}_{3#pi}^{2} #; #sigma + #it{i}#epsilon)");
    p1.add_horizontal(0);
    p1.shade_region({A,C});
    p1.set_legend(0.7, 0.6);
    p1.add_curve(bounds, [&](double s) { return std::real(pwave->basis_function(0, s+IEPS)); }, dotted(jpacColor::Blue));
    p1.add_curve(bounds, [&](double s) { return std::imag(pwave->basis_function(0, s+IEPS)); }, dotted(jpacColor::Red));

    for (int i = 1; i <= N; i++)
    {
        solver.iterate();
        if (i == 1)
        {
            p1.add_curve(bounds, [&](double s) { return std::real(pwave->basis_function(0, s+IEPS)); }, dashed(jpacColor::Blue));
            p1.add_curve(bounds, [&](double s) { return std::imag(pwave->basis_function(0, s+IEPS)); }, dashed(jpacColor::Red));
        }
        if (i == N)
        {
            p1.add_curve(bounds, [&](double s) { return std::real(pwave->basis_function(0, s+IEPS)); }, solid(jpacColor::Blue,  "Real"));
            p1.add_curve(bounds, [&](double s) { return std::imag(pwave->basis_function(0, s+IEPS)); }, solid(jpacColor::Red,  "Imaginary"));
        }
    };
    p1.save("iteration.pdf");

    timer.stop();
    timer.print_elapsed();
};