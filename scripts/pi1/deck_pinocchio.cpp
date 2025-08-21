// Plot the inhomogeneity of the KT equations when using the Deck seed term
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] - https://arxiv.org/abs/2212.11767
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

void deck_pinocchio()
{
    using namespace iterateKT;
    using iterateKT::complex;

    uint   N    = 4;    // Number of iterations
    double t    = -0.1; // Production t
    double m3pi = 1.40; // 3pi invariant mass

    kinematics kinematics = new_kinematics(m3pi, M_PION);
    solver solver(kinematics);

    // Contact piece gets just constant as driving term
    auto   contact = [&](complex sigma){return 1.;};
    // The projection function is given by our Deck loop 
    auto   Delta   = [&](complex sigma){return pi1::deck(t, m3pi*m3pi, sigma);};
    // Add isobar using the above function as our driving term
    isobar pwave   = solver.add_isobar<P_wave>({contact, Delta}, 1, id::Deck, "Deck");

    std::vector<isobar> isos = solver.get_isobars();

    // -----------------------------------------------------------------------
    plotter plotter;

    double A = kinematics->A();
    double B = kinematics->B();
    double C = kinematics->C();
    double D = kinematics->D();
    double smin = A, smax = 3.5;

    plot p1 = plotter.new_plot();
    p1.set_curve_points(1000);
    p1.set_xrange({smin, smax});
    p1.add_header("#minus#it{t}  = 0.1, #it{m}_{3#pi}^{2} = (1.4)^{2}");
    p1.set_labels("#sigma   [GeV^{2}]", "#kappa^{3} #tilde{#it{F}} (#it{t}, #it{m}^{2}_{3#pi} #; #sigma)");
    p1.add_horizontal(0);
    p1.add_vertical(D);
    p1.shade_region({A,C});
    p1.set_legend(0.225, 0.75);

    solver.timed_iterate(4);
    p1.add_curve({smin, smax}, [&](double s) { return std::real(deck->pinocchio_integral(0, s, isos)); }, solid(jpacColor::Blue,   "#alpha"));
    p1.add_curve({smin, smax}, [&](double s) { return std::imag(deck->pinocchio_integral(0, s, isos)); }, dashed(jpacColor::Blue));
    p1.add_curve({smin, smax}, [&](double s) { return std::real(deck->pinocchio_integral(1, s, isos)); }, solid(jpacColor::Red,  "#Delta"));
    p1.add_curve({smin, smax}, [&](double s) { return std::imag(deck->pinocchio_integral(1, s, isos)); }, dashed(jpacColor::Red));
   
    // Save to file
    p1.save("angular_integrals.pdf");
};