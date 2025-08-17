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

    // Contat piece gets just constant as driving term
    isobar contact = solver.add_isobar<P_wave>(1, id::Contact, "Contact");

    // The projection function is given by our Deck loop 
    auto   Delta   = [&](complex sigma){return pi1::deck(t, m3pi*m3pi, sigma);};
    // Add isobar using the above function as our driving term
    isobar deck    = solver.add_isobar<P_wave>(Delta, 1, id::Deck, "Deck");

    
    std::vector<isobar> isos = solver.get_isobars();

    // -----------------------------------------------------------------------
    timer timer;
    plotter plotter;

    timer.start();

    double A = kinematics->A();
    double B = kinematics->B();
    double C = kinematics->C();
    double D = kinematics->D();
    double smin = A, smax = 3.;

    plot p1 = plotter.new_plot();
    p1.set_curve_points(1000);
    p1.set_xrange({smin, smax});
    p1.set_labels("#sigma   [GeV^{2}]", "#tilde{#it{F}}_{#alpha} (#sigma)");
    p1.add_horizontal(0);
    p1.set_legend(0.6, 0.4);
    p1.add_header("#minus #it{t}  = 0.1, #it{m}_{3#pi}^{2} = (1.4)^{2}");
    
    plot p2 = plotter.new_plot();
    p2.set_curve_points(1000);
    p2.set_xrange({smin, smax});
    p2.set_labels("#sigma   [GeV^{2}]", "#tilde{#it{F}}_{#Delta} (#it{t}, #it{m}_{3#pi}^{2} #; #sigma)");
    p2.add_horizontal(0);
    p2.set_legend(0.6, 0.4);
    p2.add_header("#minus #it{t}  = 0.1, #it{m}_{3#pi}^{2} = (1.4)^{2}");

    p1.add_curve({smin, smax}, [&](double s) { return std::real(contact->pinocchio_integral(0, s, isos)); }, solid(jpacColor::Blue,   "Real"));
    p1.add_curve({smin, smax}, [&](double s) { return std::imag(contact->pinocchio_integral(0, s, isos)); }, solid(jpacColor::Red,    "Imaginary"));
    p2.add_curve({smin, smax}, [&](double s) { return std::real(deck   ->pinocchio_integral(1, s, isos)); }, solid(jpacColor::Green,  "Real"));
    p2.add_curve({smin, smax}, [&](double s) { return std::imag(deck   ->pinocchio_integral(1, s, isos)); }, solid(jpacColor::Orange, "Imaginary"));
   
    solver.iterate(); timer.lap();
    p1.add_curve({smin, smax}, [&](double s) { return std::real(contact->pinocchio_integral(0, s, isos)); }, dashed(jpacColor::Blue));
    p1.add_curve({smin, smax}, [&](double s) { return std::imag(contact->pinocchio_integral(0, s, isos)); }, dashed(jpacColor::Red));
    p2.add_curve({smin, smax}, [&](double s) { return std::real(deck   ->pinocchio_integral(1, s, isos)); }, dashed(jpacColor::Green));
    p2.add_curve({smin, smax}, [&](double s) { return std::imag(deck   ->pinocchio_integral(1, s, isos)); }, dashed(jpacColor::Orange));
    solver.iterate(); timer.lap();
    p1.add_curve({smin, smax}, [&](double s) { return std::real(contact->pinocchio_integral(0, s, isos)); }, dotted(jpacColor::Blue));
    p1.add_curve({smin, smax}, [&](double s) { return std::imag(contact->pinocchio_integral(0, s, isos)); }, dotted(jpacColor::Red));
    p2.add_curve({smin, smax}, [&](double s) { return std::real(deck   ->pinocchio_integral(1, s, isos)); }, dotted(jpacColor::Green));
    p2.add_curve({smin, smax}, [&](double s) { return std::imag(deck   ->pinocchio_integral(1, s, isos)); }, dotted(jpacColor::Orange));

    // Save to file
    plotter.combine({2,1}, {p1,p2}, "angular_integrals.pdf");

    timer.stop();
    timer.print_elapsed();
};