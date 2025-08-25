//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
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

void deck_iterations()
{
    using namespace iterateKT;
    using iterateKT::complex;

    uint   N    = 6;    // Number of iterations
    double t    = -0.1; // Production t
    double m3pi = 1.40; // 3pi invariant mass

    kinematics kinematics = new_kinematics(m3pi, M_PION);
    solver solver(kinematics);

    settings sets = default_settings();

    // The projection function is given by our Deck loop 
    auto   constant = [&](complex sigma){return 1.;};
    auto   deck     = [&](complex sigma){return pi1::deck(t, m3pi*m3pi, sigma);};
    isobar pwave = solver.add_isobar<P_wave>({constant, deck}, 3, id::Deck, "Deck");
    
    // -----------------------------------------------------------------------
    timer timer;
    plotter plotter;

    timer.start();

    double smin = -0., smax = 2.6;
    double A = kinematics->A();
    double B = kinematics->B();
    double C = kinematics->C();
    double D = kinematics->D();

    plot p1 = plotter.new_plot();
    p1.set_curve_points(1000);
    p1.set_xrange({smin, smax});
    p1.set_labels("#sigma   [GeV^{2}]", "#it{F}_{#Delta} (#it{t}, #it{m}_{3#pi}^{2} #; #sigma + #it{i}#epsilon)");
    p1.add_horizontal(0);
    p1.add_vertical(D);
    p1.shade_region({A,C});
    p1.add_header("#minus#it{t}  = 0.1, #it{m}_{3#pi}^{2} = (1.4)^{2}");
    p1.set_legend(0.6, 0.5);
    p1.add_curve( {smin, smax}, [&](double s) { return std::real(pwave->basis_function(1, s+IEPS)); }, "#Delta(#it{t}, #it{m}_{3#pi}^{2}; #sigma) #Omega(#sigma)");
    p1.add_dashed({smin, smax}, [&](double s) { return std::imag(pwave->basis_function(1, s+IEPS)); });
   
    std::vector<std::string> labels = {"1st", "2nd", "3rd", "4th", "5th", "6th"};
    for (int i = 1; i <= N; i++)
    {
        solver.iterate();
        p1.add_curve( {smin, smax}, [&](double s) { return std::real(pwave->basis_function(1, s+IEPS)); }, labels[i-1]);
        p1.add_dashed({smin, smax}, [&](double s) { return std::imag(pwave->basis_function(1, s+IEPS)); });
    };
    p1.save("deck_isobar.pdf");

    plot p2 = plotter.new_plot();
    p2.set_curve_points(1000);
    p2.set_xrange({smin, smax});
    p2.set_labels("#sigma   [GeV^{2}]", "#it{F}(#it{t}, #it{m}_{3#pi}^{2} #; #sigma + #it{i}#epsilon)");
    p2.add_horizontal(0);
    p2.add_vertical(D);
    p2.shade_region({A,C});
    p2.add_header("#minus#it{t}  = 0.1, #it{m}_{3#pi}^{2} = (1.4)^{2}");
    p2.set_legend(0.6, 0.6);
    p2.add_curve( {smin, smax}, [&](double s) { return std::real(pwave->basis_function(0, s+IEPS)); }, solid(jpacColor::Blue, "#alpha"));
    p2.add_curve( {smin, smax}, [&](double s) { return std::imag(pwave->basis_function(0, s+IEPS)); }, dashed(jpacColor::Blue));
    p2.add_curve( {smin, smax}, [&](double s) { return std::real(pwave->basis_function(1, s+IEPS)); }, solid(jpacColor::Red, "#Delta"));
    p2.add_curve( {smin, smax}, [&](double s) { return std::imag(pwave->basis_function(1, s+IEPS)); }, dashed(jpacColor::Red));
    p2.save("deck_comparison.pdf");

    timer.stop();
    timer.print_elapsed();
};