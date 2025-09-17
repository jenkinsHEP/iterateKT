// Basis functions for K -> 3pi
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#include <filesystem>
#include "plotter.hpp"
#include "amplitudes/kaon.hpp"
#include "isobars/pseudoscalar.hpp"

void calculate_isobars()
{
    using namespace iterateKT;

    // -----------------------------------------------------------------------
    
    // Set up general kinematics so everything knows masses
    // Assume masses are given in terms of pion mass
    kinematics kin = new_kinematics(M_KAON_AVG, M_PION_PM);

    // Set up our amplitude 
    solver solver(kin);
    
    // Add all the isobars, note the order they are added will be the order
    // the basis functions are generated
    std::vector<uint> empty = {}; // Pass empty to isobars with no sub polynomials
    isobar F0 = solver.add_isobar<I1_S0>({0, 1, 2}, 2, id::I1_S0, "F0"); 
    isobar F1 = solver.add_isobar<I1_P1>({1},       1, id::I1_P1, "F1"); 
    isobar F2 = solver.add_isobar<I1_S2>(empty,     2, id::I1_S2, "F2");
    isobar H1 = solver.add_isobar<I2_P1>({0, 1},    1, id::I2_P1, "H1"); 
    isobar H2 = solver.add_isobar<I2_S2>(empty,     2, id::I2_S2, "H2");

    // Iterate N times
    int N = 6;
    solver.timed_iterate(N);

    // Export the solution so it can be more easily recalled later
    std::string out_dir = main_dir()+"/scripts/kaon/basis_functions";
    std::filesystem::create_directory(out_dir);
    solver.export_solution(out_dir+"/basis");

    // -----------------------------------------------------------------------
    // Plot Results

    plotter plotter;
    double smin =  +0.06;
    double smax =  +0.14;

    auto plot_basis = [&](isobar isobar, int i, std::string label)
    {
        plot p = plotter.new_plot();
        p.set_curve_points(1000);
        p.set_legend(false);
        p.shade_region({kin->sth(), kin->pth()});
        p.add_vertical(kin->rth());

        label += "(#it{s} + #it{i}#epsilon)";
        p.set_labels("#it{s} [GeV^{2}]", label);

        auto rF = [&](int j){return [j,isobar,i](double s){return real(isobar->basis_function(j, i, s+IEPS));}; };
        auto iF = [&](int j){return [j,isobar,i](double s){return imag(isobar->basis_function(j, i, s+IEPS));}; };

        p.add_curve({smin, smax}, rF(1), dotted(jpacColor::Red));
        p.add_curve({smin, smax}, iF(1), dotted(jpacColor::Blue));
        p.add_curve({smin, smax}, rF(3), dashed(jpacColor::Red));
        p.add_curve({smin, smax}, iF(3), dashed(jpacColor::Blue));
        p.add_curve({smin, smax}, rF(6), solid( jpacColor::Red));
        p.add_curve({smin, smax}, iF(6), solid( jpacColor::Blue));
        return p;
    };

    std::vector<plot> Fs, Hs;
    Fs.emplace_back(plot_basis(F0, 0, "F_{0}^{#alpha}"));
    Fs.emplace_back(plot_basis(F1, 0, "F_{1}^{#alpha}"));
    Fs.emplace_back(plot_basis(F2, 0, "F_{2}^{#alpha}"));
    Fs.emplace_back(plot_basis(F0, 1, "F_{0}^{#beta}"));
    Fs.emplace_back(plot_basis(F1, 1, "F_{1}^{#beta}"));
    Fs.emplace_back(plot_basis(F2, 1, "F_{2}^{#beta}"));
    Fs.emplace_back(plot_basis(F0, 2, "F_{0}^{#gamma}"));
    Fs.emplace_back(plot_basis(F1, 2, "F_{1}^{#gamma}"));
    Fs.emplace_back(plot_basis(F2, 2, "F_{2}^{#gamma}"));
    Fs.emplace_back(plot_basis(F0, 3, "F_{0}^{#delta}"));
    Fs.emplace_back(plot_basis(F1, 3, "F_{1}^{#delta}"));
    Fs.emplace_back(plot_basis(F2, 3, "F_{2}^{#delta}"));
    plotter.combine({3,4}, Fs, out_dir+"/Fs.pdf");

    Hs.emplace_back(plot_basis(H1, 4, "H_{1}^{#mu}"));
    Hs.emplace_back(plot_basis(H2, 4, "H_{2}^{#mu}"));
    Hs.emplace_back(plot_basis(H1, 5, "H_{1}^{#nu}"));
    Hs.emplace_back(plot_basis(H2, 5, "H_{2}^{#nu}"));
    plotter.combine({2,2}, Hs, out_dir+"/Hs.pdf");
};