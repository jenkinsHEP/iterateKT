// Methods to interface K → 3π decay data with fitters
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef KAON_FITTER_HPP
#define KAON_FITTER_HPP

#include "isobars/pseudoscalar.hpp"

namespace iterateKT { namespace kaon
{
    
    // Specify the fitter interface
    struct fit
    {
        // Static identifiers for data_set types
        static const int kAll = 0, kWidth = 1, kDalitz = 2;

        // String letting us know what is being fit
        static std::string data_type(int i)
        {
            switch (i)
            {
                case kAll:     return "Width & {g, h, k}";
                case kWidth:   return "Width";
                case kDalitz:  return "{g, h, k}";
                default: return "ERROR!";
            };
        };

        // Function being minimized, sum of chi2s of individual data sets and observables
        static double fcn(std::vector<data_set> & data_vector, amplitude to_fit)
        {
            double chi2_tot = 0;
            for (auto data : data_vector) chi2_tot += chi2(data, to_fit);
            return chi2_tot;
        };

        // Indidivudal chi2 from a single data set
        static double chi2(const data_set & data, amplitude to_fit)
        {
            int type = data._type;
            to_fit->set_option(data._option);

            double chi2 = 0;

            // χ² from Γ        
            if (type == kAll || type == kWidth) chi2 += chi2_width(data, to_fit); 

            // χ² from g h k
            std::array<double,3> chi2_ghk = {0,0,0};
            if (type == kAll || type == kDalitz) chi2_ghk = chi2_dpars(data, to_fit);
            for (auto chi2_i : chi2_ghk) chi2 += chi2_i;

            return chi2;
        };

        // Compare widths
        static double chi2_width(const data_set & data, amplitude to_fit)
        {
            double gam_th = to_fit->width();
            double gam_ex = data._z[0], dgam_ex = data._dz[0];
            return norm((gam_th - gam_ex)/dgam_ex);
        };

        // Compare g, h, k
        static std::array<double,3> chi2_dpars(const data_set & data, amplitude to_fit)
        {
            double mp2  = M_PION_PM*M_PION_PM;
            auto dpars  = to_fit->get_dalitz_parameters(1E-5, {mp2, mp2});
            double g_th = dpars[0], h_th = dpars[1], k_th = dpars[3];

            bool n = (data._type == kAll);
            std::array<double,3> chi2, ghk = {dpars[0], dpars[1], dpars[3]};
            for (int i = 0; i < 3; i++) chi2[i] = norm((ghk[i]-data._z[i+n])/data._dz[i+n]);
            return chi2;
        };

        // We only fit the real parts of the parameters while the imaginary parts are
        // given by requiring Taylor invariants have vanishing imaginary parts
        static std::vector<complex> process_fitter_parameters(std::vector<complex> in_pars, amplitude amp)
        {
            double eps = 1E-5, r = amp->get_kinematics()->s0();

            //------------------------------------------------------------------------
            // First we fix the imaginary parts of the M's and N's (total 3π I=1)
            isobar F0 = amp->get_isobar(id::I1_S0);
            isobar F1 = amp->get_isobar(id::I1_P1);
            isobar F2 = amp->get_isobar(id::I1_S2);
            std::array<isobar,3> F = {F0, F1, F2};

            // Coefficients of taylor expansion up to cubic
            std::array<std::array<complex,3>,3> A, B, C;
            
            // First index is isospin, second is basis function ID
            for (uint i = 0; i < 3; i++)
            {
                for (uint n = 0; n <= 2; n++)
                {
                    A[i][n] = F[i]->basis_function(n, 0);
                    B[i][n] = F[i]->basis_derivative<1>(n, 0, eps);
                    C[i][n] = F[i]->basis_derivative<2>(n, 0, eps)/2.;
                };
            };

            // Construct the 3x3 matrix of Taylor invariants
            // First index is which invariant, and second is basis_id
            TArrayD reT_data(9), imT_data(9);
            for (uint n = 0; n < 3; n++)
            {
                std::array<complex,4> T;
                T[0] =   A[0][n] + 4./3*A[2][n] - 3*r*A[1][n] + 3*r*B[2][n] +9*r*r*C[2][n];
                T[1] =   B[0][n] - 5./3*B[2][n] + 3  *A[1][n] - 9*r*C[2][n];
                T[2] =   B[1][n] +      C[2][n];
               
                for (int j = 0; j < 3; j++)
                {
                    reT_data[3*j+n] = real(T[j]);  imT_data[3*j+n] = imag(T[j]);
                };
            };
            TMatrixD reT(3,3,reT_data.GetArray()), imT(3,3,imT_data.GetArray());

            // Now we actually solve the matrix equation relating reMu and imMu
            TMatrixD M = reT.Invert()*imT; M *= -1;

            Double_t rePars_data[3];
            for (int i = 0; i < 3; i++) rePars_data[i] = real(in_pars[i]);
            TVectorD rePars(3, rePars_data), imPars = M*rePars;
            
            // Assemble together output vector
            std::vector<complex> out_pars;
            for (int i = 0; i < 3; i++) out_pars.push_back(rePars[i]+I*imPars[i]);
            out_pars.push_back(in_pars.back()); // Last one stays real 

            return out_pars;
        };
    };
}; /* namespace iterateKT */ }; /* namespace kaon_decay */
#endif