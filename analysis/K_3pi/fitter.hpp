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

#include "K_3pi/data.hpp"

namespace iterateKT { namespace kaon
{
    // Different masses for isospin projectionss
    const double M_KAON_PM  = 0.493677;
    const double M_KAON_0   = 0.497611;
    const double M_KAON_AVG = (M_KAON_PM + M_KAON_0)/2;
    const double M_PION_PM  = 0.13957039;
    const double M_PION_0   = 0.1349768;
    const double M_PION_AVG = (M_PION_PM + M_PION_0)/2;

    // Because the kaon mass is so smal, the mass splittings between isospin projections
    // make a noticable difference in the phase-space.
    // So, even though we calculate KT isobars in isospin limit, we can integrate with the realistic 
    complex integrate_over_physical_phasespace(amplitude amp, std::array<option,2> opts, std::array<double,4> ms)
    {
        using namespace boost::math::quadrature;
        double mK = ms[0], m1 = ms[1], m2 = ms[2], m3 = ms[3];
        double sth = norm(m2+m3), pth = norm(mK-m1);
        double sigma = mK*mK + m1*m1 + m2*m2 + m3*m3;
        auto dGamma = [amp,sigma,mK,m1,m2,m3,opts](double s)
        {
            double kappa = sqrt(kallen(s, mK*mK, m1*m1))*sqrt(kallen(s,m2*m2,m3*m3))/s;
            double tmin  = (sigma - s - kappa)/2, tmax  = (sigma - s + kappa)/2;
            auto d2Gamma = [amp,s,opts](double t)
            {
                complex A, B;
                amp->set_option(opts[0]); A = amp->evaluate(s,t);
                if  (opts[0] == opts[1]) return conj(A)*A;
                amp->set_option(opts[1]); B = amp->evaluate(s,t);  
                return conj(A)*B;     
            };
            return gauss_kronrod<double,15>::integrate(d2Gamma, tmin, tmax, 0, 1.E-9, NULL);
        };
        return gauss_kronrod<double,15>::integrate(dGamma, sth, pth, 0, 1.E-9, NULL);
    };

    double physical_width(amplitude amp, option opt)
    {
        double mK, m1, m2, m3;
        switch (opt)
        {
            case (option::P_ppm): 
            {
                mK = M_KAON_PM; 
                m1 = M_PION_PM; m2 = M_PION_PM; m3 = M_PION_PM; 
                break;
            };
            case (option::L_zzz):
            {
                mK = M_KAON_0; 
                m1 = M_PION_0;  m2 = M_PION_0;  m3 = M_PION_0; 
                break;
            }; 
            case (option::P_zzp):
            {
                mK = M_KAON_PM; 
                m1 = M_PION_0; m2 = M_PION_0; m3 = M_PION_PM; 
                break;
            };
            case (option::L_pmz):
            case (option::S_pmz):
            {
                mK = M_KAON_0; 
                m1 = M_PION_PM; m2 = M_PION_PM; m3 = M_PION_0; 
                break;
            };
            default: return NaN<double>();
        };

        double prefactors = 32*pow(2*PI*mK,3)*amp->combinatorial_factor();
        complex Gam = integrate_over_physical_phasespace(amp, {opt,opt}, {mK,m1,m2,m3});
        return real(Gam) / prefactors;
    };

    inline complex interference_lambda(amplitude amp)
    {
        double mK, m1, m2, m3;
        mK = M_KAON_0; 
        m1 = M_PION_PM; m2 = M_PION_PM; m3 = M_PION_0; 
        complex num = integrate_over_physical_phasespace(amp,{option::L_pmz,option::S_pmz}, {mK,m1,m2,m3});
        complex den = integrate_over_physical_phasespace(amp,{option::L_pmz,option::L_pmz}, {mK,m1,m2,m3});
        return num/den;
    };

    // Specify the fitter interface
    struct fit
    {
        // String letting us know what is being fit
        static std::string data_type(int i)
        {
            switch (i)
            {
                case kAll:    return "Γ & {g, h, k}";
                case kHOnly:  return "Γ & h";
                case kLambda: return "λ";
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
            if (data._type == kLambda)
            {
                auto chi2s = chi2_lambda(data, to_fit);
                return chi2s[0] + chi2s[1];
            };

            to_fit->set_option(data._option);
            
            // output
            double chi2 = 0;

            // χ² from Γ
            chi2 += chi2_width(data, to_fit); 

            // χ² from g h k
            auto chi2_ghk = chi2_dpars(data, to_fit);
            for (auto chi2_i : chi2_ghk) chi2 += chi2_i;
            return chi2;
        };

        // Compare widths
        static double chi2_width(const data_set & data, amplitude to_fit)
        {
            double gam_th = physical_width(to_fit, data._option);
            double gam_ex = data._z[0], dgam_ex = data._dz[0];
            return norm((gam_th - gam_ex)/dgam_ex);
        };

        // Interference lambda parameter
        static std::array<double,2> chi2_lambda(const data_set & data, amplitude to_fit)
        {
            complex lam_th = interference_lambda(to_fit);
            complex re_ex = data._z[0], dre_ex = data._dz[0];
            double chi2_re = norm((real(lam_th) - re_ex)/dre_ex);

            complex im_ex = data._z[1], dim_ex = data._dz[1];
            double chi2_im = norm((imag(lam_th) - im_ex)/dim_ex);
            return {chi2_re, chi2_im};
        };

        // Compare g, h, k
        static std::array<double,3> chi2_dpars(const data_set & data, amplitude to_fit)
        {
            auto dpars = to_fit->get_dalitz_parameters(1.E-3);
            double g_th = dpars[0], h_th = dpars[1], k_th = dpars[3];
     
            std::array<double,3> chi2;
            if (data._type == kHOnly)
            {
                double  h_ex = data._z[1],  dh_ex = data._dz[1];
                chi2[1] = norm((h_th - h_ex)/dh_ex);
                chi2[0] = 0; chi2[2] = 0;
            }
            else
            {
                double  g_ex = data._z[1],   h_ex = data._z[2],   k_ex = data._z[3];
                double dg_ex = data._dz[1], dh_ex = data._dz[2], dk_ex = data._dz[3];
                chi2[0] = norm((g_th - g_ex)/dg_ex);
                chi2[1] = norm((h_th - h_ex)/dh_ex);
                chi2[2] = norm((k_th - k_ex)/dk_ex);
            };
            return chi2;
        };
    };
}; /* namespace iterateKT */ }; /* namespace kaon_decay */
#endif