// Amplitudes relevant for the decay of meson with JP = 1-+ into 3pi
// 
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2024)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef PI1_AMPLITUDES_HPP
#define PI1_AMPLITUDES_HPP

#include "amplitude.hpp"
#include "isobar.hpp"
#include "utilities.hpp"
#include "kinematics.hpp"
#include "settings.hpp"
#include "GKPY.hpp"

#include"isobars/pi1.hpp"
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace iterateKT
{ 
    inline settings default_settings()
    {
        settings sets;
        sets._exclusion_points        = 6;
        sets._infinitesimal           = 1E-7;
        sets._intermediate_energy     = 4;
        sets._cutoff                  = 20;
        sets._interpolation_offset    = 0.1;
        sets._interpolation_points    = {200, 10, 100};
        double xi_sth = 1E-3,  eps_sth = 1E-3;
        double xi_pth = 1E-3,  eps_pth = 1E-2;
        double xi_rth = 3E-1,  eps_rth = 3E-1;

        sets._exclusion_offsets   = {2E-1, 2E-1};
        sets._matching_intervals  = {xi_sth,  xi_pth,  xi_rth };
        sets._expansion_offsets   = {eps_sth, eps_pth, eps_rth};
        return sets;
    };

    class pi1 : public raw_amplitude
    {
        public: 
        
        // Constructor
        pi1(kinematics kin, std::string id) : raw_amplitude(kin,id)
        {};
        
        // Spin 1 decay so (2j+1) = 3
        inline double combinatorial_factor(){ return 3; };

        static constexpr double _mu  = 0.13957000;
        static constexpr double _mu2 = _mu*_mu;
        static constexpr double _eps = 1E-5;

        static inline complex deck(complex t, complex M2, complex s)
        {
            // Masses and momenta
            complex p   = csqrt(kallen(M2, t, complex(_mu2)))/2/csqrt(M2);
            complex q   = csqrt(kallen(M2, s, complex(_mu2)))/2/csqrt(M2);
            // Theres vertical cuts here, flip sign if we cross
            if (real(s-M2)>=_mu2) q *= -1;
            // phasespace
            complex rho = 2*q/csqrt(M2);
            // Kacser function
            complex k   = 4*p*q; 
            // Momentum transfer tau
            auto    tau = [&](double z)
            {
                bool above_rth = real(s) >= norm(csqrt(M2)+_mu);
                complex x = (above_rth) ? real(s) : s;
                return 2*_mu2-(M2+_mu2-t)*(M2-x+_mu2)/2/M2+z*k/2; 
            };
            // Projection of OPE
            complex Q0  = (log(_mu2-tau(-1))-log(_mu2-tau(+1)))/k;
            // Assemble the final discontinuity
            complex a   = M2-t-s+(_mu2-s)*(_mu2-t)/M2;
            complex b   = a-_mu2+tau(0);
            return PI*rho/p/p/16*((k*k-a*a)*Q0+4*b);
        };

        // Assuming a pi- pi- pi+ decay and only P-waves
        // s = (pi- + pi+)^2 
        // t = (pi- + pi+)^2
        // u = (pi- + pi-)^2
        inline complex prefactor_s(id iso_id, complex s, complex t, complex u){ return csqrt(_kinematics->kibble(s,t,u)); };
        inline complex prefactor_t(id iso_id, complex s, complex t, complex u){ return - prefactor_s(iso_id, t, s, u); };
        inline complex prefactor_u(id iso_id, complex s, complex t, complex u){ return 0.; };
    };
}; // namespace iterateKT 

#endif // PI1_AMPLITUDES_HPP