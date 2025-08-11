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
    // This defines the full amplitude, i.e. how the isobars are combined
    // Here is where we usually put the isospin combinations etc
    class pi1 : public raw_amplitude
    {
        public: 
        
        // Constructor
        pi1(kinematics kin, std::string id) : raw_amplitude(kin,id)
        {};
        
        // Spin 1 decay so (2j+1) = 3
        inline double combinatorial_factor(){ return 3; };

        static constexpr double _mu2 = 0.13957000*0.13957000;
        static constexpr double _eps = 1E-5;

        static inline complex deck(double t, double M2, complex s)
        {
            bool no_problem = (!is_zero(imag(s), _eps) || real(s) <= 4*_mu2);
            if  (no_problem) return deck_complex_plane(t, M2, s);
            else             return deck_on_real_axis (t, M2, real(s));
        };

        // If sig is complex, just evaluate all the square roots naively
        static inline complex deck_complex_plane(double t, double M2, complex s)
        {
            // Masses and momenta
            complex p2 = kallen(M2, t, _mu2)/4/M2;
            complex kappa  = csqrt(kallen(M2, t, _mu2)*kallen(M2, s, _mu2))/M2;
            complex rho    = csqrt(kallen(M2, s, _mu2))/M2;
            // Momentum transfer tau
            complex tauz = 2*_mu2-(M2+_mu2-t)*(M2-s+_mu2)/2/M2;
            complex taup = tauz+kappa/2;
            complex taum = tauz-kappa/2;
            // Assemble the final discontinuity
            complex  eta = M2-t-s+(_mu2-s)*(_mu2-t)/M2;
            complex zeta = eta-_mu2+tauz;
            complex   Q0 = log(_mu2-taum) - log(_mu2-taup)/kappa;
            
            return PI*rho/16/p2*((kappa*kappa-eta*eta)*Q0+4*zeta);
        };

        // Else we want to carefully handle the analytic continuation 
        static inline complex deck_on_real_axis(double t, double M2, double s)
        {
            // Masses and momenta
            double p = sqrt(kallen(M2, t, _mu2))/2/sqrt(M2);
            
            // Take the abs value and handle phase manually
            double aq = abs(csqrt(kallen(M2, s, _mu2))/2/sqrt(M2));
            double ak = 4*p*aq;

            int region = (s >= norm(sqrt(M2)-sqrt(_mu2)))
                       + (s >= M2 - _mu2)
                       + (s >= norm(sqrt(M2)+sqrt(_mu2)));

            complex kappa, rho;
            switch (region)
            {
                case 0: kappa = +  ak; rho = +2*  aq/sqrt(M2); break;
                case 1: kappa = +I*ak; rho = +2*I*aq/sqrt(M2); break;
                case 2: kappa = +I*ak; rho = -2*I*aq/sqrt(M2); break;
                case 3: kappa = -  ak; rho = -2*  aq/sqrt(M2); break;
                default: return NaN<complex>();
            };
            
            // Momentum transfer tau
            complex tauz = 2*_mu2-(M2+_mu2-t)*(M2-s+_mu2)/2/M2;
            complex taup = tauz+kappa/2, taum = tauz-kappa/2;

            // Assemble the final discontinuity
            complex  eta = M2-t-s+(_mu2-s)*(_mu2-t)/M2;
            complex zeta = eta-_mu2+tauz;

            bool log_problem = (region==2 && are_equal(s, M2-_mu2, 1E-2));
            complex Q0 = log_problem ? -I*PI/kappa 
                                     : log((_mu2-taum)/(_mu2-taup))/kappa;       

            return PI*rho/16/norm(p)*((norm(kappa)-eta*eta)*Q0 + 4*zeta);
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