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
#include "timer.hpp"

#include"isobars/pi1.hpp"
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace iterateKT
{ 
    inline settings default_settings()
    {
        settings sets;
        sets._exclusion_points        = 10;
        sets._infinitesimal           = 1E-8;
        sets._intermediate_energy     = 4;
        sets._cutoff                  = 20;
        sets._interpolation_offset    = 0.1;
        sets._interpolation_points    = {200, 10, 100};
        double xi_sth = 1E-3,  eps_sth = 1E-3;
        double xi_pth = 3E-3,  eps_pth = 4E-3;
        double xi_rth = 3E-1,  eps_rth = 3E-1;

        sets._exclusion_offsets   = {5E-2, 5E-2};
        sets._matching_intervals  = {xi_sth,  xi_pth,  xi_rth };
        sets._expansion_offsets   = {eps_sth, eps_pth, eps_rth};

        phase_args iso_1   = {"madrid/delta_11.dat", 1.69,  1, 1, 2};
        sets._phase_shifts = { {id::P_wave, iso_1}, {id::Contact, iso_1}, {id::Deck, iso_1} };

        return sets;
    };

    class pi1 : public raw_amplitude
    {
        public: 
        
        // Constructor
        pi1(kinematics kin) : raw_amplitude(kin) {};
        
        // Spin 1 decay so (2j+1) = 3
        inline double combinatorial_factor(){ return 3; };

        static constexpr double _mu  = 0.13957000;
        static constexpr double _mu2 = _mu*_mu;

        // ------------------------------------------------------------------------------
        // Things related to the inclusion of the bubble

        // P-wave projection of 2 particle phase space with BW factor
        // M2   -> total 3pi invariant mass
        // s    -> 2pi subsystem mass
        // lam2 -> cutoff mass squared (mass of first ignored particle)
        static inline complex bubble(complex M2, complex s, double lam2)
        {
            complex mu2 = complex(_mu2);
            complex q   = csqrt(kallen(M2, s, mu2))/2/csqrt(M2);
            complex rho = 2*q/csqrt(M2);

            // careful if we cross above the three-body cut
            // multiply by -1 to not change sign and stay on the same sheet
            bool above_3bcut = (real(s) >= real(M2)+_mu2);
            if  (above_3bcut){ q *= -1; rho *= -1; };

            complex z2 = q*q / lam2;
            complex blatt_weisskopf = 1./(1+z2);
            return z2*blatt_weisskopf*rho;
        };

        // ------------------------------------------------------------------------------
        // Things related to the inclusion of the Deck triangle

        static inline complex tau(complex t, complex M2, complex s, double z)
        {
            complex mu2 = complex(_mu2);
            complex p   = csqrt(kallen(M2, t, mu2))/2/csqrt(M2);
            complex q   = csqrt(kallen(M2, s, mu2))/2/csqrt(M2);
            complex rho = 2*q/csqrt(M2);

            // careful if we cross above the three-body cut
            // multiply by -1 to not change sign and stay on the same sheet
            bool above_3bcut = (real(s) >= real(M2)+_mu2);
            if  (above_3bcut){ q *= -1; rho *= -1; };
            return  mu2 + s - (M2+mu2-t)*(M2+s-mu2)/2/M2 + 2*p*q*z;
        };

        // 3P1 projection of the vanilla OPE
        // t    -> momentum transfer of (external) Pomeron
        // M2   -> total 3pi invariant mass
        // s    -> 2pi subsystem imvariant mass
        // mex2 -> exchanged particle mass
        static inline complex deck(complex t, complex M2, complex s, double mex2)
        {
            // Masses and momenta
            complex mu2 = complex(_mu2);
            complex p   = csqrt(kallen(M2, t, mu2))/2/csqrt(M2);
            complex q   = csqrt(kallen(M2, s, mu2))/2/csqrt(M2);
            complex rho = 2*q/csqrt(M2);

            // careful if we cross above the three-body cut
            // multiply by -1 to not change sign and stay on the same sheet
            bool above_3bcut = (real(s) >= real(M2)+_mu2);
            if  (above_3bcut){ q *= -1; rho *= -1; };
            // Momentum tranfer at costheta = 0
            complex t0  = tau(t, M2, s, 0); 
            // Angular argument
            complex z   = (mex2 - t0)/2/p/q;

            // Continution depends on the ieps used for s
            // and not that of z
            bool above_thr  = real(s) >= 4*_mu2;

            // Legendre of 2nd kind
            complex Q0;
            if (above_thr) Q0 = log(-(z+1)/(z-1))/2-I*sign(imag(s))*PI/2;
            else           Q0 = log( (z+1)/(z-1))/2;

            // Final discontinuity
            return rho*q/p*((1-z*z)*Q0+z);
        };
        static inline complex deck(complex t, complex M2, complex s)
        {
            return deck(t, M2, s, _mu2);
        };

        static inline complex deck_with_FF(complex t, complex M2, complex s, double lam2)
        {
            return deck(t, M2, s, _mu2) - deck(t, M2, s, lam2);
        };

        // Assuming a pi- pi- pi+ decay and only P-waves
        // s = (pi- + pi+)^2 
        // t = (pi- + pi+)^2
        // u = (pi- + pi-)^2
        inline complex prefactor_s(id iso_id, complex s, complex t, complex u){ return csqrt(_kinematics->kibble(s,t,u)); };
        inline complex prefactor_t(id iso_id, complex s, complex t, complex u){ return - prefactor_s(iso_id, t, s, u); };
        inline complex prefactor_u(id iso_id, complex s, complex t, complex u){ return 0.; };
    };

    // ------------------------------------------------------------------------------
    // The following are container amplitudes which holds multiple copies of the above
    // but at different bins in production t for simultaneous fits

    enum class option : unsigned int { set_m3pi2, set_lam2, set_m3pibin, set_tbin };

    class pi1_across_tbins : public raw_amplitude
    {
        public:

        pi1_across_tbins(kinematics xkin, std::array<double,4> tvals)
        : raw_amplitude(xkin), _tvals(tvals)
        {
            for (int i = 0; i < 4; i++) 
            { 
                _tbins.emplace_back(new_amplitude<pi1>(xkin));
                _tbins[i]->set_name("tbin " + to_string(i));
            };
            initialize();
        };

        inline void set_option(option opt, double x)
        {
            switch (opt)
            {
                case option::set_m3pi2: _m3pi2 = x; break;
                case option::set_lam2:  _lam2  = x; break;
                case option::set_tbin:  _current = _tbins[int(std::round(x))]; break;
                default: return;
            };
        };
        inline uint N_pars(){ return _current->N_pars(); };
        inline void set_parameters(std::vector<complex> x){ _current->set_parameters(x); };
        inline complex evaluate(complex s, complex t, complex u){ return _current->evaluate(s, t, u); };

        private:

        // Fixed m3pi and bubble cutoff sqaured
        double _m3pi2 = 1.4, _lam2 = norm(M_RHO);

        // Save each of the 4 bins as a pointer
        std::array<double,4>   _tvals;
        std::vector<amplitude> _tbins;

        // This is the one that gets called
        amplitude _current;

        inline void initialize()
        {
            int    niter = 10; // number of KT iterations to do

            // These change as t changes
            auto constant = [](complex sigma){ return complex(1.); };
            auto deck = [&](double m3pi2, double t)
            { 
                return [m3pi2,t](complex sigma){ return pi1::deck(t, m3pi2, sigma);}; 
            };

            timer timer;
            timer.start();

            // Set up all the amplitudes
            for (int i = 0; i < 4; i++)
            {
                _tbins[i]->add_isobar<P_wave>({constant, deck(_m3pi2, _tvals[i])}, 3, id::P_wave, "P-wave");
                _tbins[i]->iterate(niter);
                timer.lap("iterated tbin " + to_string(i));
            };
            timer.stop(); timer.print_elapsed();

            // At end place the first one in _current
            set_option(option::set_tbin, 0);
        };
    };
}; // namespace iterateKT 

#endif // PI1_AMPLITUDES_HPP