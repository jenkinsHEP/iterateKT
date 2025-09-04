// Isobars relevant for the decay of isospin-1/2 pseudoscalar decay to 3π
// 
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] - https://arxiv.org/abs/2403.17570
// [2] - https://arxiv.org/abs/2111.02417
// ------------------------------------------------------------------------------

#ifndef KAON_AMPLITUDES_HPP
#define KAON_AMPLITUDES_HPP

#include "isobar.hpp"
#include "utilities.hpp"
#include "kinematics.hpp"
#include "settings.hpp"
#include "phase_shift.hpp"
#include "isobars/pseudoscalar.hpp"

// For a general K_Pi1Pi2Pi3, we have:
// s = (K - π1)^2 = (π2 + π3)^2 = s1
// t = (K - π2)^2 = (π1 + π3)^2 = s2
// u = (K - π3)^2 = (π1 + π2)^2 = s3

// We assume the symmetric channels are t <-> u

namespace iterateKT
{
    inline settings default_settings()
    {
        settings sets;
        sets._exclusion_points        = 20;
        sets._exclusion_offsets       = {3.E-2, 5E-2};
        sets._infinitesimal           = 1E-8;
        sets._intermediate_energy     = 1.0;
        sets._cutoff                  = 20.0;
        sets._interpolation_offset    = 1E-4;
        sets._interpolation_points    = {400, 10, 100};

        double xi_sth = 1E-3,   eps_sth = 1E-3;
        double xi_pth = 1E-4,   eps_pth = 1E-3;
        double xi_rth = 2E-2,   eps_rth = 2E-2;
        sets._matching_intervals  = {xi_sth,  xi_pth,  xi_rth };
        sets._expansion_offsets   = {eps_sth, eps_pth, eps_rth};

        phase_args iso_0 = {"madrid/delta_00.dat", 1.69,  1, 1, 2};
        phase_args iso_1 = {"madrid/delta_11.dat", 1.69,  1, 1, 2};
        phase_args iso_2 = {"madrid/delta_02.dat", 9.99,  0, 1, 2};
        sets._phase_shifts = { {id::I0_P1, iso_1}, 
                               {id::I1_S0, iso_0}, {id::I1_P1, iso_1}, {id::I1_S2, iso_2},
                               {id::I2_P1, iso_1}, {id::I2_S2, iso_2}};
        return sets;
    };
    
    // ------------------------------------------------------------------------------
    // These are the invariant amplitudes from a 3π state of total isospin

    // I_3π = 1 amplitude
    class I1 : public raw_amplitude
    {
        public:

        I1(kinematics kin, std::string id) : raw_amplitude(kin,id){};

        inline complex prefactor_s(id iso_id, complex s, complex t, complex u)
        {
            switch (iso_id)
            {
                case id::I1_S0: return 1;
                case id::I1_S2: return -2./3;
                default: return 0;
            };
        };
        inline complex prefactor_t(id iso_id, complex s, complex t, complex u)
        {
            switch (iso_id)
            {
                case id::I1_P1: return (s-u);
                case id::I1_S2: return 1;
                default: return 0;
            };
        };
        inline complex prefactor_u(id iso_id, complex s, complex t, complex u)
        {
            return prefactor_t(iso_id, s, u, t);
        };
    };
    
    // I_3π = 2 amplitude
    class I2 : public raw_amplitude
    {
        public:

        I2(kinematics kin, std::string id) : raw_amplitude(kin,id){};

        inline complex prefactor_s(id iso_id, complex s, complex t, complex u)
        {
            return (iso_id == id::I2_S2) ? 1 : 0;
        };
        inline complex prefactor_t(id iso_id, complex s, complex t, complex u)
        {
            switch (iso_id)
            {
                case id::I2_P1: return 3*(s-u)/2;
                case id::I2_S2: return -1./2;
                default: return 0;
            };
        };
        inline complex prefactor_u(id iso_id, complex s, complex t, complex u)
        {
            return prefactor_t(iso_id, s, u, t);
        };
    };

    // ------------------------------------------------------------------------------
    // These are the physical amplitudes in the charged basis

    class KtoPPM : public raw_amplitude
    {
        public: 
        
        KtoPPM(kinematics xkin, std::string id) : raw_amplitude(xkin, id)
        {
            F = new_amplitude<I1>(xkin);
            H = new_amplitude<I2>(xkin);
        };

        inline complex evaluate(complex s, complex t, complex u)
        {
            return F->evaluate(t,s,u) + F->evaluate(u,t,s) + H->evaluate(s,t,u);
        };

        private: 
        amplitude F, H;
    };

    class KtoZZP : public raw_amplitude
    {
        public: 
        
        KtoZZP(kinematics xkin, std::string id) : raw_amplitude(xkin, id)
        {
            F = new_amplitude<I1>(xkin);
            H = new_amplitude<I2>(xkin);
        };

        inline complex evaluate(complex s, complex t, complex u)
        {
            return F->evaluate(s,t,u) + H->evaluate(s,t,u);
        };

        private: 
        amplitude F, H;
    };
};

#endif