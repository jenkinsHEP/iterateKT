// Methods to interface K → 3π decay data with fitters
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef KAON_DATA_HPP
#define KAON_DATA_HPP

#include "constants.hpp"
#include "amplitudes/kaon.hpp"
#include "data_set.hpp"

namespace iterateKT { namespace kaon
{
    // Static identifiers for data_set types
    // Have two 'types' one which contains only Gamma and h
    // and another with Gamma, g, h, k
    static const uint kAll = 0, kHOnly = 1, kLambda = 2;

    inline data_set get_lambda_data()
    {
        data_set out;
        out._id   = "KL - KS interference";
        out._type = kLambda;
        out._N    = 2;
        out._z    = {0.0334, -0.0108};
        out._dz   = {7E-4, 48E-4};
        return out;
    };

    inline data_set get_dalitz_data(option opt)
    {
        data_set out;
        switch (opt)
        {
            case option::P_ppm:
            {
                out._id = "K+ -> pi+ pi+ pi-";
                out._z  = {
                    2.9590,   // Γ
                    -0.21134, // g
                    0.0185,   // h
                    -0.00463  // k
                };
                out._dz = {
                    218E-4,   // δΓ 
                    17E-5,    // δg
                    4E-4,     // δh
                    14E-5     // δk
                };
                out._option = option::P_ppm;
                out._type = kAll;
                out._N    = 4;
                break;
            };
            case option::P_zzp:
            {
                out._id = "K+ -> pi0 pi0 pi+";
                out._z  = {
                    0.9438, // Γ
                    0.626,  // g
                    0.052,  // h
                    0.0054  // k
                };
                out._dz = {
                    150E-4, // δΓ 
                    7E-3,   // δg
                    8E-3,   // δh
                    35E-4   // δk
                };
                out._option = option::P_zzp;
                out._type = kAll;
                out._N    = 4;
                break;
            };
            case option::L_pmz:
            {
                out._id = "KL -> pi+ pi- pi0";
                out._z  = {
                    1.6200, // Γ
                    0.678,  // g
                    0.076,  // h
                    0.0099  // k
                };
                out._dz = {
                    102E-4, // δΓ 
                    8E-3,   // δg
                    6E-3,   // δh
                    15E-4   // δk
                };
                out._option = option::L_pmz;
                out._type = kAll;
                out._N    = 4;
                break;
            };
            case option::L_zzz:
            {
                out._id = "KL -> pi0 pi0 pi0";
                out._z = {
                    2.5417,  // Γ
                    -0.0061  // h
                };
                out._dz = {
                    352E-4,  // δΓ
                    10E-4    // δh
                };
                out._option = option::L_zzz;
                out._type = kHOnly;
                out._N    = 2;
                break;
            };
            default: return out;
        };
        return out;
    };
}; /* namespace iterateKT */ }; /* namespace kaon_decay */

#endif