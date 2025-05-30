// Methods to interface iterateKT::fitter
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef COMPASS_FITTER_HPP
#define COMPASS_FITTER_HPP

#include "COMPASS_pi1/data.hpp"

namespace iterateKT { namespace COMPASS
{
    struct fit
    {
        static std::string data_type(int i)
        {
            switch (i)
            {
                case kReal: return "Re (M)";
                case kImag: return "Im (M)";
                case kAbs:  return "Abs (M)";
                default: return "ERROR!";
            };
        };

        // Function to minimize
        // Filters whether we're looking at the real or imaginary parts 
        static double fcn(std::vector<data_set> & data_vector, amplitude to_fit)
        {
            double chi2 = 0;
            for (auto data : data_vector)
            {
                for (int i = 0; i < data._N; i++)
                {
                    double from_data  = data._z[i];
                    iterateKT::complex from_model = to_fit->evaluate(data._x[i], data._y[i]);  

                    switch (data._type)
                    {
                        // These two use difference of squares
                        case kReal: chi2 += norm(from_data - real(from_model));        break;
                        case kImag: chi2 += norm(from_data - imag(from_model));        break;
                        // This is a true chi2
                        case kAbs:
                        {
                            if (is_zero(data._dz[i])) continue;
                            chi2  += norm((from_data - abs(from_model)) / data._dz[i]); 
                            break;
                        };
                        default: break;
                    };
                };
            };
            return chi2;
        };
    };
}; /* namespace COMPASS */ }; /* namespace iterateKT */

#endif /* COMPASS_FITTER_HPP */