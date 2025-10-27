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
    // This fitter takes a single data set and fits to it
    struct fit_single_tbin
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

        // Dont need any additional processing, just save
        static std::vector<complex> process_parameters(std::vector<complex> pars, amplitude to_fit)
        {
            to_fit->set_parameters(pars); 
            return pars; 
        };
            
        // Function to minimize
        // Filters whether we're looking at the real or imaginary parts 
        static double fcn(std::vector<data_set> & data_vector, amplitude to_fit)
        {
            double chi2 = 0;
            for (auto data : data_vector)
            {
                to_fit->set_option(data._option);
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

    // This one on the other hand assumes we are looking at multiple tbins and have couplings with explicit t-dependence
    struct fit_all_tbins
    {
        // None of these change
        static std::string data_type(int i){ return fit_single_tbin::data_type(i); };
        static double fcn(std::vector<data_set> & data_vector, amplitude to_fit)
        {
            return fit_single_tbin::fcn(data_vector, to_fit);
        };

        // Take in the parameters with the extra form factor slopes at the end
        static std::vector<complex> process_parameters(std::vector<complex> pars, amplitude to_fit)
        {
            // We assume pars.size() = 6
            // first three are the subtraction coefficients
            // last  three are the t-slopes
            std::array<double,4> t    = {-0.12,         -0.17,         -0.26,         -0.66};
            std::array<option,4> opts = {option::tbin0, option::tbin1, option::tbin2, option::tbin3};

            for (int i = 0; i < 4; i++)
            {
                std::vector<complex> new_pars;
                // g -> g * t exp{b(t-t_0)}
                for (int j = 0; j < 2; j++) new_pars.push_back(pars[j]*t[i]*exp(pars[2+j]*(t[i]-t[0])));
                to_fit->set_option(opts[i]);
                to_fit->set_parameters(new_pars);                
            };  
            return pars;
        };
    };
}; /* namespace COMPASS */ }; /* namespace iterateKT */

#endif /* COMPASS_FITTER_HPP */