/*!
 * \file tracking_Kalman_filter.h
 * \brief Interface of a 2nd order PLL filter for carrier tracking loop
 * \author Javier Arribas, 2011. jarribas(at)cttc.es
 *
 * Class that implements 2 order PLL filter for tracking carrier loop.
 * The algorithm is described in
 * K.Borre, D.M.Akos, N.Bertelsen, P.Rinder, and S.H. Jensen,
 * A Software-Defined GPS and Galileo Receiver. A Single-Frequency Approach,
 * Birkhauser, 2007, Applied and Numerical Harmonic Analysis.
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2018  (see AUTHORS file for a list of contributors)
 *
 * GNSS-SDR is a software defined Global Navigation
 *          Satellite Systems receiver
 *
 * This file is part of GNSS-SDR.
 *
 * GNSS-SDR is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GNSS-SDR is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNSS-SDR. If not, see <https://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#ifndef GNSS_SDR_TRACKING_KALMAN_FILTER_H_
#define GNSS_SDR_TRACKING_KALMAN_FILTER_H_

#include <armadillo>
#include <gnuradio/gr_complex.h>
#include "MATH_CONSTANTS.h"

typedef std::pair<float, float> carr_nco;

/*!
 * \brief This class implements a standard Kalman filter for carrier tracking.
 *
 * The algorithm is described in:
 * TODO
 */
class Tracking_Kalman_filter
{
private:

    float d_pdi_carr = 0.0;
    uint32_t d_kf_order = 2U;
    bool d_kf_extended = false;

    // Kalman filter parameters
    arma::mat kf_P_x_ini;  // initial state error covariance matrix
    arma::mat kf_P_x;      // state error covariance matrix
    arma::mat kf_P_x_pre;  // Predicted state error covariance matrix
    arma::mat kf_P_y;      // innovation covariance matrix

    arma::mat kf_F;  // state transition matrix
    arma::mat kf_H;  // system matrix
    arma::mat kf_R;  // measurement error covariance matrix
    arma::mat kf_Q;  // system error covariance matrix

    arma::colvec kf_x_ini;  // initial state vector
    arma::colvec kf_x;      // state vector
    arma::colvec kf_x_pre;  // predicted state vector
    arma::colvec kf_z;      // measurement vector
    arma::colvec kf_y;      // innovation vector
    arma::mat kf_K;         // Kalman gain matrix

public:
    void initialize();

    void set_pdi(float pdi_carr);

    void set_KF_R(float err_variance);
    void set_KF_R(float ev_1, float ev_2);

    void set_KF_x_ini(float x_1, float x_2);
    void set_KF_x_ini(float x_1, float x_2, float x_3);

    void set_KF_P_x_ini(float P_x_1, float P_x_2);
    void set_KF_P_x_ini(float P_x_1, float P_x_2, float P_x_3);

    carr_nco get_carrier_nco(float PLL_discriminator, float pll_bw_hz);
    carr_nco get_carrier_nco(gr_complex Prompt, float d_carrier_power_snv);

    void setup_filter(float pdi_carr, uint32_t kf_order, bool kf_extended);

    Tracking_Kalman_filter();
    Tracking_Kalman_filter(float pdi_carr, uint32_t kf_order, bool kf_extended);
    ~Tracking_Kalman_filter();
};

#endif
