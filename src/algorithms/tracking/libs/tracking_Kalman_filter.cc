/*!
 * \file tracking_Kalman_filter.cc
 * \brief Implementation of a 2nd order PLL filter for tracking carrier loop.
 * \author Javier Arribas, 2011. jarribas(at)cttc.es
 *
 * Class that implements 2 order PLL filter for tracking carrier loop. The algorithm
 * is described in:
 * K.Borre, D.M.Akos, N.Bertelsen, P.Rinder, and S.~H.~Jensen, A Software-Defined
 * GPS and Galileo Receiver. A Single-Frequency Approach,
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

#include "tracking_Kalman_filter.h"

void Tracking_Kalman_filter::set_pdi(float pdi_carr)
{
    d_pdi_carr = pdi_carr;  // Summation interval for code
    setup_filter(pdi_carr, d_kf_order, d_kf_extended);
}

void Tracking_Kalman_filter::set_KF_R(float CN0_SNV_dB_Hz)
{
    double CN_lin = std::pow(10, CN0_SNV_dB_Hz / 10.0);
    double sigma2_phase_detector_cycles2 = (1.0 / (2.0 * CN_lin * d_pdi_carr)) * (1.0 + 1.0 / (2.0 * CN_lin * d_pdi_carr));

    kf_R = arma::zeros(1, 1);
    kf_R(0) = sigma2_phase_detector_cycles2;
}

void Tracking_Kalman_filter::set_KF_R(float ev_1, float ev_2)
{
    kf_R = arma::zeros(2, 2);
    kf_R(0, 0) = ev_1;
    kf_R(1, 1) = ev_2;
}

void Tracking_Kalman_filter::set_KF_x_ini(float x_1, float x_2)
{
    kf_x_ini = arma::zeros(2, 1);
    kf_x_ini(0) = x_1;
    kf_x_ini(1) = x_2;
}

void Tracking_Kalman_filter::set_KF_x_ini(float x_1, float x_2, float x_3)
{
    kf_x_ini = arma::zeros(3, 1);
    kf_x_ini(0) = x_1;
    kf_x_ini(1) = x_2;
    kf_x_ini(2) = x_3;
}

void Tracking_Kalman_filter::set_KF_P_x_ini(float P_x_1, float P_x_2)
{
    kf_P_x_ini = arma::zeros(2, 2);
    kf_P_x_ini(0, 0) = P_x_1;
    kf_P_x_ini(1, 1) = P_x_2;
}

void Tracking_Kalman_filter::set_KF_P_x_ini(float P_x_1, float P_x_2, float P_x_3)
{
    kf_P_x_ini = arma::zeros(3, 3);
    kf_P_x_ini(0, 0) = P_x_1;
    kf_P_x_ini(1, 1) = P_x_2;
    kf_P_x_ini(2, 2) = P_x_3;
}

void Tracking_Kalman_filter::initialize()
{
    // Initial Kalman Paramters
    kf_P_x = kf_P_x_ini;
    kf_x = kf_x_ini;
}

/*
 * Kalman filter using PLL arctangent discriminator
 * Req Input in [Hz/Ti]
 * The output is in [Hz/s].
 */
carr_nco Tracking_Kalman_filter::get_carrier_nco(float PLL_discriminator, float pll_bw_hz)
{
    // Kalman state prediction (time update)
    kf_x_pre = kf_F * kf_x;                        //state prediction
    kf_P_x_pre = kf_F * kf_P_x * kf_F.t() + kf_Q;  //state error covariance prediction

    // Kalman estimation (measurement update)
    kf_y = arma::zeros(1, 1);
    kf_y(0) = PLL_discriminator;

    // Kalman filter update step
    kf_P_y = kf_H * kf_P_x_pre * kf_H.t() + kf_R;                       // innovation covariance matrix
#if 0
    float zeta = 0.7;
    float Wn = pll_bw_hz * 8.0 * zeta / (4.0 * zeta * zeta + 1.0);
    float tau1 = 0.25 / (Wn * Wn);
    float tau2 = 2.0 * zeta * Wn;

    kf_K = arma::zeros(2, 1);
    kf_K(0) = tau1;
    kf_K(1) = tau2 / d_pdi_carr;
#else
    kf_K = (kf_P_x_pre * kf_H.t()) * arma::inv(kf_P_y);                 // Kalman gain
#endif
    kf_x = kf_x_pre + kf_K * kf_y;                                      // updated state estimation
    kf_P_x = (arma::eye(size(kf_P_x_pre)) - kf_K * kf_H) * kf_P_x_pre;  // update state estimation error covariance matrix

    // Set a new carrier estimation to the NCO
    return std::make_pair(static_cast<float>(kf_x(0)), static_cast<float>(kf_x(1)));
}

/*
 * Kalman filter operating on correlator outputs
 * Req Input in [Hz/Ti]
 * The output is in [Hz/s].
 */
carr_nco Tracking_Kalman_filter::get_carrier_nco(gr_complex Prompt, float d_carrier_power_snv)
{
    // Kalman state prediction (time update)
    kf_x_pre = kf_F * kf_x;                        //state prediction
    kf_P_x_pre = kf_F * kf_P_x * kf_F.t() + kf_Q;  //state error covariance prediction

    // Kalman estimation (measurement update)
    kf_z = arma::zeros(2, 1);
    kf_z(0) = Prompt.real();
    kf_z(1) = Prompt.imag();

    kf_y = arma::zeros(2, 1);
    kf_y(0) = kf_z(0) - std::sqrt(d_carrier_power_snv) * cos(kf_x_pre(0));
    kf_y(1) = kf_z(1) - std::sqrt(d_carrier_power_snv) * sin(kf_x_pre(0));

    // Observation matrix update (Jacobian)
    kf_H = arma::zeros(2, 3);
    kf_H(0, 0) = std::sqrt(d_carrier_power_snv) * -sin(kf_x_pre(0));
    kf_H(1, 0) = std::sqrt(d_carrier_power_snv) * cos(kf_x_pre(0));

    // Kalman filter update step
    kf_P_y = kf_H * kf_P_x_pre * kf_H.t() + kf_R;  // innovation covariance matrix
    kf_K = (kf_P_x_pre * kf_H.t()) * arma::inv(kf_P_y);                 // Kalman gain
    kf_x = kf_x_pre + kf_K * kf_y;                                      // updated state estimation
    kf_P_x = (arma::eye(size(kf_P_x_pre)) - kf_K * kf_H) * kf_P_x_pre;  // update state estimation error covariance matrix

    // Set a new carrier estimation to the NCO
    return std::make_pair(static_cast<float>(kf_x(0)), static_cast<float>(kf_x(1)));
}

Tracking_Kalman_filter::Tracking_Kalman_filter()
{
    d_pdi_carr = 0.001;
    d_kf_order = 2U;
    d_kf_extended = false;
    setup_filter(d_pdi_carr, d_kf_order, d_kf_extended);
}


Tracking_Kalman_filter::Tracking_Kalman_filter(float pdi_carr, uint32_t kf_order, bool kf_extended)
{
    setup_filter(pdi_carr, kf_order, kf_extended);
}


void Tracking_Kalman_filter::setup_filter(float pdi_carr, uint32_t kf_order, bool kf_extended)
{
    d_pdi_carr = pdi_carr;
    d_kf_order = kf_order;
    d_kf_extended = kf_extended;

    // covariances (static)
    double sigma2_carrier_phase = std::pow(0.5 * PI_2 / 3.0, 2);
    double sigma2_doppler = std::pow(450 / 3.0, 2);
    double sigma2_doppler_rate = std::pow(4.0 * PI_2, 2) / 12.0;

    switch (kf_order)
        {
            case (3U):
                kf_P_x_ini = arma::zeros(3, 3);
                kf_P_x_ini(0, 0) = sigma2_carrier_phase;
                kf_P_x_ini(1, 1) = sigma2_doppler;
                kf_P_x_ini(2, 2) = sigma2_doppler_rate;

                kf_Q = arma::zeros(3, 3);
                kf_Q(0, 0) = std::pow(pdi_carr, 6);
                kf_Q(1, 1) = std::pow(pdi_carr, 4);
                kf_Q(2, 2) = std::pow(pdi_carr, 2);

                kf_F = arma::zeros(3, 3);
                kf_F(0, 0) = 1.0;
                kf_F(0, 1) = PI_2 * pdi_carr;
                kf_F(0, 2) = 0.5 * PI_2 * std::pow(pdi_carr, 2);
                kf_F(1, 0) = 0.0;
                kf_F(1, 1) = 1.0;
                kf_F(1, 2) = pdi_carr;
                kf_F(2, 0) = 0.0;
                kf_F(2, 1) = 0.0;
                kf_F(2, 2) = 1.0;

                kf_x = arma::zeros(3, 1);

                break;

            default:
                kf_P_x_ini = arma::zeros(2, 2);
                kf_P_x_ini(0, 0) = sigma2_carrier_phase;
                kf_P_x_ini(1, 1) = sigma2_doppler;

                kf_Q = arma::zeros(2, 2);
                kf_Q(0, 0) = std::pow(pdi_carr, 4);
                kf_Q(1, 1) = std::pow(pdi_carr, 2);

                kf_F = arma::zeros(2, 2);
                kf_F(0, 0) = 1.0;
                kf_F(0, 1) = PI_2 * pdi_carr;
                kf_F(1, 0) = 0.0;
                kf_F(1, 1) = 1.0;

                kf_x = arma::zeros(2, 1);
        }

    if (kf_extended)
        {
            kf_R = arma::eye(2, 2);
            kf_H = arma::zeros(2, kf_order);

            kf_y = arma::zeros(2, 1);
            kf_P_y = arma::eye(2, 2);
        }
    else
        {
            double CN_dB_Hz = 30;
            kf_R = arma::eye(1, 1);
            set_KF_R(CN_dB_Hz);

            kf_H = arma::zeros(1, kf_order);
            kf_H(0, 0) = 1.0;

            kf_y = arma::zeros(1, 1);
            kf_P_y = arma::eye(1, 1);
        }
}

Tracking_Kalman_filter::~Tracking_Kalman_filter()
{
}
