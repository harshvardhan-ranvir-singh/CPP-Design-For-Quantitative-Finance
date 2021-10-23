//
// Created by Harshvardhan Singh on 23/10/2021.
//

#include "GBM_process.h"
#include "rv_library.h"
#include <cmath>

GBM_process::GBM_process(double S_0, double sig,  double r, double dt)
        :   S_0_(S_0)
        ,   sig_(sig)
        ,   r_(r)
        ,   dt_(dt)
{
    drift_ = (r_ - 0.5*sig_*sig_)*dt_;
    sgrt_ = sig_*std::sqrt(dt_);
}

double GBM_process::Next_S(double S)
{
    double w = rv::GetNormalVariate();
    return S*std::exp(drift_ + sgrt_*w);
}