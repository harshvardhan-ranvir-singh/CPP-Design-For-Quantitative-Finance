//
// Created by Harshvardhan Singh on 23/10/2021.
//

#include "Accumulator.h"
#include <cmath>

Accumulator::Accumulator(double r, double T)
        :   acc_vals_(0.0)
        ,   acc_squs_(0.0)
        ,   M_(0)
{
    discount_ = std::exp(-r*T);
}

void Accumulator::AddValue(double payoff)
{
    acc_vals_ += payoff;
    acc_squs_ += payoff*payoff;
    ++M_;
}

double Accumulator::GetOptionValue()
{
    return discount_*acc_vals_/M_;
}

double Accumulator::GetSE()
{
    return discount_*std::sqrt(acc_squs_ - acc_vals_*acc_vals_/M_)/M_;
}