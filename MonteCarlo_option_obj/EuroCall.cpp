//
// Created by Harshvardhan Singh on 23/10/2021.
//

#include "EuroCall.h"

#include <algorithm> // for std::max


EuroCall::EuroCall(double X, double T) : X_(X), T_(T) {}

double EuroCall::ComputePO(double S)
{
    return std::max(0.0, S - X_);
}