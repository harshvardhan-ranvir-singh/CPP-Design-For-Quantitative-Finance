#include <iostream>
#include <algorithm>    // std::max(library function) lives under <algorithm>
#include <cmath>        // std::sqrt() and std::exp() are library functions
#include "utility.h"
#include "rv_library.h"


//  prototypes


// next_S() is the function to compute the next value along the sample path.
double next_S(const double & S, const double & drift, const double & sgrt);
double ComputePayoff(const double & S, const double &X);



// main()

int main() {

    double S_0 = 100;   // Initial stock price
    double r = 0.05;    // Riskless rate
    double sig = 0.2;   // volatility

    double X =  100;    // Option strike
    double T = 1;       // Time to maturity

    long N = 100;       // Time steps
    long M = 1000000;     // Sample paths

    double dt = T / N;      // time steps
    double drift = (r - 0.5 * sig * sig) * dt;       // drift
    double sgrt = sig * std::sqrt(dt);      // vol

    double acc_vals = 0;    // accumulates values
    double acc_squs = 0;    // accumulates squared values

    // Generating sample paths.
    // By the time we come out of the inner loop, we have got a
    // value of a variable at time T.
    for (long j = 1; j <= M; ++j)               // for each path
    {
        double path_S = S_0;
        for (long i = 1; i <= N; ++i)           // for each time step
        {
            path_S = next_S(path_S, drift, sgrt);
        }

        double payoff = ComputePayoff(path_S, X);

        acc_vals += payoff;
        acc_squs += payoff * payoff;
    }



    double c;           // Option value
    double se;          // Standard error

    c = acc_vals / M;
    se = std::sqrt(acc_squs - acc_vals * acc_vals / M) / M;

    double discount = std::exp(-r * T);

    c *= discount;
    se *= discount;

    // we have individual payoffs. p(j) will be the payoff(max(S - X)) along the jth path

    ut::OutputLine("Option value ", c);
    ut::OutputLine("Standard error ", se);

    return ut::PauseAndReturn();
}



// next_S()
double next_S(const double & S, const double & drift, const double & sgrt)
{
    // GetNormalVariate() returns a standard normal iid variate.
    double w = rv::GetNormalVariate();
    return S * std::exp(drift + sgrt * w);
}


// ComputePayoff()
double ComputePayoff(const double & S, const double &X)
{
    return std::max(0.0, S - X);
}




