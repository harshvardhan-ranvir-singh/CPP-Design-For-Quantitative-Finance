#include "rv_library.h"
#include "utility.h"

#include "EuroCall.h"
#include "GBM_process.h"
#include "Accumulator.h"

#include <algorithm>
#include <cmath>
#include <ctime>

int main() {
    ut::OutputLine("This version of a Monte Carlo has an option object");

    double S_0 = 100;
    double r = 0.05;
    double sig = 0.20;

    double X = 100;
    double T = 1;

    // number of sample paths
    long M = 100000;
    // number of time steps
    long N = 100;

    // an option object
    EuroCall opt(X,T);
    GBM_process prc(S_0, sig, r,T/N);
    Accumulator acc(r, T);

    double st_time = clock();


    for(long j = 1; j <= M; ++j)
    {
        double S = prc.GetS0();

        for(long i = 1; i <= N; ++i)
        {
            S = prc.Next_S(S);
        }

        // call the option's ComputePO method.
        // Do not need to pass X to the method.
        double payoff = opt.ComputePO(S);
        acc.AddValue(payoff);

    }

    double c = acc.GetOptionValue();
    double se = acc.GetSE();

    double el_time = (clock() - st_time)/CLOCKS_PER_SEC;

    ut::OutputLine("Option value ", c);
    ut::OutputLine("se           ", se);
    ut::OutputLine("time taken   ", el_time);

    return ut::PauseAndReturn();

}
