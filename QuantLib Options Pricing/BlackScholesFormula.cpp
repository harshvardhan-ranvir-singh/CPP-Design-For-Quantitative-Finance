# include<iostream>
#include <ql/math/distributions/normaldistribution.hpp>
#include <cmath>


struct blackscholes {
    double value;
    double delta;
    double gamma;
    double theta;


};

blackscholes black_scholes_formula(double S, double K, double T,
                                   double r, double sigma, bool call);




using std::sqrt;
using std::exp;
using std::log;

namespace {
    QuantLib::CumulativeNormalDistribution N;
    QuantLib::NormalDistribution n;
}

blackscholes black_scholes_formula(double S, double K, double T,
                                   double r, double sigma, bool call) {
    double d1 = (1/(sigma*sqrt(T))) * (log(S/K) + (r+sigma*sigma/2)*T);
    double d2 = d1 - sigma*sqrt(T);

    blackscholes results;
    if (call) {
        results.value = N(d1)*S - N(d2)*K*exp(-r*T);
        results.delta = N(d1);
        results.theta = -(S*n(d1)*sigma)/(2*sqrt(T)) - r*K*exp(-r*T)*N(d2);
    } else {
        results.value = N(-d2)*K*exp(-r*T) -N(-d1)*S;
        results.delta = -N(-d1);
        results.theta = -(S*n(d1)*sigma)/(2*sqrt(T)) + r*K*exp(-r*T)*N(-d2);
    }
    results.gamma = n(d1)/(S*sigma*sqrt(T));

    return results;
}

int main(int, char* []) {
    blackscholes res = black_scholes_formula(100, 100, 1, 0.05, 0.2, true);
    std::cout << res.value << std::endl;
    std::cout << res.delta<< std::endl;
    std::cout << res.theta<< std::endl;
    std::cout << res.gamma<< std::endl;
    return 0;
}
