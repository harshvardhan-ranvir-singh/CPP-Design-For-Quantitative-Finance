#include <ql/option.hpp>
#include <ql/handle.hpp>
#include <ql/quote.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/yield/zerocurve.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>
#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <boost/make_shared.hpp>
#include <iostream>
#include <ql/math/distributions/normaldistribution.hpp>

using std::cout;
using std::endl;
using std::vector;
using boost::shared_ptr;
using boost::make_shared;
using namespace QuantLib;

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

class EuropeanOption : public Instrument {

  public:
    EuropeanOption(Real strike, Option::Type, const Date& exerciseDate,
                   const Handle<Quote>& u,
                   const Handle<YieldTermStructure>& r,
                   const Handle<BlackVolTermStructure>& sigma);
    bool isExpired() const;
    void performCalculations() const;
    Real delta() const;
    Real gamma() const;
    Real theta() const;
  private:
    Real strike_;
    Option::Type type_;
    const Date& exerciseDate_;
    const Handle<Quote>& u_;
    const Handle<YieldTermStructure>& r_;
    const Handle<BlackVolTermStructure>& sigma_;
    mutable Real delta_;
    mutable Real gamma_;
    mutable Real theta_;

};

EuropeanOption ::EuropeanOption(Real strike, Option::Type type, const Date& exerciseDate, const Handle<Quote>& u,
                                const Handle<YieldTermStructure>& r, const Handle<BlackVolTermStructure>& sigma)
: strike_(strike), type_(type), exerciseDate_(exerciseDate), u_(u), r_(r), sigma_(sigma){
    registerWith(u_);
    registerWith(r_);
    registerWith(sigma_);
}

bool EuropeanOption::isExpired() const {
    Date today = Settings::instance().evaluationDate();
    return today >= exerciseDate_;
}

void EuropeanOption::performCalculations() const {
    DayCounter dayCounter = r_->dayCounter();
    Date today = Settings::instance().evaluationDate();
    Time T = dayCounter.yearFraction(today, exerciseDate_);
    Rate r = r_->zeroRate(T, Continuous);
    Volatility sigma = sigma_->blackVol(T, strike_);
    blackscholes results = black_scholes_formula(u_->value(), strike_, T, r, sigma,
                                             type_ == Option::Call);

    NPV_ = results.value;
    delta_ = results.delta;
    gamma_ = results.gamma;
    theta_ = results.theta;
}

Real EuropeanOption::delta() const {

    // before returning delta_, we ensure that it is calculated: checking calculated flag (whether instrument is
    // upto date) or it's going to trigger the calculation in case its not.
    calculate();
    return delta_;
}

Real EuropeanOption::gamma() const {

    // before returning gamma_, we ensure that it is calculated: checking calculated flag (whether instrument is
    // upto date) or it's going to trigger the calculation in case its not.
    calculate();
    return gamma_;
}

Real EuropeanOption::theta() const {

    // before returning theta_, we ensure that it is calculated: checking calculated flag (whether instrument is
    // upto date) or it's going to trigger the calculation in case its not.
    calculate();
    return theta_;
}

int main() {
    Date today(1, May, 2022);
    Settings::instance().evaluationDate() = today;

    Real strike = 100.0;
    Option::Type type = Option::Put;
    Date exerciseDate = today + 360*3;

    shared_ptr<SimpleQuote> u = make_shared<SimpleQuote>(120.0);
    Handle<Quote> U(u);

    shared_ptr<YieldTermStructure> r =
        make_shared<FlatForward>(today, 0.01, Actual360());
    RelinkableHandle<YieldTermStructure> R(r);

    shared_ptr<BlackVolTermStructure> sigma =
        make_shared<BlackConstantVol>(today, TARGET(), 0.20, Actual360());
    Handle<BlackVolTermStructure> Sigma(sigma);

    EuropeanOption option(strike, type, exerciseDate, U, R, Sigma);

    cout << "value: " << option.NPV() << endl;
    cout << "delta: " << option.delta() << endl;
    cout << "gamma: " << option.gamma() << endl;
    cout << "theta: " << option.theta() << endl;

    u->setValue(115.0);
    cout << "value after u decreases: " << option.NPV() << endl;

    /*
        value: 6.12402
        delta: -0.215897
        gamma: 0.00704601
        theta: -1.70893
        value after u decreases: 7.29565
     */
    return 0;
}
