/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*!
 Copyright (C) 2005, 2006, 2007, 2009 StatPro Italia srl

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*
#include <ql/qldefines.hpp>
#if !defined(BOOST_ALL_NO_LIB) && defined(BOOST_MSVC)
#  include <ql/auto_link.hpp>
#endif
#include <ql/instruments/vanillaoption.hpp>
#include <ql/pricingengines/vanilla/binomialengine.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanengine.hpp>
#include <ql/pricingengines/vanilla/analytichestonengine.hpp>
#include <ql/pricingengines/vanilla/baroneadesiwhaleyengine.hpp>
#include <ql/pricingengines/vanilla/bjerksundstenslandengine.hpp>
#include <ql/pricingengines/vanilla/batesengine.hpp>
#include <ql/pricingengines/vanilla/integralengine.hpp>
#include <ql/pricingengines/vanilla/fdblackscholesvanillaengine.hpp>
#include <ql/pricingengines/vanilla/mceuropeanengine.hpp>
#include <ql/pricingengines/vanilla/mcamericanengine.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanvasicekengine.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/utilities/dataformatters.hpp>
#include <ql/models/shortrate/onefactormodels/vasicek.hpp>

#include <iostream>
#include <iomanip>

using namespace QuantLib;

#if defined(QL_ENABLE_SESSIONS)
namespace QuantLib {

    ThreadKey sessionId() { return {}; }

}
#endif


int main(int, char* []) {

    try {

        std::cout << std::endl;

        // set up dates
        Calendar calendar = TARGET();
        Date todaysDate(15, May, 1998);
        Date settlementDate(17, May, 1998);
        Settings::instance().evaluationDate() = todaysDate;

        // our options
        Option::Type type(Option::Put);
        Real underlying = 36;
        Real strike = 40;
        Spread dividendYield = 0.00;
        Rate riskFreeRate = 0.06;
        Volatility volatility = 0.20;
        Date maturity(17, May, 1999);
        DayCounter dayCounter = Actual365Fixed();

        std::cout << "Option type = "  << type << std::endl;
        std::cout << "Maturity = "        << maturity << std::endl;
        std::cout << "Underlying price = "        << underlying << std::endl;
        std::cout << "Strike = "                  << strike << std::endl;
        std::cout << "Risk-free interest rate = " << io::rate(riskFreeRate)
                  << std::endl;
        std::cout << "Dividend yield = " << io::rate(dividendYield)
                  << std::endl;
        std::cout << "Volatility = " << io::volatility(volatility)
                  << std::endl;
        std::cout << std::endl;
        std::string method;
        std::cout << std::endl ;

        // write column headings
        Size widths[] = { 35, 14, 14, 14 };
        std::cout << std::setw(widths[0]) << std::left << "Method"
                  << std::setw(widths[1]) << std::left << "European"
                  << std::setw(widths[2]) << std::left << "Bermudan"
                  << std::setw(widths[3]) << std::left << "American"
                  << std::endl;

        std::vector<Date> exerciseDates;
        for (Integer i=1; i<=4; i++)
            exerciseDates.push_back(settlementDate + 3*i*Months);

        ext::shared_ptr<Exercise> europeanExercise(
                                         new EuropeanExercise(maturity));

        ext::shared_ptr<Exercise> bermudanExercise(
                                         new BermudanExercise(exerciseDates));

        ext::shared_ptr<Exercise> americanExercise(
                                         new AmericanExercise(settlementDate,
                                                              maturity));

        Handle<Quote> underlyingH(
            ext::shared_ptr<Quote>(new SimpleQuote(underlying)));

        // bootstrap the yield/dividend/vol curves
        Handle<YieldTermStructure> flatTermStructure(
            ext::shared_ptr<YieldTermStructure>(
                new FlatForward(settlementDate, riskFreeRate, dayCounter)));
        Handle<YieldTermStructure> flatDividendTS(
            ext::shared_ptr<YieldTermStructure>(
                new FlatForward(settlementDate, dividendYield, dayCounter)));
        Handle<BlackVolTermStructure> flatVolTS(
            ext::shared_ptr<BlackVolTermStructure>(
                new BlackConstantVol(settlementDate, calendar, volatility,
                                     dayCounter)));
        ext::shared_ptr<StrikedTypePayoff> payoff(
                                        new PlainVanillaPayoff(type, strike));
        ext::shared_ptr<BlackScholesMertonProcess> bsmProcess(
                 new BlackScholesMertonProcess(underlyingH, flatDividendTS,
                                               flatTermStructure, flatVolTS));

        // options
        VanillaOption europeanOption(payoff, europeanExercise);
        VanillaOption bermudanOption(payoff, bermudanExercise);
        VanillaOption americanOption(payoff, americanExercise);

        // Analytic formulas:

        // Black-Scholes for European
        method = "Black-Scholes";
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                                     new AnalyticEuropeanEngine(bsmProcess)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << "N/A"
                  << std::endl;

        //Vasicek rates model for European
        method = "Black Vasicek Model";
        Real r0 = riskFreeRate;
        Real a = 0.3;
        Real b = 0.3;
        Real sigma_r = 0.15;
        Real riskPremium = 0.0;
        Real correlation = 0.5;
        ext::shared_ptr<Vasicek> vasicekProcess(new Vasicek(r0, a, b, sigma_r, riskPremium));
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                new AnalyticBlackVasicekEngine(bsmProcess, vasicekProcess, correlation)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << "N/A"
                  << std::endl;

        // semi-analytic Heston for European
        method = "Heston semi-analytic";
        ext::shared_ptr<HestonProcess> hestonProcess(
            new HestonProcess(flatTermStructure, flatDividendTS,
                              underlyingH, volatility*volatility,
                              1.0, volatility*volatility, 0.001, 0.0));
        ext::shared_ptr<HestonModel> hestonModel(
                                              new HestonModel(hestonProcess));
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                                     new AnalyticHestonEngine(hestonModel)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << "N/A"
                  << std::endl;

        // semi-analytic Bates for European
        method = "Bates semi-analytic";
        ext::shared_ptr<BatesProcess> batesProcess(
            new BatesProcess(flatTermStructure, flatDividendTS,
                             underlyingH, volatility*volatility,
                             1.0, volatility*volatility, 0.001, 0.0,
                             1e-14, 1e-14, 1e-14));
        ext::shared_ptr<BatesModel> batesModel(new BatesModel(batesProcess));
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                                                new BatesEngine(batesModel)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << "N/A"
                  << std::endl;

        // Barone-Adesi and Whaley approximation for American
        method = "Barone-Adesi/Whaley";
        americanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                       new BaroneAdesiWhaleyApproximationEngine(bsmProcess)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << "N/A"
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;

        // Bjerksund and Stensland approximation for American
        method = "Bjerksund/Stensland";
        americanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                      new BjerksundStenslandApproximationEngine(bsmProcess)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << "N/A"
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;

        // Integral
        method = "Integral";
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                                             new IntegralEngine(bsmProcess)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << "N/A"
                  << std::endl;

        // Finite differences
        Size timeSteps = 801;
        method = "Finite differences";
        ext::shared_ptr<PricingEngine> fdengine =
            ext::make_shared<FdBlackScholesVanillaEngine>(bsmProcess,
                                                          timeSteps,
                                                          timeSteps-1);
        europeanOption.setPricingEngine(fdengine);
        bermudanOption.setPricingEngine(fdengine);
        americanOption.setPricingEngine(fdengine);
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << bermudanOption.NPV()
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;

        // Binomial method: Jarrow-Rudd
        method = "Binomial Jarrow-Rudd";
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                new BinomialVanillaEngine<JarrowRudd>(bsmProcess,timeSteps)));
        bermudanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                new BinomialVanillaEngine<JarrowRudd>(bsmProcess,timeSteps)));
        americanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                new BinomialVanillaEngine<JarrowRudd>(bsmProcess,timeSteps)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << bermudanOption.NPV()
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;
        method = "Binomial Cox-Ross-Rubinstein";
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                      new BinomialVanillaEngine<CoxRossRubinstein>(bsmProcess,
                                                                   timeSteps)));
        bermudanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                      new BinomialVanillaEngine<CoxRossRubinstein>(bsmProcess,
                                                                   timeSteps)));
        americanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                      new BinomialVanillaEngine<CoxRossRubinstein>(bsmProcess,
                                                                   timeSteps)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << bermudanOption.NPV()
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;

        // Binomial method: Additive equiprobabilities
        method = "Additive equiprobabilities";
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                new BinomialVanillaEngine<AdditiveEQPBinomialTree>(bsmProcess,
                                                                   timeSteps)));
        bermudanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                new BinomialVanillaEngine<AdditiveEQPBinomialTree>(bsmProcess,
                                                                   timeSteps)));
        americanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                new BinomialVanillaEngine<AdditiveEQPBinomialTree>(bsmProcess,
                                                                   timeSteps)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << bermudanOption.NPV()
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;

        // Binomial method: Binomial Trigeorgis
        method = "Binomial Trigeorgis";
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                new BinomialVanillaEngine<Trigeorgis>(bsmProcess,timeSteps)));
        bermudanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                new BinomialVanillaEngine<Trigeorgis>(bsmProcess,timeSteps)));
        americanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                new BinomialVanillaEngine<Trigeorgis>(bsmProcess,timeSteps)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << bermudanOption.NPV()
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;

        // Binomial method: Binomial Tian
        method = "Binomial Tian";
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                      new BinomialVanillaEngine<Tian>(bsmProcess,timeSteps)));
        bermudanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                      new BinomialVanillaEngine<Tian>(bsmProcess,timeSteps)));
        americanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                      new BinomialVanillaEngine<Tian>(bsmProcess,timeSteps)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << bermudanOption.NPV()
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;

        // Binomial method: Binomial Leisen-Reimer
        method = "Binomial Leisen-Reimer";
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
              new BinomialVanillaEngine<LeisenReimer>(bsmProcess,timeSteps)));
        bermudanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
              new BinomialVanillaEngine<LeisenReimer>(bsmProcess,timeSteps)));
        americanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
              new BinomialVanillaEngine<LeisenReimer>(bsmProcess,timeSteps)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << bermudanOption.NPV()
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;

        // Binomial method: Binomial Joshi
        method = "Binomial Joshi";
        europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                    new BinomialVanillaEngine<Joshi4>(bsmProcess,timeSteps)));
        bermudanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                    new BinomialVanillaEngine<Joshi4>(bsmProcess,timeSteps)));
        americanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
                    new BinomialVanillaEngine<Joshi4>(bsmProcess,timeSteps)));
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << bermudanOption.NPV()
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;

        // Monte Carlo Method: MC (crude)
        timeSteps = 1;
        method = "MC (crude)";
        Size mcSeed = 42;
        ext::shared_ptr<PricingEngine> mcengine1;
        mcengine1 = MakeMCEuropeanEngine<PseudoRandom>(bsmProcess)
            .withSteps(timeSteps)
            .withAbsoluteTolerance(0.02)
            .withSeed(mcSeed);
        europeanOption.setPricingEngine(mcengine1);
        // Real errorEstimate = europeanOption.errorEstimate();
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << "N/A"
                  << std::endl;

        // Monte Carlo Method: QMC (Sobol)
        method = "QMC (Sobol)";
        Size nSamples = 32768;  // 2^15

        ext::shared_ptr<PricingEngine> mcengine2;
        mcengine2 = MakeMCEuropeanEngine<LowDiscrepancy>(bsmProcess)
            .withSteps(timeSteps)
            .withSamples(nSamples);
        europeanOption.setPricingEngine(mcengine2);
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << "N/A"
                  << std::endl;

        // Monte Carlo Method: MC (Longstaff Schwartz)
        method = "MC (Longstaff Schwartz)";
        ext::shared_ptr<PricingEngine> mcengine3;
        mcengine3 = MakeMCAmericanEngine<PseudoRandom>(bsmProcess)
            .withSteps(100)
            .withAntitheticVariate()
            .withCalibrationSamples(4096)
            .withAbsoluteTolerance(0.02)
            .withSeed(mcSeed);
        americanOption.setPricingEngine(mcengine3);
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << "N/A"
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << americanOption.NPV()
                  << std::endl;

        // End test
        return 0;

    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "unknown error" << std::endl;
        return 1;
    }
}
*/

/*
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
*/


//#include "europeanoption.hpp"
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

