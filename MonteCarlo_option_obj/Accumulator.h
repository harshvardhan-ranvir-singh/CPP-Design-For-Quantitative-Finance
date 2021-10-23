//
// Created by Harshvardhan Singh on 23/10/2021.
//

#ifndef MONTECARLO_OPTION_OBJ_ACCUMULATOR_H
#define MONTECARLO_OPTION_OBJ_ACCUMULATOR_H


class Accumulator
{
public:
    Accumulator(double r, double T);
    void AddValue(double payoff);

    double GetOptionValue();
    double GetSE();


private:
    double acc_vals_;
    double acc_squs_;

    double discount_;

    long M_;
};

#endif //MONTECARLO_OPTION_OBJ_ACCUMULATOR_H
