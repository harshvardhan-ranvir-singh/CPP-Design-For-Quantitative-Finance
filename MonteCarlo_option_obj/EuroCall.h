//
// Created by Harshvardhan Singh on 23/10/2021.
//

#ifndef MONTECARLO_OPTION_OBJ_EUROCALL_H
#define MONTECARLO_OPTION_OBJ_EUROCALL_H

class EuroCall
{
public:
    EuroCall(double X, double T);
    double ComputePO(double S);

private:
    double X_;
    double T_;
};

#endif //MONTECARLO_OPTION_OBJ_EUROCALL_H
