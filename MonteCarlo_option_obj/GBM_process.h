//
// Created by Harshvardhan Singh on 23/10/2021.
//

#ifndef MONTECARLO_OPTION_OBJ_GBM_PROCESS_H
#define MONTECARLO_OPTION_OBJ_GBM_PROCESS_H

class GBM_process
{
public:
    GBM_process(double S_0, double sig,  double r, double dt);
    double Next_S(double S);

    double GetS0() const {return S_0_;}

private:
    double S_0_;
    double sig_;
    double r_;

    double dt_;
    double drift_;
    double sgrt_;
};

#endif //MONTECARLO_OPTION_OBJ_GBM_PROCESS_H
