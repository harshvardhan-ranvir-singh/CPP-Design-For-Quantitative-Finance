//
// Created by Harshvardhan Singh on 12/06/2021.
//


//   rv_library.cpp
//   namespace RandomVariableStatisticalFunctions implementation file


#include "rv_library.h"


#include <cmath>
#include <cstdlib>
#include <limits>
#include <ctime>



//  double GetNormalVariate() Polar rejection.

double RandomVariableStatisticalFunctions::GetNormalVariate()
{
    static bool Spare_normal_flag = false;
    static double Spare_normal;

    double rsq, v1, v2;

    if (Spare_normal_flag == false)
    {
        do
        {
            v1 = 2.0*my_ran2() - 1.0;			//Put in range (-1,1)
            v2 = 2.0*my_ran2() - 1.0;

            rsq = v1*v1 + v2*v2;

        } while (rsq >= 1.0 || rsq == 0.0);		//Reject if outside unit circle

        double fac = std::sqrt(-2.0*std::log(rsq)/rsq);

        Spare_normal = v1*fac; 					//Generates two normals
        Spare_normal_flag = true;  				//Store one,  return the other

        return v2*fac;

    }
    else
    {
        Spare_normal_flag = false;
        return Spare_normal;
    }
}


//  double my_ran2()   modified from Numerical recipes to use statics
//  Generates a U[0,1] with period ~2.3*10^18.
//	(L'Ecuyer 1988 p747)


double RandomVariableStatisticalFunctions::my_ran2()   //modified from Numerical recipes
{
    static const long NTAB = 32;

    static const long IM1 = 2147483563;	static const long IA1 = 40014;
    static const long IQ1 = 53668;		static const long IR1 = 12211;

    static const long IM2 = 2147483399;	static const long IA2 = 40692;
    static const long IQ2 = 52774;		static const long IR2 = 3791;

    static const long IMM1 = IM1 - 1;
    static const long NDIV = 1 + IMM1/NTAB;
    static const double EPS = 3.0e-16;
    static const double RNMX = 1.0 - EPS;
    static const double AM = 1.0/double(IM1);

    static long idum2;
    static long iy;
    static long iv[NTAB];
    static long idum;

    static bool NOT_FIRST_TIME = false;

    if (NOT_FIRST_TIME == false)		//run for first time
    {
        srand(time(NULL));	//initialise ran2
        idum = -rand();

        idum = (idum == 0)? 1 : -idum;
        idum2 = idum;

        for (long j = NTAB + 7; j >= 0; --j)
        {
            long k = idum/IQ1;
            idum = IA1*(idum - k*IQ1) - k*IR1;
            if (idum < 0) idum += IM1;

            if (j < NTAB) iv[j] = idum;
        }

        iy = iv[0];

        NOT_FIRST_TIME = true;
    }

    long k = idum/IQ1;
    idum = IA1*(idum - k*IQ1) - k*IR1;
    if (idum < 0) idum += IM1;

    k = idum2/IQ2;
    idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
    if (idum2 < 0) idum2 += IM2;

    long j = iy/NDIV;
    iy = iv[j] - idum2;
    iv[j] = idum;

    if (iy < 1) iy += IMM1;

    double temp = AM*iy;

    return (temp > RNMX) ? RNMX : temp;
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//  End
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
