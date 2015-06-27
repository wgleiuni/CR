#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "DefClass.h"

int main (int argc, char *argv[])
{
    int I,Mode,numdt,num;
    double Kx,Ky,F,Omega,h;
    if (argc>1)
    {
        I=(int)atof(argv[1]);
        Mode=(int)atof(argv[2]);
        Kx=(double)atof(argv[3]);
        Ky=(double)atof(argv[4]);
        F=(double)atof(argv[5]);
        Omega=(double)atof(argv[6]);
        h=(double)atof(argv[7]);
        numdt=(int)atof(argv[8]);
        num=(int)atof(argv[9]);
    }

    CR one(I,Mode,Kx,Ky,F,Omega,h,numdt,num);
    one.go();

    return 0;
}
