/*
 * =====================================================================================
 *
 *       Filename:  DefClass.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/05/2015 07:56:39 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dr. Fritz Mehner (fgm), mehner.fritz@fh-swf.de
 *   Organization:  FH Südwestfalen, Iserlohn
 *
 * =====================================================================================
 */

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include "DefClass.h"

Particle::Particle()
{

}

double Particle::E()
{
    double out;
    switch (Mode_)
    {
        case 1:
            out=sqrt(pow(Kx_,4.0)+pow(Ky_,2.0)); break;
        case 2:
            out=pow(Kx_,2.0)+pow(Ky_,2.0); break;
        case 3:
            out=sqrt(pow(Kx_,2.0)+pow(Ky_,2.0)); break;
    }
    return out;
}

RK4::RK4()
{
    t_=0.0;
}

double RK4::dK(int N,double tKx,double tKy,double t)
{
    double out,temp;
    switch (Mode_) {
        case 1:
            temp=sqrt(pow(tKx,4.0)+pow(tKy,2.0));
            switch (N) {
                case 1:
                    out=-tKy/2.0/temp-F_*cos(Omega_*t);
                    break;
                case 2:
                    out=pow(tKx,3.0)/temp-F_*sin(Omega_*t);
                    break;
            }
            break;
        case 2:
            switch (N) {
                case 1:
                    out=-tKy-F_*cos(Omega_*t);
                    break;
                case 2:
                    out=tKx-F_*sin(Omega_*t);
                    break;
            }
            break;
        case 3:
            temp=sqrt(pow(tKx,2.0)+pow(tKy,2.0));
            switch (N) {
                case 1:
                    out=-tKy/temp-F_*cos(Omega_*t);
                    break;
                case 2:
                    out=tKx/temp-F_*sin(Omega_*t);
                    break;
            }
            break;
    }
    return out;
}

void RK4::onestep()
{
    int i;
    for (i=0;i<4;i++)
    {
        if (i==0) 
        {
            k_[i]=dK(1,Kx_,Ky_,t_);
            l_[i]=dK(2,Kx_,Ky_,t_);
        }
        else if (i==1 || i==2)
        {
            k_[i]=dK(1,Kx_+h_/2.0*k_[i-1],Ky_+h_/2.0*l_[i-1],t_+h_/2.0);
            l_[i]=dK(2,Kx_+h_/2.0*k_[i-1],Ky_+h_/2.0*l_[i-1],t_+h_/2.0);
        }
        else if (i==3)
        {
            k_[i]=dK(1,Kx_+h_*k_[i-1],Ky_+h_*l_[i-1],t_+h_);
            l_[i]=dK(2,Kx_+h_*k_[i-1],Ky_+h_*l_[i-1],t_+h_);
        }
    }
    Kx_=Kx_+h_/6.0*(k_[0]+2.0*k_[1]+2.0*k_[2]+k_[3]);
    Ky_=Ky_+h_/6.0*(l_[0]+2.0*l_[1]+2.0*l_[2]+l_[3]);
    t_=t_+h_;
}

CR::CR(int I,int Mode,double Kx,double Ky,double F,double Omega,double h,int numdt,int num)
{
    Mode_=Mode;
    Kx_=Kx;
    Ky_=Ky;
    F_=F;
    Omega_=Omega;
    numdt_=numdt;
    num_=num;
    h_=h;

    char filename[20];
    sprintf(filename,"D%d.txt",I);

    out_.open(filename,std::ostream::out);
}

void CR::record()
{
    out_ << t_/(2.0*M_PI) << "\t" << Kx_ << "\t" << Ky_ << "\t" << RK4::Particle::E() << std::endl;
}

void CR::go()
{
    int i;
    for (i=0;i<numdt_*num_;i++)
    {
        if (i%numdt_==0)
        {
            record();
        }
        RK4::onestep();
//        std::cout << t_ << "\t" << Kx_ << "\t" << Ky_ << std::endl;
    }
}
