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
    lya1_=0.0;
    lya2_=0.0;
    lKx_[0]=(double)rand();
    lKy_[0]=(double)rand();
    double norm;
    norm=sqrt(pow(lKx_[0],2.0)+pow(lKy_[0],2.0));
    lKx_[0]=lKx_[0]/norm;
    lKy_[0]=lKy_[0]/norm;
    lKx_[1]=(double)rand();
    lKy_[1]=(double)rand();
    orth();
    norm=sqrt(pow(lKx_[1],2.0)+pow(lKy_[1],2.0));
    lKx_[1]=lKx_[1]/norm;
    lKy_[1]=lKy_[1]/norm;
    lN_=0;
}

void RK4::orth()
{
    double temp;
    temp=(lKx_[0]*lKx_[1]+lKy_[0]*lKy_[1])/(lKx_[0]*lKx_[0]+lKy_[0]*lKy_[0]);
    lKx_[1]=lKx_[1]-temp*lKx_[0];
    lKy_[1]=lKy_[1]-temp*lKy_[0];
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

double RK4::dD(int N,double dKx,double dKy)
{
    double out,temp;
    switch (Mode_) {
        case 1:
            temp=sqrt(pow(Kx_,4.0)+pow(Ky_,2.0));
            switch (N) {
                case 1:
                    out=pow(Kx_,3.0)*Ky_/pow(temp,3.0)*dKx+0.5*(pow(Ky_,2.0)/pow(temp,3.0)-1/temp)*dKy;
                    break;
                case 2:
                    out=(-2.0*pow(Kx_,6.0)/pow(temp,3.0)+3.0*pow(Kx_,2.0)/temp)*dKx-pow(Kx_,3.0)*Ky_/pow(temp,3.0)*dKy;
                    break;
            }
            break;
        case 2:
            switch (N) {
                case 1:
                    out=-dKy;
                    break;
                case 2:
                    out=dKx;
                    break;
            }
            break;
        case 3:
            temp=sqrt(pow(Kx_,2.0)+pow(Ky_,2.0));
            switch (N) {
                case 1:
                    out=Kx_*Ky_/pow(temp,3.0)*dKx+(pow(Ky_,2.0)/pow(temp,3.0)-1/temp)*dKy;
                    break;
                case 2:
                    out=(1/temp-pow(Kx_,2.0)/pow(temp,3.0))*dKx-Kx_*Ky_/pow(temp,3.0)*dKy;
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

void RK4::lya()
{
    int i,j;
    for (j=0;j<2;j++)
    {
        for (i=0;i<4;i++)
        {
            if (i==0)
            {
                k_[i]=dD(1,lKx_[j],lKy_[j]);
                l_[i]=dD(2,lKx_[j],lKy_[j]);
            }
            else if (i==1 || i==2)
            {
                k_[i]=dD(1,lKx_[j]+h_/2.0*k_[i-1],lKy_[j]+h_/2.0*l_[i-1]);
                l_[i]=dD(2,lKx_[j]+h_/2.0*k_[i-1],lKy_[j]+h_/2.0*l_[i-1]);
            }
            else if (i==3)
            {
                k_[i]=dD(1,lKx_[j]+h_*k_[i-1],lKy_[j]+h_*l_[i-1]);
                l_[i]=dD(2,lKx_[j]+h_*k_[i-1],lKy_[j]+h_*l_[i-1]);
            }
        }
        lKx_[j]=lKx_[j]+h_/6.0*(k_[0]+2.0*k_[1]+2.0*k_[2]+k_[3]);
        lKy_[j]=lKy_[j]+h_/6.0*(l_[0]+2.0*l_[1]+2.0*l_[2]+l_[3]);
    }
    lya1_=lya1_+log(sqrt(pow(lKx_[0],2.0)+pow(lKy_[0],2.0)));
    double norm;
    norm=sqrt(pow(lKx_[0],2.0)+pow(lKy_[0],2.0));
    lKx_[0]=lKx_[0]/norm;
    lKy_[0]=lKy_[0]/norm;

    orth();
    lya2_=lya2_+log(sqrt(pow(lKx_[1],2.0)+pow(lKy_[1],2.0)));
    norm=sqrt(pow(lKx_[1],2.0)+pow(lKy_[1],2.0));
    lKx_[1]=lKx_[1]/norm;
    lKy_[1]=lKy_[1]/norm;
    lN_++;
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

    char filename[20],filename1[20];
    sprintf(filename,"D%d.txt",I);
    sprintf(filename1,"L%d.txt",I);

    out_.open(filename,std::ostream::out);
    out_lya_.open(filename1,std::ostream::out);
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
        if (i>(numdt_*num_/2)-1)
        {
            RK4::lya();
//            if (i%numdt_==0)
//            {
//                out_lya_ << lya1_ << "\t" << lya2_ << std::endl;
//            }
        }
    }
    out_lya_ << lya1_/(h_*lN_) << "\t" << lya2_/(h_*lN_) << std::endl;
}
