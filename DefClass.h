/*
 * =====================================================================================
 *
 *       Filename:  DefClass.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/05/2015 07:43:50 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dr. Fritz Mehner (fgm), mehner.fritz@fh-swf.de
 *   Organization:  FH Südwestfalen, Iserlohn
 *
 * =====================================================================================
 */

#ifndef DEFCLASS_H
#define DEFCLASS_H
#include <fstream>
class Particle
{
    public:
        Particle();
        double Kx_,Ky_;
    protected:
        int Mode_;
        double F_,Omega_,t_;
        double E();
};

class RK4 : protected Particle
{
    public:
        RK4();
    protected:
        void onestep();
        void lya();
        double t_,h_,lya1_,lya2_;
        int lN_;
    private:
        void orth();
        double dK(int,double,double,double);
        double dD(int,double,double);
        double k_[4],l_[4],lKx_[2],lKy_[2];
};

class CR : protected RK4
{
    public:
        CR(int,int,double,double,double,double,double,int,int);
        void go();
        void record();
    private:
        std::ofstream out_,out_lya_;
        int numdt_,num_;
};

#endif
