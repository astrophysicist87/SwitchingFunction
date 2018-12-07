#ifndef THERMO_CROSSOVER_H
#define THERMO_CROSSOVER_H

#include <vector>
#include <array>
#include "constants.hpp"
#include "particle.hpp"



struct thermparams {
    double psic;
    double a;
    double b;
    int InterpExp;
    double Escl;
    double Soft;
    double eps0_4;
    std::vector<double> Ks;
    MODEL_TYPE model;
    std::vector<Particle>  ParticleList;
};
//==========================================================================

std::vector<double> Get_Thermo(double Temp, double mu_B, double R,
	double dRdpsi, thermparams parameters);

double Get_c2(double Temp, double mu_B, double Rc, thermparams parameters);

#endif
