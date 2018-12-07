#ifndef CROSSOVER_LINE_H
#define CROSSOVER_LINE_H

#include <array>
#include <vector>
#include "thermo_crossover.hpp"


double Delta_P(double R, void * params);

double Delta_P_muc(double mu_B, void * params);

double Get_R(double psi, thermparams parameters);

double Get_k2(double Rc, thermparams parameters);

double Get_muc(double Tc, thermparams parameters);

std::vector<double> Get_Rs_quartic(double psi, double Rc, thermparams parameters);

#endif