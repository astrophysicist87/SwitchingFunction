#ifndef THERMO_PQCD_H
#define THERMO_PQCD_H

#include "constants.hpp"

double alphas_3loop(double Lambda, double SoftenConstant = 0.0);



double Total_Pressure_qgp_PQCD(double Temp, double mu_B, 
 double EnergyScaleConstant = 2*pi, double SoftenConstant = 0.0);

double Total_Energy_Dens_qgp_PQCD(double Temp, double mu_B, 
 double EnergyScaleConstant = 2*pi, double SoftenConstant = 0.0);

double Total_Baryon_Dens_qgp_PQCD(double Temp, double mu_B, 
 double EnergyScaleConstant = 2*pi, double SoftenConstant = 0.0);

double Total_Number_Dens_qgp_PQCD(double Temp, double mu_B, 
 double EnergyScaleConstant = 2*pi, double SoftenConstant = 0.0);

double Total_Entropy_Dens_qgp_PQCD(double Temp, double mu_B, 
 double EnergyScaleConstant = 2*pi, double SoftenConstant = 0.0);




#endif