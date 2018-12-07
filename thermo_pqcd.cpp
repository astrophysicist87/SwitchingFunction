#include <cmath>
#include "constants.hpp"
#include "exceptions.hpp"
#include "thermo_pqcd.hpp"

using namespace std;


/*Lambda is the energy scale.  This is formula 9.5 from the pdg review
  of QCD.  (My Lambda is pdg's mu_R.)  3-loop, so set b3=0.
  SoftenConstant == 0 is standard QCD--non-zero values soften divergences
  at small temperature.  */
double alphas_3loop(double Lambda, double SoftenConstant)
{
	double LambdaMS = 290.0;   //MeV

	double Nf = 3.0;

	double t = log(SoftenConstant + pow(Lambda/LambdaMS, 2.0) );

	double b0 = (33.0 - 2.0*Nf)/(12.0*pi);
	double b1 = (153.0 - 19.0*Nf)/(24.0*pi*pi);
	double b2 = (2857.0 - 5033.0/9.0*Nf + 325.0/27.0*Nf*Nf)/(128.0*pi*pi*pi);

//To improve:  check for t == 0.0, because 1/t == infinity == trouble
//I could throw an exception (which would break execution), or return 0,
//or return a large fixed value.  Think further
	double alphas = 1.0/(b0*t)*(1.0 - b1*log(t)/(b0*b0*t) 
	 + (b1*b1*(log(t)*log(t) - log(t) - 1.0 ) + b0*b2)/( pow(b0, 4.0)*t*t ) 
	 - ( b1*b1*b1*( pow(log(t),3.0) - 5.0/2.0*log(t)*log(t) - 2.0*log(t)
	 + 1.0/2.0  ) + 3.0*b0*b1*b2*log(t) )/( pow(b0*b0*t, 3.0) )    );

	return alphas;
}


//-----------------------------------------------------------------------------


/*Computes the pressure of QCD at finite Temp, mu_B using order
	(alpha_s)^3*log(alpha_s) perturbative results.  Formula is from eqs 1-7 
	from arXiv:1212.1797v2 by Strickland which is simply quoting the result 
	calulated by Vuorinen.  For Vuorinen's derivation, see ref 17 of 
	arXiv:1212.1797v2.  (That's arXiv:hep-ph/0212283)  
	I used the result quoted by Strickland, instead of the 
	original work by Vuorinen, because Vuorinen did not finish dimensional 
	regularization in his expression.*/

double Total_Pressure_qgp_PQCD(double Temp, double mu_B, 
 double EnergyScaleConstant, double SoftenConstant)
{
	if(Temp <=0.0)
		throw TemperatureError("Total_Pressure_qgp_PQCD expects positive"
		 " temperature.");

	if(EnergyScaleConstant <=0.0)
		throw TemperatureError("Total_Pressure_qgp_PQCD expects positive"
		 " EnergyScaleConstant.");

	double mu_q = mu_B / 3.0;           //quark chemical potential
	double mu_qh = mu_q/(2*pi*Temp);    //chem potential hat
	double mu_qh2 = pow(mu_qh, 2.0);    
	double mu_qh4 = pow(mu_qh, 4.0);

	double Lambda = EnergyScaleConstant*sqrt( Temp*Temp + pow(mu_q/pi,2.0) );
	double Lambdah = Lambda/(2.0*pi*Temp);  //hat

	double Nf = 3.0;
	double Nc = 3.0;
	double CA = Nc;

	double F0, F2, F3, F4, F5, F6, Pressure;

	double alphas = alphas_3loop(Lambda, SoftenConstant);

	F0 = 1.0 + 21.0/32.0*Nf*(1.0 + 120.0/7.0*mu_qh2 + 240.0/7.0*mu_qh4);

	F2 = -15.0/4.0*(1.0 + 5.0/12.0*Nf*(1.0 + 72.0/5.0*mu_qh2 
		   + 144.0/5.0*mu_qh4));

	F3 = 30.0 * pow((1.0 + 1.0/6.0*(1.0 + 12.0*mu_qh2)*Nf), 1.5);

	F4 = 237.223 + (15.963 + 124.773*mu_qh2 - 319.849*mu_qh4)*Nf
	     -(0.415 + 15.926*mu_qh2 + 106.719*mu_qh4)*Nf*Nf
	     +135.0/2.0*(1.0 + 1.0/6.0*(1.0 + 12.0*mu_qh2)*Nf )*
	     log(alphas/pi*(1.0 + 1.0/6.0*(1.0 + 12.0*mu_qh2)*Nf))
	     -165.0/8.0*(1.0 + 5.0/12.0*
	     	(1.0 + 72.0/5.0*mu_qh2 + 144.0/5.0*mu_qh4)*Nf)*(1.0 - 2.0/33.0*Nf)
	     *log(Lambdah);

	F5 = -sqrt(1.0 + (1.0 + 12.0*mu_qh2)/6.0*Nf)*(799.149 + 
		   (21.963 - 136.33*mu_qh2 + 482.171*mu_qh4)*Nf 
		   + (1.926 + 2.0749*mu_qh2 - 172.07*mu_qh4)*Nf*Nf )
	     + 495.0/2.0*(1.0 + (1.0 + 12.0*mu_qh2)/6.0*Nf)*(1.0 - 2.0/33.0*Nf)*
	     log(Lambdah);

	F6 = -(659.175 + (65.888 - 341.489*mu_qh2 + 1446.514*mu_qh4)*Nf
			 + (7.653 + 16.225*mu_qh2 - 516.210*mu_qh4)*Nf*Nf
			 - 1485.0/2.0*(1.0 + (1.0 + 12.0*mu_qh2)/6.0*Nf)*
			 (1.0 - 2.0/33.0*Nf)*log(Lambdah)
			 )*log(alphas/pi*(1.0 + (1.0 + 12.0*mu_qh2)/6.0*Nf)*4.0*pi*pi )
		   - 475.587*log(alphas/pi*4.0*pi*pi*CA);

  Pressure = 8.0*pi*pi/45.0*pow(Temp, 4.0)*(F0 + F2*alphas/pi + 
             F3*pow(alphas/pi, 1.5) + F4*pow(alphas/pi, 2.0) + 
             F5*pow(alphas/pi, 2.5) + F6*pow(alphas/pi, 3.0)  );

  return Pressure;
}

//finite difference:  q_B = dP/dmu_B
double Total_Baryon_Dens_qgp_PQCD(double Temp, double mu_B, 
 double EnergyScaleConstant, double SoftenConstant)
{
	if(Temp <=0.0)
		throw TemperatureError("Total_Baryon_Dens_qgp_PQCD expects positive"
		 " temperature.");

	if(EnergyScaleConstant <=0.0)
		throw TemperatureError("Total_Baryon_Dens_qgp_PQCD expects positive"
		 " EnergyScaleConstant.");

	double dmu_B = 0.01;  //reasonable step size

	double Baryon = ( Total_Pressure_qgp_PQCD(Temp, mu_B + dmu_B, 
    EnergyScaleConstant, SoftenConstant) -
	 Total_Pressure_qgp_PQCD(Temp, mu_B - dmu_B, 
	  EnergyScaleConstant, SoftenConstant) ) / (2.0*dmu_B);

	 return Baryon;
}

//This is a  bit misleading, use with caution.  It's density of 
//quarks-anti-quarks.
//Quark number density = 3 * baryon density
double Total_Number_Dens_qgp_PQCD(double Temp, double mu_B, 
 double EnergyScaleConstant, double SoftenConstant)
{
	if(Temp <=0.0)
		throw TemperatureError("Total_Number_Dens_qgp_PQCD expects positive"
		 " temperature.");

	if(EnergyScaleConstant <=0.0)
		throw TemperatureError("Total_Number_Dens_qgp_PQCD expects positive"
		 " EnergyScaleConstant.");

	return 3.0 * Total_Baryon_Dens_qgp_PQCD(Temp, mu_B, 
   EnergyScaleConstant, SoftenConstant);
}


//Finite difference: s = dP/dT
double Total_Entropy_Dens_qgp_PQCD(double Temp, double mu_B, 
 double EnergyScaleConstant, double SoftenConstant)
{
	if(Temp <=0.0)
		throw TemperatureError("Total_Entropy_Dens_qgp_PQCD expects positive"
		 " temperature.");

	if(EnergyScaleConstant <=0.0)
		throw TemperatureError("Total_Entropy_Dens_qgp_PQCD expects positive"
		 " EnergyScaleConstant.");

	double T = Temp;
	double dT = 0.01;   //reasonable step size

	double Entropy = ( Total_Pressure_qgp_PQCD(T + dT, mu_B, 
    EnergyScaleConstant, SoftenConstant) -
	 Total_Pressure_qgp_PQCD(T - dT, mu_B, 
	 	EnergyScaleConstant, SoftenConstant) ) / (2.0*dT);

	 return Entropy;
}


double Total_Energy_Dens_qgp_PQCD(double Temp, double mu_B, 
 double EnergyScaleConstant, double SoftenConstant)
{
	if(Temp <=0.0)
		throw TemperatureError("Total_Energy_Dens_qgp_PQCD expects positive"
		 " temperature.");

	if(EnergyScaleConstant <=0.0)
		throw TemperatureError("Total_Energy_Dens_qgp_PQCD expects positive"
		 " EnergyScaleConstant.");

	double Press = Total_Pressure_qgp_PQCD(Temp, mu_B, 
              EnergyScaleConstant, SoftenConstant);

	double Baryon = Total_Baryon_Dens_qgp_PQCD(Temp, mu_B, 
              EnergyScaleConstant, SoftenConstant);

	double Entropy = Total_Entropy_Dens_qgp_PQCD(Temp, mu_B, 
              EnergyScaleConstant, SoftenConstant);

	return (Temp*Entropy - Press + mu_B*Baryon); //from thermo identity
}
