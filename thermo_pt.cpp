#include <string>
#include <stdexcept>
#include <cmath>
#include <vector>

#include <gsl/gsl_errno.h>          //for modifying gsl's error handling
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_result.h>      //contains a struct for returning results 
                                    //of gsl functions

#include "constants.hpp"
#include "exceptions.hpp"
#include "particle.hpp"
#include "particlelist.hpp"
#include "thermo_pt.hpp"


//ignore errors if numerical integration has problems
#define SUPPRESS_INTEGRATION_ERRORS   

using namespace std;

//-----------------------------------------------------------------------------
//Code to compute (partial) thermal properties of one species of gas:


double equilib_distfcn(double Energy, double Temp, double mu_particle, 
 QUANTUM_TYPE quantumtype)
{
	//Returns the value of Bose-Einstein, Fermi, or Boltzmann distribution 
	//function. Depends on the type of particle, the energy, temp, and particle
	//chemical potential.

	if(Temp <=0.0)
		throw TemperatureError("equilib_distfcn expects positive temperature.");

	double a;
	if(quantumtype == BOSON)
		a = -1.0;
	else if(quantumtype == FERMION)
		a = 1.0;
	else if(quantumtype == CLASSICAL)
		a = 0.0;
	else
		throw QuantumTypeError("equilib_distfcn expects valid QUANTUM_TYPE"
		 " object.");

	double EXP_ARG_MAX = 700.0;   //needed to avoid overflow in gsl_sf_exp
	double EXP_ARG_MIN = -700.0;
	
	double argument = (Energy-mu_particle)/Temp;

	//Evaluate distribution fcn
	if(argument > EXP_ARG_MAX)  //avoid overflow
		return 0.0;
	else if(argument < EXP_ARG_MIN)  //avoid underflow
		return (1.0/a);   //a==0 could cause problems, but can avoid this by not
		                  //using classical stats at low Temp->use quantum stats
	else
		return 1.0/(gsl_sf_exp(argument) + a);
}


struct partial_Pressure_Params
{
	//store parameters passed to Partial_Pressure_Integrand() via void *p
	double Temp;
	double mass;
	double mu_particle;
	QUANTUM_TYPE quantumtype;
};

double Partial_Pressure_Integrand(double x, void *p)
{
	//Returns dimensionless integrand as a fcn of x = E/mass
	//and other parameters (like Temp, mass, particle chemical potential,...)
	//When integrated over x, this integral will compute the 
	//dimensionless component of partial pressure.

	struct partial_Pressure_Params * params = 
	 (struct partial_Pressure_Params *) p;

	double Temp = params->Temp;
	double mass = params->mass;
	double mu_particle = params->mu_particle;
	QUANTUM_TYPE quantumtype = params->quantumtype;

	//x = Energy/mass
	double E = mass*x;

	if(E < mass)
		throw ExVolValueError("Partial_Pressure_Integrand() requires"
		 "x = E/mass >= 1.");

	return pow(x*x - 1.0, 1.5)*
	 equilib_distfcn(E, Temp, mu_particle, quantumtype);
}



double Partial_Pressure_pt(double Temp, double mu_particle,
 const Particle& particle)
{
	//Compute partial pressure of a particle species in PT model, 
	//identified by particle object, at temperature Temp and 
	//particle chemical potential mu_particle.

	if(Temp <=0.0)
		throw TemperatureError("Partial_Pressure_pt expects positive"
		 " temperature.");
  
	double mass = particle.GetMass();
	double degen = particle.GetDegen();
	QUANTUM_TYPE quantumtype = particle.GetQuantumType();

	if(quantumtype == CLASSICAL)
	{
		double mu = mu_particle;
		return 1.0/(2.0*pi*pi)*degen*(mass*mass)*(Temp*Temp)*
		  gsl_sf_bessel_Kn(2, mass/Temp)*gsl_sf_exp(mu/Temp);
	}
	else
	{
		double lower_limit = 1.0;  //integrate x=E/mass from 1 to infinity

		//do some setup for gsl
		const int MaxNumIntervals = 1000;
		gsl_integration_workspace *work_ptr = 
		 gsl_integration_workspace_alloc(MaxNumIntervals);
		double abs_error = 1.0e-8;	
		double rel_error = 1.0e-8;	

		double result, error;		//where we will store results

		//create structure holding parameters to pass to integrand function:
		struct partial_Pressure_Params params;
		params.Temp = Temp;
		params.mass = mass;
		params.mu_particle = mu_particle;
		params.quantumtype = quantumtype;
		void *params_ptr = &params;

		//set up reference to the integrand function to be integrated over:
		gsl_function My_function;
		My_function.function = &Partial_Pressure_Integrand;
		My_function.params = params_ptr;


		//Disable gsl's error handling so it won't crash the code if
		//(for instance) an integral converges slowly. We'll check status after.
		gsl_set_error_handler_off();

		//invoke gls's quadrature-adaptive routine, integrate from lower limit 
		//to infinity
		int return_code;
		return_code=0;
		return_code = gsl_integration_qagiu(&My_function, lower_limit,
		 abs_error, rel_error, MaxNumIntervals, work_ptr, &result, &error);
		//Turn gsl's default error handler back on
		gsl_set_error_handler(NULL);

		gsl_integration_workspace_free(work_ptr);


		/*
		The commented-out lines below could be uncommented to
		allow some non-critical errors:

			For instance, the integral may converge slowly 
			(return_code == GSL_EDIVERGE), especially when the partial pressure
			is very tiny and negligible (result < 1.0).  Who cares if convergence
			is poor on a result that contributes negliglibly? 
			So we might ignore those errors.

			Also, there may be rounding errors (return_code == GSL_EROUND).
			Those may be problematic, but we may decide to ignore them and use
			the best integration result without crashing.

		Otherwise, throw an exception to let user know an unacceptable
		gsl error occured.*/

//#ifndef SUPPRESS_INTEGRATION_ERRORS 
 
		if(return_code != 0) //error occured
		 //if(return_code != GSL_EROUND)  
			//if( !( (return_code == GSL_EDIVERGE) && (result < 1.0) ))
			{
				std::string msg = "gsl fcn gsl_integration_qagiu() experienced "
				 "unexpected error in fcn Partial_Pressure_pt(). gsl error code was ";
				msg += std::to_string(return_code);
				msg += ". Temp = " + std::to_string(Temp);
				msg += ". mu_particle = " + std::to_string(mu_particle);
				msg += ". particle = " + particle.GetName();
				throw GenericError(msg);
			}
			
//#endif

		if(result < 0.0)  //don't allow negative pressures! That's a fail!
		{
			std::string msg = "Error: Partial_Pressure_pt() computed a "
			                  "negative pressure.";
			msg += " Temp = " + std::to_string(Temp);
			msg += " mu_particle = " + std::to_string(mu_particle);
			msg += " particle = " + particle.GetName();
			throw GenericError(msg);
		}

		//multiply integral by degeneracy, phase space factors, and units 
		//which convert integral into pressure 
		result *= degen*pow(mass, 4.0)/(6.0*pi*pi);

		return result;
	}
}

//------------------------------------------------------


double Partial_Number_Dens_Integrand(double x, void *p)
{
	//Returns dimensionless integrand as a fcn of x = E/mass
	//and other parameters (like Temp, mass, particle chemical potential,...)
	//When integrated over x, this integral will compute the 
	//dimensionless component of partial number density.

	struct partial_Pressure_Params * params =
	 (struct partial_Pressure_Params *) p;

	double Temp = params->Temp;
	double mass = params->mass;
	double mu_particle = params->mu_particle;
	QUANTUM_TYPE quantumtype = params->quantumtype;

	//x = Energy/mass
	double E = mass*x;

	if(E < mass)
		throw ExVolValueError("Partial_Number_Dens_Integrand() requires"
		 "x = E/mass >= 1.");

	return x*sqrt(x*x - 1.0)*
	 equilib_distfcn(E, Temp, mu_particle, quantumtype);
}


double Partial_Number_Dens_pt(double Temp, double mu_particle,
 const Particle& particle)
{
	//Compute partial number density of a particle species, identified by 
	//particle object at temperature Temp and particle chemical potential 
	//mu_particle.

	if(Temp <=0.0)
		throw TemperatureError("Partial_Number_Dens_pt expects positive" 
		 "temperature.");

	double mass = particle.GetMass();
	double degen = particle.GetDegen();
	QUANTUM_TYPE quantumtype = particle.GetQuantumType();

	if(quantumtype == CLASSICAL)
	{
		//valid for classical stats only:
		return Partial_Pressure_pt(Temp, mu_particle, particle)/Temp;
	}
	else
	{
		double lower_limit = 1.0;  //integrate x=E/mass from 1 to infinity

		//do some setup for gsl
		const int MaxNumIntervals = 1000;
		gsl_integration_workspace *work_ptr = 
		 gsl_integration_workspace_alloc(MaxNumIntervals);
		double abs_error = 1.0e-8;	
		double rel_error = 1.0e-8;	

		double result, error;		//where we will store results

		//create structure holding parameters to pass to integrand function:
		struct partial_Pressure_Params params;
		params.Temp = Temp;
		params.mass = mass;
		params.mu_particle = mu_particle;
		params.quantumtype = quantumtype;
		void *params_ptr = &params;

		//set up reference to the integrand function to be integrated over:
		gsl_function My_function;
		My_function.function = &Partial_Number_Dens_Integrand;
		My_function.params = params_ptr;


		//Disable gsl's error handling so it won't crash the code if
		//(for instance) an integral converges slowly. We'll check status after.
		gsl_set_error_handler_off();

		//invoke gls's quadrature-adaptive routine, integrate from lower limit to
		//infinity
		int return_code;
		return_code=0;
		return_code = gsl_integration_qagiu (&My_function, lower_limit, 
		 abs_error, rel_error, MaxNumIntervals, work_ptr, &result, &error);

		//Turn gsl's default error handler back on
		gsl_set_error_handler(NULL);

		gsl_integration_workspace_free(work_ptr);


		/*
		The commented-out lines below could be uncommented to
		allow some non-critical errors:

			For instance, the integral may converge slowly 
			(return_code == GSL_EDIVERGE), especially when the partial pressure
			is very tiny and negligible (result < 1.0).  Who cares if convergence
			is poor on a result that contributes negliglibly? 
			So we might ignore those errors.

			Also, there may be rounding errors (return_code == GSL_EROUND).
			Those may be problematic, but we may decide to ignore them and use
			the best integration result without crashing.

		Otherwise, throw an exception to let user know an unacceptable
		gsl error occured.*/

//#ifndef SUPPRESS_INTEGRATION_ERRORS 
		
		if(return_code != 0) //error occured
		 //if(return_code != GSL_EROUND)  
			//if( !( (return_code == GSL_EDIVERGE) && (result < 1.0) ))
			{
				std::string msg = "gsl fcn gsl_integration_qagiu() experienced "
				 "unexpected error in fcn Partial_Number_Dens_pt(). "
				 "gsl error code was ";
				msg += std::to_string(return_code);
				msg += ". Temp = " + std::to_string(Temp);
				msg += ". mu_particle = " + std::to_string(mu_particle);
				msg += ". particle = " + particle.GetName();
				throw GenericError(msg);
			}
			
//#endif


		if(result < 0.0) //don't allow negative pressures! That's a fail!
		{
			std::string msg = "Error: Partial_Number_Dens_pt() computed a "
			                  "negative value.";
			msg += " Temp = " + std::to_string(Temp);
			msg += " mu_particle = " + std::to_string(mu_particle);
			msg += " particle = " + particle.GetName();
			throw GenericError(msg);
		}

		//multiply integral by degeneracy, phase space factors, and units 
		//which convert integral into number density 
		result *= degen*pow(mass, 3.0)/(2.0*pi*pi);

		return result;
	}
}


//--------------------------------------------------

double Partial_Baryon_Dens_pt(double Temp, double mu_particle,
 const Particle& particle)
{
	//Return partial baryon density of this PT model gas.
	//That's just partial numbder density * baryon charge of the particle

	if(Temp <=0.0)
		throw TemperatureError("Partial_Baryon_Dens_pt expects positive"
		 " temperature.");
	
	double baryonCharge = particle.GetBaryonCharge();
	double baryonDens = baryonCharge*Partial_Number_Dens_pt(Temp, mu_particle,
 	 particle);

	return baryonDens;
}

//--------------------------------------------------


double Partial_Energy_Dens_Integrand(double x, void *p)
{
	//Returns dimensionless integrand as a fcn of x = E/mass
	//and other parameters (like Temp, mass, particle chemical potential,...)
	//When integrated over x, this integral will compute the 
	//dimensionless component of partial energy density.

	struct partial_Pressure_Params * params = 
	 (struct partial_Pressure_Params *) p;

	double Temp = params->Temp;
	double mass = params->mass;
	double mu_particle = params->mu_particle;
	QUANTUM_TYPE quantumtype = params->quantumtype;

	//x = Energy/mass
	double E = mass*x;

	if(E < mass)
		throw ExVolValueError("Partial_Energy_Dens_Integrand() requires"
		 "x = E/mass >= 1.");

	return x*x*sqrt(x*x - 1.0)*
	 equilib_distfcn(E, Temp, mu_particle, quantumtype);
}


double Partial_Energy_Dens_pt(double Temp, double mu_particle,
 const Particle& particle)
{
	//Compute partial energy density of a particle species, identified by 
	//particle object at temperature Temp and particle chemical potential 
	//mu_particle.

	if(Temp <=0.0)
		throw TemperatureError("Partial_Energy_Dens_pt expects" 
		 "positive temperature.");
  
	double mass = particle.GetMass();
	double degen = particle.GetDegen();
	QUANTUM_TYPE quantumtype = particle.GetQuantumType();

	if(quantumtype == CLASSICAL)
	{
		return 1.0/(2.0*pi*pi)*degen*pow(mass,4.0)*
					(3.0*gsl_sf_bessel_Kn(2,mass/Temp)/(pow(mass/Temp,2)) +
				gsl_sf_bessel_Kn(1,mass/Temp)/(mass/Temp))
				*gsl_sf_exp(mu_particle/Temp);
	}
	else
	{
		double lower_limit = 1.0;  //integrate x=E/mass from 1 to infinity

		//do some setup for gsl
		const int MaxNumIntervals = 1000;
		gsl_integration_workspace *work_ptr = 
		 gsl_integration_workspace_alloc(MaxNumIntervals);
		double abs_error = 1.0e-8;	
		double rel_error = 1.0e-8;	

		double result, error;		//where we will store results

		//create structure holding parameters to pass to integrand function:
		struct partial_Pressure_Params params;
		params.Temp = Temp;
		params.mass = mass;
		params.mu_particle = mu_particle;
		params.quantumtype = quantumtype;
		void *params_ptr = &params;

		//set up reference to the integrand function to be integrated over:
		gsl_function My_function;
		My_function.function = &Partial_Energy_Dens_Integrand;
		My_function.params = params_ptr;


		//Disable gsl's error handling so it won't crash the code if
		//(for instance) an integral converges slowly. We'll check status after.
		gsl_set_error_handler_off();

		//invoke gls's quadrature-adaptive routine, integrate from lower limit to
		// infinity
		int return_code;
		return_code=0;
		return_code = gsl_integration_qagiu (&My_function, lower_limit,
		 abs_error, rel_error, MaxNumIntervals, work_ptr, &result, &error);

		//Turn gsl's default error handler back on
		gsl_set_error_handler(NULL);

		gsl_integration_workspace_free(work_ptr);


		/*
		The commented-out lines below could be uncommented to
		allow some non-critical errors:

			For instance, the integral may converge slowly 
			(return_code == GSL_EDIVERGE), especially when the partial pressure
			is very tiny and negligible (result < 1.0).  Who cares if convergence
			is poor on a result that contributes negliglibly? 
			So we might ignore those errors.

			Also, there may be rounding errors (return_code == GSL_EROUND).
			Those may be problematic, but we may decide to ignore them and use
			the best integration result without crashing.

		Otherwise, throw an exception to let user know an unacceptable
		gsl error occured.*/

//#ifndef SUPPRESS_INTEGRATION_ERRORS
	  
		if(return_code != 0) //error occured
		 //if(return_code != GSL_EROUND)  
			//if( !( (return_code == GSL_EDIVERGE) && (result < 1.0) ))
			{
				std::string msg = "gsl fcn gsl_integration_qagiu() experienced "
				 "unexpected error in fcn Partial_Energy_Dens_pt(). "
				 "gsl error code was ";
				msg += std::to_string(return_code);
				msg += ". Temp = " + std::to_string(Temp);
				msg += ". mu_particle = " + std::to_string(mu_particle);
				msg += ". particle = " + particle.GetName();
				throw GenericError(msg);
			}
			
//#endif	


		if(result < 0.0)
		{
			std::string msg = "Error: Partial_Energy_Dens_pt() computed a "
			                  "negative value.";
			msg += " Temp = " + std::to_string(Temp);
			msg += " mu_particle = " + std::to_string(mu_particle);
			msg += " particle = " + particle.GetName();
			throw GenericError(msg);
		}	

		//multiply integral by degeneracy, phase space factors, and units 
		//which convert integral into energy density
		result *= degen*pow(mass, 4.0)/(2.0*pi*pi);

		return result;
	}
}

double Partial_Entropy_Dens_pt(double Temp, double mu_particle,
 const Particle& particle)
{
	if(Temp <=0.0)
		throw TemperatureError("Partial_Entropy_Dens_pt expects positive" 
		 "temperature.");

	double energydens = Partial_Energy_Dens_pt(Temp, mu_particle, particle);
	double pressure = Partial_Pressure_pt(Temp, mu_particle, particle);
	double numberdens = Partial_Number_Dens_pt(Temp, mu_particle, particle);
	//use thermodynamic relation (can be shown from kinetic theory integrals 
	//for pressure, energydens,numberdens):
	return (energydens + pressure - mu_particle*numberdens)/Temp;
}


//=============================================================================
//Code to compute thermal properties of a point particle gas of many particle 
//species

double Total_Pressure_pt(double Temp, double mu_B, 
 const std::vector<Particle>& ParticleList)
{
	if(Temp <=0.0)
		throw TemperatureError("Total_Pressure_pt expects positive temperature.");

	double mu_particle;
	double Pressure=0.0;
	//C++11 for loop means: for each particle in ParticleList:
	//auto picks type automatically; const means particle is unaltered
	for(auto const &particle: ParticleList) 
	{
		mu_particle = particle.GetParticleChemPotential(mu_B);
		Pressure += Partial_Pressure_pt(Temp, mu_particle, particle);
	}
	return Pressure;
}

double Total_Energy_Dens_pt(double Temp, double mu_B, 
 const std::vector<Particle>& ParticleList)
{
	if(Temp <=0.0)
		throw TemperatureError("Total_Energy_Dens_pt expects positive"
		 " temperature.");

	double mu_particle;
	double EnergyDens=0.0;
	//C++11 for loop means: for each particle in ParticleList:
	//auto picks type automatically; const means particle is unaltered
	for(auto const &particle: ParticleList)
 	{
		mu_particle = particle.GetParticleChemPotential(mu_B);
		EnergyDens += Partial_Energy_Dens_pt(Temp, mu_particle, particle);
	}
	return EnergyDens;
}

double Total_Number_Dens_pt(double Temp, double mu_B, 
 const std::vector<Particle>& ParticleList)
{
	if(Temp <=0.0)
		throw TemperatureError("Total_Number_Dens_pt expects positive"
		 " temperature.");

	double mu_particle;
	double NumbDens=0.0;
	//C++11 for loop means: for each particle in ParticleList:
	//auto picks type automatically; const means particle is unaltered
	for(auto const &particle: ParticleList) 
	{
		mu_particle = particle.GetParticleChemPotential(mu_B); 
		NumbDens += Partial_Number_Dens_pt(Temp, mu_particle, particle);
	}	
	return NumbDens;
}

double Total_Baryon_Dens_pt(double Temp, double mu_B, 
 const std::vector<Particle>& ParticleList)
{
	if(Temp <=0.0)
		throw TemperatureError("Total_Baryon_Dens_pt expects positive"
		 " temperature.");

	double mu_particle;
	double baryonDens=0.0;
	//C++11 for loop means: for each particle in ParticleList:
	//auto picks type automatically; const means particle is unaltered
	for(auto const &particle: ParticleList) 
	{
		mu_particle = particle.GetParticleChemPotential(mu_B);
		baryonDens += Partial_Baryon_Dens_pt(Temp, mu_particle, particle);
	}
	return baryonDens;
}

double Total_Entropy_Dens_pt(double Temp, double mu_B, 
 const std::vector<Particle>& ParticleList)
{
	if(Temp <=0.0)
		throw TemperatureError("Total_Entropy_Dens_pt expects positive" 
		 "temperature.");

	double energyDens = Total_Energy_Dens_pt(Temp, mu_B, ParticleList);
  double pressure = Total_Pressure_pt(Temp, mu_B, ParticleList);
	double baryonDens = Total_Baryon_Dens_pt(Temp, mu_B, ParticleList);
	//use standard thermodynamic relation 
	return (energyDens + pressure - mu_B*baryonDens)/Temp;
}

