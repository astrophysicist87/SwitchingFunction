#include <cmath>
#include <vector>
#include <stdio.h>
#include "constants.hpp"
#include "exceptions.hpp"
#include "particle.hpp"
#include "thermo_pt.hpp"
#include "thermo_exI.hpp"
#include "thermo_exII.hpp"

using namespace std;

//-----------------------------------------------------------------------------
//Code to find exII pressure by solving a non-linear equation.


//Find the exII excluded volume of a particle.
//If epsilon0 is specified, use that to compute
//vol = mass/epsilon0.  Else use particle's volume variable.
double PickVolume(const Particle& particle, double epsilon0 = -1.0)
{
	double volume;            //particle exII excluded volume
	if(epsilon0 > 0.0)        //use the specified epsilon0 to compute volume
		volume = particle.GetMass() / epsilon0;
	else
		volume = particle.GetExcludedVolume_exII();//use particle's preset volume
	return volume;
}


//if epsilon0 == -1.0, use particle's excluded volume,
//else use volume = mass/epsilon0
double GetEffectiveParticleChemPotential(double mu_B, double Press_exII,
const Particle& particle, double epsilon0 = -1.0)
{
	double mu_particle = particle.GetParticleChemPotential(mu_B);

	double volume = PickVolume(particle, epsilon0);

	//compute exII effective particle chemical potential
	double mu_particle_effective = mu_particle - volume*Press_exII;
	return mu_particle_effective;
}




//-----------------------------------------------
//function value:

//When we find a correct value Press_exII, this function will equal zero.
//epsilon0 is optional.  If epsilon0 is not specified, excluded volume
//is retrieved from particle object.  If epsilon0 is specified,
//then compute a particle's excluded volume: volume = mass/epsilon0.
double PexII_fcn_to_zero(double Press_exII, double Temp, double mu_B,
 const std::vector<Particle>& ParticleList, double epsilon0 = -1.0)
{
	double mu_particle_effective;
	double Sum = 0;
	for(auto const &particle: ParticleList)
	{
		mu_particle_effective = GetEffectiveParticleChemPotential(mu_B, Press_exII,
		 particle, epsilon0);

		Sum += Partial_Pressure_pt(Temp, mu_particle_effective, particle);
	}

	return (Press_exII - Sum);  //this = 0 when find solution
}


//default value for Press_exII_guess is -1.0 (means ignore)
//default value for epsilon0 = -1.0 (means ignore)
double Total_Pressure_exII(double Temp, double mu_B,
 const std::vector<Particle>& ParticleList, double epsilon0,
 double Press_exII_guess)
{
	const double Tol = 0.1;  //solver tolerance

	if(Temp <=0.0)
		throw TemperatureError("Total_Pressure_exII expects positive"
		 " temperature.");


	if(Press_exII_guess < 0.0) //no exII pressure guess provided, so generate one
	{
		if(epsilon0 > 0.0)
		{
			//user specified an epsilon0, so can use exI to generate a guess
			Press_exII_guess = Total_Pressure_exI(Temp, mu_B, epsilon0,
 			 ParticleList);
		}
		else
		{
			//no epsilon0 provided, so can't use exI: use pt to make a guess
			Press_exII_guess = Total_Pressure_pt(Temp, mu_B, ParticleList);
		}
	}

	//intialize guesses for rootfinding:
	//P is for pressure, f is for fcn to zero:  f(P) == 0 when find root P
	double P0 = Press_exII_guess; //1st guess
	double dP = 100.0;
	double P1 = P0 + dP;          //2nd guess: pick a nearby starting point

	double f0 = PexII_fcn_to_zero(P0, Temp, mu_B, ParticleList, epsilon0);
	double f1 = PexII_fcn_to_zero(P1, Temp, mu_B, ParticleList, epsilon0);

	double P2, f2, divisor, error;

	//loop and perform search, updating guesses, until find sln.
	//Uses secant method.
	for(int i = 0; i < 10*MAX_NUMBER_ITER; i++)
	{
		divisor = (f1 - f0);
		if(divisor != 0.0)
			P2 = P1 - f1*(P1 - P0)/divisor;  //secant method
		else
			P2 = Press_exII_guess*(1.0 + 0.2*i);  //fail, so shift starting point

		f2 = PexII_fcn_to_zero(P2, Temp, mu_B, ParticleList, epsilon0);
		error = std::abs(f2);
		if(error <= Tol){
			break;   //found solution within tolerance
            }
		//update variables for next iteration
		P0 = P1;
		f0 = f1;
		P1 = P2;
		f1 = f2;
	}
	if(error > Tol){printf("%8.8e",f2);
		printf("error: %f %f %f %f\n",P0,error,Temp,mu_B);
		throw ConvergenceError("Function Total_Pressure_exII failed to find"
			" a solution.");
                }
	return P2;  //this is the exII pressure which is a solution
}

//-----------------------------------------------------------------------------

double Total_Energy_Dens_exII(double Temp, double mu_B,
 const std::vector<Particle>& ParticleList, double epsilon0,
 double Press_exII_guess)
{
	if(Temp <=0.0)
		throw TemperatureError("Total_Energy_Dens_exII expects positive"
		 " temperature.");

	double mu_particle_effective, volume;
	double numerator = 0.0, denominator = 1.0;

	//need exII pressure to get effective chem potential.Use guess to speed up.
	double Press_exII = Total_Pressure_exII(Temp, mu_B, ParticleList, epsilon0,
	 Press_exII_guess);

	for(auto const &particle: ParticleList)
	{
		volume = PickVolume(particle, epsilon0);
		mu_particle_effective = GetEffectiveParticleChemPotential(mu_B,
	 	 Press_exII, particle, epsilon0);

		numerator += Partial_Energy_Dens_pt(Temp, mu_particle_effective,
		 particle);

		denominator += volume * Partial_Number_Dens_pt(Temp,
		 mu_particle_effective, particle);
	}

	return numerator/denominator;
}

//----------------------------------------------------------------------------

double Total_Baryon_Dens_exII(double Temp, double mu_B,
 const std::vector<Particle>& ParticleList, double epsilon0,
 double Press_exII_guess)
{
	if(Temp <=0.0)
		throw TemperatureError("Total_Baryon_Dens_exII expects positive"
		 " temperature.");

	double mu_particle_effective, volume;
	double numerator = 0.0, denominator = 1.0;

	//need exII pressure to get effective chem potential.Use guess to speed up.
	double Press_exII = Total_Pressure_exII(Temp, mu_B, ParticleList, epsilon0,
	 Press_exII_guess);

	for(auto const &particle: ParticleList)
	{
		volume = PickVolume(particle, epsilon0);
		mu_particle_effective = GetEffectiveParticleChemPotential(mu_B,
	 	 Press_exII, particle, epsilon0);

		numerator += Partial_Baryon_Dens_pt(Temp, mu_particle_effective,
		 particle);

		denominator += volume * Partial_Number_Dens_pt(Temp,
		 mu_particle_effective, particle);
	}

	return numerator/denominator;
}

 //----------------------------------------------------------------------------

double Total_Number_Dens_exII(double Temp, double mu_B,
 const std::vector<Particle>& ParticleList, double epsilon0,
 double Press_exII_guess)
{
	if(Temp <=0.0)
		throw TemperatureError("Total_Number_Dens_exII expects positive"
		 " temperature.");

	double mu_particle_effective, volume;
	double numerator = 0.0, denominator = 1.0;

	//need exII pressure to get effective chem potential.Use guess to speed up.
	double Press_exII = Total_Pressure_exII(Temp, mu_B, ParticleList, epsilon0,
	 Press_exII_guess);

	for(auto const &particle: ParticleList)
	{
		volume = PickVolume(particle, epsilon0);
		mu_particle_effective = GetEffectiveParticleChemPotential(mu_B,
	 	 Press_exII, particle, epsilon0);

		numerator += Partial_Number_Dens_pt(Temp, mu_particle_effective,
		 particle);

		denominator += volume * Partial_Number_Dens_pt(Temp,
		 mu_particle_effective, particle);
		}

	return numerator/denominator;
}


//----------------------------------------------------------------------------

double Total_Entropy_Dens_exII(double Temp, double mu_B,
 const std::vector<Particle>& ParticleList, double epsilon0,
 double Press_exII_guess)
{
	if(Temp <=0.0)
		throw TemperatureError("Total_Entropy_Dens_exII expects positive"
		 " temperature.");

	//find pressure, using guess to speed it up
	double Pressure = Total_Pressure_exII(Temp, mu_B, ParticleList, epsilon0,
 	 Press_exII_guess);

	//find other quantities, using pressure as a pressure guess to speed
	//the solution
	double EnergyDens = Total_Energy_Dens_exII(Temp, mu_B, ParticleList,
	 epsilon0, Pressure);

	double BaryonDens = Total_Baryon_Dens_exII(Temp, mu_B, ParticleList,
	 epsilon0, Pressure);

	double EntropyDens = (Pressure + EnergyDens - mu_B*BaryonDens)/Temp;
	return EntropyDens;
}
