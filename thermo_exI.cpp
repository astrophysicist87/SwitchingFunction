#include <cmath>
#include <string>
#include <stdexcept>
#include <vector>

#include "constants.hpp"
#include "exceptions.hpp"
#include "particle.hpp"
#include "particlelist.hpp"
#include "thermo_pt.hpp"

using namespace std;

//using namespace std;

//-----------------------------------------------------------------------------
//Functions for computing Temp, mu_B in model exI given Temp* = T_star and
//mu_B* = mu_B_star.

double Temp_exI(double T_star, double mu_B_star, double epsilon0, 
 const std::vector<Particle>& ParticleList)
{
	if(T_star <=0.0)
		throw TemperatureError("Temp_exI expects positive T_star.");

	if(epsilon0 <= 0.0)
		throw epsilon0Error("Temp_exI expects positive"
		 " epsilon0.");

	double Temp = T_star/(1.0 - 
	 Total_Pressure_pt(T_star, mu_B_star, ParticleList)/epsilon0);
	return Temp;
}

double mu_B_exI(double T_star, double mu_B_star, double epsilon0,
 const std::vector<Particle>& ParticleList)
{
	if(T_star <=0.0)
		throw TemperatureError("mu_B_exI expects positive T_star.");

	if(epsilon0 <= 0.0)
		throw epsilon0Error("mu_B_exI expects positive"
		 " epsilon0.");

	double mu_B = mu_B_star/(1.0 - 
	 Total_Pressure_pt(T_star, mu_B_star, ParticleList)/epsilon0);
	return mu_B;
}


//-----------------------------------------------------------------------------
//Compute Temp* (Tstar) and mu_B* (mu_Bstar) in model exI given Temp, mu_B.
//returns values in reference variables T_star, mu_B_star.
//Based on 3rd generation solver developed in python code.

void Find_Temp_star_mu_B_star_exI(double Temp, double mu_B, double epsilon0,
 const std::vector<Particle>& ParticleList, double& Tstar, double& mu_Bstar)
{
	if(Temp <= 0.0)
		throw TemperatureError("Find_Temp_star_mu_B_star_exI expects positive" 
			" temperature.");

	if(epsilon0 <= 0.0)
		throw epsilon0Error("Find_Temp_star_mu_B_star_exI expects positive"
		 " epsilon0.");


	double Tol_Temp = 0.00000001;       //Tolerance to which we try to know Temp.

	double TstarMin = 0.0;           //Bounds on Tstar search window.
	double TstarMax = 400.0;         //TstarMax found empirically, seems to work.
	double TempLinearCutoff = 200.0; //Below here, Temp~Tstar. Found empirically.


  double Tol_Tstar = 1e-10;  //If (TstarMax-TstarMin) < Tol_star, we say that
                             //TstarMax == TstarMin == Tstar2, and we're done


	//-----------------------------------------------------------------
	//Simple initialization.  Could improve in future.


	double Tstar0, Tstar1;     //Temp* guesses
	double f0, f1;             //fcn values at Tstar0, Tstar1 

	//generate 1st guess = Tstar0:

	if(Temp < TempLinearCutoff)  //for small Temps, Tstar ~ Temp, but Tstar <Temp
	{
		Tstar0 = Temp - 1.0;   //so guess Tstar <~ Temp
		if(Tstar0 < TstarMin)   //revise guess to avoid values below TstarMin
			Tstar0 = (TstarMin + TstarMax)/2.0;
	}
	else    //Simply start in the middle of the Tstar search window
		Tstar0 = (TstarMin + TstarMax)/2.0;

	//evaluate fcn at Tstar0 guess.  fcn==0 when you find sln:
	f0 = Temp_exI(Tstar0, (mu_B/Temp)*Tstar0, epsilon0, ParticleList) - Temp;
	if(f0 > 0.0 || f0 < -Temp) //Tstar0 is too big
	{
		TstarMax = Tstar0;
		Tstar1 = Tstar0 - 1.0;   //get a second nearby guess
	}
	else 
	{                          //Tstar0 is too small
		TstarMin = Tstar0;
		Tstar1 = Tstar0 + 1.0;
	}

	if(Tstar1 < TstarMin || Tstar1 > TstarMax)  //don't go out of bounds
		Tstar1 = (TstarMin + TstarMax)/2.0;

	//eval fcn at Tstar1.  fcn==0 when you find sln
	f1 = Temp_exI(Tstar1, (mu_B/Temp)*Tstar1, epsilon0, ParticleList) - Temp;


	//-----------------------------------------------------------------
	//Now, iterate, using previous 2 guesses in each step to pick next guess.
	//If a smart secant-method based guess fails, fall back to bisection method.


	bool UseSecantMethod = true; //use secant method to accelerate solving


	double Tstar2, divisor, f2, error;
	for(int iteration = 0; iteration < MAX_NUMBER_ITER; iteration++)
	{

		if(iteration == MAX_NUMBER_ITER / 2)
			UseSecantMethod = false;  //Secant method is not finding the solution
			                          //quickly, turn it off and use bisection method
		                            //which is guaranteed to converge.

		//make new guess Tstar2 
		if(UseSecantMethod == true)
		{
			divisor = f1 - f0;
 			Tstar2 = (Tstar0*f1 - Tstar1*f0)/divisor; //use secant method

			if(!(Tstar2 >= TstarMin && Tstar2 <= TstarMax)) //if new guess invalid,
				Tstar2 = (TstarMin + TstarMax)/2.0;   //fall back to bisection method
		}
		else
			Tstar2 = (TstarMin + TstarMax)/2.0; //use bisection method

		//find new f value at Tstar2:
		f2 = Temp_exI(Tstar2, (mu_B/Temp)*Tstar2, epsilon0, ParticleList) - Temp;
		error = fabs(f2);


		//Note: if fabs(TstarMax-TstarMin) < Tol_Tstar, you cannot
	  //improve any further because to available accuracy
		//TstarMax == TstarMin, so we should quit.

		if( (error < Tol_Temp) || ( fabs(TstarMax-TstarMin) < Tol_Tstar))   
		{
			Tstar = Tstar2;                  //found solution
			mu_Bstar = (mu_B/Temp)*Tstar;
			break;
		}	

		//update variables for next iteration:
		if(f2 > 0.0 || f2 < -Temp) //Tstar2 is too big
			TstarMax = Tstar2;
		else                       //Tstar2 is too small
			TstarMin = Tstar2;
		Tstar0 = Tstar1;
		Tstar1 = Tstar2;
		f0 = f1;
		f1 = f2;
	}

	if( (error >= Tol_Temp) && (fabs(TstarMax-TstarMin) >= Tol_Tstar) )
		throw ConvergenceError("Function Find_Temp_star_mu_B_star_exI() failed to "
			" find a solution.");
}



//-----------------------------------------------------------------------------
//Code to compute thermal properties of a exI gas of many particle 
//species


double Total_Pressure_exI(double Temp, double mu_B, double epsilon0,
 const std::vector<Particle>& ParticleList)
{
	if(Temp <=0.0)
		throw TemperatureError("Total_Pressure_exI expects positive temperature.");

	if(epsilon0 <= 0.0)
		throw epsilon0Error("Total_Pressure_exI expects positive epsilon0.");

	//get T_star, mu_B_star
	double T_star, mu_B_star;
	Find_Temp_star_mu_B_star_exI(Temp, mu_B, epsilon0, ParticleList, T_star,
	 mu_B_star);

	double Press_pt = Total_Pressure_pt(T_star, mu_B_star, ParticleList);

	double Press_exI = Press_pt/(1.0 - Press_pt/epsilon0);

	return Press_exI;
}

double Total_Energy_Dens_exI(double Temp, double mu_B, double epsilon0,
 const std::vector<Particle>& ParticleList)
{
	if(Temp <=0.0)
		throw TemperatureError("Total_Energy_Dens_exI expects positive"
		 " temperature.");

	if(epsilon0 <= 0.0)
		throw epsilon0Error("Total_Energy_Dens_exI expects positive epsilon0.");

	//get T_star, mu_B_star
	double T_star, mu_B_star;
	Find_Temp_star_mu_B_star_exI(Temp, mu_B, epsilon0, ParticleList, T_star,
	 mu_B_star);

	double EnergyDens_pt = Total_Energy_Dens_pt(T_star, mu_B_star, ParticleList);
	double EnergyDens_exI = EnergyDens_pt/(1.0 + EnergyDens_pt/epsilon0);
	return EnergyDens_exI;
}

double Total_Number_Dens_exI(double Temp, double mu_B, double epsilon0,
 const std::vector<Particle>& ParticleList)
{
	if(Temp <=0.0)
		throw TemperatureError("Total_Number_Dens_exI expects positive"
		 " temperature.");

	if(epsilon0 <= 0.0)
		throw epsilon0Error("Total_Number_Dens_exI expects positive epsilon0.");

	//get T_star, mu_B_star
	double T_star, mu_B_star;
	Find_Temp_star_mu_B_star_exI(Temp, mu_B, epsilon0, ParticleList, T_star,
	 mu_B_star);

	double EnergyDens_pt = Total_Energy_Dens_pt(T_star, mu_B_star, ParticleList);
	double NumbDens_pt = Total_Number_Dens_pt(T_star, mu_B_star, ParticleList);
	double NumberDens_exI = NumbDens_pt/(1.0 + EnergyDens_pt/epsilon0);
	return NumberDens_exI;
}

double Total_Baryon_Dens_exI(double Temp, double mu_B, double epsilon0, 
 const std::vector<Particle>& ParticleList)
{
	if(Temp <=0.0)
		throw TemperatureError("Total_Baryon_Dens_exI expects positive"
		 " temperature.");

	if(epsilon0 <= 0.0)
		throw epsilon0Error("Total_Baryon_Dens_exI expects positive epsilon0.");	

	//get T_star, mu_B_star
	double T_star, mu_B_star;
	Find_Temp_star_mu_B_star_exI(Temp, mu_B, epsilon0, ParticleList, T_star,
	 mu_B_star);

	double EnergyDens_pt = Total_Energy_Dens_pt(T_star, mu_B_star, ParticleList);
	double BaryonDens_pt = Total_Baryon_Dens_pt(T_star, mu_B_star, ParticleList);
	double BaryonDens_exI = BaryonDens_pt/(1.0 + EnergyDens_pt/epsilon0);
	return BaryonDens_exI;
}

double Total_Entropy_Dens_exI(double Temp, double mu_B, double epsilon0, 
 const std::vector<Particle>& ParticleList)
{
	if(Temp <=0.0)
		throw TemperatureError("Total_Entropy_Dens_exI expects positive" 
		 "temperature.");

	if(epsilon0 <= 0.0)
		throw epsilon0Error("Total_Entropy_Dens_exI expects positive epsilon0.");

	double EnergyDens_exI = Total_Energy_Dens_exI(Temp, mu_B, epsilon0,
	 ParticleList);
  double Press_exI = Total_Pressure_exI(Temp, mu_B, epsilon0, ParticleList);
	double BaryonDens_exI = Total_Baryon_Dens_exI(Temp, mu_B, epsilon0,
	 ParticleList);
	//use standard thermodynamic relation 
	return (EnergyDens_exI + Press_exI - mu_B*BaryonDens_exI)/Temp;
}

