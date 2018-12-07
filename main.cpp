#include <vector>
#include <algorithm>
#include <array>
#include <iostream>
#include <cmath>
#include <omp.h>

#include "constants.hpp"
#include "crossover_line.hpp"
#include "particle.hpp"
#include "particlelist.hpp"
#include "thermo_pt.hpp"
#include "thermo_exI.hpp"
#include "thermo_exII.hpp"
#include "thermo_pqcd.hpp"
#include "thermo_crossover.hpp"
#include "crossover_line.hpp"

#include "util.hpp"

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_gegenbauer.h>

#include "solver.h"

using namespace std;

void get_function_evals(vector<double> & pts, vector<double> & fevals,
						double Radius, double correction_factor,
						double R, double dRdpsi, thermparams parameters);
vector<double> get_derivatives(vector<double> & f, double stepsize, int direction);
//vector<double> get_derivatives(double stepsize, int direction);


int main(int argc, char *argv[])
{

	vector<Particle> ParticleList = GetPDGParticleList();

	/*
	List of parameters used by Switching function and models, given in
	this order:
	psic : critical angle
	a : switching function parameter
	b : switching function parameter
	r : switching function parameter
	Escl : energy scale constant
	Soft : softening constant
	eps0_4 : epsilon0^(1/4)
	Ks={k2,k4} : coefficients used to determine R(psi)
	model : model type EXI, EXII, or PT
	ParticleList : list of hadrons
	*/

	//Optimal parameters for Tc = 130 MeV, a=1, b=.8
	//for model = EXI
	thermparams parameters={1.448174,1.000000,0.800000,4,10.728936,1.000000,
		499.314608,{0.000015,0.386452},EXI,ParticleList};
	//for model = EXII
	//thermparams parameters={1.449514,1.000000,0.800000,4,10.647927,1.000000,
	//	472.080882,{0.000093,0.387725},EXII,ParticleList};

    std::vector<double> Rs;
    double R, dRdpsi, Rc, psi, psic;

	// set shape of contour
	//const double Radius = atof(argv[1]);  //MeV
	const double Radius = 200.0;  //MeV
	const double correctionFactor = atof(argv[1]);  //MeV

	// compute angle corresponding to T_min
	const double T_min = 20.0;	// MeV
	const double mu_at_T_min = correctionFactor*sqrt(Radius*Radius - T_min*T_min);	// MeV

	double left_psi_endpoint = 0.0;
	double right_psi_endpoint = atan2(mu_at_T_min / correctionFactor, T_min);
	double hw_psi = 0.5 * ( right_psi_endpoint - left_psi_endpoint );
	double cen_psi = 0.5 * ( right_psi_endpoint + left_psi_endpoint );

	// calculate derivatives at endpoints
	const int nOrder = 10, nL = 3, nR = 3;
	vector<double> fevalsL(nOrder), fevalsR(nOrder);

	// set psi points near endpoints (for evaluating derivatives)
	const double stepsize = 0.00025;
	vector<double> left_endpoints = linspace(left_psi_endpoint, left_psi_endpoint+double(nOrder-1)*stepsize, nOrder);
	vector<double> right_endpoints = linspace(right_psi_endpoint - double(nOrder-1)*stepsize, right_psi_endpoint, nOrder);

	// evaluate function near endpoints
	get_function_evals(left_endpoints, fevalsL, Radius, correctionFactor, R, dRdpsi, parameters);
	get_function_evals(right_endpoints, fevalsR, Radius, correctionFactor, R, dRdpsi, parameters);

	///*
	cout << "Left end:" << endl;
	for (int i = 0; i < left_endpoints.size(); ++i)
		cout << i << "   " << left_endpoints[i] << "   " << fevalsL[i] << endl;

	cout << "Right end:" << endl;
	for (int i = 0; i < right_endpoints.size(); ++i)
		cout << i << "   " << right_endpoints[i] << "   " << fevalsR[i] << endl;
	//*/


	// now get derivatives
	vector<double> left_derivatives = get_derivatives(fevalsL, stepsize, 1);
	vector<double> right_derivatives = get_derivatives(fevalsR, stepsize, -1);

	///*
	cout << "Left end:" << endl;
	for (int i = 0; i < left_derivatives.size(); ++i)
		cout << i << "   " << left_derivatives[i] << endl;

	cout << "Right end:" << endl;
	for (int i = 0; i < right_derivatives.size(); ++i)
		cout << i << "   " << right_derivatives[i] << endl;
	//*/

	if (true) exit (8);

	// max derivative of Chebyshev to compute
	const int kmax = max(nL, nR);

	// max order of Chebyshev to compute
	const int nD = nL + nR + 2;
	const int Nmax = nD - 1;

	double Chebyshev_array_L[ (Nmax+1) * (kmax+1) ], Chebyshev_array_R[ (Nmax+1) * (kmax+1) ];
	double Gamma [kmax+1];

	for (int i = 0; i < (Nmax+1) * (kmax+1); ++i)
	{
		Chebyshev_array_L[i] = 0.0;
		Chebyshev_array_R[i] = 0.0;
	}

	// ik: order of derivative
	for (int ik = 0; ik <= kmax; ++ik)
	{
		Gamma[ik] = (ik==0)? 0: gsl_sf_gamma(ik);
		double polyL[Nmax+1], polyR[Nmax+1];
		if (ik == 0)
		{
			for (int iN = 0; iN <= Nmax; ++iN)
			{
				Chebyshev_array_L[iN * (kmax+1) + ik] = pow(-1.0,iN);
				Chebyshev_array_R[iN * (kmax+1) + ik] = 1.0;
			}
		}
		else
		{
			//==============================
			// Do left and right endpoints
			int status = gsl_sf_gegenpoly_array (Nmax, double(ik), -1.0, polyL);
			status += gsl_sf_gegenpoly_array (Nmax, double(ik), 1.0, polyR);

			if (status > 0) cerr << "Warning: status > 0 in GSL Gegenbauer function!" << endl;

			//==============================
			// Fill arrays
			for (int iN = ik; iN <= Nmax; ++iN)	// automatically zero for k>n
			{
				double resultL = pow(2.0, ik-1.0)*iN*Gamma[ik]*polyL[iN-ik];
				double resultR = pow(2.0, ik-1.0)*iN*Gamma[ik]*polyR[iN-ik];
				Chebyshev_array_L[iN * (kmax+1) + ik] = (abs(resultL) < 1.e-100)? 0.0: resultL;
				Chebyshev_array_R[iN * (kmax+1) + ik] = (abs(resultR) < 1.e-100)? 0.0: resultR;
			}
		}

	}

	/*
	// check Chebyshev arrays
	cout << "Left check:" << endl;
	for (int iN = 0; iN <= Nmax; ++iN)
	for (int ik = 0; ik <= kmax; ++ik)
		cout << iN << "   " << ik << "   " << Chebyshev_array_L[iN * (kmax+1) + ik] << endl;

	cout << "Right check:" << endl;
	for (int iN = 0; iN <= Nmax; ++iN)
	for (int ik = 0; ik <= kmax; ++ik)
		cout << iN << "   " << ik << "   " << Chebyshev_array_R[iN * (kmax+1) + ik] << endl;
	*/


	
	double coeff_matrix [nD*nD], constraint_vector[nD];

	// fill constrain_vector
	{
		int idx = 0;
		for (int i = 0; i <= nL; ++i)
			constraint_vector[idx++] = pow(hw_psi, double(i))*left_derivatives[i];
		for (int i = 0; i <= nR; ++i)
			constraint_vector[idx++] = pow(hw_psi, double(i))*right_derivatives[i];

		// NOTE: prefactors rescale derivatives to appropriate
		// interval (-1,+1) for using Chebyshev polynomials!!!
	}

	// fill coeff_matrix
	{
		int idx = 0;
		for (int i = 0; i <= nL; ++i)
		for (int n = 0; n < nD; ++n)
			coeff_matrix[idx++] = Chebyshev_array_L[n * (kmax+1) + i];
		for (int i = 0; i <= nR; ++i)
		for (int n = 0; n < nD; ++n)
			coeff_matrix[idx++] = Chebyshev_array_R[n * (kmax+1) + i];
	}

	/*
	{
		int idx = 0;

		// check constraint_vector and coeff_matrix
		cout << "constraint_vector =" << endl;
		for (int i = 0; i < nD; ++i)
			cout << i << "   " << constraint_vector[i] << endl;

		cout << "coeff_matrix =" << endl;
		for (int i = 0; i < nD; ++i)
		{
			for (int j = 0; j < nD; ++j)
				cout << coeff_matrix[idx++] << "   ";
			cout << endl;
		}
	}
	*/


	double the_solution[nD];

	solve(coeff_matrix, constraint_vector, the_solution, nD);

	/*
	cout << "the_solution =" << endl;
	for (int i = 0; i < nD; ++i)
		cout << i << "   " << the_solution[i] << endl;
	*/

	const int nPsiPts = 1000;
	vector<double> PsiPts = linspace(left_psi_endpoint, right_psi_endpoint, nPsiPts);

	// final check
	//cout << "Final check:" << endl;
	/*
	for (int iPsi = 0; iPsi < nPsiPts; ++iPsi)
	{
		double Psi_loc = PsiPts[iPsi];
		double Psi_trans = ( left_psi_endpoint + right_psi_endpoint - 2.0 * Psi_loc ) / ( left_psi_endpoint - right_psi_endpoint );

		double exact_result = (1.0+Psi_loc)*exp(-Psi_loc*Psi_loc);
		double approx_result = 0.0;
		for (int i = 0; i < nD; ++i)
			approx_result += the_solution[i] * cos(i*acos(Psi_trans));	//cheat way to evaluate just the Chebyshev functions

		cout << Psi_loc << "   " << approx_result << "   " << exact_result << endl;
	}
	*/
	for (int iPsi = 0; iPsi < nPsiPts; ++iPsi)
	{
		double Psi_loc = PsiPts[iPsi];
		double Psi_trans = ( left_psi_endpoint + right_psi_endpoint - 2.0 * Psi_loc ) / ( left_psi_endpoint - right_psi_endpoint );

		double approx_result = 0.0;
		for (int i = 0; i < nD; ++i)
			approx_result += the_solution[i] * cos(i*acos(Psi_trans));	//cheat way to evaluate just the Chebyshev functions

		cout << atan(correctionFactor*tan(Psi_loc)) << "   " << approx_result << endl;
	}



	/*
	//const int nPsiPts = 1000;
	//vector<double> PsiPts = linspace(0.0, atan2(mu_at_T_min, T_min), nPsiPts);


	// calculate pressure along contour
	for (int iPsi = 0; iPsi < nPsiPts; ++iPsi)
	{
		double Psi_loc = PsiPts[iPsi];

		double Temp = Radius * cos(Psi_loc), mu_B = correctionFactor * Radius * sin(Psi_loc);

		// To calculate S, and the thermodynamic quantities
		psi=atan2(mu_B,Temp);
		psic=parameters.psic;

		//Determine R at critical angle
		Rc=Get_R(psic,parameters);

		//Determines R and dRdpsi at the angle of the point being calculated
		if(psi>psic)
		{
			//if below the critical line, calculates where P_H=PQCD
			R=Get_R(psi,parameters);
			dRdpsi=0.0; 
		}
		else
		{
			//if above critical line, calculates where critical line meets the 
			//quartic function T = T0 ( 1- k2*(mu/muc)^2 - k4*(mu/muc)^4 - ...)
			Rs=Get_Rs_quartic(psi,Rc,parameters);
			R=Rs[0]; dRdpsi=Rs[1];
		}

		// Returns the vector therms={Pressure,Entropy,BaryonDensity,
		// EnergyDensity,SwitchingFuncton};
	   
		std::vector<double> therms;
		therms= Get_Thermo(Temp, mu_B, R, dRdpsi, parameters);

		//printf( "%e %e %e %e %e %e %e %f %f\n",
		//		Temp, mu_B, Psi_loc,
		//		therms[0], therms[1], therms[2],
		//		therms[3], therms[4],
		//		therms[0]/(Temp*Temp*Temp*Temp) );
		std::cout << Temp << "   " << mu_B << "   " << R << "   " << psi << "   " << Psi_loc << "   "
					<< therms[0] << "   " << therms[1] << "   "
					<< therms[2] << "   " << therms[3] << "   "
					<< therms[4] << "   " << therms[0]/(Temp*Temp*Temp*Temp)
					<< std::endl;
		
	}
	*/


	return 0;	
}

//========================================================================
void get_function_evals(vector<double> & pts, vector<double> & fevals,
						double Radius, double correctionFactor,
						double R, double dRdpsi, thermparams parameters)
{
	double Tc = 130.0;
	for (int i = 0; i < (int)fevals.size(); ++i)
	{

		double R, dRdpsi;
		vector<double> therms, Rs;
		double local_psi = pts[i];

		double Temp = Radius * cos(local_psi), mu_B = correctionFactor * Radius * sin(local_psi);


		// To calculate S, and the thermodynamic quantities
		double psi=atan2(mu_B,Temp);
		double psic=parameters.psic;


		//Determine R at critical angle
		double Rc=Get_R(psic,parameters);


		//Determines R and dRdpsi at the angle of the point being calculated
		if(psi>psic)
		{

			//if below the critical line, calculates where P_H=PQCD
			R=Get_R(psi,parameters);
			dRdpsi=0.0; 

		}
		else
		{

			//if above critical line, calculates where critical line meets the 
			//quartic function T = T0 ( 1- k2*(mu/muc)^2 - k4*(mu/muc)^4 - ...)
			Rs=Get_Rs_quartic(psi,Rc,parameters);
			R=Rs[0]; dRdpsi=Rs[1];

		}


		therms = Get_Thermo( Temp, mu_B, R, dRdpsi, parameters );


		fevals[i] = therms[0] / (Tc*Tc*Tc*Tc);

		//fevals[i] = (1.0+local_psi)*exp(-local_psi*local_psi);

	}

	return;
}




