#include <vector>
#include <algorithm>
#include <array>
#include <iostream>
#include <cmath>

#include <gsl/gsl_sf_gegenbauer.h>
#include "chebyshev_library.h"

using namespace std;

namespace cheb_lib
{
	void set_nodes(vector<double> & nodes)
	{
		int number_of_points = nodes.size();
		for (int i = 1; i <= number_of_points; ++i)
		    nodes[i-1] = -cos(M_PI * (2.*i-1.) / (2.*number_of_points));

		return;
	}

	void get_nodes(double a, double b, vector<double> & nodes, vector<double> & adjnodes)
	{
		double halfwidth = 0.5*(b-a);
		int mode = 0;

		int number_of_points = nodes.size();

		set_nodes(nodes);

		switch (mode)
		{
			case 0:
				for (int i = 1; i <= number_of_points; ++i)
				   adjnodes[i-1] = (nodes[i-1] + 1.0)*halfwidth + a;
				break;
			case 1:
				break;
			case 2:
				break;
			default:
				cerr << "chebyshev(): mode = " << mode << " not supported!" << endl;
				break;
		}

		return;
	}

	void get_coefficients(vector<double> & coeffs_array, vector<double> & x_pts, vector<double> & fevals)
	{
		int npts = fevals.size();

		double dens[npts];

		for (int j = 0; j < npts; ++j)
		{
			dens[j] 			= 0.0;
			coeffs_array[j] 	= 0.0;
			for (int k = 0; k < npts; ++k)
			{
				double Tjk 		 = csf::Tfun(j, x_pts[k]);
				dens[j] 	    += Tjk*Tjk;
				coeffs_array[j] += fevals[k] * Tjk;
			}
			coeffs_array[j]    /= dens[j];
		}

	}

}

vector<double> get_derivatives(vector<double> & f, double stepsize, int direction)
{
	double adjustment = 1.0;	//heuristic
	double h = stepsize;

	double mult_factor = 1.0;

	if (direction == -1)	// backward
	{
		reverse( f.begin(), f.end() );
		mult_factor = -1.0;
	}

	vector<double> results;

	//=======================================
	// derivatives
	// function value
	results.push_back( f[0] );

	//===================
	// first derivative
	//results.push_back( mult_factor * ( -7129*f[0]+22680*f[1]-45360*f[2]+70560*f[3]-79380*f[4]
	//						+63504*f[5]-35280*f[6]+12960*f[7]-2835*f[8]+280*f[9] ) / ( 2520.0*h ) );
	results.push_back( (-147*f[0]+360*f[1]-450*f[2]+400*f[3]-225*f[4]+72*f[5]-10*f[6])/(60.0*h) );
	h *= adjustment;

	//===================
	// second derivative
	//results.push_back( ( 32575*f[0]-165924*f[1]+422568*f[2]-704368*f[3]+818874*f[4]
	//					-667800*f[5]+375704*f[6]-139248*f[7]+30663*f[8]-3044*f[9] ) / ( 5040.0*h*h ) );
	results.push_back( (812*f[0]-3132*f[1]+5265*f[2]-5080*f[3]+2970*f[4]-972*f[5]+137*f[6])/(180.0*h*h) );
	h *= adjustment;

	//===================
	// third derivative
	//results.push_back( mult_factor * ( -180920*f[0]+1145259*f[1]-3375594*f[2]+6095796*f[3]-7392546*f[4]
	//					+6185970*f[5]-3540894*f[6]+1328724*f[7]-295326*f[8]+29531*f[9] ) / ( 15120.0*h*h*h ) );
	results.push_back( (-49*f[0]+232*f[1]-461*f[2]+496*f[3]-307*f[4]+104*f[5]-15*f[6])/(8.0*h*h*h) );
	//=======================================

	if (direction == -1)	// make it face forward again
		reverse( f.begin(), f.end() );

	return ( results );
}




vector<double> get_central_derivatives(vector<double> & f, double stepsize)
{
	double adjustment = 1.0;	//heuristic
	volatile double h = stepsize;

	double mult_factor = 1.0;

	vector<double> results;

	const int shift = ( f.size() - 1 ) / 2;

	//=======================================
	// derivatives
	// function value
	results.push_back( f[shift+0] );

	//===================
	// first derivative
	results.push_back( (-2*f[shift-5]+25*f[shift-4]-150*f[shift-3]
						+600*f[shift-2]-2100*f[shift-1]+0*f[shift+0]
						+2100*f[shift+1]-600*f[shift+2]+150*f[shift+3]
						-25*f[shift+4]+2*f[shift+5])/(2520.0*h) );
	h *= adjustment;

	//===================
	// second derivative
	results.push_back( (8*f[shift-5]-125*f[shift-4]+1000*f[shift-3]
						-6000*f[shift-2]+42000*f[shift-1]-73766*f[shift+0]
						+42000*f[shift+1]-6000*f[shift+2]+1000*f[shift+3]
						-125*f[shift+4]+8*f[shift+5])/(25200.0*h*h) );
	h *= adjustment;

	//===================
	// third derivative
	results.push_back( (205*f[shift-5]-2522*f[shift-4]+14607*f[shift-3]
						-52428*f[shift-2]+70098*f[shift-1]+0*f[shift+0]
						-70098*f[shift+1]+52428*f[shift+2]-14607*f[shift+3]
						+2522*f[shift+4]-205*f[shift+5])/(30240.0*h*h*h) );




	//===================
	// fourth derivative
	results.push_back( (-82*f[shift-5]+1261*f[shift-4]-9738*f[shift-3]
						+52428*f[shift-2]-140196*f[shift-1]+192654*f[shift+0]
						-140196*f[shift+1]+52428*f[shift+2]-9738*f[shift+3]
						+1261*f[shift+4]-82*f[shift+5])/(15120.0*h*h*h*h) );


	//===================
	// fifth derivative
	results.push_back( (-3903119677054429*f[shift-5]+45636476224021016*f[shift-4]-235087900548739840*f[shift-3]
						+562049233495837630*f[shift-2]-581865071856267900*f[shift-1]+144*f[shift+0]
						+581865071856268000*f[shift+1]-562049233495837800*f[shift+2]+235087900548739840*f[shift+3]
						-45636476224021016*f[shift+4]+3903119677054429*f[shift+5])/(86469112845513500.0*h*h*h*h*h) );


	//===================
	// sixth derivative
	results.push_back( (13*f[shift-5]-190*f[shift-4]+1305*f[shift-3]
						-4680*f[shift-2]+9690*f[shift-1]-12276*f[shift+0]
						+9690*f[shift+1]-4680*f[shift+2]+1305*f[shift+3]
						-190*f[shift+4]+13*f[shift+5])/(240.0*h*h*h*h*h*h) );


	//===================
	// seventh derivative
	results.push_back( (71485708370960*f[shift-5]-743451367057984*f[shift-4]+2959508326557744*f[shift-3]
						-5833233803070335*f[shift-2]+5404319552844575*f[shift-1]-3*f[shift+0]
						-5404319552844576*f[shift+1]+5833233803070337*f[shift+2]-2959508326557744*f[shift+3]
						+743451367057984*f[shift+4]-71485708370960*f[shift+5])/(343131400180608.0*h*h*h*h*h*h*h) );


	//===================
	// eighth derivative
	results.push_back( (-1*f[shift-5]+13*f[shift-4]-69*f[shift-3]
						+204*f[shift-2]-378*f[shift-1]+462*f[shift+0]
						-378*f[shift+1]+204*f[shift+2]-69*f[shift+3]
						+13*f[shift+4]-1*f[shift+5])/(3.0*h*h*h*h*h*h*h*h) );


	//=======================================

	return ( results );
}


vector<double> get_Chebyshev_nodes( double a, double b, int number_of_points )
{
	vector<double> nodes (number_of_points), adjnodes (number_of_points);
	cheb_lib::get_nodes(a, b, nodes, adjnodes);
	vector<double> results = adjnodes;
	//results.assign(adjnodes, adjnodes + number_of_points);

	return (results);
}

vector<double> get_Chebyshev_derivatives( double a, double b, vector<double> & f, int nL_or_nR )
{
	int npts = f.size();
	int Nmax = npts - 1;
	vector<double> nodes (npts), adjnodes (npts);
	cheb_lib::get_nodes(a, b, nodes, adjnodes);

	vector<double> coeffs_array(npts);
	//double * f_array = &f[0];
	//cheb_lib::get_coefficients(coeffs_array, nodes, f_array, npts);
	cheb_lib::get_coefficients(coeffs_array, nodes, f);

	vector<double> Chebyshev_array_C(npts*(nL_or_nR+1));
	double Gamma [nL_or_nR+1];

	// interval spaced symmetrically, so
	// derivatives evaluated exactly at zero
	// ik: order of derivative
	for (int ik = 0; ik <= nL_or_nR; ++ik)
	{
		Gamma[ik] = (ik==0)? 0: gsl_sf_gamma(ik);
		double poly[Nmax+1];
		if (ik == 0)
		{
			for (int iN = 0; iN <= Nmax; ++iN)	// Nmax == npts - 1 here
				Chebyshev_array_C[iN * (nL_or_nR+1) + ik] = csf::Tfun(iN, 0.0);
				//Chebyshev_array_C[iN * (nL_or_nR+1) + ik] = cos(0.5*iN*M_PI);
		}
		else
		{
			//==============================
			// Evaluate Gegenbauer polynomials
			int status = gsl_sf_gegenpoly_array (Nmax, double(ik), 0.0, poly);

			if (status > 0) cerr << "Warning: status > 0 in GSL Gegenbauer function!" << endl;

			//==============================
			// Fill array
			for (int iN = ik; iN <= Nmax; ++iN)	// automatically zero for k>n
			{
				double resultL = pow(2.0, ik-1.0)*iN*Gamma[ik]*poly[iN-ik];
				Chebyshev_array_C[iN * (nL_or_nR+1) + ik] = resultL;
			}
		}
	}

	// store value and derivatives at (a+b)/2
	vector<double> results;
	for (int ik = 0; ik <= nL_or_nR; ++ik)
	{
		double sum = 0.0;
		for (int iN = 0; iN <= Nmax; ++iN)
			sum += coeffs_array[iN] * pow(0.5*(b-a), -double(ik))
					* Chebyshev_array_C[iN * (nL_or_nR+1) + ik];
		results.push_back( sum );
	}

	return (results);
}

