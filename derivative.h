#include <vector>
#include <algorithm>
#include <array>
#include <iostream>
#include <cmath>

using namespace std;

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




vector<double> get_central_derivatives(vector<double> & f, double stepsize, int direction)
{
	double adjustment = 1.0;	//heuristic
	double h = stepsize;

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
	//=======================================

	return ( results );
}


