
#include <vector>
#include "exceptions.hpp"
#include "util.hpp"

using namespace std;

std::vector<double> linspace(double start, double stop, int len)
{
	if(len <= 0)
		throw ExVolValueError("linspace() needs positive len.");

	std::vector<double> array(len);

	double dx = (stop - start)/( (double)len - 1.0);
	for(int i = 0; i < len; i++)
		array[i] = start + i*dx;

	return array;
}
