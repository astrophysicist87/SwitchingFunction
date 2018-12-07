#ifndef UTIL_H
#define UTIL_H

#include <vector>

//cloning numpy's linspace.  Returns vector of length len,
//evenly spaced from value start to stop.
std::vector<double> linspace(double start, double stop, int len);


#endif