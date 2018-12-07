
#include <string>

#include "constants.hpp"
#include "exceptions.hpp"

std::string Model_to_String(MODEL_TYPE model)
{
	std::string str;

  if(model == PT)
		str = "pt";
	else if(model == EXI)
		str = "exI";
	else if(model == EXII)
		str = "exII";
	else
		throw ModelError("Model_to_string expects model to be"
		 "PT, EXI, or EXII.  Unrecognized model type.");

	return str;
}