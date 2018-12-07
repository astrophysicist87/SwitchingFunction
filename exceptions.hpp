#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <string>
#include <stdexcept>


class ConvergenceError : public std::runtime_error
{
	//Thrown when a solver fails to converge.
	public:
		ConvergenceError(const std::string& msg = "")
		 : std::runtime_error(msg) {}
};

class TemperatureError : public std::domain_error
{
	//Thrown when a function receives a non-positive temperature.
	public:
		TemperatureError(const std::string& msg = "")
		 : std::domain_error(msg) {}
};

class epsilon0Error : public std::domain_error
{
	//Thrown when a function receives a non-positive epsilon0.
	public:
		epsilon0Error(const std::string& msg = "") 
		 : std::domain_error(msg) {}
};

class QuantumTypeError : public std::invalid_argument
{
	//Thrown when a function receives an unrecognized quantum type.  Expects
	//boson, fermion, or classical.
	public:
		QuantumTypeError(const std::string& msg = "")
		 : std::invalid_argument(msg) {}
};

class ModelError : public std::invalid_argument
{
	//Thrown when a function receives an unrecognized (hadronic) model.  Expects
	//MODEL_TYPE = PT, EXI, or EXII.
	public:
		ModelError(const std::string& msg = "")
		 : std::invalid_argument(msg) {}
};

class ExVolValueError : public std::invalid_argument
{
	//Thrown when a function receives an invalid argument not covered by
	//the other derived exception classes.
	public:
		ExVolValueError(const std::string& msg = "")
		 : std::invalid_argument(msg) {}
};


class GenericError: public std::runtime_error
{
	//Some kind of error occured whose type is unknown or unspecified
	public:
		GenericError(const std::string& msg = "") 
		 : std::runtime_error(msg) {}
};

class FolderError: public std::runtime_error
{
	//problem working with folder, such as failure to create folder
	public:
		FolderError(const std::string& msg = "") 
		 : std::runtime_error(msg) {}
};

class FileError: public std::runtime_error
{
	//problem working with file, such as failure to create file or write to it.
	public:
		FileError(const std::string& msg = "") 
		 : std::runtime_error(msg) {}
};


#endif
