#ifndef CONSTANTS_H
#define CONSTANTS_H


#include <string>
#include "exceptions.hpp"


const int MAX_NUMBER_ITER = 200;  //max number of iterations for solvers I coded

const double pi = 3.14159265359;


enum QUANTUM_TYPE {FERMION, BOSON, CLASSICAL};  //specifies the quantum 
                                                //statistics of a particle

enum MODEL_TYPE{PT, EXI, EXII};  //specifies a hadronic model to use



//returns a string of "pt", "exI", "exII" given model
std::string Model_to_String(MODEL_TYPE model);



#endif

