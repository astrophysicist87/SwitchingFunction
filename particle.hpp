#ifndef PARTICLE_H
#define PARTICLE_H

/*
	Contains Particle class which stores the properties of a species of
	particle.
*/


#include <string>
#include "constants.hpp"

//-----------------------------------------------------------------------------

class Particle
{
	//stores properties of one species of particle
	public:

		Particle(std::string name, double mass, double degeneracy, 
		 double baryonCharge, QUANTUM_TYPE quantumType, 
		 double excludedVolume_exII);

		

		double GetParticleChemPotential(double mu_B) const;

		std::string GetName() const;
		double GetMass() const;
		double GetDegen() const;
		double GetBaryonCharge() const;
		double GetExcludedVolume_exII() const;
		QUANTUM_TYPE GetQuantumType() const;

		void SetExIIVolProportionalToMass(double epsilon0);
		void SetQuantumType(QUANTUM_TYPE);

	private:
		std::string itsName;
		double itsMass;
		double itsDegen;                //degeneracy of the particle
		double itsBaryonCharge;
		double itsExcludedVolume_exII;  //excluded volume of particle in model EXII
		QUANTUM_TYPE itsQuantumType;
};




#endif
