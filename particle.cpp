#include <string>
#include "constants.hpp"
#include "particle.hpp"

using namespace std;

/*
	Contains Particle class which stores the properties of a species of
	particle.
*/


Particle::Particle(std::string name, double mass, double degeneracy, 
	double baryonCharge, QUANTUM_TYPE quantumType, double excludedVolume_exII)
{
	itsName = name;
	itsMass = mass;
	itsDegen = degeneracy;
	itsBaryonCharge = baryonCharge;
	itsExcludedVolume_exII = excludedVolume_exII;

	itsQuantumType = quantumType;
}


void Particle::SetExIIVolProportionalToMass(double epsilon0)
{
	//used to make exII model comparable to model exI.  Sets particle's excluded
	//volume equal to mass/epsilon0.
	itsExcludedVolume_exII = itsMass/(epsilon0);
}

void Particle::SetQuantumType(QUANTUM_TYPE qt)
{
	itsQuantumType = qt;
}



double Particle::GetParticleChemPotential(double mu_B) const
{
	//computes particle chemical potential--different from baryon chemical 
	//potential mu_B
	double mu_particle = itsBaryonCharge*mu_B;
	return mu_particle;
}


std::string Particle::GetName() const
{
	return itsName;
}

double Particle::GetMass() const
{
	return itsMass;
}

double Particle::GetDegen() const
{
	return itsDegen;
}

double Particle::GetBaryonCharge() const
{
	return itsBaryonCharge;
}

double Particle::GetExcludedVolume_exII() const
{
	return itsExcludedVolume_exII;
}

QUANTUM_TYPE Particle::GetQuantumType() const
{
	return itsQuantumType;
}



