#ifndef PARTICLELIST_H
#define PARTICLELIST_H

#include <vector>
#include "particle.hpp"


/*
	A function GetPDGParticleList() will return a vector (list) of
	confirmed resonances (from the PDG 2012 edition) made of up, down,
	and strange quarks.  Resonances with charm and heavier are excluded.

	ParticleList_SetExIIVolsPropToMass will set the exII excluded volume
	of each particle in the list to be mass/epsilon0.

	ParticleList_SetClassicalQuantum_Type will set the the QUANTUM_TYPE
	of each particle in the list to CLASSICAL.
*/

std::vector<Particle> GetPDGParticleList();


void ParticleList_SetExIIVolsPropToMass(double, std::vector<Particle>& );

void ParticleList_SetClassicalQuantum_Type(
 std::vector<Particle>& ParticleList);






#endif

