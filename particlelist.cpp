#include <vector>
#include "constants.hpp"
#include "particle.hpp"
#include "particlelist.hpp"

using namespace std;


/*
	A function GetPDGParticleList() will return a vector (list) of
	confirmed resonances (from the PDG 2012 edition) made of up, down,
	and strange quarks.  Resonances with charm and heavier are excluded.
*/




//Note: can invoke classical statistics if set particles'
// QUANTUM_TYPE = CLASSICAL






/*-----------------------------------------------------------------------------
 Global data used in this program.  Listed are confirmed resonances from
 2012 pdg.  (The mesons are listed in pdg meson summary starting pg 34.
 A concise list of mesons used are mesons on pg 77 in the table with a dot by
 them--indicating existence is confirmed.  The baryons used are those from the
 table on page 78 with 3 or 4 stars--indicating confirmation.  Baryons with 2 
 or fewer stars were not used.  Baryon properties were taken from page 79 and
 on of pdg baryon summary.) This selection strategy was used in 
 arXiv:0901.1430v1.

 When masses are listed as a range, like 1200-1500, I picked the midpoint of 
 the range.

 When degeneracies are not listed, I make up conservative (ie, small) values 
 so I underestimate instead of overestimate the resonance contribution.  

 Note: There is no need to list (explicitly) the anti-particles: 
  For mesons, I account for anti-particles via the degeneracies.
  For baryons, I create an anti-baryon from each listed baryon using
   the opposite baryon number.

 Entries below are of the form 
  (particlename, mass (MeV), degeneracy, baryonNumber, QUANTUM_TYPE, 
   excludedVolume_exII)
*/


std::vector<Particle> GetPDGParticleList()
{

	std::vector<Particle> ParticleList;

	//-----------------------------------------------------------------------------
	//Light Unflavored Mesons:

	ParticleList.push_back(Particle("pi+", 139.6, 2.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("pi0", 135.0, 1.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("eta", 547.8, 1.0, 0.0, BOSON, 0.0));

	//Many sources omit this resonance due to its huge width
	//ParticleList.push_back(Particle("f0(500)", 475.0, 1.0, 0.0, BOSON, 0.0));

	ParticleList.push_back(Particle("rho(770)", 775.5, 9.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("omega(782)", 782.6, 3.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("eta'(958)", 957.8, 1.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("f0(980)", 990.0, 1.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("a0(980)", 980.0, 3.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("phi(1020)", 1019.5, 3.0, 0.0, BOSON, 0.0));

	ParticleList.push_back(Particle("h1(1170)", 1170.0, 3.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("b1(1235)", 1229.5, 9.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("a1(1260)", 1230.0, 9.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("f2(1270)", 1275.1, 5.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("f1(1285)", 1282.1, 3.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("eta(1295)", 1294.0, 1.0, 0.0, BOSON, 0.0));

	//This resonance was omitted in arXiv:0901.1430v1.  Why?  because of the
	// broad width?  Typo? The width is no worse than the a1(1260).
	ParticleList.push_back(Particle("pi(1300)", 1300.0, 3.0, 0.0, BOSON, 0.0));

	ParticleList.push_back(Particle("a2(1310)", 1318.3, 15.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("f0(1370)", 1350.0, 1.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("pi1(1400)", 1354.0, 9.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("eta(1405)", 1408.9, 1.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("f1(1420)", 1426.4, 3.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("omega(1420)", 1425.0, 3.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("a0(1450)", 1474.0, 3.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("rho(1450)", 1465.0, 9.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("eta(1475)", 1476, 1.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("f0(1500)", 1505.0, 1.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("f2'(1525)", 1525.0, 5.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("pi1(1600)", 1662.0, 9.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("eta2(1645)", 1617.0, 5.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("omega(1650)", 1670.0, 3.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("omega3(1670)", 1667.0, 7.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("pi2(1670)", 1672.2, 15.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("phi(1680)", 1680, 3.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("rho3(1690)", 1688.8, 21.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("rho(1700)", 1720.0, 9.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("f0(1710)", 1720.0, 1.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("pi(1800)", 1812.0, 3.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("phi3(1850)", 1854.0, 7.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("pi2(1880)", 1895.0, 15.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("f2(1950)", 1944.0, 5.0, 0.0, BOSON, 0.0));

	//Passing 2 GeV threshold--some papers end here
	ParticleList.push_back(Particle("f2(2010)", 2011.0, 5.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("a4(2040)", 1996.0, 27.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("f4(2050)", 2018.0, 9.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("phi(2170)", 2175.0, 3.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("f2(2300)", 2297.0, 5.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("f2(2340)", 2339.0, 5.0, 0.0, BOSON, 0.0));


	//-----------------------------------------------------------------------------
	//Strange Mesons (pdg starting on pg 39):

	ParticleList.push_back(Particle("K+", 493.7, 2.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("K0", 497.6, 2.0, 0.0, BOSON, 0.0));

	//I do not include KL or KS since I think they are not independent of K0 and 
	//anti-K0

	//Here I write the 4 K* states analogously to the 4 K states
	ParticleList.push_back(Particle("K*(892)+", 891.7, 6.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("K*(892)0", 895.9, 6.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("K*2(1430)+", 1425.6, 10.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("K*2(1430)0", 1432.4, 10.0, 0.0, BOSON, 0.0));

	/*Here I double degeneracy from pdg value to include particle and
		anti-particle in a single listing.  (e.g., K1 should be K1 and K1_bar each 
		with degen 6; I combined them into K1 with degen 12.0)*/
	ParticleList.push_back(Particle("K1(1270)", 1272.0, 12.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("K1(1400)", 1403.0, 12.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("K*(1410)", 1414.0, 12.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("K*0(1430)", 1425.0, 4.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("K*(1680)", 1717.0, 12.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("K2(1770)", 1773.0, 20.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("K3*(1780)", 1776.0, 28.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("K2(1820)", 1816.0, 20.0, 0.0, BOSON, 0.0));

	//Passing 2 GeV threshold--some papers end here
	ParticleList.push_back(Particle("K4*(2045)", 2045.0, 36.0, 0.0, BOSON, 0.0));



	/*  Exclude particles with charm:

	//-----------------------------------------------------------------------------
	//Charmed Mesons (pdg starting on pg 43):

	ParticleList.push_back(Particle("D+", 1869.6, 2.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("D0", 1864.9, 2.0, 0.0, BOSON, 0.0));

	//CHECK THESE:

	ParticleList.push_back(Particle("D*(2007)0", 2007.0, 6.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("D*(2010)+", 2010.3, 6.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("D*(2400)0", 2318.0, 2.0, 0.0, BOSON, 0.0));

	// D*(2400)+  not confirmed, not listed here

	ParticleList.push_back(Particle("D1(2420)0", 2421.3, 6.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("D2*(2460)0", 2462.6, 10.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("D2*(2460)+", 2464.4, 10.0, 0.0, BOSON, 0.0));

	//-----------------------------------------------------------------------------
	//Charmed, Strange Mesons (pdg starting on pg 48):

	//careful with the degeneracy: must double pdg value to include anti-particle

	ParticleList.push_back(Particle("Ds+", 1968.5, 2.0, 0.0, BOSON, 0.0));

	//This is not certain: J is not listed, so I use value J=1 as suggested when 
	//finding degeneracy
	ParticleList.push_back(Particle("Ds*+", 2112.3, 6.0, 0.0, BOSON, 0.0));

	ParticleList.push_back(Particle("Ds0*(2317)+", 2317.8, 2.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("Ds1*(2460)+", 2459.6, 6.0, 0.0, BOSON, 0.0));
	ParticleList.push_back(Particle("Ds1(2536)+", 2535, 6.0, 0.0, BOSON, 0.0));

	//this is not certain: using suggested J=2 value
	ParticleList.push_back(Particle("Ds2*(2573)", 2571.9, 10.0, 0.0, BOSON, 0.0));

	//-----------------------------------------------------------------------------
	//c-cbar Mesons:

	ParticleList.push_back(Particle("etac(1S)", 2981.0, 1.0, 0.0, BOSON, 0.0));

	*/




	//-----------------------------------------------------------------------------
	//N Baryons (pgd starting on page 79):

	//Baryons have baryon number, so anti-baryons must be explicitly listed.
	//To save typing, I will enter baryons and automate the creation of 
	//anti-baryons

	std::vector<Particle> BaryonList;


	BaryonList.push_back(Particle("p", 938.8, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("n", 939.6, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("N(1440)", 1440.0, 4.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("N(1520)", 1520.0, 8.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("N(1535)", 1535.0, 4.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("N(1650)", 1655.0, 4.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("N(1675)", 1675.0, 12.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("N(1680)", 1685.0, 12.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("N(1700)", 1700.0, 8.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("N(1710)", 1710.0, 4.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("N(1720)", 1720.0, 8.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("N(1875)", 1875.0, 8.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("N(1900)", 1900.0, 8.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("N(2190)", 2190.0, 16.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("N(2200)", 2250.0, 20.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("N(2250)", 2275.0, 20.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("N(2600)", 2600.0, 24.0, 1.0, FERMION, 0.0));


	//-----------------------------------------------------------------------------
	//Delta Baryons (pgd starting on pg 82)

	BaryonList.push_back(Particle("Delta(1232)", 1232.0, 16.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Delta(1600)", 1600.0, 16.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Delta(1620)", 1630.0, 8.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Delta(1700)", 1700.0, 16.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Delta(1905)", 1880.0, 24.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Delta(1910)", 1890.0, 8.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Delta(1920)", 1920.0, 16.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Delta(1930)", 1950.0, 24.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Delta(1950)", 1930.0, 32.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Delta(2420)", 2420.0, 48.0, 1.0, FERMION, 0.0));

	//-----------------------------------------------------------------------------
	//Lambda Baryons (pdg starting on pg 83)

	BaryonList.push_back(Particle("Lambda", 1115.7, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Lambda(1405)", 1405.1, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Lambda(1520)", 1519.5, 4.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Lambda(1600)", 1600.0, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Lambda(1670)", 1670.0, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Lambda(1690)", 1690.0, 4.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Lambda(1800)", 1800.0, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Lambda(1810)", 1810.0, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Lambda(1820)", 1820.0, 6.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Lambda(1830)", 1830.0, 6.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Lambda(1890)", 1890.0, 4.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Lambda(2100)", 2100.0, 8.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Lambda(2110)", 2110.0, 6.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Lambda(2350)", 2350.0, 10.0, 1.0,FERMION, 0.0));

	//-----------------------------------------------------------------------------
	//Sigma Baryons (pdg starting on pg 84)

	BaryonList.push_back(Particle("Sigma+", 1189.4, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Sigma0", 1192.6, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Sigma-", 1197.4, 2.0, 1.0, FERMION, 0.0));

	BaryonList.push_back(Particle("Sigma(1385)", 1385.0, 12.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Sigma(1660)", 1660.0, 6.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Sigma(1670)", 1670.0, 12.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Sigma(1750)", 1750.0, 6.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Sigma(1775)", 1775.0, 18.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Sigma(1915)", 1915.0, 18.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Sigma(1940)", 1940.0, 12.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Sigma(2030)", 2030.0, 24.0, 1.0, FERMION, 0.0));

	//J not listed, so I will be conservative and assume J=1/2 which implies 
	//degen=6.0
	BaryonList.push_back(Particle("Sigma(2250)", 2250.0, 6.0, 1.0, FERMION, 0.0));

	//-----------------------------------------------------------------------------
	// Xi baryons (cascades) (pdg starting on pg 86)

	BaryonList.push_back(Particle("Xi0", 1314.9, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Xi-", 1321.7, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Xi(1530)0", 1531.8, 4.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Xi(1530)-", 1535.0, 4.0, 1.0, FERMION, 0.0));

	//J not listed, suggested value is J=1/2
	BaryonList.push_back(Particle("Xi(1690)", 1690.0, 4.0, 1.0, FERMION, 0.0));

	BaryonList.push_back(Particle("Xi(1820)", 1823, 8.0, 1.0, FERMION, 0.0));

	//J is not listed, I will guess a conservative J=1/2 to estimate 
	//degeneracy = 4.0
	BaryonList.push_back(Particle("Xi(1950)", 1950, 4.0, 1.0, FERMION, 0.0));

	//J is not certain
	BaryonList.push_back(Particle("Xi(2030)", 2025.0, 12.0, 1.0, FERMION, 0.0));

	//-----------------------------------------------------------------------------
	//Omega baryons

	BaryonList.push_back(Particle("Omega-", 1672.5, 4.0, 1.0, FERMION, 0.0));

	//J not listed, I will guess the conservative J=1/2
	BaryonList.push_back(Particle("Omega(2250)-", 2252.0, 2.0, 1.0, FERMION, 0.0));


	/* Exclude charmed baryons
	//-----------------------------------------------------------------------------
	//Charmed baryons

	BaryonList.push_back(Particle("Lambdac+", 2286.5, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Lambdac(2595)+", 2592.3, 2.0, 1.0,FERMION,0.0));
	BaryonList.push_back(Particle("Lambdac(2625)+", 2628.1, 4.0, 1.0,FERMION,0.0));
	BaryonList.push_back(Particle("Lambdac(2880)+", 2881.5, 6.0, 1.0,FERMION,0.0));

	//J not listed, I guess the conservative J=1/2
	BaryonList.push_back(Particle("Lambdac(2940)+", 2939.3, 2.0, 1.0,FERMION,0.0));

	BaryonList.push_back(Particle("Sigmac(2455)++", 2454.0, 2.0, 1.0,FERMION,0.0));
	BaryonList.push_back(Particle("Sigmac(2455)+", 2452.9, 2.0, 1.0, FERMION,0.0));
	BaryonList.push_back(Particle("Sigmac(2455)0", 2453.7, 2.0, 1.0, FERMION,0.0));

	BaryonList.push_back(Particle("Sigmac(2520)++", 2517.9, 2.0, 1.0,FERMION,0.0));
	BaryonList.push_back(Particle("Sigmac(2520)+", 2517.5, 2.0, 1.0, FERMION,0.0));
	BaryonList.push_back(Particle("Sigmac(2520)0", 2518.8, 2.0, 1.0, FERMION,0.0));

	//J not listed, I guess the conservative J=1/2
	BaryonList.push_back(Particle("Sigmac(2800)++", 2801.0, 2.0, 1.0,FERMION,0.0));
	BaryonList.push_back(Particle("Sigmac(2800)+", 2792.0, 2.0, 1.0, FERMION,0.0));
	BaryonList.push_back(Particle("Sigmac(2800)0", 2806.0, 2.0, 1.0, FERMION,0.0));

	//Check degeneracy: I think Xic+ and Xic0 are in isospin doublet, so their 
	//only degen comes from J
	BaryonList.push_back(Particle("Xic+", 2467.8, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Xic0", 2470.9, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Xic'+", 2575.6, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Xic'0", 2577.9, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Xic(2645)+", 2645.9, 4.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Xic(2645)0", 2645.9, 4.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Xic(2790)+", 2789.1, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Xic(2790)0", 2791.8, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Xic(2815)+", 2816.6, 4.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Xic(2815)0", 2819.6, 4.0, 1.0, FERMION, 0.0));

	//These have no J listed, I use J=1/2 to estimage degeneracy
	BaryonList.push_back(Particle("Xic(2980)+", 2971.4, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Xic(2980)0", 2968.0, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Xic(3080)+", 3077.0, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Xic(3080)0", 3079.9, 2.0, 1.0, FERMION, 0.0));

	BaryonList.push_back(Particle("Omegac0", 2695.2, 2.0, 1.0, FERMION, 0.0));
	BaryonList.push_back(Particle("Omegac(2770)0", 2770.0, 4.0, 1.0,FERMION, 0.0));
	*/


	//--------------------------------------------------------------
	//Now include baryons and anti-baryons in ParticleList:

	//C++11 for loop means: for each baryon b in BaryonList:
	for(Particle &b: BaryonList) 
	{
		ParticleList.push_back(b);
		ParticleList.push_back(Particle("anti-" + b.GetName(), b.GetMass(),
		 b.GetDegen(), -b.GetBaryonCharge(), b.GetQuantumType(), 
		 b.GetExcludedVolume_exII() ));
	}

	return ParticleList;
}


//-----------------------------------------------------------------------------


void ParticleList_SetExIIVolsPropToMass(double epsilon0,
 std::vector<Particle>& ParticleList)
{
	for(Particle &p: ParticleList) 
	{
		p.SetExIIVolProportionalToMass(epsilon0);
	}
}


void ParticleList_SetClassicalQuantum_Type(std::vector<Particle>& ParticleList)
{
	for(Particle &p: ParticleList) 
	{
		p.SetQuantumType(CLASSICAL);
	}
}







