
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <array>
#include <limits>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>          //for modifying gsl's error handling



#include "constants.hpp"
#include "exceptions.hpp"
#include "particle.hpp"
#include "thermo_pt.hpp"
#include "thermo_exI.hpp"
#include "thermo_exII.hpp"
#include "thermo_pqcd.hpp"
#include "thermo_crossover.hpp"
#include "crossover_line.hpp"
#define SUPPRESS_INTEGRATION_ERRORS 

using namespace std;


std::vector<double> Get_Thermo(double Temp, double mu_B, double R, 
    double dRdpsi, thermparams parameters){

    double psic=parameters.psic; double a=parameters.a; double b=parameters.b;
    double EnergyScaleConstant = parameters.Escl; double SoftenConstant = parameters.Soft;
    double epsilon0 = parameters.eps0_4;
    epsilon0=pow(epsilon0,4);
    double Press_exII_guess=-1;
    double InterpExp;
    InterpExp= parameters.InterpExp;
    MODEL_TYPE model=parameters.model;
    const std::vector<Particle>&  ParticleList=parameters.ParticleList;
    std::vector<double> therm; 
    if(R<0.0){
        therm = {NAN,NAN,NAN,NAN,NAN};
        return therm;
    }
    //printf("T,mu,R,dRdpsi: %f %f %f %f \n",Temp,mu_B,R,dRdpsi);

    double psi =atan2(mu_B,Temp);
  
	if(Temp <= 0.0)
		throw TemperatureError("Get_Thermo expects positive Temp.");

	if( InterpExp <= 0)
		throw ExVolValueError("Get_Thermo expects positive InterpolationExponent");

	if(EnergyScaleConstant <= 0.0)
		throw ExVolValueError("Get_Thermo expects positive EnergyScaleConstant.");

    if(model == EXI){
        if(epsilon0 <= 0.0){
            printf("%f \n",epsilon0);
			throw ExVolValueError("Get_Thermo expects positive epsilon0 for model EXI");
        }
    }


	double  Press_qcd, Press_h, Press, Entropy_qcd, Entropy_h, Entropy,
	 Baryon_qcd, Baryon_h, Baryon, Energy;


    Press_qcd = Total_Pressure_qgp_PQCD(Temp, mu_B,
   EnergyScaleConstant, SoftenConstant);
    Entropy_qcd = Total_Entropy_Dens_qgp_PQCD(Temp, mu_B,
   EnergyScaleConstant, SoftenConstant);
    Baryon_qcd = Total_Baryon_Dens_qgp_PQCD(Temp, mu_B,
   EnergyScaleConstant, SoftenConstant);


    if(model == PT){
        Press_h = Total_Pressure_pt(Temp, mu_B, ParticleList);
        Entropy_h = Total_Entropy_Dens_pt(Temp, mu_B, ParticleList);
        Baryon_h = Total_Baryon_Dens_pt(Temp, mu_B, ParticleList);
        }
    else if(model == EXI){
        Press_h = Total_Pressure_exI(Temp, mu_B, epsilon0, ParticleList);
        Entropy_h = Total_Entropy_Dens_exI(Temp, mu_B, epsilon0, ParticleList);
        Baryon_h = Total_Baryon_Dens_exI(Temp, mu_B, epsilon0, ParticleList);
        }
    else if(model == EXII){
        Press_h = Total_Pressure_exII(Temp, mu_B,
         ParticleList, epsilon0, Press_exII_guess);
        Entropy_h = Total_Entropy_Dens_exII(Temp, mu_B,
         ParticleList, epsilon0, Press_exII_guess);
        Baryon_h = Total_Baryon_Dens_exII(Temp, mu_B,
         ParticleList, epsilon0, Press_exII_guess);
         }
    else{
        throw ModelError("Unrecognized MODEL_TYPE in Get_Thermo");
    }


    double theta, x,targ,tharg, rr;
    double  dSdT, dSdmu, eta1, eta2;
    double dpsidT,dpsidmu, dthetadT,dthetadmu,deta2dtheta, deta1dpsi;
    double S;
    rr=Temp*Temp + mu_B*mu_B;
    if(psi>=psic && rr<=R*R){
        Press=Press_h;
        Entropy=Entropy_h;
        Baryon=Baryon_h;
        S=0.0;


    }else if(psi>=psic &&  rr>=R*R){
        Press=Press_qcd;
        Entropy=Entropy_qcd;
        Baryon=Baryon_qcd;
        S=1.0;

    }else{
        theta=pow(rr/(R*R),-InterpExp/2.0);

        //x=abs(psi/psic);
        x=fabs(psi/psic);
        if(psi==0.0){
            tharg=0.0;
            eta1=1.0;
            deta1dpsi = 0.0;
        }else{
            tharg=a*(b-x)/(x-x*x);
            eta1=0.5*(1.0+tanh(tharg));
            deta1dpsi = (1.0/(cosh(tharg)*cosh(tharg)))/(psi*(1.0-psi/psic));
            deta1dpsi = deta1dpsi*(tharg*(0.5-psi/psic)-a/2.0);
			/*printf("%f\n", eta1);
			printf("%f\n", tharg);
			printf("%f\n", x);
			printf("%f\n", psi);
			printf("%f\n", psic);
			printf("%f\n", psi/psic);
			printf("%f\n", std::abs(psi/psic));
			printf("%f\n", fabs(psi/psic));
			std::cout << "Check: "
						<< psi << "   " << psic << "   " << psi/psic << "   " 
						<< abs((psi/psic)) << "   " << fabs(psi/psic) << std::endl;
			*/
       }   


        targ=pi/pow(2.0,theta)-0.5*pi;
        eta2=tan(targ);
        if(targ<=0.0 && eta2>0.0) eta2=-eta2;
        S= 0.5+atan2(eta2,eta1)/(pi);


        dpsidT=-mu_B/rr; dpsidmu=Temp/rr;
        deta2dtheta=(-pi/(cos(targ)*cos(targ)))*log(2)*pow(2,-theta);
        dthetadT=InterpExp*theta/rr;
        dthetadmu=dthetadT*(Temp*dRdpsi/R - mu_B);
        dthetadT=dthetadT*(-mu_B*dRdpsi/R - Temp);
		//printf("%f   %f   %f   %f   %f   %f   %f   %f\n", eta1, eta2, deta2dtheta, dthetadmu, deta1dpsi, dpsidmu, dthetadT, dpsidT);
        dSdmu=(eta1*deta2dtheta*dthetadmu-eta2*deta1dpsi*dpsidmu)/(eta1*eta1 + eta2*eta2)/pi; 
        dSdT=(eta1*deta2dtheta*dthetadT-eta2*deta1dpsi*dpsidT)/(eta1*eta1 + eta2*eta2)/pi; 
        Press =Press_qcd*S + Press_h*(1.0 - S);
        Entropy=Entropy_qcd*S + Entropy_h*(1.0-S) + (Press_qcd - Press_h)*dSdT;
        Baryon=Baryon_qcd*S + Baryon_h*(1.0-S) + (Press_qcd - Press_h)*dSdmu;

		//printf("%f   %f   %f   %f   %f   %f\n", Entropy_qcd, Entropy_h, dSdT, Baryon_qcd, Baryon_h, dSdmu);
    }

	Energy = Temp*Entropy - Press + mu_B*Baryon;
	therm = {Press,Entropy,Baryon,Energy,S};
   

    return therm;
}

double Get_c2(double Temp, double mu_B, double Rc, thermparams parameters){

    double dx=1.0;
    std::vector<double> Rs;
    double R,dRdpsi,psi;
    double s,q,sl,sr,su,sd,ql,qr,qu,qd;
    std::vector<double> therm;
    double psic=parameters.psic;
    psi=atan2(mu_B,Temp);
    if(psi>psic){ 
        R=Get_R(psi,parameters);
        dRdpsi=0.0;
    }else{
        Rs=Get_Rs_quartic(psi,Rc,parameters);
        R=Rs[0]; dRdpsi=Rs[1];
    }
    therm=Get_Thermo(Temp, mu_B, R, dRdpsi, parameters);
    s=therm[1]; q=therm[2];
    psi=atan2(mu_B,Temp+dx);
    if(psi>psic){ 
        R=Get_R(psi,parameters);
        dRdpsi=0.0;
    }else{
        Rs=Get_Rs_quartic(psi,Rc,parameters);
        R=Rs[0]; dRdpsi=Rs[1];
    }
    therm=Get_Thermo(Temp+dx, mu_B, R, dRdpsi, parameters);
    su=therm[1]; qu=therm[2];


    psi=atan2(mu_B,Temp-dx);
    if(psi>psic){ 
        R=Get_R(psi,parameters);
        dRdpsi=0.0;
    }else{
        Rs=Get_Rs_quartic(psi,Rc,parameters);
        R=Rs[0]; dRdpsi=Rs[1];
    }
    therm=Get_Thermo(Temp-dx, mu_B, R, dRdpsi, parameters);
    sd=therm[1]; qd=therm[2];
    double qT,qmu,sT,smu;
    qT=(qu-qd)/(2*dx); sT=(su-sd)/(2*dx);




    if(mu_B==0.0){

        psi=atan2(mu_B+dx,Temp);
        if(psi>psic){ 
            R=Get_R(psi,parameters);
            dRdpsi=0.0;
        }else{
            Rs=Get_Rs_quartic(psi,Rc,parameters);
            R=Rs[0]; dRdpsi=Rs[1];
        }
        therm=Get_Thermo(Temp, mu_B+dx, R, dRdpsi, parameters);
        sr=therm[1]; qr=therm[2];


        qmu=(qr-q)/(dx); smu=(sr-s)/(dx);

    }else{


        psi=atan2(mu_B+dx,Temp);
        if(psi>psic){ 
            R=Get_R(psi,parameters);
            dRdpsi=0.0;
        }else{
            Rs=Get_Rs_quartic(psi,Rc,parameters);
            R=Rs[0]; dRdpsi=Rs[1];
        }
        therm=Get_Thermo(Temp, mu_B+dx, R, dRdpsi, parameters);
        sr=therm[1]; qr=therm[2];

        psi=atan2(mu_B-dx,Temp);
        if(psi>psic){ 
            R=Get_R(psi,parameters);
            dRdpsi=0.0;
        }else{
            Rs=Get_Rs_quartic(psi,Rc,parameters);
            R=Rs[0]; dRdpsi=Rs[1];
        }
        therm=Get_Thermo(Temp, mu_B-dx, R, dRdpsi, parameters);
        sl=therm[1]; ql=therm[2];

        qmu=(qr-ql)/(2*dx); smu=(sr-sl)/(2*dx);


    }

    
    double z = (s*qT -q*sT)/(q*smu -s*qmu);
    if(mu_B==0) z=0;
    double c2=s/(Temp*s - mu_B*q);
    c2=c2*(s+q*z)/(sT+smu*z);
    if(!std::isfinite(c2)) printf("!!! %f | %e | %e \n",(q*smu -s*qmu),s,(s*qT -q*sT));
    return c2;

}
