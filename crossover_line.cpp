#include <cmath>
#include <vector>
#include <array>
#include <omp.h>

#include <gsl/gsl_roots.h> 
#include <gsl/gsl_errno.h>          //for modifying gsl's error handling

#include "exceptions.hpp"
#include "particle.hpp"
#include "particlelist.hpp"
#include "thermo_pt.hpp"
#include "thermo_exI.hpp"
#include "thermo_exII.hpp"
#include "thermo_pqcd.hpp"
#include "thermo_crossover.hpp"
#include "util.hpp"

using namespace std;

//=============================================================================
//Fcns to return best-fit values of crossover EOS.  Best fit values
//were found using MinimizeChiSquare_TCrossover_gsl() in chisquare.cpp
//model specifies the MODEL_TYPE (PT, EXI, EXII) and
//UsingQuantumStatistics=true means use parameters that were computed using
//quantum statistics.

//=============================================================================
//Generate R as a function of psi.
//=============================================================================

struct delPparams{
  MODEL_TYPE model;
  const std::vector<Particle>&  ParticleList;
  double Press_exII_guess;
  double epsilon0;
  double SoftenConstant;
  double EnergyScaleConstant;
  double variant;
  
};
//Objective functions for root finding:
double Delta_P(double R, void * params){
  
  delPparams* p = (delPparams* ) params;
  double psi=p->variant;
  double Temp=R*cos(psi);
  double mu_B=R*sin(psi);
  double Press_h,Press_qcd;
  if(Temp<=0.0) Temp=3.0;
  Press_qcd = Total_Pressure_qgp_PQCD(Temp, mu_B, 
   p->EnergyScaleConstant, p->SoftenConstant);
  if(p->model == PT)
    Press_h = Total_Pressure_pt(Temp, mu_B, p->ParticleList);
  else if(p->model == EXI)
    Press_h = Total_Pressure_exI(Temp, mu_B,p->epsilon0, p->ParticleList);
  else if(p->model == EXII)
    Press_h = Total_Pressure_exII(Temp, mu_B, 
     p->ParticleList, p->epsilon0, p->Press_exII_guess);
  else
    throw ModelError("Unrecognized MODEL_TYPE in "
     "Total_Pressure_crossover_Tempcrossover");
  return Press_qcd-Press_h;
  
} 


double Delta_P_muc(double mu_B, void * params){
  
  delPparams* p = (delPparams* ) params;
  double Tc=p->variant;
  double Press_h,Press_qcd;
  if(Tc<=0.0) Tc=3.0;
  Press_qcd = Total_Pressure_qgp_PQCD(Tc, mu_B, 
   p->EnergyScaleConstant, p->SoftenConstant);
  if(p->model == PT)
    Press_h = Total_Pressure_pt(Tc, mu_B, p->ParticleList);
  else if(p->model == EXI)
    Press_h = Total_Pressure_exI(Tc, mu_B,p->epsilon0, p->ParticleList);
  else if(p->model == EXII)
    Press_h = Total_Pressure_exII(Tc, mu_B, 
     p->ParticleList, p->epsilon0, p->Press_exII_guess);
  else
    throw ModelError("Unrecognized MODEL_TYPE in "
     "Total_Pressure_crossover_Tempcrossover");
  
  return Press_qcd-Press_h;
} 

/////////////////////////////////////////////////////////////////////////

//Finds the distance R = sqrt(T^2 + mu^2) to the crossover line from the
//origin in the T-mu plane at the specified angle psi
double Get_R(double psi, thermparams parameters) { 

  const gsl_root_fsolver_type * T = gsl_root_fsolver_brent;
  gsl_root_fsolver * s = gsl_root_fsolver_alloc (T);
  int status;
  int iter = 0, max_iter = 50;
  double tol=.000000001;

  double R_min; double R_max;
  R_min=300.0 +(300.0/1.6)*psi;
  R_max=250.0*pow(cos(psi)*cos(psi) + sin(psi)*sin(psi)/36.0 ,-.5); //Estimated upper bound for crossover line
  
  double Soft=parameters.Soft; double Escl=parameters.Escl; 
  double epsilon0=parameters.eps0_4;
  epsilon0=pow(epsilon0,4);
  MODEL_TYPE model=parameters.model;
  std::vector<Particle> ParticleList=parameters.ParticleList;

  double R, dR;
  R=1.0;
  gsl_function F;
  struct delPparams params = {model,ParticleList,-1,
    epsilon0, Soft, Escl, psi};
  while(Delta_P(R_max,  &params)>0) {R_max=R_max-40;}

  while(Delta_P(R_min, &params)<0 && R_min<R_max) {R_min=R_min+25;}
  //printf("middle R \n");
  if(R_min>=R_max) {
    return -1.0;
  }else{
    if(!gsl_finite(Delta_P(R_min,  &params))){
      printf("!!!!Infinite DeltaP\n");
      return -1.0;
    }
    F.function = &Delta_P;
    F.params = &params;
    gsl_set_error_handler_off();
    status=gsl_root_fsolver_set(s, &F, R_min, R_max);
    for(iter=0;iter<=max_iter;iter++){

      status = gsl_root_fsolver_iterate (s);
      dR=abs(R-gsl_root_fsolver_root (s));
      if(dR < tol) break;
      R = gsl_root_fsolver_root (s);
      R_min = gsl_root_fsolver_x_lower (s);
      R_max = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (R_min, R_max, 0, 0.001);
      //printf("%4.7f %4.7f | %f %f\n", R_min,R_max,Delta_P(R_min,  &params),Delta_P(R_max,  &params));
    }
    return R;
  }
 }

////////////////////////////////////////////////////////////////////////////////
//Quartic Crossover
////////////////////////////////////////////////////////////////////////////////







//Finds the crossover mu_c for the given value of T_c
double Get_muc(double Tc, thermparams parameters) { 

  const gsl_root_fsolver_type * T = gsl_root_fsolver_brent;
  gsl_root_fsolver * s = gsl_root_fsolver_alloc (T);
 
  int iter = 0, max_iter = 50;
  double tol=.0000001;

  double mu_min; double mu_max;
  
  double Soft=parameters.Soft; double Escl=parameters.Escl; 
  double epsilon0=parameters.eps0_4;
  epsilon0=pow(epsilon0,4);

  MODEL_TYPE model=parameters.model;
  std::vector<Particle> ParticleList=parameters.ParticleList;

  double muc, dmu;
  muc=1.0;
  gsl_function F;
  struct delPparams params = {model,ParticleList,-1,
    epsilon0, Soft, Escl, Tc};

  mu_min=450;
  while(Delta_P_muc(mu_min, &params)<0 && mu_max<4000){
    mu_min=mu_min+100;
  }
  mu_max=mu_min+600;
  while(Delta_P_muc(mu_max, &params)>0 && mu_max<4000){
    mu_max=mu_max+400;
    //if(mu_max>2000) mu_max+600;
    if(mu_max>2000) mu_max=mu_max+600;
    //printf("%f %f ||%f %f \n", mu_min,mu_max,Delta_P_muc(mu_min, &params),Delta_P_muc(mu_max, &params));
  }
  if(mu_max>=4000) {
    return -1.0;
  }else{
    if(!gsl_finite(Delta_P_muc(mu_min,  &params))){
      printf("!!!!Infinite DeltaP\n");
      return -1.0;
    }
    F.function = &Delta_P_muc;
    F.params = &params;
    gsl_set_error_handler_off();
    gsl_root_fsolver_set(s, &F, mu_min, mu_max);

    for(iter=0;iter<=max_iter;iter++){
      gsl_root_fsolver_iterate (s);
      dmu=abs(muc-gsl_root_fsolver_root (s));
      if(dmu < tol) break;
      muc = gsl_root_fsolver_root (s);
      mu_min = gsl_root_fsolver_x_lower (s);
      mu_max = gsl_root_fsolver_x_upper (s);
      //status = gsl_root_test_interval (mu_min, mu_max, 0, 0.001);
    }
    gsl_root_fsolver_free (s);
    return muc;
  }
}



// Gets the appropriate value for k2 assuming it is 
// constrained to have continuous derivative at critical point.
double Get_k2(double Rc, thermparams parameters){


  double Soft=parameters.Soft; double Escl=parameters.Escl;
  double epsilon0=parameters.eps0_4; double psic=parameters.psic;
  epsilon0=pow(epsilon0,4);
  MODEL_TYPE model=parameters.model;
  std::vector<double> Ks=parameters.Ks;
  std::vector<Particle> ParticleList=parameters.ParticleList;
  if(Rc==-1.0) return -100.0;
  double Entropy_qcd, Baryon_qcd, Entropy_h=0, Baryon_h=0,dTdmu,k2,Temp,mu_B,B;
  Temp=Rc*cos(psic); mu_B= Rc*sin(psic);
  Entropy_qcd = Total_Entropy_Dens_qgp_PQCD(Temp, mu_B,
    Escl, Soft);
  Baryon_qcd = Total_Baryon_Dens_qgp_PQCD(Temp, mu_B,
    Escl, Soft);

  if(model == PT){
    Entropy_h = Total_Entropy_Dens_pt(Temp, mu_B, ParticleList);
    Baryon_h = Total_Baryon_Dens_pt(Temp, mu_B, ParticleList);
  }
  else if(model == EXI){
    Entropy_h = Total_Entropy_Dens_exI(Temp, mu_B, epsilon0, ParticleList);
    Baryon_h = Total_Baryon_Dens_exI(Temp, mu_B, epsilon0, ParticleList);
  }
  else if(model == EXII){
    Entropy_h = Total_Entropy_Dens_exII(Temp, mu_B,
     ParticleList, epsilon0, -1);
    Baryon_h = Total_Baryon_Dens_exII(Temp, mu_B,
     ParticleList, epsilon0, -1);
  }
  dTdmu=-(Baryon_h-Baryon_qcd)/(Entropy_h-Entropy_qcd);
  

  B=(mu_B/Temp)*dTdmu;
  k2=-B;
  int len=Ks.size();
  for(int iter=1;iter<len;iter++){
    k2=k2 + (B-(2+2*iter))*Ks[iter];
  }
  k2=k2/(2-B);

  
  return k2;
}






//Finds the distance R = sqrt(T^2 + mu^2) to the crossover line from the
//origin in the T-mu plane at the specified angle psi. Also determines the
//derivative dRdpsi at this point. Returns as vector Rs={R,dRdpsi}.
std::vector<double> Get_Rs_quartic(double psi, double Rc, thermparams parameters) {

  double psic=parameters.psic;
  double Tc=Rc*cos(psic); double muc=Rc*sin(psic);
  std::vector<double> Ks=parameters.Ks;
  double T0=1; double tol=.00001;
  int len=Ks.size();
  for(int iter=0;iter<len;iter++){
    T0=T0 - Ks[iter];
  }
  T0=Tc/T0;
  if(psi==0) return {T0,0.0};

 
  //Manually implement root finding for 0=1 - mu / T0tan(psi) - k2 mu^2 - k4 mu^4 - ...
  double m, m_min, m_max,fmin,fmax,f;
  double tpsi=tan(psi);

  m_max=T0*tpsi/muc;
  if(Ks[0]<0) m_max=m_max-Ks[0]; 
  m_min=T0/(T0 +(muc/tpsi) -Tc);
  fmin=1 - m_min*muc/(T0*tpsi);
  fmax=1 - m_max*muc/(T0*tpsi);  
  for(int iter=0;iter<len;iter++){
    fmin=fmin - Ks[iter]*pow(m_min,2+2*iter);
    fmax=fmax - Ks[iter]*pow(m_max,2+2*iter);
  }
  //printf(":::%f %f | %f %f\n", psi,Tc,muc,k2);
  //printf("%f %f | %f %f\n", m_min,m_max,fmin,fmax);
  if(fmin<fmax) throw TemperatureError("Get_Rs_quartic couldn't determine crossover mu.");
  f=.5*(fmax+fmin);
  while(abs(f)>tol){
    m=(m_min+m_max)/2.0;
    f=1 - m*muc / (T0*tpsi);
    for(int iter=0;iter<len;iter++){
      f=f - Ks[iter]*pow(m,2+2*iter);
    }
    if(f<0){
      fmax=f; m_max=m;
    } else {
      fmin=f; m_min=m;
    }
    //printf("%f %f | %f %f\n", m_min,m_max,fmin,fmax);
  }

  m=m_min*fmax +m_max*fmin;
  m=m/(fmin+fmax);



  double R,dRdpsi;
  R=m*muc/sin(psi);

  dRdpsi=muc/(T0*tpsi);
  for(int iter=0;iter<len;iter++){
    dRdpsi=dRdpsi + (2+2*iter)*Ks[iter]*pow(m,1+2*iter);
  }
  dRdpsi=((muc*muc*m)/(T0*pow(sin(psi),3)))/dRdpsi;
  dRdpsi= dRdpsi-(R/tpsi);

  if(Rc==-1.0) R=-1.0;

  return {R,dRdpsi};
}
