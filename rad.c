//KORAL - rad.c
//radiation-related routines

#include "ko.h"


//*********************************************************************
//******* calculates total opacity over dx[] ***************************
//**********************************************************************
int
calc_tautot(ldouble *pp, ldouble *xx, ldouble *dx, ldouble *tautot)
{
  ldouble rho=pp[0];
  ldouble u=pp[1];  
  ldouble pr=(GAMMA-1.)*(u);
  ldouble T=pr*MU_GAS*M_PROTON/K_BOLTZ/rho;

  //xx[0] holds time
  ldouble kappa=calc_kappa(rho,T,xx[1],xx[2],xx[3]);
  ldouble chi=kappa+calc_kappaes(rho,T,xx[1],xx[2],xx[3]);

  tautot[0]=chi*dx[0];
  tautot[1]=chi*dx[1];
  tautot[2]=chi*dx[2];

  return 0;
}

//**********************************************************************
//******* calculates abs opacity over dx[] ***************************
//**********************************************************************
int
calc_tauabs(ldouble *pp, ldouble *xx, ldouble *dx, ldouble *tauabs)
{
  ldouble rho=pp[0];
  ldouble u=pp[1];  
  ldouble pr=(GAMMA-1.)*(u);
  ldouble T=pr*MU_GAS*M_PROTON/K_BOLTZ/rho;

  //xx[0] holds time
  ldouble kappa=calc_kappa(rho,T,xx[1],xx[2],xx[3]);

  tauabs[0]=kappa*dx[0];
  tauabs[1]=kappa*dx[1];
  tauabs[2]=kappa*dx[2];

  return 0;
}

      
//**********************************************************************
//******* solves implicitidly four-force source terms *********************
//******* in the lab frame, returns ultimate deltas ***********************
//******* the fiducial approach *****************************************
//**********************************************************************

int f_implicit_lab(ldouble *uu0,ldouble *uu,ldouble *pp,ldouble dt,ldouble gg[][5], ldouble GG[][5],ldouble tup[][4], ldouble tlo[][4],ldouble *f)
{
  ldouble Rij[4][4];
  ldouble pp2[NV];
  int iv;
  for(iv=0;iv<NV;iv++)
    pp2[iv]=pp[iv];

  //opposite changes in gas quantities
  uu[1] = uu0[1] - (uu[6]-uu0[6]);
  uu[2] = uu0[2] - (uu[7]-uu0[7]);
  uu[3] = uu0[3] - (uu[8]-uu0[8]);
  uu[4] = uu0[4] - (uu[9]-uu0[9]);

  //calculating primitives  
  int corr,fixup[2];
  if(u2p(uu,pp2,gg,GG,tup,tlo,&corr,fixup)<0) return -1;

  //radiative four-force
  ldouble Gi[4];
  calc_Gi(pp2,gg,GG,Gi); 
  indices_21(Gi,Gi,gg);
 
  f[0] = uu[6] - uu0[6] + dt * Gi[0];
  f[1] = uu[7] - uu0[7] + dt * Gi[1];
  f[2] = uu[8] - uu0[8] + dt * Gi[2];
  f[3] = uu[9] - uu0[9] + dt * Gi[3];
 
  return 0;
} 

int
print_state_implicit_lab (int iter, ldouble *x, ldouble *f)
{
  printf ("iter = %3d x = % .3e % .3e % .3e % .3e "
	  "f(x) = % .3e % .3e % .3e % .3e\n",
	  iter,
	  x[0],x[1]/x[0],x[2]/x[0],x[3]/x[0],f[0],f[1],f[2],f[3]);
}

int
solve_implicit_lab(int ix,int iy,int iz,ldouble dt,ldouble* deltas)
{
  int i1,i2,i3,iv,i,j;
  ldouble J[4][4],iJ[4][4];
  ldouble pp[NV],uu[NV],uu0[NV],uup[NV]; 
  ldouble f1[4],f2[4],f3[4],xxx[4];
  ldouble gg[4][5],GG[4][5], tlo[4][4],tup[4][4];

  pick_g(ix,iy,iz,gg);
  pick_G(ix,iy,iz,GG);
  pick_T(tmuup,ix,iy,iz,tup);
  pick_T(tmulo,ix,iy,iz,tlo);

  for(iv=0;iv<NV;iv++)
    {
      pp[iv]=get_u(p,iv,ix,iy,iz);      
      uu[iv]=get_u(u,iv,ix,iy,iz);  
      uu0[iv]=uu[iv];
    }

  ldouble EPS = 1.e-8;
  ldouble CONV = 1.e-6 ;

  int verbose=0;
  int iter=0;

  if(verbose) 
    {
      ldouble xx,yy,zz;
      xx=get_x(ix,0);
      yy=get_x(iy,1);
      zz=get_x(iz,2);
      printf("=== i: %d %d %d\n=== x: %e %e %e\n",ix,iy,iz,xx,yy,zz);
      print_Nvector(pp,NV);
      print_metric(gg);
    }
  do
    {
      iter++;
      
      for(i=0;i<NV;i++)
	{
	  uup[i]=uu[i];
	}

      //values at zero state
      if(f_implicit_lab(uu0,uu,pp,dt,gg,GG,tup,tlo,f1)<0) return -1;
 
      //calculating approximate Jacobian
      for(i=0;i<4;i++)
	{
	  for(j=0;j<4;j++)
	    {
	      ldouble del;
	      if(uup[j+6]==0.) del=EPS*uup[6];
	      else del=EPS*uup[j+6];
	      uu[j+6]=uup[j+6]-del;

	      if(f_implicit_lab(uu0,uu,pp,dt,gg,GG,tup,tlo,f2)<0) return -1;
     
	      J[i][j]=(f2[i] - f1[i])/(uu[j+6]-uup[j+6]);

	      uu[j+6]=uup[j+6];
	    }
	}

      //inversion
      inverse_44matrix(J,iJ);

      //updating x
      for(i=0;i<4;i++)
	{
	  xxx[i]=uup[i+6];
	}

      for(i=0;i<4;i++)
	{
	  for(j=0;j<4;j++)
	    {
	      xxx[i]-=iJ[i][j]*f1[j];
	    }
	}

      if(verbose>0)    print_state_implicit_lab (iter,xxx,f1); 

      for(i=0;i<4;i++)
	{
	  uu[i+6]=xxx[i];
	}
  
      //test convergence
      for(i=0;i<4;i++)
	{
	  f3[i]=(uu[i+6]-uup[i+6]);
	  f3[i]=fabs(f3[i]/uup[6]);
	}

      if(f3[0]<CONV && f3[1]<CONV && f3[2]<CONV && f3[3]<CONV)
	{
	  if(verbose) printf("success ===\n");
	  break;
	}

      if(iter>50)
	{
	  printf("iter exceeded in solve_implicit_lab()\n");	  
	  return -1;
	}
     
    }
  while(1);

  deltas[0]=uu[6]-uu0[6];
  deltas[1]=uu[7]-uu0[7];
  deltas[2]=uu[8]-uu0[8];
  deltas[3]=uu[9]-uu0[9];
  
  return 0;
}


//**********************************************************************
//******* solves implicitly gas - radiation interaction ****************
//******* in the fluid frame, returns vector of deltas *****************
//******* the explicit-implicit approximate approach *******************
//******* used as the fail-safe backup method **************************
//**********************************************************************
int
solve_implicit_ff(int ix,int iy,int iz,ldouble dt,ldouble* deltas)
{
  int i1,i2,i3,iv;
  ldouble pp[NV];
  ldouble tup[4][4],tlo[4][4];
  pick_T(tmuup,ix,iy,iz,tup);
  pick_T(tmulo,ix,iy,iz,tlo);
  ldouble gg[4][5];
  pick_g(ix,iy,iz,gg);
  ldouble GG[4][5];
  pick_G(ix,iy,iz,GG);

  for(iv=0;iv<NV;iv++)
    {
      pp[iv]=get_u(p,iv,ix,iy,iz);      
    }

  //transforming radiative primitives to ortonormal fluid frame
  prad_lab2ff(pp,pp,gg,GG,tup);
  
  //four-force in the fluid frame
  ldouble Gi[4];
  calc_Gi_ff(pp,Gi);
  
  //implicit flux:
  ldouble rho=pp[0];
  ldouble u=pp[1];  
  ldouble E=pp[6];  
  ldouble pr=(GAMMA-1.)*(u);
  ldouble T=pr*MU_GAS*M_PROTON/K_BOLTZ/rho;
  ldouble xx=get_x(ix,0);
  ldouble yy=get_x(iy,1);
  ldouble zz=get_x(iz,2);
  ldouble kappa=calc_kappa(rho,T,xx,yy,zz);
  ldouble chi=kappa+calc_kappaes(rho,T,xx,yy,zz);  
  ldouble B = SIGMA_RAD*pow(T,4.)/Pi;

  ldouble Fold[3]={pp[7],pp[8],pp[9]};
  ldouble Fnew[3];
  Fnew[0]=Fold[0]/(1.+dt*chi); 
  Fnew[1]=Fold[1]/(1.+dt*chi);
  Fnew[2]=Fold[2]/(1.+dt*chi);
 
  deltas[1]=Fnew[0]-Fold[0];
  deltas[2]=Fnew[1]-Fold[1];
  deltas[3]=Fnew[2]-Fold[2];

  //solving in parallel for E and u
  calc_LTE_ff(rho,&u,&E,dt,0); 

  deltas[0]=E-pp[6];

  return 0;

}

//**********************************************************************
//******* solves explicitly gas - radiation interaction  *******************
//******* in the lab frame, returns vector of deltas **********************
//**********************************************************************
int
solve_explicit_lab(int ix,int iy,int iz,ldouble dt,ldouble* deltas)
{
  int i1,i2,i3,iv;
  ldouble pp[NV];
  ldouble eup[4][4],elo[4][4];
  pick_T(emuup,ix,iy,iz,eup);
  pick_T(emulo,ix,iy,iz,elo);
  ldouble gg[4][5];
  pick_g(ix,iy,iz,gg);
  ldouble GG[4][5];
  pick_G(ix,iy,iz,GG);
  ldouble gdet=gg[3][4];

  for(iv=0;iv<NV;iv++)
    {
      pp[iv]=get_u(p,iv,ix,iy,iz);      
    }

  ldouble Gi[4];
  calc_Gi(pp,gg,GG,Gi);
  
  deltas[0]=-Gi[0]*dt;
  deltas[1]=-Gi[1]*dt;
  deltas[2]=-Gi[2]*dt;
  deltas[3]=-Gi[3]*dt;


  return 0;

}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//numerical LTE solver used only by solve_implicit_ff()
//TODO: to be replaced with something faster or abandoned
//**********************************************************************
//**********************************************************************
//**********************************************************************

struct calc_LTE_ff_parameters
{
  ldouble rho,E,u,kappa,dt;
  int verbose;
};

double
f_calc_LTE_ff (double u, void *params)
{
  struct calc_LTE_ff_parameters *p 
    = (struct calc_LTE_ff_parameters *) params;
  
  ldouble rho=p->rho;
  ldouble Eold=p->E;
  ldouble uold=p->u;
  ldouble kappa=p->kappa;
  ldouble dt=p->dt;
  int verbose=p->verbose;

  ldouble pr= (GAMMA-1.)*(ldouble)u;
  ldouble T = pr*MU_GAS*M_PROTON/K_BOLTZ/rho;
  ldouble B = SIGMA_RAD*pow(T,4.)/Pi;
  ldouble Enew=(Eold +4.*Pi*kappa*B*dt)/(1.+kappa*dt);

  return u - uold - dt*(-4.*Pi*kappa*B + kappa*Enew);
}
     
double
df_calc_LTE_ff (double u, void *params)
{  
  struct calc_LTE_ff_parameters *p 
    = (struct calc_LTE_ff_parameters *) params;

  ldouble rho=p->rho;
  ldouble Eold=p->E;
  ldouble uold=p->u;
  ldouble kappa=p->kappa;
  ldouble dt=p->dt;
  int verbose=p->verbose;

  ldouble pr= (GAMMA-1.)*(ldouble)u;
  ldouble T = pr*MU_GAS*M_PROTON/K_BOLTZ/rho;
  ldouble B = SIGMA_RAD*pow(T,4.)/Pi;
  ldouble Enew=(Eold +4.*Pi*kappa*B*dt)/(1.+kappa*dt);

  ldouble dBdu=4.*B/(ldouble)u;
  ldouble dEdu=4.*Pi*kappa*dt*dBdu/(1.+dt*kappa);

  ldouble dfdu=1.+dt*kappa*4.*Pi*dBdu-dt*kappa*dEdu;

  return (double)dfdu;
}
     
void
fdf_calc_LTE_ff (double u, void *params, 
	       double *y, double *dy)
{ 
   struct calc_LTE_ff_parameters *p 
    = (struct calc_LTE_ff_parameters *) params;

  ldouble rho=p->rho;
  ldouble Eold=p->E;
  ldouble uold=p->u;
  ldouble kappa=p->kappa;
  ldouble dt=p->dt;
  int verbose=p->verbose;

  ldouble pr= (GAMMA-1.)*(ldouble)u;
  ldouble T = pr*MU_GAS*M_PROTON/K_BOLTZ/rho;
  ldouble B = SIGMA_RAD*pow(T,4.)/Pi;
  ldouble Enew=(Eold +4.*Pi*kappa*B*dt)/(1.+kappa*dt);

  ldouble dBdu=4.*B/(ldouble)u;
  ldouble dEdu=4.*Pi*kappa*dt*dBdu/(1.+dt*kappa);

  ldouble dfdu=1.+dt*kappa*4.*Pi*dBdu-dt*kappa*dEdu;
  ldouble f=u - uold - dt*(-4.*Pi*kappa*B + kappa*Enew);
 
  *y = f;
  *dy = dfdu;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
int
calc_LTE_ff(ldouble rho,ldouble *uint, ldouble *E,ldouble dt, int verbose)
{
  struct calc_LTE_ff_parameters cltep;
  cltep.rho=rho;
  cltep.u=*uint;
  cltep.E=*E;
  cltep.verbose=verbose;
  
  if(cltep.E<EFLOOR && 0)
    {
      printf("imposing EFLOOR\n");
      cltep.E=EFLOOR;
    }
  
  ldouble p=(GAMMA-1.)*cltep.u;
  ldouble Tgas=p*MU_GAS*M_PROTON/K_BOLTZ/rho;
  cltep.kappa=calc_kappa(rho,Tgas,-1.,-1.,-1.);

  ldouble tlte = 1./cltep.kappa;

  cltep.dt=(double)dt;

  //solving for E
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fdfsolver_type *TT;
  gsl_root_fdfsolver *s;
  double x0, x = cltep.u;
  gsl_function_fdf FDF;
     
  FDF.f = &f_calc_LTE_ff;
  FDF.df = &df_calc_LTE_ff;
  FDF.fdf = &fdf_calc_LTE_ff;
  FDF.params = &cltep;
  
  TT = gsl_root_fdfsolver_newton;
  s = gsl_root_fdfsolver_alloc (TT);
  gsl_root_fdfsolver_set (s, &FDF, x);

  do
    {
      iter++;
      status = gsl_root_fdfsolver_iterate (s);
      x0 = x;
      x = gsl_root_fdfsolver_root (s);
      status = gsl_root_test_delta (x, x0, 0, U2PRADPREC);
    }
  while (status == GSL_CONTINUE && iter < max_iter);
     
  if(iter>=max_iter) 
    {
      printf("lte in: %e %e %e %e\n",rho,*uint,*E,dt);
      my_err("iter lte did not work\n");
      return -1;
    }

  gsl_root_fdfsolver_free (s);

  *uint=x;
  ldouble pp1 = (GAMMA-1.)*(x);
  ldouble Ttu = pp1*MU_GAS*M_PROTON/K_BOLTZ/cltep.rho;
  ldouble Bp1 = SIGMA_RAD*pow(Ttu,4.)/Pi;
  *E=(cltep.E+4.*Pi*cltep.kappa*Bp1*dt)/(1.+cltep.kappa*dt);
  
  if(*uint<0 && 0)
    {
      *uint=1.e-30;
      printf("imposing ufloor in lte\n");
      getch();
    }
  return 0;
}
//**********************************************************************
//end of the numerical LTE solver
//**********************************************************************


//**********************************************************************
//******* opacities ****************************************************
//**********************************************************************
//absorption
ldouble calc_kappa(ldouble rho, ldouble T,ldouble x,ldouble y,ldouble z)
{
#include PR_KAPPA //PROBLEMS/XXX/kappa.c
}

//scattering
ldouble calc_kappaes(ldouble rho, ldouble T,ldouble x,ldouble y,ldouble z)
{  
#include PR_KAPPAES //PROBLEMS/XXX/kappaes.c
}

//**********************************************************************
//****** takes radiative stress tensor and gas primitives **************
//****** and calculates contravariant four-force ***********************
//**********************************************************************
int
calc_Gi(ldouble *pp, ldouble gg[][5], ldouble GG[][5], ldouble Gi[4])
{
  int i,j,k;

  //radiative stress tensor in the lab frame
  ldouble Rij[4][4];
  calc_Rij(pp,gg,GG,Rij);

  //the four-velocity of fluid in lab frame
  ldouble ucon[4],ucov[4],vpr[3];
  ucon[1]=pp[2];
  ucon[2]=pp[3];
  ucon[3]=pp[4];
  conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);

  //covariant four-velocity
  indices_21(ucon,ucov,gg);  

  //gas properties
  ldouble rho=pp[0];
  ldouble u=pp[1];
  ldouble p= (GAMMA-1.)*(ldouble)u;
  ldouble Tgas=p*MU_GAS*M_PROTON/K_BOLTZ/rho;
  ldouble B = SIGMA_RAD*pow(Tgas,4.)/Pi;
  ldouble kappa=calc_kappa(rho,Tgas,-1.,-1.,-1.);
  ldouble kappaes=calc_kappaes(rho,Tgas,-1.,-1.,-1.);
  ldouble chi=kappa+kappaes;

  //contravariant four-force in the lab frame

  //R^ab u_a u_b
  ldouble Ruu=0.;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      Ruu+=Rij[i][j]*ucov[i]*ucov[j];

  ldouble Ru;
  for(i=0;i<4;i++)
    {
      Ru=0.;
      for(j=0;j<4;j++)
	Ru+=Rij[i][j]*ucov[j];
      Gi[i]=-chi*Ru - (kappaes*Ruu + kappa*4.*Pi*B)*ucon[i];
    }

  //  print_4vector(Gi);getchar();

  return 0;
}

//**********************************************************************
//******* takes fluid frame E,F^i in place of radiative primitives   *******
//******* and calculates force G^\mu in the fluid frame ****************
//**********************************************************************
int
calc_Gi_ff(ldouble *pp, ldouble Gi[4])
{
  ldouble rho=pp[0];
  ldouble u=pp[1];
  ldouble E=pp[6];
  ldouble F[3]={pp[7],pp[8],pp[9]};

  ldouble p= (GAMMA-1.)*(ldouble)u;
  ldouble T = p*MU_GAS*M_PROTON/K_BOLTZ/rho;
  ldouble B = SIGMA_RAD*pow(T,4.)/Pi;
  ldouble Tgas=p*MU_GAS*M_PROTON/K_BOLTZ/rho;
  ldouble kappa=calc_kappa(rho,Tgas,-1.,-1.,-1.);
  ldouble kappaes=calc_kappaes(rho,Tgas,-1.,-1.,-1.);
  ldouble chi=kappa+kappaes;

  Gi[0]=kappa*(E-4.*Pi*B);
  Gi[1]=chi*F[0];
  Gi[2]=chi*F[1];
  Gi[3]=chi*F[2];

  return 0;
}

//***********************************************************************************
//******* takes primitives and closes Rij in arbitrary frame ****************************
//***********************************************************************************
int
calc_Rij(ldouble *pp0, ldouble gg[][5], ldouble GG[][5], ldouble Rij[][4])
{
  ldouble pp[NV],Erf;
  int verbose=0;
  int i,j;
 //relative velocity
  ldouble urfcon[4];
  //covariant formulation

#ifdef LABRADFLUXES

#if(0)
  //artificially puts pp=uu and converts them to urf and Erf using the regular converter
  u2p_rad_urf(pp0,pp,gg,GG,&i);

#else //from labfluxes

  for(i=0;i<NV;i++)
    pp[i]=pp0[i];
  //R^0_mu
  ldouble A[4]={pp[6],pp[7],pp[8],pp[9]};

  //indices up
  indices_12(A,A,GG);

  //covariant formulation
  
  //g_munu R^0mu R^0nu
  ldouble gRR=gg[0][0]*A[0]*A[0]+gg[0][1]*A[0]*A[1]+gg[0][2]*A[0]*A[2]+gg[0][3]*A[0]*A[3]+
    gg[1][0]*A[1]*A[0]+gg[1][1]*A[1]*A[1]+gg[1][2]*A[1]*A[2]+gg[1][3]*A[1]*A[3]+
    gg[2][0]*A[2]*A[0]+gg[2][1]*A[2]*A[1]+gg[2][2]*A[2]*A[2]+gg[2][3]*A[2]*A[3]+
    gg[3][0]*A[3]*A[0]+gg[3][1]*A[3]*A[1]+gg[3][2]*A[3]*A[2]+gg[3][3]*A[3]*A[3];
 
  //the quadratic equation for u^t of the radiation rest frame (urf[0])
  ldouble a,b,c;
  a=16.*gRR;
  b=8.*(gRR*GG[0][0]+A[0]*A[0]);
  c=gRR*GG[0][0]*GG[0][0]-A[0]*A[0]*GG[0][0];
  ldouble delta=b*b-4.*a*c;
  urfcon[0]=sqrtl((-b-sqrtl(delta))/2./a);
  if(isnan(urfcon[0])) 
    {
      //      printf("err\n");
      //      print_4vector(Av);
      
      //TODO: gtph
      ldouble utaim=10.;
      ldouble Afac = sqrtl((-1.-utaim*utaim*gg[0][0])/(A[1]*A[1]*gg[1][1]+A[2]*A[2]*gg[2][2]+A[3]*A[3]*gg[3][3]));
      
      urfcon[0]=utaim;
      urfcon[1]=Afac*A[1];
      urfcon[2]=Afac*A[2];
      urfcon[3]=Afac*A[3];
    }

  //radiative energy density in the radiation rest frame
  Erf=3.*A[0]/(4.*urfcon[0]*urfcon[0]+GG[0][0]);

  //four-velocity of the rest frame
  urfcon[1]=3./(4.*Erf*urfcon[0])*(A[1]-1./3.*Erf*GG[0][1]);
  urfcon[2]=3./(4.*Erf*urfcon[0])*(A[2]-1./3.*Erf*GG[0][2]);
  urfcon[3]=3./(4.*Erf*urfcon[0])*(A[3]-1./3.*Erf*GG[0][3]);

  //lab frame:
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      Rij[i][j]=4./3.*Erf*urfcon[i]*urfcon[j]+1./3.*Erf*GG[i][j];

  return 0;

#endif



#else
  for(i=0;i<NV;i++)
    pp[i]=pp0[i];
#endif

  //radiative energy density in the radiation rest frame
  Erf=pp[6];
 

  urfcon[0]=0.;
  urfcon[1]=pp[7];
  urfcon[2]=pp[8];
  urfcon[3]=pp[9];
  //converting to lab four-velocity
  conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,gg,GG);
  //lab frame stress energy tensor:
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      Rij[i][j]=4./3.*Erf*urfcon[i]*urfcon[j]+1./3.*Erf*GG[i][j];

  //  if(pp[7]!=0.) {print_tensor(Rij); getchar();}

  return 0;
}

//**********************************************************************
//******* takes E and F^i from primitives (artificial) **********************
//******* and calculates radiation stress ******************************
//******* tensor R^ij in fluid frame using M1 closure scheme ***********
//**********************************************************************
int
calc_Rij_ff(ldouble *pp, ldouble Rij[][4])
{
  ldouble E=pp[6];
  ldouble F[3]={pp[7],pp[8],pp[9]};

  ldouble nx,ny,nz,nlen,f;

  nx=F[0]/E;
  ny=F[1]/E;
  nz=F[2]/E;

  nlen=sqrt(nx*nx+ny*ny+nz*nz);
 
#ifdef EDDINGTON_APR
  f=1./3.;
#else  
  if(nlen>=1.)
    {
      f=1.;
    }
  else //M1
    f=(3.+4.*(nx*nx+ny*ny+nz*nz))/(5.+2.*sqrt(4.-3.*(nx*nx+ny*ny+nz*nz)));  
#endif
  
  if(nlen>0) 
    {
      nx/=nlen;
      ny/=nlen;
      nz/=nlen;
    }
  else
    {
      ;
    }
 
  Rij[0][0]=E;
  Rij[0][1]=Rij[1][0]=F[0];
  Rij[0][2]=Rij[2][0]=F[1];
  Rij[0][3]=Rij[3][0]=F[2];

  Rij[1][1]=E*(.5*(1.-f) + .5*(3.*f - 1.)*nx*nx);
  Rij[1][2]=E*(.5*(3.*f - 1.)*nx*ny);
  Rij[1][3]=E*(.5*(3.*f - 1.)*nx*nz);

  Rij[2][1]=E*(.5*(3.*f - 1.)*ny*nx);
  Rij[2][2]=E*(.5*(1.-f) + .5*(3.*f - 1.)*ny*ny);
  Rij[2][3]=E*(.5*(3.*f - 1.)*ny*nz);

  Rij[3][1]=E*(.5*(3.*f - 1.)*nz*nx);
  Rij[3][2]=E*(.5*(3.*f - 1.)*nz*ny);
  Rij[3][3]=E*(.5*(1.-f) + .5*(3.*f - 1.)*nz*nz);

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//returns rad primitives for an atmosphere
int
set_radatmosphere(ldouble *pp,ldouble *xx,ldouble gg[][5],ldouble GG[][5],int atmtype)
{
#ifdef RADIATION  
  if(atmtype==0) //fixed Erf, urf of normal observer
    {
      pp[6]=ERADATMMIN; 
      ldouble ucon[4];
      calc_normalobs_4vel(GG,ucon);
      conv_vels(ucon,ucon,VEL4,VELPRIMRAD,gg,GG);
      pp[7]=ucon[1]; 
      pp[8]=ucon[2];
      pp[9]=ucon[3];
    }
  if(atmtype==1) //optically thin atmosphere, photon flying radially 
    {
      ldouble ucon[4];
      calc_photonrad_4vel(gg,GG,ucon);
      conv_vels(ucon,ucon,VEL4,VELPRIMRAD,gg,GG);

      pp[6]=ERADATMMIN; 
      pp[7]=ucon[1]; 
      pp[8]=ucon[2];
      pp[9]=ucon[3];

      //printf("%g > ",xx[1]);  print_4vector(&pp[6]); getchar();
    }
  if(atmtype==2) //optically thin atmosphere, scalings from numerical solution of radiall influx
    {
      ldouble gammamax=10.;
      //ldouble xxout[4]={0.,get_x(NX-1,0),0.,0.};
      //      coco_N(xxout,xxout,MYCOORDS,BLCOORDS);
      //      ldouble rout=xxout[1];

      ldouble rout=2.; //normalize at r_BL=2

      ldouble xxBL[4];
      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
      ldouble r=xxBL[1];
     
      pp[6]=ERADATMMIN*(rout/r)*(rout/r)*(rout/r)*(rout/r);

      ldouble ut[4]={0.,-gammamax*pow(r/rout,1.),0.,0.};

      ldouble ggBL[4][5],GGBL[4][5];
      calc_g_arb(xxBL,ggBL,KERRCOORDS);
      calc_G_arb(xxBL,GGBL,KERRCOORDS);

      conv_vels(ut,ut,VELR,VEL4,ggBL,GGBL);

      trans2_coco(xxBL,ut,ut,KERRCOORDS,MYCOORDS);

      conv_vels(ut,ut,VEL4,VELPRIM,gg,GG);
      
      pp[7]=ut[1];      
      pp[8]=ut[2];      
      pp[9]=ut[3];

    }
#endif
  return 0;
}

//**********************************************************************
//suplementary routines for conversions
//**********************************************************************
ldouble calc_PEQ_ufromTrho(ldouble T,ldouble rho)
{
  ldouble p=K_BOLTZ*rho*T/MU_GAS/M_PROTON;	      
  ldouble u=p/(GAMMA-1.);
  return u;
}

ldouble calc_PEQ_Tfromurho(ldouble u,ldouble rho)
{
  ldouble p=u*(GAMMA-1.);
  ldouble T=p/(K_BOLTZ*rho/MU_GAS/M_PROTON);
  return T;
}

ldouble calc_LTE_EfromT(ldouble T)
{
  return 4.*SIGMA_RAD*T*T*T*T;
}

ldouble calc_LTE_TfromE(ldouble E )
{
  return sqrt(sqrt((E/4./SIGMA_RAD)));
}


ldouble calc_LTE_Efromurho(ldouble u,ldouble rho)
{
  ldouble p=(GAMMA-1.)*(u);
  ldouble T=p*MU_GAS*M_PROTON/K_BOLTZ/rho;

  return calc_LTE_EfromT(T);
}

/************************************************************************/
/******* calculates wavespeeds in the lab frame takin 1/@3 in ************/
/******* radiative rest frame and boosting it to lab frame *****************/
/******* using the HARM algorithm **************************************/
/************************************************************************/
int
calc_rad_wavespeeds(ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble tautot[3],ldouble *aval,int verbose)
{
  int i,j;
  
  //metric
  ldouble g00=gg[0][0];
  ldouble g03=gg[0][3];
  ldouble g30=g03;
  ldouble g11=gg[1][1];
  ldouble g22=gg[2][2];
  ldouble g33=gg[3][3];

  //inversed metric
  ldouble G00=GG[0][0];
  ldouble G03=GG[0][3];
  ldouble G11=GG[1][1];
  ldouble G22=GG[2][2];
  ldouble G33=GG[3][3];
  ldouble G30=G03;

#ifdef LABRADFLUXES
  //artificially puts pp=uu and converts them to urf and Erf using the regular converter
  u2p_rad_urf(pp,pp,gg,GG,&i);
#endif
  
  //radiative energy density in the radiation rest frame
  ldouble Erf=pp[6];
  //relative four-velocity
  ldouble urfcon[4];
  urfcon[0]=0.;
  urfcon[1]=pp[7];
  urfcon[2]=pp[8];
  urfcon[3]=pp[9];

  //converting to lab four-velocity
  conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,gg,GG);

  //square of radiative wavespeed in radiative rest frame
  ldouble rv2rad = 1./3.;
  ldouble rv2,rv2tau;

  //**********************************************************************
  //algorithm from HARM to transform the fluid frame wavespeed into lab frame
  //**********************************************************************

  ldouble Acov[4],Acon[4],Bcov[4],Bcon[4],Asq,Bsq,Au,Bu,AB,Au2,Bu2,AuBu,A,B,discr,wspeed2;
  ldouble axl,axr,ayl,ayr,azl,azr;
  axl=axr=ayl=ayr=azl=azr=1.;
   
  //**********************************************************************
  //**********************************************************************
  int dim;
  for(dim=0;dim<3;dim++)
    {
      //characterisitic limiter based on the optical depth
      //TODO: validate against opt.thick tests
      if(tautot[dim]>0.) 
	{
	  rv2tau=4./3./tautot[dim]*4./3./tautot[dim];
	  rv2=my_min(rv2rad,rv2tau);		     
	}
      else
	rv2=rv2rad;
      
      Acov[0]=0.;
      Acov[1]=0.;
      Acov[2]=0.;
      Acov[3]=0.;
      Acov[dim+1]=1.;
      indices_12(Acov,Acon,GG);
  
      Bcov[0]=1.;
      Bcov[1]=0.;
      Bcov[2]=0.;
      Bcov[3]=0.;
      indices_12(Bcov,Bcon,GG);

      Asq = dot(Acon,Acov);
      Bsq = dot(Bcon,Bcov);
      Au = dot(Acov, urfcon);
      Bu = dot(Bcov, urfcon);
      AB = dot(Acon, Bcov);
      Au2 = Au * Au;
      Bu2 = Bu * Bu;
      AuBu = Au * Bu;

      wspeed2=rv2;
      B = 2. * (AuBu * (1.0-wspeed2)  - AB*wspeed2);
      A = Bu2 * (1.0 - wspeed2) - Bsq * wspeed2;
      discr = 4.0 * wspeed2 * ((AB * AB - Asq * Bsq) * wspeed2 + (2.0 * AB * Au * Bu - Asq * Bu2 - Bsq * Au2) * (wspeed2 - 1.0));
      if(discr<0.) {printf("x1discr in ravespeeds lt 0\n"); discr=0.;}
      discr = sqrt(discr);
      ldouble cst1 = -(-B + discr) / (2. * A);
      ldouble cst2 = -(-B - discr) / (2. * A);  

      axl = my_min(cst1,cst2);
      axr = my_max(cst1,cst2);

      aval[dim*2+0]=axl;
      aval[dim*2+1]=axr;
    }

  /*
  if(fabs(pp[7])>1.e-3)
    {
      print_Nvector(pp,NV);
      print_Nvector(&aval[0],6);
      getchar();
    }
  */

  return 0;
}


/************************************************************************/
/******* aplying radiative four velocity-related changes in conserved ******/
/************************************************************************/
int
apply_rad_source_del4(int ix,int iy,int iz,ldouble *del4)
{
  ldouble delapl[NV];

  delapl[0]=0.;
  delapl[1]=-del4[0];
  delapl[2]=-del4[1];
  delapl[3]=-del4[2];
  delapl[4]=-del4[3];
  delapl[5]=0.;
  delapl[6]=del4[0];
  delapl[7]=del4[1];
  delapl[8]=del4[2];
  delapl[9]=del4[3];

  int iv;
  for(iv=0;iv<NV;iv++)
    {
      set_u(u,iv,ix,iy,iz, get_u(u,iv,ix,iy,iz)+delapl[iv] );
    }

  return 0;
}


/************************************************************************/
/******* explicit radiative source term  ***********************************/
/************************************************************************/
int explicit_rad_source_term(int ix,int iy, int iz,ldouble dt, ldouble gg[][5], ldouble GG[][5])
{
  set_cflag(RADSOURCETYPEFLAG,ix,iy,iz,RADSOURCETYPEEXPLICIT);   

  ldouble del4[4],delapl[NV];
  int iv;

  //new primitives before the source operator
  calc_primitives(ix,iy,iz);
  //applied explicitly directly in lab frame
  solve_explicit_lab(ix,iy,iz,dt,del4);
  indices_21(del4,del4,gg);

  apply_rad_source_del4(ix,iy,iz,del4);

  set_cflag(RADSOURCEWORKEDFLAG,ix,iy,iz,0); 

  return 0;
}

/************************************************************************/
/******* implicit radiative source term in fluid frame and transported to lab  - backup method */
/************************************************************************/
int implicit_ff_rad_source_term(int ix,int iy, int iz,ldouble dt, ldouble gg[][5], ldouble GG[][5],ldouble tlo[][4], ldouble tup[][4],ldouble *pp)
{
  set_cflag(RADSOURCETYPEFLAG,ix,iy,iz,RADSOURCETYPEIMPLICITFF); 

  ldouble del4[4],delapl[NV];
  int iv;

  calc_primitives(ix,iy,iz);
  
  if(solve_implicit_ff(ix,iy,iz,dt,del4)<0) 
    {
      //failure, keeping u[] intact, reporting
      set_cflag(RADSOURCEWORKEDFLAG,ix,iy,iz,-1); 
      return -1;
    }

  trans2_on2cc(del4,del4,tlo);
  boost2_ff2lab(del4,del4,pp,gg,GG);
  indices_21(del4,del4,gg);

  apply_rad_source_del4(ix,iy,iz,del4);

  set_cflag(RADSOURCEWORKEDFLAG,ix,iy,iz,0); 

  return 0;
}

/************************************************************************/
/******* calculated Gi, compares with conserved and determines **********/
/******* if implicit step necessary ****************************************/
/************************************************************************/
int test_if_rad_implicit(int ix,int iy, int iz,ldouble dt, ldouble gg[][5], ldouble GG[][5], ldouble *del4)
{
  int iv;
  ldouble delapl[NV],uu[NV],pp[NV];

  //calculating radforce
  //new primitives
  calc_primitives(ix,iy,iz);
  //rad-for-force
  solve_explicit_lab(ix,iy,iz,dt,del4);
  //del4[] will be passed up
  indices_21(del4,del4,gg); 
  
#if(0)

  //basing on maximal tautot
  ldouble TAULIMIT=0.1;
  ldouble dx[3],xx[4];
  get_xx(ix,iy,iz,xx);
  dx[0]=get_size_x(ix,0)*sqrt(gg[1][1]);
  dx[1]=get_size_x(iy,1)*sqrt(gg[2][2]);
  dx[2]=get_size_x(iz,2)*sqrt(gg[3][3]);
  ldouble tautot[3],taumax;
  calc_tautot(pp,xx,dx,tautot);
  taumax=my_max(tautot[0],my_max(tautot[1],tautot[2]));
  if(taumax<TAULIMIT) 
    return 0; //can do explicit
  else
    return 1; //must do implicitxs

#else 

  //based on comparison of Gi[i] with uu[i]
  ldouble DULIMIT=0.1;

  //gettin' pp & uu
  for(iv=0;iv<NV;iv++)
    {
      uu[iv]=get_u(u,iv,ix,iy,iz);
      pp[iv]=get_u(p,iv,ix,iy,iz);
    }

  //changes to conserved
  delapl[0]=0.;
  delapl[1]=-del4[0];
  delapl[2]=-del4[1];
  delapl[3]=-del4[2];
  delapl[4]=-del4[3];
  delapl[5]=0.;
  delapl[6]=del4[0];
  delapl[7]=del4[1];
  delapl[8]=del4[2];
  delapl[9]=del4[3];
  
  //comparing with conserved to get the largest change
  ldouble maxdu=-1., uval;

  for(iv=1;iv<NV;iv++)
    {
 
      if(iv==5) continue; //skip entropy
      if(iv==3 && NY==1) continue; //skip y-momentum
      if(iv==8 && NY==1) continue; //skip y-momentum
      if(iv==4 && NZ==1) continue; //skip y-momentum
      if(iv==9 && NZ==1) continue; //skip y-momentum

      uval=uu[iv];

      if(fabs(uval)<SMALL) //to avoid dividing by 0
	maxdu=my_max(maxdu,1./SMALL);
      else
	maxdu=my_max(maxdu,delapl[iv]/uval);
    }

  //debug
  //printf("%d %d %d: if implicit = %e | %d\n",ix,iy,iz,maxdu,maxdu < DULIMIT ? 0 : 1); getchar();

  if(maxdu<0.1)
    return 0; //can do explicit
  else
    return 1; //must do implicit

#endif

}
  

/************************************************************************/
/******* explicit radiative source term with sub-stepping - inefficient *****/
/************************************************************************/
int explicit_substep_rad_source_term(int ix,int iy, int iz,ldouble dt, ldouble gg[][5], ldouble GG[][5])
{
  set_cflag(RADSOURCETYPEFLAG,ix,iy,iz,RADSOURCETYPEEXPLICITSUBSTEP); 

  int iv;
  ldouble del4[4],delapl[NV],uu[NV],pp[NV];
  double fdt, fdta, maxfu=-1., fu, uval,futau;
  //calculating reference time step basing on maximal tautot
  ldouble dx[3],xx[4];
  get_xx(ix,iy,iz,xx);
  dx[0]=get_size_x(ix,0)*sqrt(gg[1][1]);
  dx[1]=get_size_x(iy,1)*sqrt(gg[2][2]);
  dx[2]=get_size_x(iz,2)*sqrt(gg[3][3]);
  ldouble tautot[3],taumax;
  calc_tautot(pp,xx,dx,tautot);
  taumax=my_max(tautot[0],my_max(tautot[1],tautot[2]));
  futau=1./taumax;
  //reference time step only approximate
  fdt=fdta=0.;
  do
    {
      //new primitives
      //TODO: check here and there if worked
      calc_primitives(ix,iy,iz);
      //conserved
      for(iv=0;iv<NV;iv++)
	{
	  uu[iv]=get_u(u,iv,ix,iy,iz);
	  pp[iv]=get_u(p,iv,ix,iy,iz);
	}

      //vector of changes of conserved assuming original dt which only multiplies source terms
      solve_explicit_lab(ix,iy,iz,dt,del4);
      indices_21(del4,del4,gg);
      //changes to conserved
      delapl[0]=0.;
      delapl[1]=-del4[0];
      delapl[2]=-del4[1];
      delapl[3]=-del4[2];
      delapl[4]=-del4[3];
      delapl[5]=0.;
      delapl[6]=del4[0];
      delapl[7]=del4[1];
      delapl[8]=del4[2];
      delapl[9]=del4[3];
		
#if(0) //my old dtsub esitmation based on single dimension
      //comparing with conserved to get the largest change
      maxfu=-1.;

      //fluxes can be zero and their relative change large
      //so far considering only energy densities
      for(iv=1;iv<NV;iv++)
	{
	  if(iv==5) continue;
	  uval=get_u(u,iv,ix,iy,iz);
	  if(fabs(uval)<SMALL)  //to avoid dividing by 0
	    fu=futau;
	  else
	    fu=fabs(delapl[iv]/uval);

	  if(ix==NX/2 ) printf("> %d %e %e %e\n",iv,fu,uval,delapl[iv]);
	  if(fu>maxfu) maxfu=fu;
	}
		  
      if(maxfu<MAXEXPLICITSUBSTEPCHANGE)
	fdt=1.;
      else
	fdt=MAXEXPLICITSUBSTEPCHANGE/maxfu;

      if(ix==NX/2 ) printf("----\n%e\n",maxfu);
		  
#else //Jon's spacetime
      //substep
      ldouble Umhd,Urad,Gtot,iUmhd,iUrad,idtsub,dtsub;
      Umhd=Urad=Gtot=0.;
      for(iv=0;iv<4;iv++)
	{
	  Umhd+=uu[1+iv]*uu[1+iv]*GG[iv][iv]; //GG?
	  Urad+=uu[6+iv]*uu[6+iv]*GG[iv][iv]; //GG?
	  Gtot+=del4[iv]*del4[iv]*GG[iv][iv]; //GG?
	}

      iUmhd=1.0/(fabs(Umhd)+SMALL);
      iUrad=1.0/(fabs(Urad)+SMALL);
      idtsub=SMALL+fabs(Gtot*my_max(iUmhd,iUrad));
      dtsub=1./idtsub;
      /*
	if(ix==42 && pp[7]!=0 && 0)
	{
	print_Nvector(pp,NV);
	print_4vector(del4);
	printf("----\n%e %e %e > %e\n",Gtot,Umhd,Urad,dtsub);getchar();
	}
      */
      if(dtsub>1.)
	fdt=1.;
      else
	fdt=dtsub;

      if(fdt<1.e-6 && 0) 
	{
	  printf("----\n%e %e %e > %e\n",Gtot,Umhd,Urad,dtsub);
	  printf("%d %d %d %f %f\n",ix,iy,iz,fdt,fdta);
	  getchar();
	}
#endif		 
		  
      if(fdta+fdt>1.) fdt=1.-fdta;
	
      for(iv=0;iv<NV;iv++)
	{
	  uu[iv]+=delapl[iv]*fdt;
	  set_u(u,iv,ix,iy,iz, uu[iv]);
	}

      fdta+=fdt;

    }
  while(fdta<1.);
	    
  set_cflag(RADSOURCEWORKEDFLAG,ix,iy,iz,0); 

  return 0;
}

/************************************************************************/
/******* implicit radiative source term in lab frame *********************/
/******* solves numerically in 4D **************************************/
/************************************************************************/
int implicit_lab_rad_source_term(int ix,int iy, int iz,ldouble dt, ldouble gg[][5], ldouble GG[][5],ldouble tlo[][4], ldouble tup[][4],ldouble *pp)
{
  ldouble del4[4],delapl[NV];
  int iv;
  
  set_cflag(RADSOURCETYPEFLAG,ix,iy,iz,RADSOURCETYPEIMPLICITLAB); 

  if(solve_implicit_lab(ix,iy,iz,dt,del4)<0)
    {
      //numerical implicit in 4D did not work
      //use the explicit-implicit backup method
      if(implicit_ff_rad_source_term(ix,iy,iz,dt,gg,GG,tlo,tup,pp)<0)
	{
	  //this one failed too - failure
	  set_cflag(RADSOURCEWORKEDFLAG,ix,iy,iz,-1); 
	  return -1;	  
	}
    }
  else
    {
      //success in lab frame
      apply_rad_source_del4(ix,iy,iz,del4);
    }

  set_cflag(RADSOURCEWORKEDFLAG,ix,iy,iz,0); 

  return 0;
}
