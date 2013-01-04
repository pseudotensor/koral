//KORAL - rad.c
//radiation-related routines

#include "ko.h"

//**********************************************************************
//******* calculates total opacity as in the fluid frame ***************
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
//******* calculates total opacity as in the fluid frame ***************
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
//******* solves implicitidly four-force source terms ******************
//******* in the lab frame, returns ultimate deltas ********************
//******* the fiducial approach ****************************************
//**********************************************************************

//TODO: clean arguments
int f_implicit_lab(ldouble *uu0,ldouble *uu,ldouble *pp,ldouble dt,ldouble gg[][5], ldouble GG[][5],ldouble *f)
{
  ldouble Rij[4][4];
  ldouble ppp[NV];

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
  u2p(uu,pp2,gg);

  //radiative four-force
  ldouble Gi[4];
  //covariant calculation
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
  printf ("iter = %3d x = % .3Le % .3Le % .3Le % .3Le "
	  "f(x) = % .3Le % .3Le % .3Le % .3Le\n",
	  iter,
	  x[0],x[1]/x[0],x[2]/x[0],x[3]/x[0],f[0],f[1],f[2],f[3]);
}

int
solve_implicit_lab(int ix,int iy,int iz,ldouble dt,ldouble* deltas)
{
  int i1,i2,i3,iv,i,j;
  ldouble J[4][4],iJ[4][4];
  ldouble pp[NV],uu[NV],uu0[NV],uup[NV]; 
  ldouble f1[4],f2[4],f3[4],x[4];
  ldouble gg[4][5];
  ldouble GG[4][5];
  ldouble tup[4][4],tlo[4][4];
  pick_T(tmuup,ix,iy,iz,tup);
  pick_T(tmulo,ix,iy,iz,tlo);
  pick_g(ix,iy,iz,gg);
  pick_G(ix,iy,iz,GG);

  for(iv=0;iv<NV;iv++)
    {
      pp[iv]=get_u(p,iv,ix,iy,iz);      
      uu[iv]=get_u(u,iv,ix,iy,iz);  
      uu0[iv]=uu[iv];
    }

  ldouble EPS = 1.e-6;
  ldouble CONV = 1.e-6 ;

  int verbose=0;
  int iter=0;

  do
    {
      iter++;
      
      for(i=0;i<NV;i++)
	{
	  uup[i]=uu[i];
	}

      //values at zero state
      f_implicit_lab(uu0,uu,pp,dt,gg,GG,f1);
 
      //calculating approximate Jacobian
      for(i=0;i<4;i++)
	{
	  for(j=0;j<4;j++)
	    {
	      ldouble del;
	      if(uup[j+6]==0.) del=EPS*uup[6];
	      else del=EPS*uup[j+6];
	      uu[j+6]=uup[j+6]-del;

	      f_implicit_lab(uu0,uu,pp,dt,gg,GG,f2);
     
	      J[i][j]=(f2[i] - f1[i])/(uu[j+6]-uup[j+6]);

	      uu[j+6]=uup[j+6];
	    }
	}

      //inversion
      inverse_44matrix(J,iJ);

      //updating x
      for(i=0;i<4;i++)
	{
	  x[i]=uup[i+6];
	}

      for(i=0;i<4;i++)
	{
	  for(j=0;j<4;j++)
	    {
	      x[i]-=iJ[i][j]*f1[j];
	    }
	}

      if(verbose>0)    print_state_implicit_lab (iter,x,f1); 

      for(i=0;i<4;i++)
	{
	  uu[i+6]=x[i];
	}
  
      //test convergence
      for(i=0;i<4;i++)
	{
	  f3[i]=(uu[i+6]-uup[i+6]);
	  f3[i]=fabs(f3[i]/uup[6]);
	}

      if(f3[0]<CONV && f3[1]<CONV && f3[2]<CONV && f3[3]<CONV)
	break;

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
//******* solves explicitly gas - radiation interaction  *************
//******* in the fluid frame, returns vector of deltas******************
//**********************************************************************
int
solve_explicit_ff(int ix,int iy,int iz,ldouble dt,ldouble* deltas)
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
//TODO: to be replaced with something faster
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
      printf("lte in: %Le %Le %Le %Le\n",rho,*uint,*E,dt);
      my_err("iter lte did not work\n");
      return -1;
    }

  gsl_root_fdfsolver_free (s);

  *uint=x;
  ldouble pp1 = (GAMMA-1.)*(x);
  ldouble Ttu = pp1*MU_GAS*M_PROTON/K_BOLTZ/cltep.rho;
  ldouble Bp1 = SIGMA_RAD*powl(Ttu,4.)/Pi;
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

  ldouble vr=pp[2];
  ldouble vth=pp[3];
  ldouble vph=pp[4];

  ldouble gtt=gg[0][0];
  ldouble gtph=gg[0][3];
  ldouble grr=gg[1][1];
  ldouble gthth=gg[2][2];
  ldouble gphph=gg[3][3];
  ldouble gtr=gg[1][2];
  ldouble grph=gg[1][3];

  ldouble ut2=-1./(gtt + 2.*vph*gtph + 2.*vr*gtr + 2.*vr*vph*grph + vr*vr*grr + vth*vth*gthth + vph*vph*gphph );
  if(ut2<0.) 
    {
      ut2=0.;
    }
  ldouble ut=sqrtl(ut2);

  ucon[0]=ut;
  ucon[1]=vr*ut;
  ucon[2]=vth*ut;
  ucon[3]=vph*ut;
  indices_21(ucon,ucov,gg);  

  //gas properties
  ldouble rho=pp[0];
  ldouble u=pp[1];
 
  ldouble p= (GAMMA-1.)*(ldouble)u;
  ldouble T = p*MU_GAS*M_PROTON/K_BOLTZ/rho;
  ldouble B = SIGMA_RAD*pow(T,4.)/Pi;
  ldouble Tgas=p*MU_GAS*M_PROTON/K_BOLTZ/rho;
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

  return 0;
}

//**********************************************************************
//******* takes fluid frame E,F^i in place of radiative primitives   ***
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
//******* takes primitives (=conserved) and closes Rij in arbitrary frame ************
//***********************************************************************************
int
calc_Rij(ldouble *pp, ldouble gg[][5], ldouble GG[][5], ldouble Rij[][4])
{
  int verbose=0;
 
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
  ldouble urf[4],Erf;
  urf[0]=sqrtl((-b-sqrtl(delta))/2./a);
  if(isnan(urf[0])) urf[0]=1.;

  //radiative energy density in the radiation rest frame
  Erf=3.*A[0]/(4.*urf[0]*urf[0]+GG[0][0]);

  //four-velocity of the rest frame
  urf[1]=3./(4.*Erf*urf[0])*(A[1]-1./3.*Erf*GG[0][1]);
  urf[2]=3./(4.*Erf*urf[0])*(A[2]-1./3.*Erf*GG[0][2]);
  urf[3]=3./(4.*Erf*urf[0])*(A[3]-1./3.*Erf*GG[0][3]);

  //lab frame:
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      Rij[i][j]=4./3.*Erf*urf[i]*urf[j]+1./3.*Erf*GG[i][j];

  return 0;
}

//**********************************************************************
//******* takes E and F^i from primitives (which are assumed to replace ***
//******* the original primitives R^t_mu in pp) ************************
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

  nlen=sqrtl(nx*nx+ny*ny+nz*nz);
  
 
#ifdef EDDINGTON_APR
  f=1./3.;
#else  
  if(nlen>=1.)
    {
      f=1.;
    }
  else //M1
    f=(3.+4.*(nx*nx+ny*ny+nz*nz))/(5.+2.*sqrtl(4.-3.*(nx*nx+ny*ny+nz*nz)));  
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
  return sqrtl(sqrtl((E/4./SIGMA_RAD)));
}


ldouble calc_LTE_Efromurho(ldouble u,ldouble rho)
{
  ldouble p=(GAMMA-1.)*(u);
  ldouble T=p*MU_GAS*M_PROTON/K_BOLTZ/rho;

  return calc_LTE_EfromT(T);
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//return maximal wavespeeds eigenvalue for the radiation part
int max_eigen_lr_radx(ldouble nx,ldouble ny,ldouble nz,int idim,ldouble *al, ldouble *ar)
{
  ldouble j21,j22,j23,j24,j31,j32,j33,j34,j41,j42,j43,j44;
  ldouble a,b,c,d,e;
  ldouble cccbrtb,cccrt,x1,x2,x3,x4;
  gsl_complex z1,z2,z3,z4;
  ldouble e1,e2,e3,e4;

  ldouble nlen=sqrtl(nx*nx+ny*ny+nz*nz);
  if(nlen>1.)
    {
      nx/=nlen;
      ny/=nlen;
      nz/=nlen;
    }
  
  if(fabs(nx)<1.e-20 && fabs(ny)<1.e-20 && fabs(nz)<1.e-20)
    {
      *al=-sqrtl(1./3.);
      *ar=sqrtl(1./3.);
      return 0;
    }
  
  if(idim==0)
    {
      j21=(0. + 24.*Power(nx,2) - 68.*Power(nx,4) + 28.*Power(ny,2) - 64.*Power(nx,2)*Power(ny,2) + 4.*Power(ny,4) + 28.*Power(nz,2) - 64.*Power(nx,2)*Power(nz,2) + 8.*Power(ny,2)*Power(nz,2) + 4.*Power(nz,4) + 15.*Power(nx,2)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 20.*Power(nx,4)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 13.*Power(ny,2)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 16.*Power(nx,2)*Power(ny,2)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 4.*Power(ny,4)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 13.*Power(nz,2)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 16.*Power(nx,2)*Power(nz,2)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 8.*Power(ny,2)*Power(nz,2)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 4.*Power(nz,4)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2));
      j22=nx*((-41. + 12.*Power(nx,2) + 12.*Power(ny,2) + 12.*Power(nz,2) - 20.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) + (0.5*Power(nx,2)*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - (1.*Power(nx,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2) + (1.*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/(Power(nx,2) + Power(ny,2) + Power(nz,2)));
      j23=ny*((-41. + 12.*Power(nx,2) + 12.*Power(ny,2) + 12.*Power(nz,2) - 20.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) + (0.5*Power(nx,2)*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - (1.*Power(nx,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2));
      j24=nz*((-41. + 12.*Power(nx,2) + 12.*Power(ny,2) + 12.*Power(nz,2) - 20.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) + (0.5*Power(nx,2)*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - (1.*Power(nx,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2));
      j31=(0. - 4.*nx*ny - 72.*Power(nx,3)*ny - 72.*nx*Power(ny,3) - 72.*nx*ny*Power(nz,2) + 2.*nx*ny*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*Power(nx,3)*ny*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*nx*Power(ny,3)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*nx*ny*Power(nz,2)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2));
      j32=(ny*((0.5*Power(nx,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(nx,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
      j33=(nx*((0.5*Power(ny,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(ny,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
      j34=(nx*ny*((0.5*nz*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*nz*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
      j41=(0. - 4.*nx*nz - 72.*Power(nx,3)*nz - 72.*nx*Power(ny,2)*nz - 72.*nx*Power(nz,3) + 2.*nx*nz*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*Power(nx,3)*nz*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*nx*Power(ny,2)*nz*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*nx*Power(nz,3)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2));
      j42=(nz*((0.5*Power(nx,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(nx,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
      j43=(nx*nz*((0.5*ny*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*ny*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
      j44=(nx*((0.5*Power(nz,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(nz,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);

      a=1;
      b=-j22 - j33 - j44;
      c=-j21 - j23*j32 + j22*j33 - j24*j42 - j34*j43 + j22*j44 + j33*j44;
      d=-(j24*j41) + j24*j33*j42 - j24*j32*j43 + j22*j34*j43 - j22*j33*j44 + j21*(j33 + j44) - j23*(j31 + j34*j42 - j32*j44);
      e=j24*j33*j41 - j23*j34*j41 - j24*j31*j43 + j21*j34*j43 + j23*j31*j44 - j21*j33*j44;
    }   

  if(idim==1)
    {
      j21=(0. - 4.*nx*ny - 72.*Power(nx,3)*ny - 72.*nx*Power(ny,3) - 72.*nx*ny*Power(nz,2) + 2.*nx*ny*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*Power(nx,3)*ny*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*nx*Power(ny,3)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*nx*ny*Power(nz,2)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2));
      j22=(ny*((0.5*Power(nx,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(nx,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
      j23=(nx*((0.5*Power(ny,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(ny,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
      j24=(nx*ny*((0.5*nz*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*nz*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
      j31=(0. + 28.*Power(nx,2) + 4.*Power(nx,4) + 24.*Power(ny,2) - 64.*Power(nx,2)*Power(ny,2) - 68.*Power(ny,4) + 28.*Power(nz,2) + 8.*Power(nx,2)*Power(nz,2) - 64.*Power(ny,2)*Power(nz,2) + 4.*Power(nz,4) + 13.*Power(nx,2)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 4.*Power(nx,4)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 15.*Power(ny,2)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 16.*Power(nx,2)*Power(ny,2)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 20.*Power(ny,4)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 13.*Power(nz,2)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 8.*Power(nx,2)*Power(nz,2)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 16.*Power(ny,2)*Power(nz,2)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 4.*Power(nz,4)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2));
      j32=nx*((-41. + 12.*Power(nx,2) + 12.*Power(ny,2) + 12.*Power(nz,2) - 20.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) + (0.5*Power(ny,2)*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - (1.*Power(ny,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2));
      j33=ny*((-41. + 12.*Power(nx,2) + 12.*Power(ny,2) + 12.*Power(nz,2) - 20.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) + (0.5*Power(ny,2)*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - (1.*Power(ny,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2) + (1.*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/(Power(nx,2) + Power(ny,2) + Power(nz,2)));
      j34=nz*((-41. + 12.*Power(nx,2) + 12.*Power(ny,2) + 12.*Power(nz,2) - 20.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) + (0.5*Power(ny,2)*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - (1.*Power(ny,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2));
      j41=(0. - 4.*ny*nz - 72.*Power(nx,2)*ny*nz - 72.*Power(ny,3)*nz - 72.*ny*Power(nz,3) + 2.*ny*nz*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*Power(nx,2)*ny*nz*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*Power(ny,3)*nz*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*ny*Power(nz,3)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2));
      j42=(ny*nz*((0.5*nx*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*nx*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
      j43=(nz*((0.5*Power(ny,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(ny,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
      j44=(ny*((0.5*Power(nz,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(nz,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);

      a=1;
      b=-j22 - j33 - j44;
      c=-j31 - j23*j32 + j22*j33 - j24*j42 - j34*j43 + j22*j44 + j33*j44;
      d=-(j21*j32) - j34*j41 + j24*j33*j42 - j23*j34*j42 - j24*j32*j43 + j31*j44 + j23*j32*j44 + j22*(j31 + j34*j43 - j33*j44);
      e=-(j24*j32*j41) + j22*j34*j41 + j24*j31*j42 - j21*j34*j42 - j22*j31*j44 + j21*j32*j44;

    }

  if(idim==2)
    {
      j21=(0. - 4.*nx*nz - 72.*Power(nx,3)*nz - 72.*nx*Power(ny,2)*nz - 72.*nx*Power(nz,3) + 2.*nx*nz*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*Power(nx,3)*nz*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*nx*Power(ny,2)*nz*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*nx*Power(nz,3)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2));
      j22=(nz*((0.5*Power(nx,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(nx,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
      j23=(nx*nz*((0.5*ny*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*ny*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
      j24=(nx*((0.5*Power(nz,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(nz,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
      j31=(0. - 4.*ny*nz - 72.*Power(nx,2)*ny*nz - 72.*Power(ny,3)*nz - 72.*ny*Power(nz,3) + 2.*ny*nz*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*Power(nx,2)*ny*nz*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*Power(ny,3)*nz*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*ny*Power(nz,3)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2));
      j32=(ny*nz*((0.5*nx*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*nx*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
      j33=(nz*((0.5*Power(ny,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(ny,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
      j34=(ny*((0.5*Power(nz,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(nz,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
      j41=(0. + 28.*Power(nx,2) + 4.*Power(nx,4) + 28.*Power(ny,2) + 8.*Power(nx,2)*Power(ny,2) + 4.*Power(ny,4) + 24.*Power(nz,2) - 64.*Power(nx,2)*Power(nz,2) - 64.*Power(ny,2)*Power(nz,2) - 68.*Power(nz,4) + 13.*Power(nx,2)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 4.*Power(nx,4)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 13.*Power(ny,2)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 8.*Power(nx,2)*Power(ny,2)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 4.*Power(ny,4)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 15.*Power(nz,2)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 16.*Power(nx,2)*Power(nz,2)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 16.*Power(ny,2)*Power(nz,2)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 20.*Power(nz,4)*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2));
      j42=nx*((-41. + 12.*Power(nx,2) + 12.*Power(ny,2) + 12.*Power(nz,2) - 20.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) + (0.5*Power(nz,2)*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - (1.*Power(nz,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2));
      j43=ny*((-41. + 12.*Power(nx,2) + 12.*Power(ny,2) + 12.*Power(nz,2) - 20.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) + (0.5*Power(nz,2)*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - (1.*Power(nz,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2));
      j44=nz*((-41. + 12.*Power(nx,2) + 12.*Power(ny,2) + 12.*Power(nz,2) - 20.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) + (0.5*Power(nz,2)*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - (1.*Power(nz,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2) + (1.*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrtl(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/(Power(nx,2) + Power(ny,2) + Power(nz,2)));

      a=1;
      b=-j22 - j33 - j44;
      c=-(j23*j32) - j41 - j24*j42 - j34*j43 + j33*j44 + j22*(j33 + j44);
      d=-(j21*j42) - j23*j34*j42 + j33*(j41 + j24*j42) - j31*j43 - j24*j32*j43 + j23*j32*j44 + j22*(j41 + j34*j43 - j33*j44);
      e=j23*j32*j41 - j22*j33*j41 - j23*j31*j42 + j21*j33*j42 + j22*j31*j43 - j21*j32*j43;
    }

  //attemp to solve analytically:  
  int ret=gsl_poly_complex_solve_quartic(b/a,c/a,d/a,e/a,&z1,&z2,&z3,&z4);
  
  //if didn't work - solve numerically
  if(isnan(GSL_REAL(z1)) || isnan(GSL_REAL(z1)) || isnan(GSL_REAL(z1)) || isnan(GSL_REAL(z1)))
    {
            
      double coef[5] = { (double)e,(double)d,(double)c,(double)b,(double)a };  
      double z[8];      
      gsl_poly_complex_workspace * w 
	= gsl_poly_complex_workspace_alloc (5);      
      int result=gsl_poly_complex_solve (coef, 5, w, z);      
      gsl_poly_complex_workspace_free (w);
      if(result!=GSL_SUCCESS)
	{
	  printf("padaka w rozwiazywaniu czwormianu\n");
	  printf("%10Lf %10Le %10Le %10Le %e %e %e %e\n",sqrtl(nx*nx+ny*ny+nz*nz),nx,ny,nz,z[0],z[2],z[4],z[6]);
	  getchar();
	}

      e1=z[0];
      e2=z[2];
      e3=z[4];
      e4=z[6];
    }
  else
    {
      e1=GSL_REAL(z1);
      e2=GSL_REAL(z2);
      e3=GSL_REAL(z3);
      e4=GSL_REAL(z4);
    }
  
  *al=0.;
  if(e1<*al) *al=e1;
  if(e2<*al) *al=e2;
  if(e3<*al) *al=e3;
  if(e4<*al) *al=e4;
  
  *ar=0.;
  if(e1>*ar) *ar=e1;
  if(e2>*ar) *ar=e2;
  if(e3>*ar) *ar=e3;
  if(e4>*ar) *ar=e4;

  return 0;
}

