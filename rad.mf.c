//KORAL - rad.c
//radiation-related routines adjusted for multi rad fluids

#include "ko.h"

//***********************************************************************************
//******* takes primitives and closes Rij in arbitrary frame for all fluids *********
//***********************************************************************************
int
calc_Rij_mf(ldouble *pp0, ldouble gg[][5], ldouble GG[][5], ldouble Rij[][4][4])
{
#ifdef MULTIRADFLUID
  ldouble pp[NV],Erf;
  int verbose=0;
  int i,j,irf;
 //relative velocity
  ldouble urfcon[4];
  //covariant formulation

#ifdef LABRADFLUXES
  my_err("LABRADFLUXES not implemented in calc_Rij_mf\n");
#else
  for(i=0;i<NV;i++)
    pp[i]=pp0[i];
#endif


  for(irf=0;irf<NRF;irf++)
    {
      //radiative energy density in the radiation rest frame
      Erf=pp[EE(irf)];

      urfcon[0]=0.;
      urfcon[1]=pp[FX(irf)];
      urfcon[2]=pp[FY(irf)];
      urfcon[3]=pp[FZ(irf)];
      //converting to lab four-velocity
      conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,gg,GG);
      //lab frame stress energy tensor:
      for(i=0;i<4;i++)
	for(j=0;j<4;j++)
	  Rij[irf][i][j]=4./3.*Erf*urfcon[i]*urfcon[j]+1./3.*Erf*GG[i][j];

    }

  

  return 0;
#endif
}

/************************************************************************/
/******* calculates wavespeeds in the lab frame takin 1/@3 in ************/
/******* radiative rest frame and boosting it to lab frame ***************/
/******* using the HARM algorithm for all fluids ***********************/
/************************************************************************/
/************************************************************************/
/******* currently returns one wavespeed for all fluids******************/
/************************************************************************/
int
calc_rad_wavespeeds_mf(ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble tautot[3],ldouble *aval,int verbose)
{
#ifdef MULTIRADFLUID
  int i,j,irf;
  
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
  my_err("LABRADFLUXES not implemented for MULTIRADFLUID\n");
#endif

  for(i=0;i<6;i++)
    aval[i]=0.;
  
  for(irf=0;irf<NRF;irf++)
    {
      //radiative energy density in the radiation rest frame
      ldouble Erf=pp[EE(irf)];
      //relative four-velocity
      ldouble urfcon[4];
      urfcon[0]=0.;
      urfcon[1]=pp[FX(irf)];
      urfcon[2]=pp[FY(irf)];
      urfcon[3]=pp[FZ(irf)];

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

	  aval[dim*2+0]=my_min(axl,aval[dim*2+0]);
	  aval[dim*2+1]=my_max(axr,aval[dim*2+1]);
	}
    }
 

  return 0;
#endif
}


//**********************************************************************
//******* takes E and F^i from primitives (artificial) **********************
//******* and calculates radiation stress ******************************
//******* tensor R^ij in fluid frame using M1 closure scheme ***********
//**********************************************************************
int
calc_Rij_ff_mf(ldouble *pp, ldouble Rij[][4][4])
{
#ifdef MULTIRADFLUID
  int irf;

  for(irf=0;irf<NRF;irf++)
    {

      ldouble E=pp[EE(irf)];
      ldouble F[3]={pp[FX(irf)],pp[FY(irf)],pp[FZ(irf)]};

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
 
      Rij[irf][0][0]=E;
      Rij[irf][0][1]=Rij[irf][1][0]=F[0];
      Rij[irf][0][2]=Rij[irf][2][0]=F[1];
      Rij[irf][0][3]=Rij[irf][3][0]=F[2];

      Rij[irf][1][1]=E*(.5*(1.-f) + .5*(3.*f - 1.)*nx*nx);
      Rij[irf][1][2]=E*(.5*(3.*f - 1.)*nx*ny);
      Rij[irf][1][3]=E*(.5*(3.*f - 1.)*nx*nz);

      Rij[irf][2][1]=E*(.5*(3.*f - 1.)*ny*nx);
      Rij[irf][2][2]=E*(.5*(1.-f) + .5*(3.*f - 1.)*ny*ny);
      Rij[irf][2][3]=E*(.5*(3.*f - 1.)*ny*nz);

      Rij[irf][3][1]=E*(.5*(3.*f - 1.)*nz*nx);
      Rij[irf][3][2]=E*(.5*(3.*f - 1.)*nz*ny);
      Rij[irf][3][3]=E*(.5*(1.-f) + .5*(3.*f - 1.)*nz*nz);

    }

 

  return 0;
#endif
}

