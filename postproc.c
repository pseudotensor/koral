//KORAL - misc.c
//routines for postprocessing 
//used both on the go and separately

#include "ko.h"

/*********************************************/
/* calculates radial profiles - L(r) etc. */
/*********************************************/
int calc_radialprofiles(ldouble profiles[][NX])
{
  //adjust NRADPROFILES in problem.h

  int ix,iy,iz,iv;
  ldouble xx[4],xxBL[4],dx[3],mdot,rho,ucon[4],utcon[4],ucon3[4];
  ldouble ucov[4],pp[NV],gg[4][5],GG[4][5],ggBL[4][5],GGBL[4][5];
  ldouble tautot[3],tauabs[3];

  //search for appropriate radial index
  for(ix=0;ix<NX;ix++)
    {
      //vertically integrated/averaged profiles

      for(iv=0;iv<NRADPROFILES;iv++)
	profiles[iv][ix]=0.;

      if(NZ==1) //phi-symmetry
	{
	  iz=0;
	  for(iy=0;iy<NY;iy++)
	    {
	      //calc_primitives_local(ix,iy,iz,pp);
	      //use old ones instead
	      for(iv=0;iv<NV;iv++)
		{
		  pp[iv]=get_u(p,iv,ix,iy,iz);
		}

	      get_xx(ix,iy,iz,xx);
	      dx[0]=get_size_x(ix,0);
	      dx[1]=get_size_x(iy,1);
	      dx[2]=get_size_x(iz,2);
	      pick_g(ix,iy,iz,gg);
	      pick_G(ix,iy,iz,GG);

	      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);

	      struct geometry geom;
	      fill_geometry_arb(ix,iy,iz,&geom,MYCOORDS);

	      struct geometry geomBL;
	      fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);
	      
	      trans_pmhd_coco(pp,pp,MYCOORDS,BLCOORDS,xx,&geom,&geomBL);

	      rho=pp[0];

	      utcon[1]=pp[2];
	      utcon[2]=pp[3];
	      utcon[3]=pp[4];

	      conv_vels(utcon,ucon3,VELPRIM,VEL3,geomBL.gg,geomBL.GG);
	      conv_vels_both(utcon,ucon,ucov,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
	      //conv_velscov(utcon,ucov,VELPRIM,VEL4,geomBL.gg,geomBL.GG);

	      dx[0]=dx[0]*sqrt(geomBL.gg[0][0]);
	      dx[1]=dx[1]*sqrt(geomBL.gg[2][2]);
	      dx[2]=2.*M_PI*sqrt(geomBL.gg[3][3]);

	      calc_tautot(pp,xxBL,dx,tautot);
	      calc_tauabs(pp,xxBL,dx,tauabs);

	      //#ifdef CGSOUTPUT
	      rho=rhoGU2CGS(rho);
	      dx[0]=lenGU2CGS(dx[0]); //dr
	      dx[1]=lenGU2CGS(dx[1]); //dth
	      dx[2]=lenGU2CGS(dx[2]); //dph
	      //#endif

	      //surface density (2nd column)
	      profiles[0][ix]+=rho*dx[1];
	      //rest mass flux (3)
	      profiles[1][ix]+=-rho*ucon[1]*dx[1]*dx[2];
	      //rho-weighted minus radial velocity (4)
	      profiles[2][ix]+=-ucon[1]*rho*dx[1];
	      //rho-weighted u_phi (5)
	      profiles[3][ix]+=ucov[3]*rho*dx[1];	
	      //abs optical depth (7)
	      profiles[5][ix]+=tauabs[1];	
	      //tot optical depth (8)
	      profiles[6][ix]+=tautot[1];	
	    }
	  //normalizing by sigma
	  profiles[2][ix]/=profiles[0][ix];
	  profiles[3][ix]/=profiles[0][ix];
	  //Keplerian u_phi (6)
	  ldouble r=xxBL[1];
	  profiles[4][ix]=(r*r/(sqrt(r*(r*r-3.*r))));  
	  //net accretion rate at given radius (7)
	  profiles[7][ix]=fabs(calc_mdot(xxBL[1],0)/calc_mdotEdd());
	  //inflow accretion rate at given radius (8)
	  profiles[8][ix]=fabs(calc_mdot(xxBL[1],1)/calc_mdotEdd());
	  //outflow accretion rate at given radius (9)
	  profiles[9][ix]=fabs(calc_mdot(xxBL[1],2)/calc_mdotEdd());
	  //luminosity at given radius (12)
	  profiles[10][ix]=calc_lum(xxBL[1])/calc_lumEdd();
	  //location of the photosphere (8)
	  profiles[11][ix]=calc_photloc(ix);
	}
    }

  return 0;
}

/*********************************************/
/* calculates scalars - total mass, accretion rate etc. */
/*********************************************/
int calc_scalars(ldouble *scalars,ldouble t)
{
  //adjust NSCALARS in problem.h


  //total mass inside the domain (2nd column)
  scalars[0]=calc_totalmass();

  //accretion rate through horizon (3)
  ldouble mdot=calc_mdot(r_horizon_BL(BHSPIN),0);
  scalars[1]=-mdot/calc_mdotEdd();

  //luminosity (4) at 0.5*rmax
  ldouble xx[4],xxBL[4];
  get_xx(NX-1,0,0,xx);
  coco_N(xx,xxBL,MYCOORDS,BLCOORDS);

  scalars[2]=calc_lum(xxBL[1]/2.)/calc_lumEdd();

  //magnetic flux through horizon parameter (5)
  scalars[3]=calc_Bflux(r_horizon_BL(BHSPIN),0)/2./sqrt(fabs(mdot));

  /*********************************************/
  //L1 ERRRORS for some problems
  /*********************************************/

#ifdef CALCL1_HDWAVE
  //temporarily here: L1 error for HDWAVE
  ldouble L1=0;
  int i;
  for(i=0;i<NX;i++)
    {
      ldouble xx=get_x(i,0);
      ldouble om=1./CC*2.*Pi;
      ldouble myrho=RHOZERO*(1.+AAA*cos(KK*xx-om*t));
      ldouble myuint=UINT*(1.+GAMMA*AAA*cos(KK*xx-om*t));
      ldouble mycs=1./CC;
      ldouble myvx=AAA*cos(KK*xx-om*t)*mycs;
      //L1 in rho:
      L1+=fabs(get_u(p,0,i,0,0)-myrho);
    }
  scalars[0]=L1/(ldouble)NX;
#endif

#ifdef CALCL1_HUBBLE
  //temporarily here: L1 error for HDWAVE
  ldouble L1=0;
  int i;
  for(i=0;i<NX;i++)
    {
      ldouble xx=get_x(i,0);      
      ldouble myrho=RHO0 / (1.+VPRIME*t);
      ldouble myuint=UINT0 / pow(1.+VPRIME*t,GAMMA);
      ldouble myvx=VPRIME*xx / (1.+VPRIME*t);
      //L1 in rho:
      L1+=fabs(get_u(p,0,i,0,0)-myrho);
    }
  scalars[0]=L1/(ldouble)NX;
#endif

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//integrates mass in the domain
ldouble
calc_totalmass()
{
  int ix,iy,iz;
  ldouble xx[4],dx[3],mass,rho,gdet;
  
  mass=0.;
  for(iz=0;iz<NZ;iz++)
    {
      for(iy=0;iy<NY;iy++)
	{
	  for(ix=0;ix<NX;ix++)
	    {
	      get_xx(ix,iy,iz,xx);
	      dx[0]=get_size_x(ix,0);
	      dx[1]=get_size_x(iy,1);
	      dx[2]=get_size_x(iz,2);
	      gdet=calc_gdet(xx);
	      rho=get_u(p,0,ix,iy,iz);
	      mass+=rho*dx[0]*dx[1]*dx[2]*gdet;
	    }
	}
    }
  return mass;
}
	  
//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates the Eddington mass accretion rate
ldouble
calc_mdotEdd()
{
  ldouble mcgs=2.23e18*MASS; //g/s

  //#ifdef CGSOUTPUT
  return mcgs;
  //#else
  //return 1.;
  //#endif
}

//**********************************************************************
//**********************************************************************
//*********************************************************************
//calculates the Eddington luminosity
ldouble
calc_lumEdd()
{
  ldouble Lcgs=1.25e38*MASS; //erg/s

  //#ifdef CGSOUTPUT
  return Lcgs;
  //#else
  //return 1.;
  //#endif
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates luminosity by integrating positive flux from the axis up to tau=1 surface
//normalized to total sphere, taken at radius radius
ldouble
calc_lum(ldouble radius)
{
  if(MYCOORDS != BLCOORDS && MYCOORDS != KSCOORDS && MYCOORDS != MKS1COORDS && MYCOORDS != MKS2COORDS)
    return -1.; //no BH

  int ix,iy,iz,iv;
  ldouble xx[4],xxBL[4],dx[3],pp[NV],Fr;

  //search for appropriate radial index
  for(ix=0;ix<NX;ix++)
    {
      get_xx(ix,0,0,xx);
      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
      if(xxBL[1]>radius) break;
    }

  ldouble lum=0.,tau=0.;

  if(NZ==1) //phi-symmetry
    {
      iz=0; 
      for(iy=0;iy<NY;iy++)
	{
	  for(iv=0;iv<NV;iv++)
	    pp[iv]=get_u(p,iv,ix,iy,iz);

	  ldouble tautot[3];

	  struct geometry geomBL;
	  fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
	  struct geometry geom;
	  fill_geometry(ix,iy,iz,&geom);
	  

	  get_xx(ix,iy,iz,xx);
	  dx[0]=get_size_x(ix,0);
	  dx[1]=get_size_x(iy,1);
	  dx[2]=get_size_x(iz,2);
	  dx[0]=dx[0]*sqrt(geom.gg[1][1]);
	  dx[1]=dx[1]*sqrt(geom.gg[2][2]);
	  dx[2]=2.*M_PI*sqrt(geom.gg[3][3]);

	  coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
	  calc_tautot(pp,xxBL,dx,tautot);

	  tau+=tautot[1];
	  if(tau>1.) break;
	  
	  trans_prad_coco(pp,pp,MYCOORDS,KERRCOORDS,xx,&geom,&geomBL);
	  prad_lab2on(pp,pp,&geomBL);

	  Fr=pp[FX(0)];	
	  if(Fr<0.) Fr=0.;

	  //#ifdef CGSOUTPUT
	  //always!
	  Fr=fluxGU2CGS(Fr);
	  dx[1]=lenGU2CGS(dx[1]);
	  dx[2]=lenGU2CGS(dx[2]);
	  //#endif

	  lum+=Fr*dx[1]*dx[2];
	}

      return lum;
    }
  else
    return -1;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates theta corresponding to integrated tau from the axis
ldouble
calc_photloc(int ix)
{
  if(MYCOORDS != BLCOORDS && MYCOORDS != KSCOORDS && MYCOORDS != MKS1COORDS && MYCOORDS != MKS2COORDS)
    return -1.; //no BH

  ldouble tau=0.,pp[NV],xx[4],xxBL[4],dx[3];

  int iz=0; int iy,iv; 

  if(NZ==1)
    {
      for(iy=0;iy<NY;iy++)
	{
	  for(iv=0;iv<NV;iv++)
	    pp[iv]=get_u(p,iv,ix,iy,iz);

	  ldouble tautot[3];

	  struct geometry geomBL;
	  fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
	  struct geometry geom;
	  fill_geometry(ix,iy,iz,&geom);

	  get_xx(ix,iy,iz,xx);
	  dx[0]=get_size_x(ix,0);
	  dx[1]=get_size_x(iy,1);
	  dx[2]=get_size_x(iz,2);
	  dx[0]=dx[0]*sqrt(geom.gg[1][1]);
	  dx[1]=dx[1]*sqrt(geom.gg[2][2]);
	  dx[2]=2.*M_PI*sqrt(geom.gg[3][3]);

	  coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
	  calc_tautot(pp,xxBL,dx,tautot);
	  tau+=tautot[1];
	  if(tau>1.) break;
	}
      return xxBL[2];
    }
  else
    return -1;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates rest mass flux through r=radius within range of thetas
//normalized to 2pi in phi
//type == 0 (net)
//type == 1 (inflow only)
//type == 2 (outflow only)
ldouble
calc_mdot(ldouble radius,int type)
{
  if(MYCOORDS != BLCOORDS && MYCOORDS != KSCOORDS && MYCOORDS != MKS1COORDS && MYCOORDS != MKS2COORDS)
    return -1.; //no BH

  int ix,iy,iz,iv;
  ldouble xx[4],xxBL[4],dx[3],mdot,rho,ucon[4],pp[NV],gg[4][5],GG[4][5],ggBL[4][5],GGBL[4][5];

  //search for appropriate radial index
  for(ix=0;ix<NX;ix++)
    {
      get_xx(ix,0,0,xx);
      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
      if(xxBL[1]>radius) break;
    }

  mdot=0.;

  if(NZ==1) //phi-symmetry
    {
      iz=0;
      for(iy=0;iy<NY;iy++)
	{
	  for(iv=0;iv<NVMHD;iv++)
	    pp[iv]=get_u(p,iv,ix,iy,iz);

	  get_xx(ix,iy,iz,xx);
	  dx[0]=get_size_x(ix,0);
	  dx[1]=get_size_x(iy,1);
	  dx[2]=get_size_x(iz,2);
	  pick_g(ix,iy,iz,gg);
	  pick_G(ix,iy,iz,GG);

	  coco_N(xx,xxBL,MYCOORDS,BLCOORDS);

	  struct geometry geom;
	  fill_geometry_arb(ix,iy,iz,&geom,MYCOORDS);

	  struct geometry geomBL;
	  fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

	  trans_pmhd_coco(pp,pp,MYCOORDS,BLCOORDS,xx,&geom,&geomBL);

	  rho=pp[0];

	  ucon[1]=pp[2];
	  ucon[2]=pp[3];
	  ucon[3]=pp[4];

	  conv_vels(ucon,ucon,VELPRIM,VEL4,geomBL.gg,geomBL.GG);

	  dx[1]=dx[1]*sqrt(geom.gg[2][2]);
	  dx[2]=2.*M_PI*sqrt(geom.gg[3][3]);

	  //#ifdef CGSOUTPUT
	  //always
	  rho=rhoGU2CGS(rho);
	  ucon[1]=velGU2CGS(ucon[1]);
	  dx[1]=lenGU2CGS(dx[1]);
	  dx[2]=lenGU2CGS(dx[2]);
	  //#endif

	  if(type==0 || (type==1 && ucon[1]<0.) || (type==2 && ucon[1]>0.))
	    mdot+=rho*ucon[1]*dx[1]*dx[2];
	}
    }
  else
    return -1;

  return mdot;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates magnetic flux through r=radius within range of thetas
//normalized to 2pi in phi
//type == 0 (default)
ldouble
calc_Bflux(ldouble radius,int type)
{
  if(MYCOORDS != BLCOORDS && MYCOORDS != KSCOORDS && MYCOORDS != MKS1COORDS && MYCOORDS != MKS2COORDS)
    return -1.; //no BH

  int ix,iy,iz,iv;
  ldouble xx[4],xxBL[4],dx[3],Psi,rho,ucon[4],pp[NV],gg[4][5],GG[4][5],ggBL[4][5],GGBL[4][5];

  //search for appropriate radial index
  for(ix=0;ix<NX;ix++)
    {
      get_xx(ix,0,0,xx);
      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
      if(xxBL[1]>radius) break;
    }

  Psi=0.;

  if(NZ==1) //phi-symmetry
    {
      iz=0;
      for(iy=0;iy<NY;iy++)
	{
	  for(iv=0;iv<NVMHD;iv++)
	    pp[iv]=get_u(p,iv,ix,iy,iz);

	  get_xx(ix,iy,iz,xx);
	  dx[0]=get_size_x(ix,0);
	  dx[1]=get_size_x(iy,1);
	  dx[2]=get_size_x(iz,2);
	  pick_g(ix,iy,iz,gg);
	  pick_G(ix,iy,iz,GG);

	  coco_N(xx,xxBL,MYCOORDS,BLCOORDS);

	  struct geometry geom;
	  fill_geometry_arb(ix,iy,iz,&geom,MYCOORDS);

	  struct geometry geomBL;
	  fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

	  trans_pmhd_coco(pp,pp,MYCOORDS,BLCOORDS,xx,&geom,&geomBL);

	  ldouble Br=fabs(pp[B1]);

	  dx[1]=dx[1]*sqrt(geom.gg[2][2]);
	  dx[2]=2.*M_PI*sqrt(geom.gg[3][3]);

	  //#ifdef CGSOUTPUT
	  //always
	  Br=sqrt(endenGU2CGS(1.))*Br/CCC;
	  dx[1]=lenGU2CGS(dx[1]);
	  dx[2]=lenGU2CGS(dx[2]);
	  //#endif

	  
	  if(type==0 || (type==1 && ucon[1]<0.) || (type==2 && ucon[1]>0.))
	    Psi+=Br*dx[1]*dx[2];
	  

	}
    }
  else
    return -1;

  return Psi;
}
