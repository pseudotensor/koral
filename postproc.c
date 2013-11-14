//KORAL - misc.c
//routines for postprocessing 
//used both on the go and separately

#include "ko.h"

/*********************************************/
/* calculates radial profiles - L(r) etc. */
/* uses mostly primitives, but may use avg quantities as well */
/* however, this part would work only when postprocessing with */
/* ./avg, otherwise pavg hold non-normalized numbers */
/*********************************************/
int calc_radialprofiles(ldouble profiles[][NX])
{
  //adjust NRADPROFILES in problem.h

  int ix,iy,iz,iv,i,j;
  ldouble xx[4],xxBL[4],dx[3],dxph[3],dxcgs[3],mdot,rho,ucon[4],utcon[4],ucon3[4];
  ldouble rhouconr,Tij[4][4],Rij[4][4],Trt,Rrt,bsq,bcon[4],bcov[4];
  ldouble ucov[4],pp[NV],gg[4][5],GG[4][5],ggBL[4][5],GGBL[4][5];
  ldouble tautot[3],tauabs[3];
  ldouble avgsums[NV+NAVGVARS][NX];

  //search for appropriate radial index
  for(ix=0;ix<NX;ix++)
    {
      //vertically integrated/averaged profiles

      for(iv=0;iv<NRADPROFILES;iv++)
	profiles[iv][ix]=0.;

      for(iv=0;iv<NAVGVARS;iv++)
	avgsums[iv][ix]=0.;

 #ifdef BHDISK_PROBLEMTYPE
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
	      
	      pick_g(ix,iy,iz,gg);
	      pick_G(ix,iy,iz,GG);

	      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);

	      struct geometry geom;
	      fill_geometry_arb(ix,iy,iz,&geom,MYCOORDS);

	      struct geometry geomBL;
	      fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);
	      
	      //to BL     
	      trans_pmhd_coco(pp,pp,MYCOORDS,BLCOORDS,xx,&geom,&geomBL);

	      /*
		dxph[0]=get_size_x(ix,0)*sqrt(geom.gg[1][1]);
		dxph[1]=get_size_x(iy,1)*sqrt(geom.gg[2][2]);
		dxph[2]=get_size_x(iz,2)*sqrt(geom.gg[3][3]);
		dx[0]=dxph[0]/sqrt(geomBL.gg[1][1]);
		dx[1]=dxph[1]/sqrt(geomBL.gg[2][2]);
		dx[2]=dxph[2]/sqrt(geomBL.gg[3][3]);
	      */
	      ldouble dxph[3];
	      ldouble xx1[4],xx2[4];
	      xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(ix+1,1);xx2[2]=get_xb(iy,1);xx2[3]=get_xb(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	      dx[0]=fabs(xx2[1]-xx1[1]);
	      xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(ix,1);xx2[2]=get_xb(iy+1,1);xx2[3]=get_xb(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	      dx[1]=fabs(xx2[2]-xx1[2]);
	      dx[2]=2.*M_PI;
	      dxph[0]=dx[0]*sqrt(geomBL.gg[1][1]);
	      dxph[1]=dx[1]*sqrt(geomBL.gg[2][2]);
	      dxph[2]=dx[2]*sqrt(geomBL.gg[3][3]);
	      

	      if(doingavg)
		{
		  
	
		  rho=get_uavg(pavg,RHO,ix,iy,iz);
		  bsq=get_uavg(pavg,AVGBSQ,ix,iy,iz);
		  utcon[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  utcon[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  utcon[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  rhouconr=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz);

		  
		  Trt=get_uavg(pavg,AVGRHOUCONUCOV(1,0),ix,iy,iz)
		    + GAMMA*get_uavg(pavg,AVGUUUCONUCOV(1,0),ix,iy,iz)
		    + get_uavg(pavg,AVGBSQUCONUCOV(1,0),ix,iy,iz)
		    - get_uavg(pavg,AVGBCONBCOV(1,0),ix,iy,iz); 

#ifdef RADIATION  
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz); 
		  indices_2221(Rij,Rij,geomBL.gg);

		  Rrt = Rij[1][0];
#endif

		}
	      else
		{
		  rho=pp[0];
		  utcon[1]=pp[2];
		  utcon[2]=pp[3];
		  utcon[3]=pp[4];
		  
#ifdef MAGNFIELD
		  calc_bcon_4vel(pp,ucon,ucov,bcon);
		  indices_21(bcon,bcov,geomBL.gg); 
		  bsq = dot(bcon,bcov);
#endif

		  conv_vels_both(utcon,ucon,ucov,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
		  rhouconr=rho*utcon[1];

		  calc_Tij(pp,&geomBL,Tij);
		  indices_2221(Tij,Tij,geomBL.gg);

		  Trt = Tij[1][0];

#ifdef RADIATION
		  calc_Rij(pp,&geomBL,Rij);
		  indices_2221(Rij,Rij,geomBL.gg);

		  Rrt = Rij[1][0];
#endif
		}


	      
	      calc_tautot(pp,xxBL,dx,tautot);
	      calc_tauabs(pp,xxBL,dx,tauabs);


	      //surface density (2) (column)
	      profiles[0][ix]+=rho*dxph[1];

	      //surface density in the inflow (23)
	      if(utcon[1]<0.)
		profiles[21][ix]+=rho*dxph[1];

	      
	      //rest mass flux (3)
	      profiles[1][ix]+=-rhouconr*dx[1]*dx[2]*geomBL.gdet;

	      //total mhd energy flux (14)
	      profiles[12][ix]+=(-Trt)*dx[1]*dx[2]*geomBL.gdet;
	      //outflowin mhd energy flux (15)
	      if(utcon[1]>0.)
		profiles[13][ix]+=(-Trt)*dx[1]*dx[2]*geomBL.gdet;
	      //jet mhd energy flux (16)
	      if(utcon[1]>0. && bsq>rho)
		profiles[14][ix]+=(-Trt)*dx[1]*dx[2]*geomBL.gdet;

	      //total rad energy flux (17)
	      profiles[15][ix]+=(-Rrt)*dx[1]*dx[2]*geomBL.gdet;
	      //outflowin rad energy flux (18)
	      if(utcon[1]>0.)
		profiles[16][ix]+=(-Rrt)*dx[1]*dx[2]*geomBL.gdet;
	      //jet rad energy flux (19)
	      if(utcon[1]>0. && bsq>rho)
		profiles[17][ix]+=(-Rrt)*dx[1]*dx[2]*geomBL.gdet;
	      
	      //outflowin mass flux (20)
	      if(utcon[1]>0.)
		profiles[18][ix]+=(-rhouconr)*dx[1]*dx[2]*geomBL.gdet;
	      //jet mass flux (21)
	      if(utcon[1]>0. && bsq>rho)
		profiles[19][ix]+=(-rhouconr)*dx[1]*dx[2]*geomBL.gdet;

	      //rho-weighted minus radial velocity (4)
	      profiles[2][ix]+=-utcon[1]*rho*dxph[1];

	      //rho-weighted minus radial velocity in the inflow (24)
	      if(utcon[1]<0.)
		profiles[22][ix]+=-utcon[1]*rho*dxph[1];

	      //abs optical depth (7)
	      profiles[5][ix]+=tauabs[1];	
	      //tot optical depth (8)
	      profiles[6][ix]+=tautot[1];

	      for(iv=0;iv<NV+NAVGVARS;iv++)
		avgsums[iv][ix]+=get_uavg(pavg,iv,ix,iy,iz)*dxph[1];

	      //<(rho+u+bsq/2)u^r><u_phi> (5)
	      //profiles[3][ix]+=get_uavg(pavg,AVGWUCON(1),ix,iy,iz)*get_uavg(pavg,AVGUCOV(3),ix,iy,iz)*dxph[1];
	      profiles[3][ix]+=get_uavg(pavg,AVGRHOUCOV(3),ix,iy,iz)*dxph[1];
	    }

	  //normalizing by sigma
	  profiles[2][ix]/=profiles[0][ix];
	  profiles[22][ix]/=profiles[21][ix];

	  //normalizing by <(rho+u+bsq/2)u^r>
	  //profiles[3][ix]/=avgsums[AVGWUCON(1)][ix];
	  profiles[3][ix]/=profiles[0][ix];
 
	  //Keplerian u_phi (6)
	  ldouble r=xxBL[1];
	  profiles[4][ix]=((r*r-2.*BHSPIN*sqrt(r)+BHSPIN*BHSPIN)/(sqrt(r*(r*r-3.*r+2.*BHSPIN*sqrt(r)))));  



	  //net accretion rate at given radius (9)
	  profiles[7][ix]=fabs(calc_mdot(xxBL[1],0));
	  //inflow accretion rate at given radius (10)
	  profiles[8][ix]=fabs(calc_mdot(xxBL[1],1));
	  //outflow accretion rate at given radius (11)
	  profiles[9][ix]=fabs(calc_mdot(xxBL[1],2));
	  //luminosity at given radius (12)
	  profiles[10][ix]=calc_lum(xxBL[1],0);
	  //luminosity at given radius (22)
	  profiles[20][ix]=calc_lum(xxBL[1],1);
	  //location of the photosphere (13)
	  profiles[11][ix]=calc_photloc(ix);




	}

#endif
    }

  return 0;
}

/*********************************************/
/* calculates scalar s - total mass, accretion rate etc. */
/*********************************************/
int calc_scalars(ldouble *scalars,ldouble t)
{
  //adjust NSCALARS in problem.h


  //total mass inside the domain (2nd column)
  scalars[0]=calc_totalmass();

  //accretion rate through horizon (3)
  ldouble mdotscale = rhoGU2CGS(1.)*velGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.);

  ldouble mdot=calc_mdot(r_horizon_BL(BHSPIN),0);
  scalars[1]=-mdot;
  scalars[4]=-mdot*mdotscale/calc_mdotEdd();

  //luminosity (4) at 0.5*rmax
  ldouble xx[4],xxBL[4];
  get_xx(NX-1,0,0,xx);
  coco_N(xx,xxBL,MYCOORDS,BLCOORDS);

  scalars[2]=calc_lum(xxBL[1]/2.,0)/calc_lumEdd();

  //magnetic flux through horizon parameter (5)
  ldouble Bflux=calc_Bflux(r_horizon_BL(BHSPIN),0);
  //ldouble Bfluxcgs=Bflux*sqrt(endenGU2CGS(1.))*lenGU2CGS(1.)*lenGU2CGS(1.)/CCC;
  //scalars[3]=Bfluxcgs/2./sqrt(fabs(mdot));

  scalars[3]=Bflux;///sqrt(fabs(mdot))*sqrt(4.*M_PI)/2.;

  /*********************************************/
  //L1 ERRRORS for some problems
  /*********************************************/

#ifdef CALCL1_RMHDWAVE
  //temporarily here: L1 error for RMHDWAVE
  ldouble L1=0;
  int i;
  for(i=0;i<NX;i++)
    {
      calc_primitives(i,0,0,0);
      ldouble xx=get_x(i,0);
      ldouble dx=get_size_x(i,0);
      ldouble myrho=RHOZERO+DRRE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx)-DRIM/DRRE*sin(OMRE*t-KK*xx));
      //L1 in rho:
      L1+=fabs(get_u(p,RHO,i,0,0)-myrho)*dx;
    }
  scalars[0]=L1;///(ldouble)NX;
#endif

#ifdef CALCL1_HDWAVE
  //temporarily here: L1 error for HDWAVE
  ldouble L1=0;
  int i;
  for(i=0;i<NX;i++)
    {
      calc_primitives(i,0,0,0);
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
  //temporarily here: L1 error for HUBBLE
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
  ldouble mcgs=1.09649*2.23e18*MASS; //g/s assuming eta=0.057

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
calc_lum(ldouble radius,int type)
{
#ifndef BHDISK_PROBLEMTYPE
    return -1.; //no BH
#endif

#ifdef RADIATION

    int ix,iy,iz,iv,i,j;
    ldouble xx[4],xxBL[4],dx[3],pp[NV],Fr;
    ldouble Rij[4][4],Rtt,ehat,ucongas[4];
    ldouble tautot[3];
    ldouble gdet;

 
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

	  struct geometry geomBL;
	  fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
	  struct geometry geom;
	  fill_geometry(ix,iy,iz,&geom);

	  if(doingavg)
	    {
	      

	      ldouble dxph[3];
	      ldouble xx1[4],xx2[4];
	      xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(ix+1,1);xx2[2]=get_xb(iy,1);xx2[3]=get_xb(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	      dx[0]=fabs(xx2[1]-xx1[1]);
	      xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(ix,1);xx2[2]=get_xb(iy+1,1);xx2[3]=get_xb(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	      dx[1]=fabs(xx2[2]-xx1[2]);
	      dx[2]=2.*M_PI;
	      dxph[0]=dx[0]*sqrt(geomBL.gg[1][1]);
	      dxph[1]=dx[1]*sqrt(geomBL.gg[2][2]);
	      dxph[2]=dx[2]*sqrt(geomBL.gg[3][3]);
	      gdet=geomBL.gdet;

	      PLOOP(iv)
		pp[iv]=get_uavg(pavg,iv,ix,iy,iz);

	      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
	      calc_tautot(pp,xxBL,dxph,tautot);

	      tau+=tautot[1];

	      if(type==0) //R^r_t outside photosphere
		{
		  if(tau>1.) break;
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
		  indices_2221(Rij,Rij,geomBL.gg);
		  Fr=-Rij[1][0];
		  if(Fr<0.) Fr=0.;
		}
	      else if(type==1) //R^r_t everywhere in outflow
		{
		  //ehat=get_uavg(pavg,AVGEHAT,ix,iy,iz);
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
		  indices_2221(Rij,Rij,geomBL.gg);
		  ldouble uconr=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  Fr=-Rij[1][0];// + ehat*uconr);
		  if(uconr<0. || Fr<0.) Fr=0.;
		}
	      else
		Fr=0.;
	    }
	  else
	    {
	      get_xx(ix,iy,iz,xx);
	      dx[0]=get_size_x(ix,0);
	      dx[1]=get_size_x(iy,1);
	      dx[2]=2.*M_PI;
	      gdet=geom.gdet;

	      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
	      calc_tautot(pp,xxBL,dx,tautot);

	      tau+=tautot[1];

	      if(type==0) //R^r_t outside photosphere
		{
		  if(tau>1.) break;	  
		  //trans_prad_coco(pp,pp,MYCOORDS,KERRCOORDS,xx,&geom,&geomBL);
		  //prad_lab2on(pp,pp,&geomBL);
		  calc_Rij(pp,&geom,Rij); 
		  indices_2221(Rij,Rij,geom.gg);
		  Fr=-Rij[1][0];
		  if(Fr<0.) Fr=0.;
		}
	      else if(type==1) //R^r_t in the outflow region
		{
		  //calc_ff_Rtt(pp,&Rtt,ucongas,&geom);
		  //ehat=-Rtt;
		  calc_Rij(pp,&geom,Rij); 
		  indices_2221(Rij,Rij,geom.gg);
		  ucongas[1]=pp[2];
		  ucongas[2]=pp[3];
		  ucongas[3]=pp[4];	      
		  conv_vels(ucongas,ucongas,VELPRIM,VEL4,geom.gg,geom.GG);
		  Fr=-Rij[1][0];// + ehat*ucongas[1];
		  if(ucongas[1]<0.)
		    Fr=0.;
		}
	      else
		Fr=0.;
	    }

	  //#ifdef CGSOUTPUT
	  //never!
	  //Fr=fluxGU2CGS(Fr);
	  //dx[1]=lenGU2CGS(dx[1]);
	  //dx[2]=lenGU2CGS(dx[2]);
	  //#endif		  


	  //printf("%e %e %e -> %e\n",Fr,gdet,dx[1],dx[2],lum);getch();

	  lum+=gdet*Fr*dx[1]*dx[2];
	}
      
      return lum;
    }
  else

#endif
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
#ifndef BHDISK_PROBLEMTYPE
  return -1.; //no BH
#endif

  int ix,iy,iz,iv;
  ldouble xx[4],xxBL[4],dx[3],mdot,gdet,rho,rhouconr,ucon[4],pp[NV],gg[4][5],GG[4][5],ggBL[4][5],GGBL[4][5];

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
	  
	  

	  if(doingavg)
	    {
	      rho=get_uavg(pavg,RHO,ix,iy,iz);
	      ucon[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      ucon[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      ucon[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      rhouconr=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz);	      	      
	      gdet=geomBL.gdet;
	      ldouble dxph[3];
	      ldouble xx1[4],xx2[4];
	      xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(ix+1,1);xx2[2]=get_xb(iy,1);xx2[3]=get_xb(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	      dx[0]=fabs(xx2[1]-xx1[1]);
	      xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(ix,1);xx2[2]=get_xb(iy+1,1);xx2[3]=get_xb(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	      dx[1]=fabs(xx2[2]-xx1[2]);
	      dx[2]=2.*M_PI;
	      /*
	      dxph[0]=get_size_x(ix,0)*sqrt(geom.gg[1][1]);
	      dxph[1]=get_size_x(iy,1)*sqrt(geom.gg[2][2]);
	      dxph[2]=get_size_x(iz,2)*sqrt(geom.gg[3][3]);
	      dx[0]=dxph[0]/sqrt(geomBL.gg[1][1]);
	      dx[1]=dxph[1]/sqrt(geomBL.gg[2][2]);
	      dx[2]=dxph[2]/sqrt(geomBL.gg[3][3]);
	      */
	    }
	  else
	    {
	      rho=pp[0];
	      ucon[1]=pp[2];
	      ucon[2]=pp[3];
	      ucon[3]=pp[4];
	      
	      conv_vels(ucon,ucon,VELPRIM,VEL4,geom.gg,geom.GG);
	      rhouconr=rho*ucon[1];
	      gdet=geom.gdet;	
	      dx[1]=dx[1];	     
	      dx[2]=2.*M_PI;
	      /*
	      trans_pmhd_coco(pp,pp,MYCOORDS,BLCOORDS,xx,&geom,&geomBL);
	      ldouble dxph[3];
	      conv_vels(ucon,ucon,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
	      rhouconr=rho*ucon[1];
	      ldouble xx1[4],xx2[4];
	      xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(ix+1,1);xx2[2]=get_xb(iy,1);xx2[3]=get_xb(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	      dx[0]=fabs(xx2[1]-xx1[1]);
	      xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(ix,1);xx2[2]=get_xb(iy+1,1);xx2[3]=get_xb(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	      dx[1]=fabs(xx2[2]-xx1[2]);
	      dx[2]=2.*M_PI;
	      gdet=geomBL.gdet;
	      */
	      
	    }

	 

	  

	  //#ifdef CGSOUTPUT
	  //always

	  /*
	  rho=rhoGU2CGS(rho);
	  ucon[1]=velGU2CGS(ucon[1]);
	  dx[1]=lenGU2CGS(dx[1]);
	  dx[2]=lenGU2CGS(dx[2]);
	  */
	  //#endif

	  if(type==0 || (type==1 && ucon[1]<0.) || (type==2 && ucon[1]>0.))
	    mdot+=gdet*rhouconr*dx[1]*dx[2];	     
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

  #ifdef MAGNFIELD

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

	  dx[0]=get_size_x(ix,0);
	  dx[1]=get_size_x(iy,1);
	  dx[2]=get_size_x(iz,2);

	  //coco_N(xx,xxBL,MYCOORDS,BLCOORDS);

	  struct geometry geom;
	  fill_geometry_arb(ix,iy,iz,&geom,MYCOORDS);

	  /*
	  struct geometry geomBL;
	  fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

	  trans_pmhd_coco(pp,pp,MYCOORDS,BLCOORDS,xx,&geom,&geomBL);
	  */

	  ldouble Br=fabs(pp[B1]);

	  dx[1]=dx[1];//*sqrt(geom.gg[2][2]);
	  dx[2]=2.*M_PI;//*sqrt(geom.gg[3][3]);

	  //#ifdef CGSOUTPUT
	  //always
	  //Br=sqrt(endenGU2CGS(1.))*Br/CCC;
	  //dx[1]=lenGU2CGS(dx[1]);
	  //dx[2]=lenGU2CGS(dx[2]);
	  //#endif

	  
	  if(type==0 || (type==1 && ucon[1]<0.) || (type==2 && ucon[1]>0.))
	    Psi+=geom.gdet*Br*dx[1]*dx[2];
	  

	}
    }
  else
    return -1;

  return Psi;
#else
  return -1;
#endif
}
