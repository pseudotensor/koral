//KORAL - misc.c
//misceleanous routines

#include "ko.h"

/*********************************************/
/* calculates radial profiles - L(r) etc. */
/*********************************************/
int calc_radialprofiles(ldouble profiles[][NX])
{
  int ix,iy,iz,iv;
  ldouble xx[4],xxBL[4],dx[3],mdot,rho,ucon[4],ucon3[4];
  ldouble ucov[4],pp[NV],gg[4][5],GG[4][5],ggBL[4][5],GGBL[4][5];
  ldouble tautot[3];

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
	      calc_primitives(ix,iy,iz);

	      for(iv=0;iv<NVHD;iv++)
		pp[iv]=get_u(p,iv,ix,iy,iz);

	      get_xx(ix,iy,iz,xx);
	      dx[0]=get_size_x(ix,0);
	      dx[1]=get_size_x(iy,1);
	      dx[2]=get_size_x(iz,2);
	      pick_g(ix,iy,iz,gg);
	      pick_G(ix,iy,iz,GG);

	      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);

	      calc_g_arb(xxBL,ggBL,KERRCOORDS);
	      calc_G_arb(xxBL,GGBL,KERRCOORDS);

	      trans_phd_coco(pp,pp,MYCOORDS,BLCOORDS,xx,gg,GG,ggBL,GGBL);

	      rho=pp[0];

	      ucon[1]=pp[2];
	      ucon[2]=pp[3];
	      ucon[3]=pp[4];

	      conv_vels(ucon,ucon3,VELPRIM,VEL3,ggBL,GGBL);
	      conv_vels(ucon,ucon,VELPRIM,VEL4,ggBL,GGBL);
	     
	      indices_21(ucon,ucov,ggBL);

	      dx[0]=dx[0]*sqrt(ggBL[0][0]);
	      dx[1]=dx[1]*sqrt(ggBL[2][2]);
	      dx[2]=2.*M_PI*sqrt(ggBL[3][3]);

	      calc_tautot(pp,xxBL,dx,tautot);

	      //surface density
	      profiles[0][ix]+=rho*dx[1]*dx[2];
	      //rest mass flux
	      profiles[1][ix]+=-rho*ucon[1]*dx[1]*dx[2];
	      //rho-weighted radial velocity
	      profiles[2][ix]+=ucon[1]*rho*dx[1]*dx[2];
	      //rho-weighted u_phi
	      profiles[3][ix]+=ucov[3]*rho*dx[1]*dx[2];	
	      //optical depth
	      profiles[5][ix]+=tautot[1];	
	    }
	  //normalizing by sigma
	  profiles[2][ix]/=profiles[0][ix];
	  profiles[3][ix]/=profiles[0][ix];
	  //Keplerian u_phi
	  ldouble r=xxBL[1];
	  profiles[4][ix]=(r*r/(sqrt(r*(r*r-3.*r))));	
	}
    }

  return 0;
}

/*********************************************/
/* calculates scalars - total mass, accretion rate etc. */
/*********************************************/
int calc_scalars(ldouble *scalars)
{
  //total mass inside the domain
  scalars[0]=calc_totalmass();

  //accretion rate through horizon
  scalars[1]=calc_mdot(r_horizon_BL(BHSPIN))/calc_mdotEdd();
  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//integrates mass in the domain, assumes flat spacetime and uniform mesh
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
	      rho=get_u(u,0,ix,iy,iz);
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
  ldouble mcgs=2.23e18*MASS; //g/cm

#ifdef CGSOUTPUT
  return mcgs;
#else
  return 1.;
#endif
}

//**********************************************************************
//**********************************************************************
//*********************************************************************
//calculates the Eddington luminosity
ldouble
calc_lumEdd()
{
  ldouble Lcgs=1.25e38*MASS; //erg/s

#ifdef CGSOUTPUT
  return Lcgs;
#else
  return 1.;
#endif
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates rest mass flux through r=radius within range of thetas
//normalized to 2pi in phi
ldouble
calc_mdot(ldouble radius)
{
  if(MYCOORDS != BLCOORDS && MYCOORDS != KSCOORDS && MYCOORDS != MKS1COORDS)
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
	  for(iv=0;iv<NVHD;iv++)
	    pp[iv]=get_u(p,iv,ix,iy,iz);

	  get_xx(ix,iy,iz,xx);
	  dx[0]=get_size_x(ix,0);
	  dx[1]=get_size_x(iy,1);
	  dx[2]=get_size_x(iz,2);
	  pick_g(ix,iy,iz,gg);
	  pick_G(ix,iy,iz,GG);

	  coco_N(xx,xxBL,MYCOORDS,BLCOORDS);

	  calc_g_arb(xxBL,ggBL,KERRCOORDS);
	  calc_G_arb(xxBL,GGBL,KERRCOORDS);

	  trans_phd_coco(pp,pp,MYCOORDS,BLCOORDS,xx,gg,GG,ggBL,GGBL);

	  rho=pp[0];

	  ucon[1]=pp[2];
	  ucon[2]=pp[3];
	  ucon[3]=pp[4];

	  conv_vels(ucon,ucon,VELPRIM,VEL4,ggBL,GGBL);

	  dx[1]=dx[1]*sqrt(ggBL[2][2]);
	  dx[2]=2.*M_PI*sqrt(ggBL[3][3]);

#ifdef CGSOUTPUT
	  rho=rhoGU2CGS(rho);
	  ucon[1]=velGU2CGS(ucon[1]);
	  dx[1]=lenGU2CGS(dx[1]);
	  dx[2]=lenGU2CGS(dx[2]);
#endif

	  mdot+=rho*ucon[1]*dx[1]*dx[2];
	}
    }
  else
    return -1;

  return -mdot*2.;
}
	  
//**********************************************************************
//**********************************************************************
//**********************************************************************
//calls gnuplot for 2d
int
convert_out2gif_2d(char *fname,char *fname2,int niter,ldouble t)
{
  FILE *fgnu=fopen("plot.gp","w");
  char bufor[50];

  //PROBLEMS/XXX/out2gid_2d.c
  #include PR_OUT2GIF_2D

  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  int i=system("gnuplot plot.gp ");
  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calls gnuplot for 1d
int
convert_out2gif_1d(char *fname,char *fname2,int niter,ldouble t)
{
  FILE *fgnu=fopen("plot.gp","w");
  char bufor[50];

  //PROBLEMS/XXX/out2gid_1d.c
  #include PR_OUT2GIF_1D

  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  int i=system("gnuplot plot.gp ");
  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
ldouble my_min(ldouble a, ldouble b)
{
  if(a<b) return a;
  else return b;
}

ldouble my_sign(ldouble x)
{
  if(x>0.) return 1.;
  if(x<0.) return -1.;
  if(x==0.) return 0.;
  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//initializes arrays
int
initialize_arrays()
{
  int i,j,k;

  
  //grid 
  x=(ldouble*)malloc((NX+NY+NZ+6*NG)*sizeof(ldouble));
  xb=(ldouble*)malloc((NX+1+NY+1+NZ+1+6*NG)*sizeof(ldouble));

  //primitives at cell centers
  p=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  //temporary
  pt0=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  //primitives at cell centers after reconstruction
  px=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  py=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  pz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));

  //conserved averages
  u=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  //source terms at cell centers
  s=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));

  //conserved at cell centers
  cellflag=(int*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NFLAGS*sizeof(int));
 
  //metric at cell centers
  g=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*gSIZE*sizeof(ldouble));
  //metric at cell centers
  G=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*gSIZE*sizeof(ldouble));
  //Kristofels at cell centers
  gKr=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*64*sizeof(ldouble));

  //LNRF basis one-forms
  emuup=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //LNRF basis vectors
  emulo=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //LNRF basis one-forms bx
  emuupbx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //LNRF basis vectors bx
  emulobx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //LNRF basis one-forms by
  emuupby=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*16*sizeof(ldouble));
  //LNRF basis vectors by
  emuloby=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*16*sizeof(ldouble));
  //LNRF basis one-forms bz
  emuupbz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*16*sizeof(ldouble));
  //LNRF basis vectors bz
  emulobz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*16*sizeof(ldouble));

  //ortonormal tetrad one-forms
  tmuup=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //ortonormal tetrad vectors
  tmulo=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //ortonormal tetrad one-forms bx
  tmuupbx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //ortonormal tetrad vectors bx
  tmulobx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //ortonormal tetrad one-forms by
  tmuupby=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*16*sizeof(ldouble));
  //ortonormal tetrad vectors by
  tmuloby=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*16*sizeof(ldouble));
  //ortonormal tetrad one-forms bz
  tmuupbz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*16*sizeof(ldouble));
  //ortonormal tetrad vectors bz
  tmulobz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*16*sizeof(ldouble));

  //as above but for the suplementary system of coordinates
  //LNRF2= basis one-forms
  emuup2=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //LNRF basis vectors
  emulo2=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //LNRF basis one-forms bx
  emuupbx2=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //LNRF basis vectors bx
  emulobx2=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //LNRF basis one-forms by
  emuupby2=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*16*sizeof(ldouble));
  //LNRF basis vectors by
  emuloby2=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*16*sizeof(ldouble));
  //LNRF basis one-forms bz
  emuupbz2=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*16*sizeof(ldouble));
  //LNRF basis vectors bz
  emulobz2=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*16*sizeof(ldouble));
 
  //ortonormal tetrad one-forms
  tmuup2=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //ortonormal tetrad vectors
  tmulo2=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //ortonormal tetrad one-forms bx
  tmuupbx2=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //ortonormal tetrad vectors bx
  tmulobx2=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //ortonormal tetrad one-forms by
  tmuupby2=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*16*sizeof(ldouble));
  //ortonormal tetrad vectors by
  tmuloby2=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*16*sizeof(ldouble));
  //ortonormal tetrad one-forms bz
  tmuupbz2=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*16*sizeof(ldouble));
  //ortonormal tetrad vectors bz
  tmulobz2=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*16*sizeof(ldouble));

  //left-interpolated primitives at cell x-faces
  pbLx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG )*(NZ+2*NG)*NV*sizeof(ldouble));
  //right-interpolated primitives at cell x-faces
  pbRx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  //left-interpolated primitives at cell y-faces
  pbLy=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*NV*sizeof(ldouble));
  //right-interpolated primitives at cell y-faces
  pbRy=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*NV*sizeof(ldouble));
  //left-interpolated primitives at cell z-faces
  pbLz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*NV*sizeof(ldouble));
  //right-interpolated primitives at cell z-faces
  pbRz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*NV*sizeof(ldouble));

   //left-interpolated primitives at cell x-faces
  sbLx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG )*(NZ+2*NG)*NV*sizeof(ldouble));
  //right-interpolated primitives at cell x-faces
  sbRx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  //left-interpolated primitives at cell y-faces
  sbLy=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*NV*sizeof(ldouble));
  //right-interpolated primitives at cell y-faces
  sbRy=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*NV*sizeof(ldouble));
  //left-interpolated primitives at cell z-faces
  sbLz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*NV*sizeof(ldouble));
  //right-interpolated primitives at cell z-faces
  sbRz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*NV*sizeof(ldouble));

 //left-interpolated conserved at cell x-faces - unused
  ubLx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG )*(NZ+2*NG)*NV*sizeof(ldouble));
  //right-interpolated conserved at cell x-faces
  ubRx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  //left-interpolated conserved at cell y-faces
  ubLy=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*NV*sizeof(ldouble));
  //right-interpolated conserved at cell y-faces
  ubRy=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*NV*sizeof(ldouble));
  //left-interpolated conserved at cell z-faces
  ubLz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*NV*sizeof(ldouble));
  //right-interpolated conserved at cell z-faces
  ubRz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*NV*sizeof(ldouble));

  //corrected flux at x faces
  flbx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG )*(NZ+2*NG)*NV*sizeof(ldouble));
  //corrected flux at x faces
  flby=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1 )*(NZ+2*NG)*NV*sizeof(ldouble));
  //corrected flux at x faces
  flbz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG )*(NZ+2*NG+1)*NV*sizeof(ldouble));
  //flux based on left-interpolated conserved at cell x-faces
  flLx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG )*(NZ+2*NG)*NV*sizeof(ldouble));
  //flux based on right-interpolated conserved at cell x-faces
  flRx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  //flux based on left-interpolated conserved at cell y-faces
  flLy=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*NV*sizeof(ldouble));
  //flux based on right-interpolated conserved at cell y-faces
  flRy=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*NV*sizeof(ldouble));
  //flux based on left-interpolated conserved at cell z-faces
  flLz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*NV*sizeof(ldouble));
  //flux based on right-interpolated conserved at cell z-faces
  flRz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*NV*sizeof(ldouble));


  //Krzysie at cell x-faces
  gKrbx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG )*(NZ+2*NG)*64*sizeof(ldouble));
   //Krzysie at cell x-faces
  gKrby=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1 )*(NZ+2*NG)*64*sizeof(ldouble));
 //Krzysie at cell x-faces
  gKrbz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG )*(NZ+2*NG+1)*64*sizeof(ldouble));

  //metric at cell x-faces
  gbx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG )*(NZ+2*NG)*gSIZE*sizeof(ldouble));
  //metric at cell y-faces
  gby=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*gSIZE*sizeof(ldouble));
  //metric at cell z-faces
  gbz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*gSIZE*sizeof(ldouble));
  //metric at cell x-faces
  Gbx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG )*(NZ+2*NG)*gSIZE*sizeof(ldouble));
  //metric at cell y-faces
  Gby=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*gSIZE*sizeof(ldouble));
  //metric at cell z-faces
  Gbz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*gSIZE*sizeof(ldouble));


  //indices of the ghost cells
  gcidx=(int**)malloc(3*sizeof(int*));
  gcidx[0]=(int*)malloc(2*NG*sizeof(int));
  gcidx[1]=(int*)malloc(2*NG*sizeof(int));
  gcidx[2]=(int*)malloc(2*NG*sizeof(int));

  //auxiliary primitive arrays
  ut0=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  ut1=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  ut2=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  ut3=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  ut4=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  du=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  u_bak=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  p_bak=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  u_step1=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  u_step2=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  u_step3=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  u_step4=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));

  //wavespeeds hd and rad - max(al,ar)
  ahdx=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  ahdy=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  ahdz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  aradx=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  arady=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  aradz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));

  //wavespeeds hd and rad - leftgoing
  ahdxl=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  ahdyl=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  ahdzl=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  aradxl=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  aradyl=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  aradzl=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));

  //wavespeeds hd and rad - lrightgoing
  ahdxr=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  ahdyr=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  ahdzr=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  aradxr=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  aradyr=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  aradzr=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
int
free_arrays()
{
  int i1,i,j;
  free(cellflag);
  free(x);
  free(xb);
  free(p);
  free(pt0);
  free(px);
  free(py);
  free(pz);
  free(u);
  free(g);
  free(G);
  free(gKr);
  free(gKrbx);
  free(gKrby);
  free(gKrbz);
  free(emuup);
  free(emulo);
  free(emuupbx);
  free(emulobx);
  free(emuupby);
  free(emuloby);
  free(emuupbz);
  free(emulobz);
  free(tmuup);
  free(tmulo);
  free(tmuupbx);
  free(tmulobx);
  free(tmuupby);
  free(tmuloby);
  free(tmuupbz);
  free(tmulobz);
  free(emuup2);
  free(emulo2);
  free(emuupbx2);
  free(emulobx2);
  free(emuupby2);
  free(emuloby2);
  free(emuupbz2);
  free(emulobz2);
  free(tmuup2);
  free(tmulo2);
  free(tmuupbx2);
  free(tmulobx2);
  free(tmuupby2);
  free(tmuloby2);
  free(tmuupbz2);
  free(tmulobz2);
  free(pbLx);
  free(pbRx);
  free(pbLy);
  free(pbRy);
  free(pbLz);
  free(pbRz);
  free(sbLx);
  free(sbRx);
  free(sbLy);
  free(sbRy);
  free(sbLz);
  free(sbRz);
  free(ubLx);
  free(ubRx);
  free(ubLy);
  free(ubRy);
  free(ubLz);
  free(ubRz);
  free(flbx);
  free(flby);
  free(flbz);
  free(flLx);
  free(flRx);
  free(flLy);
  free(flRy);
  free(flLz);
  free(flRz);
  free(gbx);
  free(gby);
  free(gbz);
  free(Gbx);
  free(Gby);
  free(Gbz);
 
  free(gcidx[0]);
  free(gcidx[1]);
  free(gcidx[2]);
  free(gcidx);

  free(ut0);
  free(ut1);
  free(ut2);
  free(ut3);
  free(ut4);
  free(du);
  free(u_bak);
  free(p_bak);
  free(u_step1);
  free(u_step2);
  free(u_step3);
  free(u_step4);
  free(aradx);
  free(arady);
  free(aradz);
  free(ahdx);
  free(ahdy);
  free(ahdz);
  free(aradxl);
  free(aradyl);
  free(aradzl);
  free(ahdxl);
  free(ahdyl);
  free(ahdzl);
  free(aradxr);
  free(aradyr);
  free(aradzr);
  free(ahdxr);
  free(ahdyr);
  free(ahdzr);
 
  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//inverse 4by4 matrix
int
inverse_44matrix(ldouble a[][4], ldouble ia[][4])
{

  ldouble mat[16],dst[16];
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      mat[i*4+j]=a[i][j];

  ldouble	tmp[12]; ldouble	src[16]; ldouble det;
  /* transpose matrix */
  for (i = 0; i <4; i++)
    {
      src[i]=mat[i*4];
      src[i+4]=mat[i*4+1];
      src[i+8]=mat[i*4+2];
      src[i+12]=mat[i*4+3];
    }
  /* calculate pairs for first 8 elements (cofactors) */
  tmp[0] = src[10] * src[15];
  tmp[1] = src[11] * src[14];
  tmp[2] = src[9] * src[15];
  tmp[3] = src[11] * src[13]; 
  tmp[4] = src[9] * src[14]; 
  tmp[5] = src[10] * src[13];
  tmp[6] = src[8] * src[15];
  tmp[7] = src[11] * src[12];
  tmp[8] = src[8] * src[14];
  tmp[9] = src[10] * src[12];
  tmp[10] = src[8] * src[13];
  tmp[11] = src[9] * src[12];
  /* calculate first 8 elements (cofactors) */
  dst[0] = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7]; 
  dst[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];
  dst[1] = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7]; 
  dst[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7]; 
  dst[2] = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];
  dst[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7]; 
  dst[3] = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6]; 
  dst[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6]; 
  dst[4] = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3]; 
  dst[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3]; 
  dst[5] = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3]; 
  dst[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];
  dst[6] = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3]; 
  dst[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];
  dst[7] = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];
  dst[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];
  /* calculate pairs for second 8 elements (cofactors) */
  tmp[0] = src[2]*src[7]; 
  tmp[1] = src[3]*src[6];
  tmp[2] = src[1]*src[7];
  tmp[3] = src[3]*src[5]; 
  tmp[4] = src[1]*src[6];
  tmp[5] = src[2]*src[5];
  tmp[6] = src[0]*src[7];
  tmp[7] = src[3]*src[4];
  tmp[8] = src[0]*src[6];
  tmp[9] = src[2]*src[4];
  tmp[10] = src[0]*src[5];
  tmp[11] = src[1]*src[4];
  /* calculate second 8 elements (cofactors) */
  dst[8] = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15]; 
  dst[8] -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];
  dst[9] = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15]; 
  dst[9] -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15]; 
  dst[10] = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];
  dst[10]-= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15]; 
  dst[11] = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];
  dst[11]-= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14]; 
  dst[12] = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];
  dst[12]-= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10]; 
  dst[13] = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10]; 
  dst[13]-= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8]; 
  dst[14] = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8]; 
  dst[14]-= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9]; 
  dst[15] = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9]; 
  dst[15]-= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8];
  /* calculate determinant */
  det=src[0]*dst[0]+src[1]*dst[1]+src[2]*dst[2]+src[3]*dst[3];

 
  /* calculate matrix inverse */
  det = 1/det; 

  if(isnan(det))
    //    my_err("det in inverse 4x4 zero\n");
    return -1;

  for (j = 0; j < 16; j++)
    dst[j] *= det;

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	ia[i][j]= dst[i*4+j];
	if(isnan(ia[i][j])) return -1;
      }

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//Heavyside step function around 0
//1       /-------
//       | 
//0 ____/0
//x9 determines the sharpness and says where step function equals 0.9
ldouble
step_function(ldouble x,ldouble x9)
{
  ldouble k=1.47222/x9;
  return 1./(1.+exp(-2.*k*x));
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//prints error message and gets chars
int
my_err(char *message)
{
  char bufor[200];
  sprintf(bufor,"|err| : %s\n",message);
  printf("%s",bufor);
  getchar();
  return 0;  
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
int
getch()
{
  getchar();
  return 0;
}


