//KORAL - silo.c
//routines for writing a silo file with quadratic mesh
//used both on the go and separately

#include "ko.h"
#include <silo.h>
#include <string.h>

/*********************************************/
/* writes silo file in dumps
/*********************************************/
int fprint_silofile(ldouble time, int num, char* folder, char* prefix)
{
  char bufor[50];
  sprintf(bufor,"%s/%s%04d.silo",folder,prefix,num);
 
  DBfile *file = NULL;/* The Silo file pointer */
  char *coordnames[3];/* Names of the coordinates */
  ldouble *nodex;/* The coordinate arrays */
  ldouble *nodey;
  ldouble *nodez;
  ldouble *coordinates[3];/* The array of coordinatearrays */
  int dimensions[3];/* The number of nodes */
   
  /* Create the Silo file */
  file = DBCreate(bufor, DB_CLOBBER, DB_LOCAL, NULL,DB_PDB);

  /* Name the coordinate axes ‘X’ and ‘Y’ */
  coordnames[0] = strdup("X");
  coordnames[1] = strdup("Y");
  coordnames[2] = strdup("Z");
  
  /* Give the cartesian coordinates of the mesh */
  int ix,iy,iz,iv;
  int i,j;
  ldouble pp[NV],uu[NV],xxvec[4],xxveccar[4],xxvecsph[4],xx1[4],xx2[4];
  
  int nx=NX;
  int ny=NY;
  int nz=NZ;

#ifdef FULLPHI //printing one more cell in phi to close the sphere
  nz++;
#endif

  nodex=(ldouble *)malloc(nx*ny*nz*sizeof(ldouble));
  nodey=(ldouble *)malloc(nx*ny*nz*sizeof(ldouble));
  nodez=(ldouble *)malloc(nx*ny*nz*sizeof(ldouble));

  ldouble *rho = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *uint = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *temp = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *forcebal1 = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *forcebal2 = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Omega = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *muBe = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Qtheta = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *divB = (ldouble*)malloc(nx*ny*nz*sizeof(double));

  

  #ifdef TRACER
  ldouble *tracer = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  #endif 

  ldouble *vx = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *vy = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *vz = (ldouble*)malloc(nx*ny*nz*sizeof(double));

  ldouble *Edotx = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Edoty = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Edotz = (ldouble*)malloc(nx*ny*nz*sizeof(double));


  #ifdef MAGNFIELD
  ldouble *bsq = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Bx = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *By = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Bz = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *phi = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  #endif

  #ifdef RADIATION
  ldouble *tautot = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *tauabs = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *forcebal3 = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Erad = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Ehat = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Fx = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Fy = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Fz = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *uradx = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *urady = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *uradz = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  #endif

  for(iz=0;iz<nz;iz++)
    {
      
      
      for(iy=0;iy<ny;iy++)
	{
	  
#ifdef PRINTXGC_RIGHT
	  for(ix=NG;ix<nx+NG;ix++)
#else
	  for(ix=0;ix<nx;ix++)
#endif

	    {
	      int iix,iiy,iiz;
	      iix=ix;
	      iiy=iy;
	      iiz=iz;
	      if(iiz>=NZ) iiz-=NZ;


	      struct geometry geom;
	      fill_geometry(iix,iiy,iiz,&geom);
	      struct geometry geomout;
	      fill_geometry_arb(iix,iiy,iiz,&geomout,OUTCOORDS);
	      
	      ldouble dxph[3],dx[3];
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
	      dxph[0]=dx[0]*sqrt(geomout.gg[1][1]);
	      dxph[1]=dx[1]*sqrt(geomout.gg[2][2]);
	      dxph[2]=dx[2]*sqrt(geomout.gg[3][3]);

	      get_xx(iix,iiy,iiz,xxvec);
	      coco_N(xxvec,xxvecsph,MYCOORDS,SPHCOORDS);
	      coco_N(xxvec,xxveccar,MYCOORDS,MINKCOORDS);
	      ldouble r=xxvecsph[1];
	      ldouble th=xxvecsph[2];
	      ldouble ph=xxvecsph[3];



	      //gdet and coordinates of cells +- 1 in radius
	      ldouble gdet1,gdet2,gdet,gdetu;
	      gdet=geomout.gdet;
	      gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
	      gdetu=1.;
#endif

	      get_xx(iix-1,iiy,iiz,xx1);
	      coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
	      gdet1=calc_gdet_arb(xx1,BLCOORDS);

	      get_xx(iix+1,iiy,iiz,xx2);
	      coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
	      gdet2=calc_gdet_arb(xx2,BLCOORDS);

	      if(OUTCOORDS==BLCOORDS && geomout.xx<r_horizon_BL(BHSPIN))
	      continue;

#ifdef PRINTXGC_RIGHT
	      int nodalindex=iz*(ny*nx) + iy*nx + ix-NG;
#else
	      int nodalindex=iz*(ny*nx) + iy*nx + ix;
#endif
	      for(iv=0;iv<NV;iv++)
		{
		  if(doingavg)
		    pp[iv]=get_uavg(pavg,iv,ix,iy,iz); //this should not be used later but it is
		  else
		    pp[iv]=get_u(p,iv,iix,iiy,iiz);
		}


	      //coordinates
	      nodex[nodalindex]=xxveccar[1];
	      nodey[nodalindex]=xxveccar[2];
	      nodez[nodalindex]=xxveccar[3];

	      //primitives to OUTCOORDS
              #ifdef RADIATION
	      trans_prad_coco(pp, pp, MYCOORDS,OUTCOORDS, xxvec,&geom,&geomout);
              #endif
	      trans_pmhd_coco(pp, pp, MYCOORDS,OUTCOORDS, xxvec,&geom,&geomout);

	      //magnetic fields
#ifdef MAGNFIELD
	      ldouble bcon[4],bcov[4];
	      if(doingavg==0)
		{
		  calc_bcon_prim(pp,bcon,&geomout);
		  indices_21(bcon,bcov,geomout.gg); 
		  bsq[nodalindex] = dot(bcon,bcov);
		}
	      else
		{
		  bcon[1]=get_uavg(pavg,AVGBCON(1),ix,iy,iz);
		  bcon[2]=get_uavg(pavg,AVGBCON(2),ix,iy,iz);
		  bcon[3]=get_uavg(pavg,AVGBCON(3),ix,iy,iz);
		  bsq[nodalindex]=get_uavg(pavg,AVGBSQ,ix,iy,iz);
		}
#endif
	      
	      //velocities etc
	      ldouble vel[4],vcov[4],vcon[4],velprim[4];
	      ldouble Tit[4],Tij[4][4];
	      ldouble dpdr; //d/dr (gdet * p)
	      ldouble gracen; //gdet T^k_l Gamma^l_kr
	      ldouble w;//entalphy

	      if(doingavg==0) //using snapshot date
		{
		  rho[nodalindex]=pp[RHO];
		  uint[nodalindex]=pp[UU];
		  vel[1]=pp[VX];
		  vel[2]=pp[VY];
		  vel[3]=pp[VZ];
		  
		  conv_vels(vel,vel,VELPRIM,VEL4,geomout.gg,geomout.GG);
		  for(i=0;i<4;i++) vcon[i]=vel[i];
		  indices_21(vel,vcov,geomout.gg);
 
		  Omega[nodalindex]=vel[3]/vel[0];

		  calc_Tij(pp,&geomout,Tij);
		  indices_2221(Tij,Tij,geomout.gg);

		  Tit[1]=Tij[1][0];
		  Tit[2]=Tij[2][0];
		  Tit[3]=Tij[3][0];

		  //muBe[nodalindex]=Tij[1][0]/(rho[nodalindex]*vel[1])-1.;
		  #ifdef MAGNFIELD
		  muBe[nodalindex]=-(rho[nodalindex]*vcov[0]+
				     GAMMA*uint[nodalindex]*vcov[0]+
				     bsq[nodalindex])/rho[nodalindex]-1.;


		  Qtheta[nodalindex]=2.*M_PI/Omega[nodalindex]/dx[1]*fabs(bcon[2])/sqrt(rho[nodalindex]);
		  divB[nodalindex]=calc_divB(ix,iy,iz);
                  #endif
		  

		  dpdr = (gdet2*GAMMA*get_u(p,UU,iix+1,iiy,iiz)-gdet1*GAMMA*get_u(p,UU,iix-1,iiy,iiz)) / (xx2[1]-xx1[1]);
		  gracen=0.;
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      gracen += gdet*Tij[i][j]*get_gKr(j,1,i,ix,iy,iz);

		  forcebal1[nodalindex]=-dpdr;
		  forcebal2[nodalindex]=gracen;
		}
	      else //using averaged data
		{
		  rho[nodalindex]=get_uavg(pavg,RHO,ix,iy,iz);
		  uint[nodalindex]=get_uavg(pavg,UU,ix,iy,iz);
		 
		  vel[0]=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  vel[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  vel[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  vel[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  for(i=0;i<4;i++) vcon[i]=vel[i];
		  indices_21(vel,vcov,geomout.gg);
 
		  conv_vels_ut(vel,velprim,VEL4,VELPRIM,geomout.gg,geomout.GG);
		 
		  Omega[nodalindex]=vel[3]/vel[0];

		  pp[VX]=vel[1]; //updates pp[VI] to have rho-weighted velocities there
		  pp[VY]=vel[2];
		  pp[VZ]=vel[3];

		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Tij[i][j]=get_uavg(pavg,AVGRHOUCONUCOV(1,0),ix,iy,iz)
			+ GAMMA*get_uavg(pavg,AVGUUUCONUCOV(1,0),ix,iy,iz)
			+ get_uavg(pavg,AVGBSQUCONUCOV(1,0),ix,iy,iz)
			+ delta(i,j)*(GAMMA*get_uavg(pavg,UU,ix,iy,iz) + 1./2.*get_uavg(pavg,AVGBSQ,ix,iy,iz))
			- get_uavg(pavg,AVGBCONBCOV(1,0),ix,iy,iz); 

		  
		  Tit[1]=Tij[1][0];
		  Tit[2]=Tij[2][0];
		  Tit[3]=Tij[3][0];

		  #ifdef MAGNFIELD
		  //		  muBe[nodalindex]=Tij[1][0]/(rho[nodalindex]*vel[1])-1.;
		  muBe[nodalindex]=-(rho[nodalindex]*vcov[0]+
				     GAMMA*uint[nodalindex]*vcov[0]+
				     bsq[nodalindex])/rho[nodalindex]-1.;

		  Qtheta[nodalindex]=2.*M_PI/Omega[nodalindex]/dx[1]*fabs(bcon[2])/sqrt(rho[nodalindex]);
		  divB[nodalindex]=calc_divB(ix,iy,iz);
		  #endif

		  dpdr = (gdet2*GAMMA*get_uavg(pavg,UU,iix+1,iiy,iiz)-gdet1*GAMMA*get_uavg(pavg,UU,iix-1,iiy,iiz)) / (xx2[1]-xx1[1]);
		  gracen=0.;
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      gracen += gdet*Tij[i][j]*get_gKr(j,1,i,ix,iy,iz);

		  forcebal1[nodalindex]=-dpdr;
		  forcebal2[nodalindex]=gracen;
		}

	     
	      //fdef CGSOUTPUT
	      //rho[nodalindex]=rhoGU2CGS(rho[nodalindex]);
	      //uint[nodalindex]=endenGU2CGS(uint[nodalindex]);
	      //ndif

	      ldouble temploc=calc_PEQ_Tfromurho(uint[nodalindex],rho[nodalindex]);
	      temp[nodalindex]=temploc;
	      
	      #ifdef TRACER
	      tracer[nodalindex]=pp[TRA];
	      #endif

	      //default, but can be non-ortonormal VEL4
	      vx[nodalindex]=vel[1];
	      vy[nodalindex]=vel[2];
	      vz[nodalindex]=vel[3];

	      Edotx[nodalindex]=Tit[1];
	      Edoty[nodalindex]=Tit[2];
	      Edotz[nodalindex]=Tit[3];

	      //transform to cartesian
	      if (MYCOORDS==SCHWCOORDS || MYCOORDS==KSCOORDS || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS || MYCOORDS==MKS1COORDS || MYCOORDS==MKS2COORDS)
		{
		  vel[2]*=r;
		  vel[3]*=r*sin(th);
		  
		  vx[nodalindex] = sin(th)*cos(ph)*vel[1] 
		    + cos(th)*cos(ph)*vel[2]
		    - sin(ph)*vel[3];

		  vy[nodalindex] = sin(th)*sin(ph)*vel[1] 
		    + cos(th)*sin(ph)*vel[2]
		    + cos(ph)*vel[3];

		  vz[nodalindex] = cos(th)*vel[1] 
		    - sin(th)*vel[2];

		  Tit[2]*=r;
		  Tit[3]*=r*sin(th);
		  
		  Edotx[nodalindex] = sin(th)*cos(ph)*Tit[1] 
		    + cos(th)*cos(ph)*Tit[2]
		    - sin(ph)*Tit[3];

		  Edoty[nodalindex] = sin(th)*sin(ph)*Tit[1] 
		    + cos(th)*sin(ph)*Tit[2]
		    + cos(ph)*Tit[3];

		  Edotz[nodalindex] = cos(th)*Tit[1] 
		    - sin(th)*Tit[2];
		}
		

	      #ifdef MAGNFIELD
	      //magnetic field
	      

	      	      
	      //to ortonormal	      
	      //trans2_cc2on(bcon,bcon,geomout.tup);

	      Bx[nodalindex]=bcon[1];
	      By[nodalindex]=bcon[2];
	      Bz[nodalindex]=bcon[3];

	      if(iy==0)
		{
		  //phi[nodalindex]=geomout.gdet*pp[B1]*dx[1];//*get_size_x(iy,1);
		  phi[nodalindex]=geom.gdet*get_u(p,B1,ix,iy,iz)*get_size_x(iy,1)*2.*M_PI;
		}
	      else
		{
#ifdef PRINTXGC_RIGHT
		  int idx=iz*(ny*nx) + (iy-1)*nx + ix-NG;
#else
		  int idx=iz*(ny*nx) + (iy-1)*nx + ix;
#endif
		  //phi[nodalindex]=phi[idx]+geomout.gdet*pp[B1]*dx[1];//*get_size_x(iy,1);
		  phi[nodalindex]=phi[idx]+geom.gdet*get_u(p,B1,ix,iy,iz)*get_size_x(iy,1)*2.*M_PI;
		}


	      //transform to cartesian
	      if (MYCOORDS==SCHWCOORDS || MYCOORDS==KSCOORDS || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS || MYCOORDS==MKS1COORDS || MYCOORDS==MKS2COORDS)
		{
		  bcon[2]*=r;
		  bcon[3]*=r*sin(th);

		  Bx[nodalindex] = sin(th)*cos(ph)*bcon[1] 
		    + cos(th)*cos(ph)*bcon[2]
		    - sin(ph)*bcon[3];

		  By[nodalindex] = sin(th)*sin(ph)*bcon[1] 
		    + cos(th)*sin(ph)*bcon[2]
		    + cos(ph)*bcon[3];

		  Bz[nodalindex] = cos(th)*bcon[1] 
		    - sin(th)*bcon[2];
		}
	      #endif

	      #ifdef RADIATION

	      ldouble Rtt,ehat,ugas[4],urad[4],rvel[4],Rij[4][4],Gi[4];

	      ldouble tauabsloc = vcon[0]*calc_kappa(rho[nodalindex],temploc,geomout.xx,geomout.yy,geomout.zz);
	      ldouble tautotloc = vcon[0]*calc_kappaes(rho[nodalindex],temploc,geomout.xx,geomout.yy,geomout.zz);

	      if(doingavg==0)
		{

		  calc_ff_Rtt(pp,&Rtt,ugas,&geomout);
		  ehat=-Rtt;  	      							  
		  //prad_lab2on(pp,pp,&geomout);
		  //rvel[1]=pp[FX0];
		  //rvel[2]=pp[FY0];
		  //rvel[3]=pp[FZ0];
		  //conv_vels(rvel,rvel,VELPRIM,VEL4,geomout.gg,geomout.GG);
		  calc_Rij(pp,&geomout,Rij); //calculates R^munu in OUTCOORDS
		  indices_2221(Rij,Rij,geomout.gg);

		  
		}
	      else
		{
		  ehat=get_uavg(pavg,AVGEHAT,ix,iy,iz);
		  /*
		  ldouble Rtop[4];
		  Rtop[0]=get_uavg(pavg,AVGRIJ(0,0),ix,iy,iz);
		  Rtop[1]=get_uavg(pavg,AVGRIJ(0,1),ix,iy,iz);
		  Rtop[2]=get_uavg(pavg,AVGRIJ(0,2),ix,iy,iz);
		  Rtop[3]=get_uavg(pavg,AVGRIJ(0,3),ix,iy,iz);
		  */
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
		  
		}

	      //correcting rad-velocities basing on <R^t_mu>
	      int radcorr;
	      //print_primitives(pp);
	      
	      p2u(pp,uu,&geomout);
	      //print_conserved(uu);
	      uu[EE0]=gdetu*Rij[0][0];
	      uu[FX0]=gdetu*Rij[0][1];
	      uu[FY0]=gdetu*Rij[0][2];
	      uu[FZ0]=gdetu*Rij[0][3];
	      //print_conserved(uu);
	      u2p_rad(uu,pp,&geomout,&radcorr);
	      
	      //print_primitives(pp);getchar();
	      
	      //four fource
	      calc_Gi(pp,&geomout,Gi); 
	      indices_21(Gi,Gi,geomout.gg);
	      
	      Ehat[nodalindex]=ehat;
	      Erad[nodalindex]=Rij[0][0];

	      Fx[nodalindex]=Rij[1][0];
	      Fy[nodalindex]=Rij[2][0];
	      Fz[nodalindex]=Rij[3][0];

	      urad[1]=pp[FX0];
	      urad[2]=pp[FY0];
	      urad[3]=pp[FZ0];
	      conv_vels(urad,urad,VELPRIM,VEL4,geomout.gg,geomout.GG);

	      uradx[nodalindex]=urad[1];
	      urady[nodalindex]=urad[2];
	      uradz[nodalindex]=urad[3];

	      forcebal3[nodalindex]=gdet*Gi[1];
	      	
	      if(iy==0)
		{
		  tautot[nodalindex]=tautotloc*dxph[1];//get_size_x(iy,1)*sqrt(geomout.gg[2][2]);
		  tauabs[nodalindex]=tauabsloc*dxph[1];//get_size_x(iy,1)*sqrt(geomout.gg[2][2]);
		}
	      else
		{
#ifdef PRINTXGC_RIGHT
		  int idx=iz*(ny*nx) + (iy-1)*nx + ix-NG;
#else
		  int idx=iz*(ny*nx) + (iy-1)*nx + ix;
#endif
		  if(iy<=NY/2) //proper integration only in the upper half
		    {
		      tautot[nodalindex]=tautot[idx]+tautotloc*dxph[1];
		      tauabs[nodalindex]=tauabs[idx]+tauabsloc*dxph[1];
		    }
		  else
		    {
		      idx=iz*(ny*nx) + (NY-iy-1)*nx + ix-NG;
		      tautot[nodalindex]=tautot[idx];
		      tauabs[nodalindex]=tauabs[idx];
		    }
		}

	      //transform to cartesian
	      if (MYCOORDS==SCHWCOORDS || MYCOORDS==KSCOORDS || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS || MYCOORDS==MKS1COORDS || MYCOORDS==MKS2COORDS)
		{
		  Rij[2][0]*=r;
		  Rij[3][0]*=r*sin(th);

		  Fx[nodalindex] = sin(th)*cos(ph)*Rij[1][0] 
		    + cos(th)*cos(ph)*Rij[2][0]
		    - sin(ph)*Rij[3][0];

		  Fy[nodalindex] = sin(th)*sin(ph)*Rij[1][0] 
		    + cos(th)*sin(ph)*Rij[2][0]
		    + cos(ph)*Rij[3][0];

		  Fz[nodalindex] = cos(th)*Rij[1][0] 
		    - sin(th)*Rij[2][0];

		  urad[2]*=r;
		  urad[3]*=r*sin(th);

		  uradx[nodalindex] = sin(th)*cos(ph)*urad[1]
		    + cos(th)*cos(ph)*urad[2]
		    - sin(ph)*urad[3];

		  urady[nodalindex] = sin(th)*sin(ph)*urad[1]
		    + cos(th)*sin(ph)*urad[2]
		    + cos(ph)*urad[3];

		  uradz[nodalindex] = cos(th)*urad[1]
		    - sin(th)*urad[2];
		}

	      if(iy==NY/2 && 0)
		{
		  printf("%d %d %d | %e %e %e | %e %e %e | %e %e %e | %e %e %e\n",ix,iy,iz,
			 xxvec[1],xxvec[2],xxvec[3],
			 xxveccar[1],xxveccar[2],xxveccar[3],
			 //			 urad[1],urad[2],urad[3],
			 //Rij[0][1],Rij[0][2],Rij[0][3]
			 uradx[nodalindex],urady[nodalindex],uradz[nodalindex],
			 Fx[nodalindex],Fy[nodalindex],Fz[nodalindex]
			 );getchar();
		}

	      #endif
	  
	      
	      
	    }
	}
    }
  
  /* assign grid */
  int ndim;
  if(ny==1 && nz==1) //1d
    {
      ndim=1;
      dimensions[0]=nx;
      coordinates[0]=nodex;
    }
  else if(nz==1) //2d
    {
      ndim=2;

      dimensions[0] = nx;
      dimensions[1] = ny;

      coordinates[0] = nodex;
      coordinates[1] = nodey;
#ifdef SILO2D_XZPLANE
      coordinates[1] = nodez;
#endif
    }
  else if(ny==1) //2d, switch order
    {
      ndim=2;
      
      /* How many nodes in each direction? */
      dimensions[0] = nx;
      dimensions[1] = nz;

      /* Assign coordinates to coordinates array */
      coordinates[0] = nodex;
      coordinates[1] = nodez; 
    }
  else //3d
    {
      ndim=3;
      
      /* How many nodes in each direction? */
      dimensions[0] = nx;
      dimensions[1] = ny;
      dimensions[2] = nz;

      /* Assign coordinates to coordinates array */
      coordinates[0] = nodex;
      coordinates[1] = nodey;
      coordinates[2] = nodez;
    }      
     
  /* Write out the mesh to the file */
  DBPutQuadmesh(file, "mesh1", coordnames, coordinates,
  		dimensions, ndim, DB_DOUBLE, DB_NONCOLLINEAR, NULL);

  /* Write scalars */
  DBoptlist *optList = DBMakeOptlist(1);
  DBAddOption(optList, DBOPT_DTIME, (void*)&time);
  DBPutQuadvar1(file, "rho","mesh1", rho,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);

  DBPutQuadvar1(file, "uint","mesh1", uint,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);

  DBPutQuadvar1(file, "temp","mesh1", temp,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);

  DBPutQuadvar1(file, "omega","mesh1", Omega,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);

  DBPutQuadvar1(file, "muBe","mesh1", muBe,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);

  DBPutQuadvar1(file, "Qtheta","mesh1", Qtheta,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);

  DBPutQuadvar1(file, "divB","mesh1", divB,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);

  DBPutQuadvar1(file, "forcebal1","mesh1", forcebal1,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);

  DBPutQuadvar1(file, "forcebal2","mesh1", forcebal2,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);

  #ifdef TRACER
  DBPutQuadvar1(file, "tracer","mesh1", tracer,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);
  #endif

  #ifdef RADIATION
  DBPutQuadvar1(file, "erad","mesh1", Erad,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);
  DBPutQuadvar1(file, "ehat","mesh1", Ehat,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);
  DBPutQuadvar1(file, "tautot","mesh1", tautot,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);
  DBPutQuadvar1(file, "tauabs","mesh1", tauabs,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);
  DBPutQuadvar1(file, "forcebal3","mesh1", forcebal3,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);
  #endif

  #ifdef MAGNFIELD
  DBPutQuadvar1(file, "bsq","mesh1", bsq,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);
  DBPutQuadvar1(file, "phi","mesh1", phi,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);
  #endif


  /* Write vectors */
  optList = DBMakeOptlist(1);
  DBAddOption(optList, DBOPT_DTIME, (void*)&time);
  char *names[3];  
  ldouble *vels[3];

  //velocity
  vels[0]=vx;
  vels[1]=vy;
  vels[2]=vz;
#ifdef SILO2D_XZPLANE
  vels[1]=vz;
#endif

  names[0] = strdup("vel1");
  names[1] = strdup("vel2");
  names[2] = strdup("vel3");
  DBPutQuadvar(file, "velocity","mesh1", 3, names, vels, 
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);

  //en. flux
  vels[0]=Edotx;
  vels[1]=Edoty;
  vels[2]=Edotz;
#ifdef SILO2D_XZPLANE
  vels[1]=Edotz;
#endif

  names[0] = strdup("Tit1");
  names[1] = strdup("Tit2");
  names[2] = strdup("Tit3");
  DBPutQuadvar(file, "en_flux","mesh1", 3, names, vels, 
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);

  #ifdef MAGNFIELD
  //magn field
  vels[0]=Bx;
  vels[1]=By;
  vels[2]=Bz;
#ifdef SILO2D_XZPLANE
  vels[1]=Bz;
#endif
  names[0] = strdup("B1");
  names[1] = strdup("B2");
  names[2] = strdup("B3");
  DBPutQuadvar(file, "magn_field","mesh1", 3, names, vels, 
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);
  #endif

  #ifdef RADIATION 
  //radiative flux
  vels[0]=Fx;
  vels[1]=Fy;
  vels[2]=Fz;
#ifdef SILO2D_XZPLANE
  vels[1]=Fz;
#endif
  names[0] = strdup("F1");
  names[1] = strdup("F2");
  names[2] = strdup("F3");
  DBPutQuadvar(file, "rad_flux","mesh1", 3, names, vels, 
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);

  //rad rest frame velocity
  vels[0]=uradx;
  vels[1]=urady;
  vels[2]=uradz;
#ifdef SILO2D_XZPLANE
  vels[1]=uradz;
#endif
  names[0] = strdup("urad1");
  names[1] = strdup("urad2");
  names[2] = strdup("urad3");
  DBPutQuadvar(file, "rad_vel","mesh1", 3, names, vels, 
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);
  #endif
 
  /* Close the Silo file */
  DBClose(file);

  free(nodex);
  free(nodey);
  free(nodez);

  free(rho);
  free(uint);
  free(temp);
  #ifdef TRACER
  free(tracer);
  #endif

  #ifdef RADIATION
  free(Erad);
  free(Ehat);
  free(Fx);
  free(Fy);
  free(Fz);
  #endif

  #ifdef MAGNFIELD
  free(bsq);
  free(Bx);
  free(By);
  free(Bz);
  #endif
 
  
  free(vx);
  free(vy);
  free(vz);

  return (0);
}
