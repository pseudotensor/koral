//KORAL - physics.c
//magnetic field related routines

#include "ko.h"

/* calculate magnetic field four-vector knowing gas four-velocity ucov */
void calc_bcon_4vel(double *pr, double *ucon, double *ucov, double *bcon) 
{
  int j ;

  bcon[0] = pr[B1]*ucov[1] + pr[B2]*ucov[2] + pr[B3]*ucov[3] ;
  for(j=1;j<4;j++)
    {
      bcon[j] = (pr[B1-1+j] + bcon[0]*ucon[j])/ucon[0] ;
    }

  return ;
}

/* calculate magnetic field four-vector from primitives */
void calc_bcon_prim(double *pp, double *bcon, void* ggg) 
{
  int j;
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble ucon[4],ucov[4];
  ucon[0]=0;
  ucon[1]=pp[2];
  ucon[2]=pp[3];
  ucon[3]=pp[4];

  conv_vels(ucon,ucon,VELPRIM,VEL4,geom->gg,geom->GG);
  indices_21(ucon,ucov,geom->gg);

  calc_bcon_4vel(pp,ucon,ucov,bcon);

 

  return;
}

/* calculate B^i from bcon and primitives */
void calc_Bcon_prim(double *pp, double *bcon,double *Bcon, void* ggg) 
{
  int j;
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble ucon[4],ucov[4];
  ucon[0]=0;
  ucon[1]=pp[2];
  ucon[2]=pp[3];
  ucon[3]=pp[4];

  conv_vels(ucon,ucon,VELPRIM,VEL4,geom->gg,geom->GG);
  indices_21(ucon,ucov,geom->gg);

 
  

  Bcon[0]=0.;

  for(j=1;j<4;j++)
    {
      Bcon[j] = bcon[j]*ucon[0] - bcon[0]*ucon[j];
    }
  
  
  return;
}

/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/***********************************************************************************************/
/* wrappers for missing cells / dimensions */
int fl_x(int i)
{ if(NX==1) return 0;
  /*
  if(i<0) {
    #ifdef PERIODIC_XBC
    return NX+i;
    #else
    return 0;
    #endif 
  }
  if(i>=NX) {
    #ifdef PERIODIC_XBC
    return i-NX;
    #else
    return NX-1;
    #endif 
  }
  */
  return i; }

int fl_y(int i)
{ if(NY==1) return 0;
  /*
  if(i<0) {
    #ifdef PERIODIC_YBC
    return NY+i;
    #else
    return 0;
    #endif 
  }
  if(i>=NY) {
    #ifdef PERIODIC_YBC
    return i-NY;
    #else
    return NY-1;
    #endif 
  }
  */
  return i; }

int fl_z(int i)
{ if(NZ==1) return 0;
  /*
  if(i<0) {
    #ifdef PERIODIC_ZBC
    return NZ+i;
    #else
    return 0;
    #endif 
  }
  if(i>=NZ) {
    #ifdef PERIODIC_ZBC
    return i-NZ;
    #else
    return NZ-1;
    #endif 
  }
  */
  return i; }

int flx_x(int i)
{ return i; }
int flx_y(int i)
{ return fl_y(i); }
int flx_z(int i)
{ return fl_z(i); }

int fly_x(int i)
{ return fl_x(i); }
int fly_y(int i)
{ return i; }
int fly_z(int i)
{ return fl_z(i); }

int flz_x(int i)
{ return fl_x(i); }
int flz_y(int i)
{ return fl_y(i); }
int flz_z(int i)
{ return i; }


/***********************************************************************************************/
/***********************************************************************************************/
/* flux constrained transport */
int
flux_ct()
{
#ifdef MAGNFIELD

  //requires GDETIN = 1
  if(GDETIN==0)
    {
      my_err("MAGNFIELD requires GDETIN==1\n");
      exit(0);
    }

  //TOTH algorithm from HARM's fluxct.c
  ldouble coefemf[4];

  if((NY>1)&&(NZ>1)) coefemf[1]=0.25;
  else coefemf[1]=0.5; 
  if((NX>1)&&(NZ>1)) coefemf[2]=0.25;
  else coefemf[2]=0.5; 
  if((NX>1)&&(NY>1)) coefemf[3]=0.25;
  else coefemf[3]=0.5; 
  
  int ix,iy,iz,iv,ii;
  //#pragma omp parallel for private(ix,iy,iz,iv,ii) schedule (static)
  //calculating EMF at corners
  for(ii=0;ii<Nloop_4;ii++) 
    {
      ix=loop_4[ii][0];
      iy=loop_4[ii][1];
      iz=loop_4[ii][2]; 
      
      ////////////////////
      // EMF1
      ////////////////////
      
#if((NY>1)||(NZ>1))
      set_emf(1,ix,iy,iz,
	      coefemf[1] * (
                            #if(NY>1)
			    + get_ub(flby,B3,fly_x(ix),fly_y(iy),fly_z(iz),1) 
			    + get_ub(flby,B3,fly_x(ix),fly_y(iy),fly_z(iz-1),1)
                            #endif
                            #if(NZ>1)
			    - get_ub(flbz,B2,flz_x(ix),flz_y(iy),flz_z(iz),2) 
			    - get_ub(flbz,B2,flz_x(ix),flz_y(iy-1),flz_z(iz),2)
                            #endif
			    ));
#else  // end if doing EMF1
      set_emf(1,ix,iy,iz,0.); // not really 0, but differences in emf will be 0, and that's all that matters
#endif // end if not doing EMF1
      
	////////////////////
	// EMF2
	////////////////////
#if((NX>1)||(NZ>1))
      set_emf(2,ix,iy,iz,
	      coefemf[2] * (
                            #if(NZ>1)
			    + get_ub(flbz,B1,flz_x(ix),flz_y(iy),flz_z(iz),2) 
			    + get_ub(flbz,B1,flz_x(ix-1),flz_y(iy),flz_z(iz),2)
                            #endif
                            #if(NX>1)
			    - get_ub(flbx,B3,flx_x(ix),flx_y(iy),flx_z(iz),0) 
			    - get_ub(flbx,B3,flx_x(ix),flx_y(iy),flx_z(iz-1),0)
                            #endif
			    ));
#else  
      set_emf(2,ix,iy,iz,0.);
#endif 

	////////////////////
	// EMF3
	////////////////////
#if((NX>1)||(NY>1))
      set_emf(3,ix,iy,iz,
	      coefemf[3] * (
                            #if(NX>1)
			    + get_ub(flbx,B2,flx_x(ix),flx_y(iy),flx_z(iz),0) 
			    + get_ub(flbx,B2,flx_x(ix),flx_y(iy-1),flx_z(iz),0)
                            #endif
                            #if(NY>1)
			    - get_ub(flby,B1,fly_x(ix),fly_y(iy),fly_z(iz),1) 
			    - get_ub(flby,B1,fly_x(ix-1),fly_y(iy),fly_z(iz),1)
                            #endif
			    ));
#else  
      set_emf(3,ix,iy,iz,0.);
#endif
    }
  
  #pragma omp barrier

  //reset certain emfs at the boundaries to ensure stationarity
  adjust_fluxcttoth_emfs();

  //#pragma omp parallel for private(ix,iy,iz,iv,ii) schedule (static)
  for(ii=0;ii<Nloop_4;ii++) // 0...NX
    {
      ix=loop_4[ii][0];
      iy=loop_4[ii][1];
      iz=loop_4[ii][2]; 

      /////////////////////////////////////
      // F1
      ////////////////////////////////////
#if(NX>1)
	if(iy<NY && iz<NZ) //no need to fill x-face fluxes for iy=NY etc., 
	  {
	    //	    if(iy==0 && iz==0) printf("%d > %d\n",PROCID,ix);
	    set_ubx(flbx,B1,ix,iy,iz,0.);
	    set_ubx(flbx,B2,ix,iy,iz,0.5 * (get_emf(3,ix,iy,iz) + get_emf(3,ix,iy+1,iz)));
	    set_ubx(flbx,B3,ix,iy,iz,-0.5 * (get_emf(2,ix,iy,iz) + get_emf(2,ix,iy,iz+1)));
	  }
#endif

      /////////////////////////////////////
      // F2
      ////////////////////////////////////
#if(NY>1)
      if(ix<NX && iz<NZ)	
	{
	  set_uby(flby,B1,ix,iy,iz,-0.5 * (get_emf(3,ix,iy,iz) + get_emf(3,ix+1,iy,iz)));
	  set_uby(flby,B2,ix,iy,iz,0.);
	  set_uby(flby,B3,ix,iy,iz,0.5 * (get_emf(1,ix,iy,iz) + get_emf(1,ix,iy,iz+1)));
	}
#endif
      
			    
      /////////////////////////////////////
      // F3
      ////////////////////////////////////
#if(NZ>1)
	if(ix<NX && iy<NY)	
	{
	  set_ubz(flbz,B1,ix,iy,iz,0.5 * (get_emf(2,ix,iy,iz) + get_emf(2,ix+1,iy,iz)));
	  set_ubz(flbz,B2,ix,iy,iz,-0.5 * (get_emf(1,ix,iy,iz) + get_emf(1,ix,iy+1,iz)));
	  set_ubz(flbz,B3,ix,iy,iz,0.);
	}
#endif
    }
  
#endif //MAGNFIELD
      return 0;
}


// resets emf's (E_3) at the boundaries in x1-x2 plane to zero
int adjust_fluxcttoth_emfs()
{
  //don't know what to do here
  return 0;
}

//calculates magnetic field from vector potential given in pinput[B1..B3]
int
calc_BfromA(ldouble* pinput, int ifoverwrite)
{
  #ifdef MAGNFIELD

  int ix,iy,iz,iv,ii;
  
  //#pragma omp parallel for private(ix,iy,iz,iv,ii) schedule (static)
  for(ii=0;ii<Nloop_4;ii++) //all corners of the inner domain
    {      
      ix=loop_4[ii][0];
      iy=loop_4[ii][1];
      iz=loop_4[ii][2]; 

      //calculating A_i on corners by averaging neighbouring cell centers
      ldouble A[3];

      for(iv=0;iv<3;iv++)
	{
	  if(NY==1 && NZ==1)
	    A[iv]=1./2.*(get_u(pinput,B1+iv,ix,iy,iz) + get_u(pinput,B1+iv,ix-1,iy,iz));

	  if(NY>1 && NZ==1)
	    A[iv]=1./4.*(get_u(pinput,B1+iv,ix,iy,iz) + get_u(pinput,B1+iv,ix,iy-1,iz) + 
			 get_u(pinput,B1+iv,ix-1,iy,iz) + get_u(pinput,B1+iv,ix-1,iy-1,iz));

	  if(NZ>1 && NY==1)
	    A[iv]=1./4.*(get_u(pinput,B1+iv,ix,iy,iz) + get_u(pinput,B1+iv,ix,iy,iz-1) + 
			 get_u(pinput,B1+iv,ix-1,iy,iz) + get_u(pinput,B1+iv,ix-1,iy,iz-1));

	  if(NZ>1 && NY>1)
	    A[iv]=1./8.*(get_u(pinput,B1+iv,ix,iy,iz) + get_u(pinput,B1+iv,ix,iy-1,iz) + 
			 get_u(pinput,B1+iv,ix-1,iy,iz) + get_u(pinput,B1+iv,ix-1,iy-1,iz)
			 +get_u(pinput,B1+iv,ix,iy,iz-1) + get_u(pinput,B1+iv,ix,iy-1,iz-1) 
			 + get_u(pinput,B1+iv,ix-1,iy,iz-1) + get_u(pinput,B1+iv,ix-1,iy-1,iz-1));

	  //saving to pvecpot used inside calc_BfromA_core(); 
 	  set_u(pvecpot,B1+iv,ix,iy,iz,A[iv]);
	}
            
    } //cell loop
  
  //calculating curl and B
  //new components of B^i in pvecpot[1...3]
  calc_BfromA_core();
  
  //overwriting vector potential with magnetic fields (e.g., at init)  
  if(ifoverwrite)
    {
      //#pragma omp parallel for private(ix,iy,iz,iv,ii) schedule (static)
      for(ii=0;ii<Nloop_02;ii++) //domain and ghost cells
	{
	  ix=loop_02[ii][0];
	  iy=loop_02[ii][1];
	  iz=loop_02[ii][2]; 

	  struct geometry geom;
	  fill_geometry(ix,iy,iz,&geom);
      
	  ldouble pp[NV],uu[NV];
	  PLOOP(iv)
	    pp[iv]=get_u(p,iv,ix,iy,iz);

	  pp[B1]=get_u(pvecpot,1,ix,iy,iz);
	  pp[B2]=get_u(pvecpot,2,ix,iy,iz);
	  pp[B3]=get_u(pvecpot,3,ix,iy,iz);

	  p2u(pp,uu,&geom);

	  set_u(p,B1,ix,iy,iz,pp[B1]);
	  set_u(p,B2,ix,iy,iz,pp[B2]);
	  set_u(p,B3,ix,iy,iz,pp[B3]);
	  set_u(u,B1,ix,iy,iz,uu[B1]);
	  set_u(u,B2,ix,iy,iz,uu[B2]);
	  set_u(u,B3,ix,iy,iz,uu[B3]);     


	}
    }

 
#endif //MAGNFIELD

  return 0;
}

/***********************************************************************************************/
/** calculates B-field from A given on corners in B1-B3 primitives of pvecpot *******************/
//new components of B^i in pvecpot[1...3]
/***********************************************************************************************/
int
calc_BfromA_core()
{
  #ifdef MAGNFIELD

  int ix,iy,iz,iv,ii;
  
  if(NY==1 && NZ==1)
    {
      my_err("1D calc_BfromA_core() not implemented.\n"); exit(-1);
    }

  if(NY>1 && NZ>1)
    if(PROCID==0) printf("calc_BfromA_core(): warning: assumes d/dz A_i = 0. \n");


  //#pragma omp parallel for private(ix,iy,iz,iv,ii) schedule (static)
  for(ii=0;ii<Nloop_0;ii++) //domain + one layer
    {
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2]; 

      struct geometry geom;
      fill_geometry(ix,iy,iz,&geom);

      ldouble B[4];
      ldouble dA[4][4];
      dA[1][1]=dA[2][2]=dA[3][3]=0.;

     
      /* flux-ct-compatible, axisymmetric! */

      dA[1][2] = -(get_u(pvecpot,B1,ix,iy,iz) - get_u(pvecpot,B1,ix,iy+1,iz)
		   + get_u(pvecpot,B1,ix+1,iy,iz) - get_u(pvecpot,B1,ix+1,iy+1,iz))/(2.*get_size_x(iy,1)) ;
      dA[1][1] = -(get_u(pvecpot,B1,ix,iy,iz) + get_u(pvecpot,B1,ix,iy+1,iz)
		   - get_u(pvecpot,B1,ix+1,iy,iz) - get_u(pvecpot,B1,ix+1,iy+1,iz))/(2.*get_size_x(ix,0)) ;
      dA[1][3] = 0.;

      dA[2][2] = -(get_u(pvecpot,B2,ix,iy,iz) - get_u(pvecpot,B2,ix,iy+1,iz)
		   + get_u(pvecpot,B2,ix+1,iy,iz) - get_u(pvecpot,B2,ix+1,iy+1,iz))/(2.*get_size_x(iy,1)) ;
      dA[2][1] = -(get_u(pvecpot,B2,ix,iy,iz) + get_u(pvecpot,B2,ix,iy+1,iz)
		   - get_u(pvecpot,B2,ix+1,iy,iz) - get_u(pvecpot,B2,ix+1,iy+1,iz))/(2.*get_size_x(ix,0)) ;
      dA[2][3] = 0.;

      dA[3][2] = -(get_u(pvecpot,B3,ix,iy,iz) - get_u(pvecpot,B3,ix,iy+1,iz)
		   + get_u(pvecpot,B3,ix+1,iy,iz) - get_u(pvecpot,B3,ix+1,iy+1,iz))/(2.*get_size_x(iy,1)) ;
      dA[3][1] = -(get_u(pvecpot,B3,ix,iy,iz) + get_u(pvecpot,B3,ix,iy+1,iz)
		   - get_u(pvecpot,B3,ix+1,iy,iz) - get_u(pvecpot,B3,ix+1,iy+1,iz))/(2.*get_size_x(ix,0)) ;
      dA[3][3] = 0.;


      B[1]=(dA[2][3] - dA[3][2])/geom.gdet;
      B[2]=(dA[3][1] - dA[1][3])/geom.gdet;
      B[3]=(dA[1][2] - dA[2][1])/geom.gdet;

      set_u(pvecpot,1,ix,iy,iz,B[1]);
      set_u(pvecpot,2,ix,iy,iz,B[2]);
      set_u(pvecpot,3,ix,iy,iz,B[3]);
     
    } //cell loop
  
#endif //MAGNFIELD

  return 0;
}


/***********************************************************************************************/
/** calculates div B for given cell *****************************************************/
/***********************************************************************************************/
ldouble
calc_divB(int ix,int iy,int iz)
{
  //within domain:
  if(!if_indomain(ix,iy,iz)) return 0.; //do not calculate in ghost cells
  
  ldouble divB;
  
  //TODO: so far 2d only
  //this is corner based 
  divB = (pick_gdet(ix,iy,iz)*get_u(p,B1,ix,iy,iz) + pick_gdet(ix,iy-1,iz)*get_u(p,B1,ix,iy-1,iz) 
	  - pick_gdet(ix-1,iy,iz)*get_u(p,B1,ix-1,iy,iz) - pick_gdet(ix-1,iy-1,iz)*get_u(p,B1,ix-1,iy-1,iz))/(2.*(get_x(ix+1,0)-get_x(ix,0)))
    + (pick_gdet(ix,iy,iz)*get_u(p,B2,ix,iy,iz) + pick_gdet(ix-1,iy,iz)*get_u(p,B2,ix-1,iy,iz) 
       - pick_gdet(ix,iy-1,iz)*get_u(p,B2,ix,iy-1,iz) - pick_gdet(ix-1,iy-1,iz)*get_u(p,B2,ix-1,iy-1,iz))/(2.*(get_x(iy+1,1)-get_x(iy,1)));

  divB/=pick_gdet(ix,iy,iz);

  return divB;  
}

ldouble
calc_Qtheta(int ix, int iy, int iz)
{
  if(!doingavg)
    {
      struct geometry geom;
      fill_geometry(ix,iy,iz,&geom);

      ldouble rho=get_u(p,RHO,ix,iy,iz);
      ldouble bcon[4];
      calc_bcon_prim(&get_u(p,0,ix,iy,iz),bcon,&geom);
      ldouble ucon[4];
      ucon[1]=get_u(p,VX,ix,iy,iz);
      ucon[2]=get_u(p,VY,ix,iy,iz);
      ucon[3]=get_u(p,VZ,ix,iy,iz);
      conv_vels(ucon,ucon,VELPRIM,VEL4,geom.gg,geom.GG);
      ldouble Omega = ucon[3]/ucon[0];
      ldouble dx=get_xb(iy+1,1)-get_xb(iy,1);

      return 2.*M_PI/fabs(Omega)/dx*fabs(bcon[2])/sqrt(rho);
    }
  else //doingavg, in BL
    {
      struct geometry geomBL;
      fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

      ldouble bcon2 = get_uavg(pavg,AVGBCON(2),ix,iy,iz);
      ldouble Omega = get_uavg(pavg,AVGUCON(3),ix,iy,iz)/get_uavg(pavg,AVGUCON(0),ix,iy,iz);
      ldouble rho = get_uavg(pavg,RHO,ix,iy,iz);
      ldouble dx;
      ldouble xx1[4],xx2[4];
      xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
      xx2[0]=0.;xx2[1]=get_xb(ix,0);xx2[2]=get_xb(iy+1,1);xx2[3]=get_xb(iz,2);
      coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
      coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
      dx=fabs(xx2[2]-xx1[2]);

      return 2.*M_PI/fabs(Omega)/dx*fabs(bcon2)/sqrt(rho);
    }
}

//calculates sqrt(g_rr g_phph) b^r b^phi and b^2
int
calc_angle_brbphibsq(int ix, int iy, int iz, ldouble *brbphi, ldouble *bsq, ldouble *bcon,ldouble *bcov)
{
  int i;

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  struct geometry geomBL;
  fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

  ldouble pp[NV];

  if(doingavg)
    {
      *bsq = get_uavg(pavg,AVGBSQ,ix,iy,iz);
      ldouble bcon1bcon3 = 
	get_uavg(pavg,AVGBCONBCOV(1,0),ix,iy,iz)*geomBL.GG[3][0] + 
	get_uavg(pavg,AVGBCONBCOV(1,1),ix,iy,iz)*geomBL.GG[3][1] +
	get_uavg(pavg,AVGBCONBCOV(1,2),ix,iy,iz)*geomBL.GG[3][2] +
	get_uavg(pavg,AVGBCONBCOV(1,3),ix,iy,iz)*geomBL.GG[3][3];
      *brbphi = sqrt(geomBL.gg[1][1]*geomBL.gg[3][3])*bcon1bcon3;
    }
  else
    {
      PLOOP(i)
	pp[i]=get_u(p,i,ix,iy,iz);
      
      //to BL     
      trans_pmhd_coco(pp,pp,MYCOORDS,BLCOORDS,geom.xxvec,&geom,&geomBL);
      
      //b^mu
      calc_bcon_prim(pp, bcon, &geomBL);
      
      //b_mu
      indices_21(bcon,bcov,geomBL.gg); 

      *bsq = dot(bcon,bcov);
      *brbphi = sqrt(geomBL.gg[1][1]*geomBL.gg[3][3])*bcon[1]*bcon[3];
    }

  return 0;
}

//calculates b^p b^phi and b^2
//here cannot work on averages (bpoloidal bphi is not averaged)
int
calc_angle_bpbphibsq(int ix, int iy, int iz, ldouble *bpbphi, ldouble *bsq, ldouble *bcon, ldouble *bcov)
{
  int i;

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  struct geometry geomBL;
  fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

  ldouble pp[NV];
  PLOOP(i)
    pp[i]=get_u(p,i,ix,iy,iz);
	      
  //to BL     
  trans_pmhd_coco(pp,pp,MYCOORDS,BLCOORDS,geom.xxvec,&geom,&geomBL);
  
  //b^mu
  calc_bcon_prim(pp, bcon, &geomBL);
  
  //b_mu
  indices_21(bcon,bcov,geomBL.gg); 

  *bsq = dot(bcon,bcov);

  ldouble br = sqrt(geomBL.gg[1][1])*bcon[1];
  ldouble bth = sqrt(geomBL.gg[2][2])*bcon[2];
  ldouble bp = sqrt(br*br + bth*bth)*my_sign(bcon[1]);

  *bpbphi = bp*sqrt(geomBL.gg[3][3])*bcon[3];

  return 0;
}

//calculates curl of a 3-vector field from array ptu starting from index idx
//returns to curl[1..4]
//TODO: idx not used
int
calc_curl(ldouble *ptu, ldouble idx, int ix, int iy, int iz, void* ggg, ldouble *curl)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble dA[4][4];

  dA[1][1]=dA[2][2]=dA[3][3]=0.;

  //d_1 A_2
  if(NX==1)
    dA[1][2]=0.;
  else
    {
      if(ix>-NG && ix<NX+NG-1)
	dA[1][2]=(get_u(ptu,B2,ix+1,iy,iz)-get_u(ptu,B2,ix-1,iy,iz))/(get_x(ix+1,0)-get_x(ix-1,0));
      if(ix==-NG)
	dA[1][2]=(-3.*get_u(ptu,B2,ix,iy,iz)+4.*get_u(ptu,B2,ix+1,iy,iz)-get_u(ptu,B2,ix,iy,iz+2))/(get_x(ix+2,0)-get_x(ix,0));
      if(ix==NX+NG-1)
	dA[1][2]=(3.*get_u(ptu,B2,ix,iy,iz)-4.*get_u(ptu,B2,ix-1,iy,iz)+get_u(ptu,B2,ix-2,iy,iz))/(get_x(ix,0)-get_x(ix-2,0));
    }

  //d_1 A_3
  if(NX==1)
    dA[1][3]=0.;
  else
    {
      if(ix>-NG && ix<NX+NG-1)
	dA[1][3]=(get_u(ptu,B3,ix+1,iy,iz)-get_u(ptu,B3,ix-1,iy,iz))/(get_x(ix+1,0)-get_x(ix-1,0));
      if(ix==-NG)
	dA[1][3]=(-3.*get_u(ptu,B3,ix,iy,iz)+4.*get_u(ptu,B3,ix+1,iy,iz)-get_u(ptu,B3,ix,iy,iz+2))/(get_x(ix+2,0)-get_x(ix,0));
      if(ix==NX+NG-1)
	dA[1][3]=(3.*get_u(ptu,B3,ix,iy,iz)-4.*get_u(ptu,B3,ix-1,iy,iz)+get_u(ptu,B3,ix-2,iy,iz))/(get_x(ix,0)-get_x(ix-2,0));
    }

  //d_2 A_3
  if(NY==1)
    dA[2][3]=0.;
  else
    {
      if(iy>-NG && iy<NY+NG-1)
	dA[2][3]=(get_u(ptu,B3,ix,iy+1,iz)-get_u(ptu,B3,ix,iy-1,iz))/(get_x(iy+1,1)-get_x(iy-1,1));
      if(iy==-NG)
	dA[2][3]=(-3.*get_u(ptu,B3,ix,iy,iz)+4.*get_u(ptu,B3,ix,iy+1,iz)-get_u(ptu,B3,ix,iy+2,iz))/(get_x(iy+2,1)-get_x(iy,1));
      if(iy==NY+NG-1)
	dA[2][3]=(3.*get_u(ptu,B3,ix,iy,iz)-4.*get_u(ptu,B3,ix,iy-1,iz)+get_u(ptu,B3,ix,iy-2,iz))/(get_x(iy,1)-get_x(iy-2,1));
    }

  //d_2 A_1
  if(NY==1)
    dA[2][1]=0.;
  else
    {
      if(iy>-NG && iy<NY+NG-1)
	dA[2][1]=(get_u(ptu,B1,ix,iy+1,iz)-get_u(ptu,B1,ix,iy-1,iz))/(get_x(iy+1,1)-get_x(iy-1,1));
      if(iy==-NG)
	dA[2][1]=(-3.*get_u(ptu,B1,ix,iy,iz)+4.*get_u(ptu,B1,ix,iy+1,iz)-get_u(ptu,B1,ix,iy+2,iz))/(get_x(iy+2,1)-get_x(iy,1));
      if(iy==NY+NG-1)
	dA[2][1]=(3.*get_u(ptu,B1,ix,iy,iz)-4.*get_u(ptu,B1,ix,iy-1,iz)+get_u(ptu,B1,ix,iy-2,iz))/(get_x(iy,1)-get_x(iy-2,1));
    }

  //d_3 A_2
  if(NZ==1)
    dA[3][2]=0.;
  else
    {
      if(iz>-NG && iz<NZ+NG-1)
	dA[3][2]=(get_u(ptu,B2,ix,iy,iz+1)-get_u(ptu,B2,ix,iy,iz-1))/(get_x(iz+1,2)-get_x(iz-1,2));
      if(iz==-NG)
	dA[3][2]=(-3.*get_u(ptu,B2,ix,iy,iz)+4.*get_u(ptu,B2,ix,iy,iz+1)-get_u(ptu,B2,ix,iy,iz+2))/(get_x(iz+2,2)-get_x(iz,2));
      if(iz==NZ+NG-1)
	dA[3][2]=(3.*get_u(ptu,B2,ix,iy,iz)-4.*get_u(ptu,B2,ix,iy,iz-1)+get_u(ptu,B2,ix,iy,iz-2))/(get_x(iz,2)-get_x(iz-2,2));
    }

  //d_3 A_1
  if(NZ==1)
    dA[3][1]=0.;
  else
    {
      if(iz>-NG && iz<NZ+NG-1)
	dA[3][1]=(get_u(ptu,B1,ix,iy,iz+1)-get_u(ptu,B1,ix,iy,iz-1))/(get_x(iz+1,2)-get_x(iz-1,2));
      if(iz==-NG)
	dA[3][1]=(-3.*get_u(ptu,B1,ix,iy,iz)+4.*get_u(ptu,B1,ix,iy,iz+1)-get_u(ptu,B1,ix,iy,iz+2))/(get_x(iz+2,2)-get_x(iz,2));
      if(iz==NZ+NG-1)
	dA[3][1]=(3.*get_u(ptu,B1,ix,iy,iz)-4.*get_u(ptu,B1,ix,iy,iz-1)+get_u(ptu,B1,ix,iy,iz-2))/(get_x(iz,2)-get_x(iz-2,2));
    }

  //gdet B^i = d_j A_k eps^ijk

  curl[1]=(dA[2][3] - dA[3][2])/geom->gdet;
  curl[2]=(dA[3][1] - dA[1][3])/geom->gdet;
  curl[3]=(dA[1][2] - dA[2][1])/geom->gdet;

  return 0;
}

/***********************************************************************************************/
/************************* mimics alpha-dynamo in axisymmetric sims involvin MRI ***************/
/***********************************************************************************************/
int
mimic_dynamo(ldouble dtin)
{
#ifdef MAGNFIELD
#ifdef MIMICDYNAMO
#ifdef BHDISK_PROBLEMTYPE

  int ix,iy,iz,iv,ii;
  ldouble dt;
  int verbose=1,stop=0;
  int verix=75;
  int veriy=100;
  ldouble bcon[4],bcov[4],bsqin,bsqout,ugasin,ugasout;
  ugasin=ugasout=0.;


  //#pragma omp parallel for private(ix,iy,iz,iv,ii) schedule (static)
  for(ii=0;ii<Nloop_6;ii++) //inner domain plus 1-cell layer including corners
    {      
      ix=loop_6[ii][0];
      iy=loop_6[ii][1];
      iz=loop_6[ii][2]; 

      calc_primitives(ix,iy,iz,0,0);

      struct geometry geom;
      fill_geometry(ix,iy,iz,&geom);


      if(verbose)// && ix+TOI==verix && iy+TOJ==veriy)
	{
	  /*
	  printf("before dynamo: \n");
	  print_primitives(&get_u(p,0,ix,iy,iz));
	  print_conserved(&get_u(u,0,ix,iy,iz));
	  */

	  set_u(p_bak_subzone,UU,ix,iy,iz,0.);

	  calc_bcon_prim(&get_u(p,0,ix,iy,iz), bcon, &geom);
	  indices_21(bcon,bcov,geom.gg); 
	  bsqin = dot(bcon,bcov);

	  //ugasin=get_u(p,UU,ix,iy,iz);

	  //printf("bsqin: %e\n",bsqin);
	}


      //let's use ptemp1 to hold the extra vector potential
      set_u(ptemp1,B1,ix,iy,iz,0.);
      set_u(ptemp1,B2,ix,iy,iz,0.);
      set_u(ptemp1,B3,ix,iy,iz,0.);

      ldouble Aphi,Pk,Omk,Om;
      ldouble xxBL[4],bcon[4],bcov[4],Bp;
      ldouble ucon[4],ucov[4];
      int j;
      ldouble angle,bbphi,bsq;

      //magnetic field angle      
      calc_angle_brbphibsq(ix,iy,iz,&bbphi,&bsq,bcon,bcov);
      angle=-bbphi/bsq;

      //BL radius
      coco_N(geom.xxvec,xxBL,MYCOORDS, BLCOORDS);

      //to avoid BH
      if(xxBL[1]<1.0001*rhorizonBL) continue;

      

      //timescale
      Omk = 1./(BHSPIN+sqrt(xxBL[1]*xxBL[1]*xxBL[1]));
      Pk = 2.*M_PI/Omk;

      //local omega based
      //Pk = 2.*M_PI/Om;

      //angle
      ldouble facangle=0.;
      if(isfinite(angle))
	{
	  if(angle<-1.) angle=-1.;
	  facangle = my_max(0., (THETAANGLE-angle)/THETAANGLE);
	}

      //radius
      ldouble facradius = step_function(xxBL[1]-1.*rISCOBL,.1*rISCOBL);

      //pre(gas+rad)
      ldouble prermhd = GAMMAM1*get_u(p,UU,ix,iy,iz);
#ifdef RADIATION
      ldouble Rtt,Ehat,uconr[4],prad;
      calc_ff_Rtt(&get_u(p,0,ix,iy,iz),&Rtt,uconr,&geom);
      Ehat=-Rtt; 
      prad=Ehat/3.;
      prermhd+=prad;		
#endif
      
      //magnetic beta
      ldouble beta = bsq/2./prermhd;
      //bsq/rho 
      ldouble betarho = bsq/2./get_u(p,RHO,ix,iy,iz);

      ldouble facmag1 = step_function(1.-beta,0.1);      
      ldouble facmag2 = step_function(.1-betarho,.01);
      
    
#ifdef CALCHRONTHEGO //calculate HR on the go
      int gix=ix+TOI; //global radial index
      if(gix<0) gix=0; 
      if(gix>=TNX) gix=TNX-1;
      ldouble HRDTHETA = scaleth_otg[gix];
      //HRDTHETA *= 2.; //because scale height/photosphere = 1/3 --- 0.55
      if(HRDTHETA > 0.9*M_PI/2.) HRDTHETA=0.9*M_PI/2.; //to avoid the axis
#else      
      ldouble HRDTHETA = EXPECTEDHR * M_PI/2.;
#endif

      ldouble zH = (M_PI/2. - xxBL[2])/HRDTHETA;
      ldouble zHpow = 1.;
      ldouble faczH = my_max(0.,pow(1. - zH*zH,zHpow));

      
      //ldouble facmagnetization = my_min(facmag1,facmag2);		       
      //ldouble facmagnetization = my_min(faczH,my_min(facmag1,facmag2));
      ldouble facmagnetization = faczH;

      //the extra vector potential
      ldouble effalpha=ALPHADYNAMO;

      //dynamo proportional to vertical gravity ~ z
      #ifdef ALPHAFLIPSSIGN
      effalpha = - (M_PI/2. - xxBL[2])/(HRDTHETA/2.) * ALPHADYNAMO;  //2 to get average alpha = alphadynamo
      #endif

      //timestep
      dt=dtin;
      #ifdef SELFTIMESTEP
      dt=1./get_u_scalar(cell_tsteps,ix,iy,iz); //individual time step
      #endif

      ldouble Bphi=get_u(p,B3,ix,iy,iz);
      Aphi = effalpha * (HRDTHETA/(M_PI/2.)) / 0.4 
	* dt / Pk  * xxBL[1] * geom.gg[3][3] * Bphi
	* facradius 
	* facmagnetization 
	* facangle;

      //saving vector potential to ptemp1
      set_u(ptemp1,B3,ix,iy,iz,Aphi); 

//damping azimuthal component of magnetic field if beta exceeds DAMPBETA
#ifdef DAMPBETA      
      //local four-velocity
      /*     
      ucon[0]=0;
      ucon[1]=get_u(p,VX,ix,iy,iz);
      ucon[2]=get_u(p,VY,ix,iy,iz);
      ucon[3]=get_u(p,VZ,ix,iy,iz);
      conv_vels_both(ucon,ucon,ucov,VELPRIM,VEL4,geom.gg,geom.GG);
      
      ldouble bphi,dbphi,Bcon[4];
      
      bphi = bcon[3];

      //double range for damping
      zH/=2.;
      faczH = my_max(0.,pow(1. - zH*zH,zHpow));

      dbphi = - ALPHABETA 
	* facvel
	* facradius
	* faczH
	* dt / Pk 
	* my_max(0.,beta - BETASATURATED) / BETASATURATED 
	* bphi;

      if((dbphi+bphi)*bphi<0.) dbphi=-bphi; //not to overshoot zero

     
      bphi+=dbphi;
      bcon[3]=bphi;

      Bcon[1]=get_u(p,B1,ix,iy,iz);
      Bcon[2]=get_u(p,B2,ix,iy,iz);
      Bcon[3]=get_u(p,B3,ix,iy,iz);

      Bcon[3]=(bcon[3]*ucon[0]-(Bcon[1]*ucov[1]+Bcon[2]*ucov[2])*ucon[3])/
	(1.+ucon[3]*ucov[3]);

      ldouble Bphinew=Bcon[3];
      const ldouble maxBphichange=1.e-1;
      ldouble facdamp=1.;
      ldouble Bphichange=fabs((Bphinew-Bphi)/Bphi);

      if (Bphichange>maxBphichange)
        facdamp=maxBphichange/Bphichange;

      Bcon[3]=Bphi + (Bphinew-Bphi)*facdamp;

      set_u(p,B3,ix,iy,iz,Bcon[3]);
*/
      
      
      ldouble dBphi = - ALPHABETA 
	* facradius
	* faczH
	* dt / Pk 
	* my_max(0.,beta - BETASATURATED) / BETASATURATED 
	* Bphi;



      
      if((dBphi+Bphi)*Bphi<0.) dBphi=-Bphi; //not to overshoot zero 
      
      /*
      if(verbose && ix+TOI==verix && iy+TOJ==veriy)
	{
	  printf("dBphi: %e \n",dBphi);
	  if( dBphi!=0.) stop=1;
	}
      */
                                                   
      set_u(p,B3,ix,iy,iz,Bphi+dBphi);    
      
#endif    


    }

  //once the whole array is filled with cell centered A^phi we can 
  //calculate the extra magnetic field returned through pvecpot[1..3]
  calc_BfromA(ptemp1,0);  
   
  //and superimpose it on the original one
  //#pragma omp parallel for private(ix,iy,iz,iv,ii) schedule (static)
  for(ii=0;ii<Nloop_0;ii++) //domain only!
    {
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2];

      struct geometry geom;
      ldouble B[4]; ldouble xxBL[4];

      fill_geometry(ix,iy,iz,&geom);
      //BL radius
      coco_N(geom.xxvec,xxBL,MYCOORDS, BLCOORDS);

      B[1]=get_u(pvecpot,1,ix,iy,iz);
      B[2]=get_u(pvecpot,2,ix,iy,iz);
      B[3]=get_u(pvecpot,3,ix,iy,iz);

#ifdef PRESERVEBSQ
      //abandoned at 3/18/14
#endif

      set_u(p,B1,ix,iy,iz,get_u(p,B1,ix,iy,iz)+B[1]);
      set_u(p,B2,ix,iy,iz,get_u(p,B2,ix,iy,iz)+B[2]);
      set_u(p,B3,ix,iy,iz,get_u(p,B3,ix,iy,iz)+B[3]);

      
      if(verbose)// && ix+TOI==verix && iy+TOJ==veriy)
	{
	  //printf("dBvec = %e %e %e\n",B[1],B[2],B[3]);
	
	  ldouble pptemp[NV];
	  PLOOP(iv) pptemp[iv]=get_u(p,iv,ix,iy,iz);
	  pptemp[B1]+=B[1];
	  pptemp[B2]+=B[2];
	  pptemp[B3]+=B[3];
	  

	  calc_bcon_prim(&get_u(p,0,ix,iy,iz), bcon, &geom);
	  indices_21(bcon,bcov,geom.gg); 
	  bsqout = dot(bcon,bcov);

	  //printf("bsqout: %e\n",bsqout);

	  
	}

      ldouble uutemp[NV];
      p2u(&get_u(p,0,ix,iy,iz),uutemp,&geom);
      set_u(u,B1,ix,iy,iz,uutemp[B1]);
      set_u(u,B2,ix,iy,iz,uutemp[B2]);
      set_u(u,B3,ix,iy,iz,uutemp[B3]);

     
      if(verbose)// && ix+TOI==verix && iy+TOJ==veriy)
	{
	  /*
	  printf("after dynamo: \n");
	  print_primitives(&get_u(p,0,ix,iy,iz));
	  print_conserved(&get_u(u,0,ix,iy,iz));
	  */

	  ldouble pptemp[NV];
	  PLOOP(iv) pptemp[iv]=get_u(p,0,ix,iy,iz);
	  int corr[3],fixup[3];
	  u2p(&get_u(u,0,ix,iy,iz),pptemp,&geom,corr,fixup,0);
	  /*
	    printf("after inversion: \n");
	  print_primitives(pptemp);	
	  */

	  ugasout=pptemp[UU];			
	  ugasin=get_u(p,UU,ix,iy,iz);
	  
	  /*
	  printf("dbsq/2 = %e\n",(bsqout - bsqin)/2.);
	  printf("dugas = %e\n",ugasout - ugasin);
	  */

	  //if(stop) getch();

	  set_u(p_bak_subzone,UU,ix,iy,iz,ugasout-ugasin);
	  //if(iy==100) printf("%d > %e\n",ix,ugasout-ugasin);

	}

    }

  //done consistently only inside the inner domain
  //need for set_bc / exchange information afterwards or executed as the last operator
  //done
     

#endif
#endif
#endif

return 0;
}
