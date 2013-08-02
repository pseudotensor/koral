//KORAL - physics.c
//magnetic field related routines

#include "ko.h"

/* calculate magnetic field four-vector */
void bcon_calc(double *pr, double *ucon, double *ucov, double *bcon) 
{
  int j ;

  bcon[0] = pr[B1]*ucov[1] + pr[B2]*ucov[2] + pr[B3]*ucov[3] ;
  for(j=1;j<4;j++)
    {
      bcon[j] = (pr[B1-1+j] + bcon[0]*ucon[j])/ucon[0] ;
      //printf("> %d %e %e %e %e\n",j,bcon[j],pr[B1-1+j],bcon[0]*ucon[j],ucon[j]);
    }
  return ;
}

/***********************************************************************************************/
/***********************************************************************************************
/***********************************************************************************************/
/***********************************************************************************************
/* wrappers for missing cells / dimensions */
int fl_x(int i)
{ if(NX==1) return 0;
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
  return i; }

int fl_y(int i)
{ if(NY==1) return 0;
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
  return i; }

int fl_z(int i)
{ if(NZ==1) return 0;
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
/***********************************************************************************************
/* flux constrained transport */
int
flux_ct()
{
#ifdef MAGNFIELD
  //TOTH algorithm from HARM's fluxct.c
  ldouble coefemf[4];

  if((NY>1)&&(NZ>1)) coefemf[1]=0.25;
  else coefemf[1]=0.5; 
  if((NX>1)&&(NZ>1)) coefemf[2]=0.25;
  else coefemf[2]=0.5; 
  if((NX>1)&&(NY>1)) coefemf[3]=0.25;
  else coefemf[3]=0.5; 
  
  int ix,iy,iz,iv,ii;
#pragma omp parallel for private(ix,iy,iz,iv) schedule (static)
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
      //TODO: debug here
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
  
  //reset certain emfs at the boundaries to ensure stationarity
  adjust_fluxcttoth_emfs();

#pragma omp parallel for private(ix,iy,iz,iv) schedule (static)
  for(ii=0;ii<Nloop_4;ii++) // 0...NX
    {
      ix=loop_4[ii][0];
      iy=loop_4[ii][1];
      iz=loop_4[ii][2]; 

      /////////////////////////////////////
      // F1
      ////////////////////////////////////
#if(NX>1)
      if(iy<NY && iz<NZ) //no need to fill x-face fluxes for iy=NY etc.
	{
	  set_ubx(flbx,B1,ix,iy,iz,0.);
	  set_ubx(flbx,B2,ix,iy,iz,0.5 * (get_emf(3,ix,iy,iz) + get_emf(3,ix,iy+1,iz)));
	  set_ubx(flbx,B3,ix,iy,iz,-0.5 * (get_emf(2,ix,iy,iz) + get_emf(2,ix,iy,iz+1)));
	}
#endif

      /////////////////////////////////////
      // F2
      ////////////////////////////////////
#if(NY>1)
      if(ix<NX && iz<NZ && 1)	
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

/***********************************************************************************************/
/** calculates B-field from A given in B1-B3 primitives ****************************************/
/** uses pinit as a temporary holder ***********************************************************/
/** assumes set_bc already done ****************************************************************/
/***********************************************************************************************/
int
calc_BfromA()
{
  #ifdef MAGNFIELD

  int ix,iy,iz,iv,ii;
  //A_mu converted to code coordinates in init.c
  
  #pragma omp parallel for private(ix,iy,iz,iv) schedule (static)
  for(ii=0;ii<Nloop_02;ii++) //domain and ghost cells
    {
      ix=loop_02[ii][0];
      iy=loop_02[ii][1];
      iz=loop_02[ii][2]; 

      ldouble B[4];
      ldouble dA[4][4];

      dA[1][1]=dA[2][2]=dA[3][3]=0.;

      //d_1 A_2
      if(NX==1)
	dA[1][2]=0.;
      else
	{
	  if(ix>-NG && ix<NX+NG-1)
	    dA[1][2]=(get_u(p,B2,ix+1,iy,iz)-get_u(p,B2,ix-1,iy,iz))/(get_x(ix+1,0)-get_x(ix-1,0));
	  if(ix==-NG)
	    dA[1][2]=(-3.*get_u(p,B2,ix,iy,iz)+4.*get_u(p,B2,ix+1,iy,iz)-get_u(p,B2,ix,iy,iz+2))/(get_x(ix+2,0)-get_x(ix,0));
	  if(ix==NX+NG-1)
	    dA[1][2]=(3.*get_u(p,B2,ix,iy,iz)-4.*get_u(p,B2,ix-1,iy,iz)+get_u(p,B2,ix-2,iy,iz))/(get_x(ix,0)-get_x(ix-2,0));
	}

      //d_1 A_3
      if(NX==1)
	dA[1][3]=0.;
      else
	{
	  if(ix>-NG && ix<NX+NG-1)
	    dA[1][3]=(get_u(p,B3,ix+1,iy,iz)-get_u(p,B3,ix-1,iy,iz))/(get_x(ix+1,0)-get_x(ix-1,0));
	  if(ix==-NG)
	    dA[1][3]=(-3.*get_u(p,B3,ix,iy,iz)+4.*get_u(p,B3,ix+1,iy,iz)-get_u(p,B3,ix,iy,iz+2))/(get_x(ix+2,0)-get_x(ix,0));
	  if(ix==NX+NG-1)
	    dA[1][3]=(3.*get_u(p,B3,ix,iy,iz)-4.*get_u(p,B3,ix-1,iy,iz)+get_u(p,B3,ix-2,iy,iz))/(get_x(ix,0)-get_x(ix-2,0));
	}

      //d_2 A_3
      if(NY==1)
	dA[2][3]=0.;
      else
	{
	  if(iy>-NG && iy<NY+NG-1)
	    dA[2][3]=(get_u(p,B3,ix,iy+1,iz)-get_u(p,B3,ix,iy-1,iz))/(get_x(iy+1,1)-get_x(iy-1,1));
	  if(iy==-NG)
	    dA[2][3]=(-3.*get_u(p,B3,ix,iy,iz)+4.*get_u(p,B3,ix,iy+1,iz)-get_u(p,B3,ix,iy+2,iz))/(get_x(iy+2,1)-get_x(iy,1));
	  if(iy==NY+NG-1)
	    dA[2][3]=(3.*get_u(p,B3,ix,iy,iz)-4.*get_u(p,B3,ix,iy-1,iz)+get_u(p,B3,ix,iy-2,iz))/(get_x(iy,1)-get_x(iy-2,1));
	}

      //d_2 A_1
      if(NY==1)
	dA[2][1]=0.;
      else
	{
	  if(iy>-NG && iy<NY+NG-1)
	    dA[2][1]=(get_u(p,B1,ix,iy+1,iz)-get_u(p,B1,ix,iy-1,iz))/(get_x(iy+1,1)-get_x(iy-1,1));
	  if(iy==-NG)
	    dA[2][1]=(-3.*get_u(p,B1,ix,iy,iz)+4.*get_u(p,B1,ix,iy+1,iz)-get_u(p,B1,ix,iy+2,iz))/(get_x(iy+2,1)-get_x(iy,1));
	  if(iy==NY+NG-1)
	    dA[2][1]=(3.*get_u(p,B1,ix,iy,iz)-4.*get_u(p,B1,ix,iy-1,iz)+get_u(p,B1,ix,iy-2,iz))/(get_x(iy,1)-get_x(iy-2,1));
	}

      //d_3 A_2
      if(NZ==1)
	dA[3][2]=0.;
      else
	{
	  if(iz>-NG && iz<NZ+NG-1)
	    dA[3][2]=(get_u(p,B2,ix,iy,iz+1)-get_u(p,B2,ix,iy,iz-1))/(get_x(iz+1,2)-get_x(iz-1,2));
	  if(iz==-NG)
	    dA[3][2]=(-3.*get_u(p,B2,ix,iy,iz)+4.*get_u(p,B2,ix,iy,iz+1)-get_u(p,B2,ix,iy,iz+2))/(get_x(iz+2,2)-get_x(iz,2));
	  if(iz==NZ+NG-1)
	    dA[3][2]=(3.*get_u(p,B2,ix,iy,iz)-4.*get_u(p,B2,ix,iy,iz-1)+get_u(p,B2,ix,iy,iz-2))/(get_x(iz,2)-get_x(iz-2,2));
	}

       //d_3 A_1
      if(NZ==1)
	dA[3][1]=0.;
      else
	{
	  if(iz>-NG && iz<NZ+NG-1)
	    dA[3][1]=(get_u(p,B1,ix,iy,iz+1)-get_u(p,B1,ix,iy,iz-1))/(get_x(iz+1,2)-get_x(iz-1,2));
	  if(iz==-NG)
	    dA[3][1]=(-3.*get_u(p,B1,ix,iy,iz)+4.*get_u(p,B1,ix,iy,iz+1)-get_u(p,B1,ix,iy,iz+2))/(get_x(iz+2,2)-get_x(iz,2));
	  if(iz==NZ+NG-1)
	    dA[3][1]=(3.*get_u(p,B1,ix,iy,iz)-4.*get_u(p,B1,ix,iy,iz-1)+get_u(p,B1,ix,iy,iz-2))/(get_x(iz,2)-get_x(iz-2,2));
	}

      //B^i = d_j A_k eps^ijk
      //saving temporarily to pinit

      B[1]=dA[2][3] - dA[3][2];
      B[2]=dA[3][1] - dA[1][3];
      B[3]=dA[1][2] - dA[2][1];

      set_u(pinit,B1,ix,iy,iz,B[1]);
      set_u(pinit,B2,ix,iy,iz,B[2]);
      set_u(pinit,B3,ix,iy,iz,B[3]);
     
    } //cell loop

  //overwriting vector potential with magnetic fields
  #pragma omp parallel for private(ix,iy,iz,iv) schedule (static)
  for(ii=0;ii<Nloop_02;ii++) //domain and ghost cells
    {
      ix=loop_02[ii][0];
      iy=loop_02[ii][1];
      iz=loop_02[ii][2]; 

      struct geometry geom;
      fill_geometry(ix,iy,iz,&geom);
      
      ldouble pp[NV],uu[NV];
      PLOOP(iv)
	pp[iv]=get_u(p,B1,ix,iy,iz);
      pp[B1]=get_u(pinit,B1,ix,iy,iz);
      pp[B2]=get_u(pinit,B2,ix,iy,iz);
      pp[B3]=get_u(pinit,B3,ix,iy,iz);

      p2u(pp,uu,&geom);

      set_u(p,B1,ix,iy,iz,pp[B1]);
      set_u(p,B2,ix,iy,iz,pp[B2]);
      set_u(p,B3,ix,iy,iz,pp[B3]);
      set_u(u,B1,ix,iy,iz,uu[B1]);
      set_u(u,B2,ix,iy,iz,uu[B2]);
      set_u(u,B3,ix,iy,iz,uu[B3]);     
    }
  
#endif //MAGNFIELD

  return 0;
}
