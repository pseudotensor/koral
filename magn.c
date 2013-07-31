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
int flx_x(int i)
{ return i;}
int flx_y(int i)
{ if(NY==1) return 0;
  if(i<0) return 0;
  if(i>=NY) return NY-1; }
int flx_z(int i)
{ if(NZ==1) return 0;
  if(i<0) return 0;
  if(i>=NZ) return NZ-1; }

int fly_x(int i)
{ if(NX==1) return 0;
  if(i<0) return 0;
  if(i>=NX) return NX-1; }
int fly_y(int i)
{ return i;}
int fly_z(int i)
{ if(NZ==1) return 0;
  if(i<0) return 0;
  if(i>=NZ) return NZ-1; }

int flz_x(int i)
{ if(NX==1) return 0;
  if(i<0) return 0;
  if(i>=NX) return NX-1; }
int flz_y(int i)
{ if(NY==1) return 0;
  if(i<0) return 0;
  if(i>=NY) return NY-1; }
int flz_z(int i)
{ return i;}

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
      set_emf(1,ix,iy,iz,
	      coefemf[1] * (
                            #if(NY>1)
			    + get_ub(flby,fly_x(ix),fly_y(iy),fly_z(iz),B3,1) 
			    + get_ub(flby,fly_x(ix),fly_y(iy),fly_z(iz-1),B3,1)
                            #endif
                            #if(NZ>1)
			    - get_ub(flbz,flz_x(ix),flz_y(iy),flz_z(iz),B2,2) 
			    - get_ub(flbz,flz_x(ix),flz_y(iy-1),flz_z(iz),B2,2)
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
			    + get_ub(flbz,flz_x(ix),flz_y(iy),flz_z(iz),B1,2) 
			    + get_ub(flbz,flz_x(ix-1),flz_y(iy),flz_z(iz),B1,2)
                            #endif
                            #if(NX>1)
			    - get_ub(flbx,flx_x(ix),flx_y(iy),flx_z(iz),B3,0) 
			    - get_ub(flbx,flx_x(ix),flx_y(iy),flx_z(iz-1),B3,0)
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
			    + get_ub(flbx,flx_x(ix),flx_y(iy),flx_z(iz),B2,0) 
			    + get_ub(flbx,flx_x(ix),flx_y(iy-1),flx_z(iz),B2,0)
                            #endif
                            #if(NY>1)
			    - get_ub(flby,fly_x(ix),fly_y(iy),fly_z(iz),B1,1) 
			    - get_ub(flby,fly_x(ix-1),fly_y(iy),fly_z(iz),B1,1)
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
      //TODO: check if this is the right loop

      ix=loop_4[ii][0];
      iy=loop_4[ii][1];
      iz=loop_4[ii][2]; 

      /////////////////////////////////////
      // F1
      ////////////////////////////////////
#if(NX>1)
      if(iy<NY && iz<NZ) //no need to fill x-face fluxes for iy=NY
	{
	  set_ubx(flbx,ix,iy,iz,B1,0.);
	  set_ubx(flbx,ix,iy,iz,B2,0.5 * (get_emf(3,ix,iy,iz) + get_emf(3,ix,iy+1,iz)));
	  set_ubx(flbx,ix,iy,iz,B3,-0.5 * (get_emf(2,ix,iy,iz) + get_emf(2,ix,iy,iz+1)));
	}
#endif

      /////////////////////////////////////
      // F2
      ////////////////////////////////////
#if(NY>1)
      if(ix<NX && iz<NZ)	
	{
	  set_uby(flby,ix,iy,iz,B1,-0.5 * (get_emf(3,ix,iy,iz) + get_emf(3,ix+1,iy,iz)));
	  set_uby(flby,ix,iy,iz,B2,0.);
	  set_uby(flby,ix,iy,iz,B3,0.5 * (get_emf(1,ix,iy,iz) + get_emf(1,ix,iy,iz+1)));
	}
#endif

      /////////////////////////////////////
      // F3
      ////////////////////////////////////
#if(NZ>1)
      if(ix<NX && iy<NY)	
	{
	  set_ubz(flbz,ix,iy,iz,B1,0.5 * (get_emf(2,ix,iy,iz) + get_emf(2,ix+1,iy,iz)));
	  set_ubz(flbz,ix,iy,iz,B2,-0.5 * (get_emf(1,ix,iy,iz) + get_emf(1,ix,iy+1,iz)));
	  set_ubz(flbz,ix,iy,iz,B3,0.);
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
