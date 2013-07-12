//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

int iix,iiy,iiz,iv;

/***********************************************/
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
/***********************************************/

/***********************************************/
/***********************************************/
if(ix>=NX) //outflow
  {
    //outflow BC
    for(iv=0;iv<NV;iv++)
       {
	 pp[iv]=get_u(p,iv,NX-1,iy,iz);
       }

    //zero radial velocity
    if(pp[VX]<0.) 
      {
	pp[VX]=0.;		
      }

    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,&geom);
    //end of floor section
     
    p2u(pp,uu,&geom);
 
    return 0;
  }
 
if(ix<0) //bulk inflow
  {
    //hydro atmosphere
    pp[RHO]=RHOZERO;
    pp[UU]=calc_PEQ_ufromTrho(TEMPZERO,RHOZERO);
    pp[VX]=VELBULK;
    pp[VY]=0.;
    pp[VZ]=0.;

    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,&geom);
    //end of floor section
     
    p2u(pp,uu,&geom);
 
    return 0;
  }
 
if(iy<0) //beam
  {
    if(geom.xx<BEAMX1 || geom.xx>BEAMX2)
      {
	//hydro atmosphere
	pp[RHO]=RHOZERO;
	pp[UU]=calc_PEQ_ufromTrho(TEMPZERO,RHOZERO);
	pp[VX]=VELBULK;
	pp[VY]=0.;
	pp[VZ]=0.;
      }
    else
      {
	//hydro atmosphere
	pp[RHO]=RHOBEAM;
	pp[UU]=calc_PEQ_ufromTrho(TEMPBEAM,RHOBEAM);
	pp[VX]=0;
	pp[VY]=VELBEAM;
	pp[VZ]=0.;
      }

    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,&geom);
    //end of floor section
     
    p2u(pp,uu,&geom);
 
    return 0;
  }
 
if(iy>=NY) //bulk motion
  {
    //hydro atmosphere
    pp[RHO]=RHOZERO;
    pp[UU]=calc_PEQ_ufromTrho(TEMPZERO,RHOZERO);
    pp[VX]=VELBULK;
    pp[VY]=0.;
    pp[VZ]=0.;

    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,&geom);
    //end of floor section
     
    p2u(pp,uu,&geom);
 
    return 0;
  }

/***********************************************/
/***********************************************/
//periodic in z:
iiz=iz;
iiy=iy;
iix=ix;
if(iz<0) iiz=iz+NZ;
if(iz>NZ-1) iiz=iz-NZ;

for(iv=0;iv<NV;iv++)
  {
    uu[iv]=get_u(u,iv,iix,iiy,iiz);
    pp[iv]=get_u(p,iv,iix,iiy,iiz);      
  }

//and that is all
 
return 0;

