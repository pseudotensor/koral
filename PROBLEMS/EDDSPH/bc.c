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
if(ix>=NX) //analytical solution within the torus and atmosphere outside
  {
    //hydro atmosphere
    set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,2);
    //rad atmosphere
#ifdef RADIATION
    set_radatmosphere(pp,geom.xxvec,geom.gg,geom.GG,1);
#endif

    //outflow of radiation
    for(iv=NVHD;iv<NV;iv++)
       {
	 pp[iv]=get_u(p,iv,NX-1,iy,iz);
       }

    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,&geom);
    //end of floor section

    p2u(pp,uu,&geom);
 
    return 0;
  }
 
/***********************************************/
/***********************************************/
if(ix<0) //outflow near 
   {
     for(iv=0;iv<NVHD;iv++)
       {
	 pp[iv]=get_u(p,iv,0,iy,iz);
       }

     pp[6]=10000.*ERADATMMIN;
     pp[7]=pp[6]/2.;
     pp[8]=pp[9]=0.;
     
     prad_on2lab(pp,pp,&geom);
 
     //testing if interpolated primitives make sense
     check_floors_hd(pp,VELPRIM,&geom);
     //end of floor section

     p2u(pp,uu,&geom);
     return 0;
   }

/***********************************************/
/***********************************************/
//periodic in phi:

if(iy<0 || iy>NY-1)
  {
    iiz=iz;
    iiy=iy;
    iix=ix;
    if(iy<0) iiy=iy+NY;
    if(iy>NY-1) iiy=iy-NY;

    for(iv=0;iv<NV;iv++)
      {
	uu[iv]=get_u(u,iv,iix,iiy,iiz);
	pp[iv]=get_u(p,iv,iix,iiy,iiz);      
      }
  }

/***********************************************/
/***********************************************/
//periodic in phi:
if(iz<0 || iz>NZ-1)
  {
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
  }

//and that is all
 
return 0;

