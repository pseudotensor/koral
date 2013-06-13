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
    set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);

    //rad atmosphere
#ifdef RADIATION
    set_radatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
#endif
    p2u(pp,uu,&geom);
    return 0;
  }
 
/***********************************************/
/***********************************************/
if(ix<0) //outflow near BH
   {
     //extrapolating along MYCOORDS
     //gc radialcoordinate
     ldouble r=get_x(ix,0);
     //iix=0
     ldouble r0=get_x(0,0);
     ldouble r1=get_x(1,0);     

     for(iv=0;iv<NV;iv++)
       {
	 //linear extrapolation
	 pp[iv]=get_u(p,iv,0,iy,iz)+(get_u(p,iv,1,iy,iz)-get_u(p,iv,0,iy,iz))*(r-r0)/(r1-r0);
	 //constant
	 //pp[iv]=get_u(p,iv,0,iiy,iiz);
       }
 
      //testing if interpolated primitives make sense
      check_floors_hd(pp,VELPRIM,&geom);
      //end of floor section

      p2u(pp,uu,&geom);
      return 0;
   }

/***********************************************/
/***********************************************/
//reflections/outflow in theta 
if(iy<0.) //spin axis 
  {      
    iiy=-iy-1;
    iiz=iz;
    iix=ix;
    for(iv=0;iv<NV;iv++)
      {
	//theta vector components
	if(iv==3 || iv==8)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }
    
    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,&geom);
    //end of floor section

    p2u(pp,uu,&geom);
    return 0;
  }

/***********************************************/
/***********************************************/
if(iy>=NY) //equatorial plane
  {
    iiy=NY-(iy-NY)-1;
    iiz=iz;
    iix=ix;
  	  
    for(iv=0;iv<NV;iv++)
      {
	if(iv==3 || iv==8)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,&geom);
    //end of floor section

    p2u(pp,uu,&geom); 
    return 0; 
  }
   
/***********************************************/
/***********************************************/
//periodic in phi:
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

