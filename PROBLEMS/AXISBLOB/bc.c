/**********************/
//geometries
int iix,iiy,iiz,iv;  	  
ldouble xxvec[4], xxcart[4];

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomcart;
fill_geometry_arb(ix,iy,iz,&geomcart,MINKCOORDS);

ldouble r=geom.xx;
ldouble th=geom.yy;

//radius
if(BCtype==XBCHI)
  {
    iix=NX-1;
    iiy=iy;
    iiz=iz;

    //copying everything
    for(iv=0;iv<=NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

    if(pp[VX]<0.) pp[VX]=0.;

    p2u(pp,uu,&geom);
    return 0;  
  }
 else if(BCtype==XBCLO) //outflow near BH
   {
     iix=0;
     iiy=iy;
     iiz=iz;

     
    //copying everything
    for(iv=0;iv<=NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

    if(pp[VX]>0.) pp[VX]=0.;

    p2u(pp,uu,&geom);
    return 0;
   }

//reflections in theta 
if(BCtype==YBCLO) //the axis 
  {      
    //simple reflection
    /*
    iiy=-iy-1;
    iiz=iz;
    iix=ix;

    for(iv=0;iv<NV;iv++)
      {
	if(iv==VY || iv==B2 || iv==VZ || iv==B3)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }
    */

    //from the other side - only for openMP! for MPI use CORRECTAXIS!    
    iiy=-iy-1;
    iiz=(iz+NZ/2)%NZ;
    iix=ix;

    for(iv=0;iv<NV;iv++)
      {
	if(iv==VY || iv==B2 || iv==VZ || iv==B3)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }
    
    
    p2u(pp,uu,&geom);
    return 0;
  }

if(BCtype==YBCHI) //~equatorial plane
  {
    iiy=NY-1;
    iiz=iz;
    iix=ix;
  
    //copying everything
    for(iv=0;iv<=NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

    if(pp[VY]<0.) pp[VY]=0.;

    p2u(pp,uu,&geom);
    return 0;  
  }
   
//periodic in phi:
iiz=iz;
iiy=iy;
iix=ix;
if(BCtype==ZBCLO) iiz=iz+NZ;
if(BCtype==ZBCHI) iiz=iz-NZ;

for(iv=0;iv<NV;iv++)
  {
    uu[iv]=get_u(u,iv,iix,iiy,iiz);
    pp[iv]=get_u(p,iv,iix,iiy,iiz);      
  }

//testing if interpolated primitives make sense
//check_floors_hd(pp,VELPRIM,gg,GG);
//end of floor section
 
return 0;

