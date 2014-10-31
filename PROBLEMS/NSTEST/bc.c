/**********************/
//geometries
int iix,iiy,iiz,iv;  	  
ldouble xxvec[4], xxcart[4];

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

ldouble r=geomBL.xx;
ldouble th=geomBL.yy;

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
 else if(BCtype==XBCLO) //NS surface
   {
     iix=0;
     iiy=iy;
     iiz=iz;
     
     
    //copying everything
    for(iv=0;iv<=NV;iv++)
      {
	pp[iv]=get_u(p,iv,ix,iy,iz);
      }

    struct geometry geomBL0;
    fill_geometry_arb(0,iy,iz,&geomBL0,BLCOORDS);
    ldouble r0=geomBL0.xx;
    ldouble scfac = r0*r0*r0/(r*r*r);

     if(ifinit) //magnetic field fixed
       {
	 //then use initial dipolar magnetic field from the first domain cell scaled down in both
	 /*
	 pp[B1]=get_u(pproblem1,B1,0,iy,iz)*scfac;
	 pp[B2]=get_u(pproblem1,B2,0,iy,iz)*scfac;
	 pp[B3]=get_u(pproblem1,B3,0,iy,iz)*scfac;
	 */
	 pp[B1]=get_u(p,B1,0,iy,iz)*scfac;
	 pp[B2]=get_u(p,B2,0,iy,iz)*scfac;
	 pp[B3]=get_u(p,B3,0,iy,iz)*scfac;
       }

     //pp[B2]=get_u(p,B2,0,iy,iz)*scfac;
     //pp[B3]=get_u(p,B3,0,iy,iz)*scfac;

    //density fixed, temparture too
     ldouble rfac=pow(r/RMIN,-6.);

     pp[RHO]=RHO_AMB*rfac;
     pp[UU]=U_AMB;//*rfac;
    pp[5]=calc_Sfromu(pp[0],pp[1]);

    //velocities zero - not rotating, settling down
    pp[VX]=pp[VY]=pp[VZ]=0.;
    



    p2u(pp,uu,&geom);
    return 0;
   }

//reflections in theta 
if(BCtype==YBCLO) //the axis 
  {      
    //simple reflection
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
    
    p2u(pp,uu,&geom);
    return 0;
  }

if(BCtype==YBCHI) //equatorial plane
  {
    //simple reflection
    iiy=NY-(iy-NY)-1;
    iiz=iz;
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

 
return 0;

