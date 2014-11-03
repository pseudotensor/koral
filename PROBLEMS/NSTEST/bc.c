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

    struct geometry geomBL0,geom0;
    fill_geometry_arb(0,iy,iz,&geomBL0,BLCOORDS);
    fill_geometry(0,iy,iz,&geom0);
    ldouble r0=geomBL0.xx;
    ldouble scfac = 1.;
    #ifdef MONOPOLE
    scfac=r0*r0/r/r;
    #endif
    #ifdef DIPOLE
    scfac=r0*r0*r0/r/r/r;
    #endif

     if(ifinit) //magnetic field fixed
       {
	 //then use initial dipolar magnetic field from the first domain cell scaled down in both
	 /*
	 pp[B1]=get_u(pproblem1,B1,0,iy,iz)*scfac;
	 pp[B2]=get_u(pproblem1,B2,0,iy,iz)*scfac;
	 pp[B3]=get_u(pproblem1,B3,0,iy,iz)*scfac;
	 */
	 pp[B1]=get_u(p,B1,0,iy,iz)*sqrt(geom0.gg[1][1])*scfac/sqrt(geom.gg[1][1]);
	 pp[B2]=get_u(p,B2,0,iy,iz)*sqrt(geom0.gg[2][2])*scfac/sqrt(geom.gg[2][2]);
	 pp[B3]=get_u(p,B3,0,iy,iz)*sqrt(geom0.gg[3][3])*scfac/sqrt(geom.gg[3][3]);
       }

     pp[B2]=get_u(p,B2,0,iy,iz)*sqrt(geom0.gg[2][2])*scfac/sqrt(geom.gg[2][2]);
     pp[B3]=get_u(p,B3,0,iy,iz)*sqrt(geom0.gg[3][3])*scfac/sqrt(geom.gg[3][3]);

    //density fixed, temparture too
     ldouble rfac=pow(r/RMIN,-6.);
     rfac=1.;

   

    //velocities zero - not rotating, settling down
     pp[VX]=pp[VY]=pp[VZ]=0.;

     //pure rotation
     //pp[VZ]=OMEGA;
     
    
     ldouble Bcon0[4]={0.,get_u(p,B1,0,iy,iz),get_u(p,B2,0,iy,iz),get_u(p,B3,0,iy,iz)},Bcov0[4];
     indices_21(Bcon0,Bcov0,geom0.gg);
     ldouble vcon0[4]={0.,get_u(p,VX,0,iy,iz),get_u(p,VY,0,iy,iz),get_u(p,VZ,0,iy,iz)};
     ldouble vpar0 = dot(vcon0,Bcov0) / sqrt(dot(Bcon0,Bcov0)) * my_sign(Bcon0[1]);

     ldouble vpar;
     if(vpar0>0.) vpar=0.5;
     else vpar=vpar0;

     ldouble Bcon[4]={0.,pp[B1],pp[B2],pp[B3]},Bcov[4];
     indices_21(Bcon,Bcov,geom.gg);
    
     pp[VX]=pp[B1]*vpar/sqrt(dot(Bcon,Bcov))* my_sign(Bcon0[1]);;
     pp[VY]=pp[B2]*vpar/sqrt(dot(Bcon,Bcov))* my_sign(Bcon0[1]);;
     pp[VZ]=OMEGA+pp[B3]*vpar/sqrt(dot(Bcon,Bcov))* my_sign(Bcon0[1]);;
     
     if(vpar0>0.)
       {
	 pp[RHO]=RHO_AMB*rfac;
	 pp[UU]=U_AMB*rfac;
       }
     else
       {
	 pp[RHO]=get_u(p,RHO,0,iy,iz);
	 pp[UU]=get_u(p,UU,0,iy,iz);
       }


     pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU]);


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
