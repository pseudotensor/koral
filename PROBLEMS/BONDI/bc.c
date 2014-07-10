//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  


struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geoml;
fill_geometry(NX-1,iy,iz,&geoml);

struct geometry geomr;
fill_geometry(0,iy,iz,&geomr);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

struct geometry geomBLl;
fill_geometry_arb(NX-1,iy,iz,&geomBLl,KERRCOORDS);


/**********************/

//radius
if(ix>=NX) //analytical solution at rout only
  {
    ldouble rho,rho0,uint,uint0,ur,url,rhol;

    rho0=get_u(pproblem1,RHO,ix,iy,iz);
    uint0=get_u(pproblem1,UU,ix,iy,iz);

    //last but one cell
    url=get_u(p,VX,NX-1,iy,iz);
    rhol=get_u(p,RHO,NX-1,iy,iz);

    rho = rhol;

    uint = uint0; //to keep pressure fixed
    if(calc_PEQ_Tfromurho(uint,rho) < TAMB) //too cold
      uint = calc_PEQ_ufromTrho(TAMB,rho);

    ldouble uconl[4]={0.,url,0.,0.};
    conv_vels(uconl,uconl,VELPRIM,VEL4,geoml.gg,geoml.GG);

    //velocity in the ghost cell to keep mdot const
    ldouble mdot=rhol*uconl[1]*geoml.gdet;
    ldouble ucon[4]={0.,0.,0.,0.};
    ucon[1]=mdot/rho/geom.gdet;
    conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);
    
    pp[0]=rho;
    pp[1]=uint;
    pp[2]=ucon[1];
    pp[3]=ucon[2];
    pp[4]=ucon[3];
    pp[5]=calc_Sfromu(rho,uint);

#ifdef RADIATION
    //outflow / no inflow
    ldouble El,E;
    
    //last but one cell
    url=get_u(p,FX0,NX-1,iy,iz);
    El=get_u(p,EE0,NX-1,iy,iz);
    uconl[0]=uconl[2]=uconl[3]=0.;
    uconl[1]=url;
    conv_vels(uconl,uconl,VELPRIM,VEL4,geoml.gg,geoml.GG);
    trans2_coco(geoml.xxvec,uconl,uconl,MYCOORDS,BLCOORDS);

    //ghost cell
    ucon[0]=ucon[2]=ucon[3]=0.;
    if(uconl[1]>0.)
      {
	E=El;
	ucon[1]=uconl[1];
      }
    else
      {
	ucon[1]=0.;
	E=El;
      }

    conv_vels(ucon,ucon,VEL4,VELPRIM,geomBL.gg,geomBL.GG);
    pp[6]=E;
    pp[7]=ucon[1];
    pp[8]=ucon[2];
    pp[9]=ucon[3]; 

#ifdef LIKEINFRAGILE
    pp[6]=get_u(pproblem1,EE0,ix,iy,iz);
    pp[7]=get_u(pproblem1,FX0,ix,iy,iz);
    pp[8]=get_u(pproblem1,FY0,ix,iy,iz);
    pp[9]=get_u(pproblem1,FZ0,ix,iy,iz);
#endif

    //transforming rad primitives from BL to MYCOORDS
    trans_prad_coco(pp, pp, KERRCOORDS, MYCOORDS,geomBL.xxvec,&geomBL,&geom);
    
#endif	 
    
    p2u(pp,uu,&geom);
 
    //printf("%d > %e %e %e %e\n",ix,url,uconl[1],ucon[1],pp[VX]);getch();
    return 0.;
  }
 else if(ix<0) //outflow near BH
   {
     iix=0;
     iiy=iy;
     iiz=iz;
   
     //copying primitives with gdet taken into account
     for(iv=0;iv<NV;iv++)
       { 
	 if(iv==VX)
	   {
	     //first cell
	     ldouble urr=get_u(p,VX,0,iy,iz);
	     ldouble rhor=get_u(p,RHO,0,iy,iz);
	     ldouble rho=get_u(p,RHO,ix,iy,iz);
	     ldouble uconr[4]={0.,urr,0.,0.};
	     conv_vels(uconr,uconr,VELPRIM,VEL4,geomr.gg,geomr.GG);
	     ldouble mdot=rhor*uconr[1]*geomr.gdet;
	     uconr[1]=mdot/rho/geom.gdet;
	     conv_vels(uconr,uconr,VEL4,VELPRIM,geom.gg,geom.GG);
	     pp[VX]=uconr[1];
	   }
	 else
	   {
	
	     pp[iv]=get_u(p,iv,iix,iiy,iiz);
	   }

       }

     
     p2u(pp,uu,&geom);
     return 0;
   }

iix=ix;
iiz=iz;
iiy=iy;

//periodic
while(iiz<0)    iiz+=NZ;
while(iiz>=NZ)    iiz-=NZ; 
//periodic
while(iiy<0)    iiy+=NY;
while(iiy>=NY)    iiy-=NY; 


for(iv=0;iv<NV;iv++)
  {
    uu[iv]=get_u(u,iv,iix,iiy,iiz);
    pp[iv]=get_u(p,iv,iix,iiy,iiz);      
  }
  
return 0;
  
