//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geoml;
fill_geometry(global_ix2-1,iy,iz,&geoml);

struct geometry geomr;
fill_geometry(global_ix1,iy,iz,&geomr);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

struct geometry geomBLl;
fill_geometry_arb(global_ix2-1,iy,iz,&geomBLl,KERRCOORDS);


/**********************/

//radius
if(ix>=NX) //total boundary, properties of the galaxy
  {
    ldouble rho,rho0,uint,uintl,uint0,ur,url,rhol;

    //calculating Bondi-related values at the boundary
    ldouble RMAXout=geomBL.xx;
    ldouble mdotscale = rhoGU2CGS(1.)*velGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.);
    ldouble mdotout = MDOT * calc_mdotEdd() / mdotscale;
    ldouble urout = -sqrt(1./2./RMAXout);
    ldouble rhoout = -mdotout / (4.*M_PI *urout* RMAXout * RMAXout);
    ldouble csout = sqrt(1./2./RMAXout);   //cs2 = GM/2R for gamma=5/3
    ldouble uintout = csout * csout * rhoout / GAMMA / GAMMAM1;
    ldouble Eout=PRADGASINIT * GAMMAM1*uintout*3.;

    rho0=rhoout;
    uint0=uintout;
    rho=rho0; 
    uint = uint0; 

    //last but one cell
    url=get_u(p,VX,NX-1,iy,iz);
    rhol=get_u(p,RHO,NX-1,iy,iz);
    uintl=get_u(p,UU,NX-1,iy,iz);
    
    #ifdef FIX_PRESSURE
    //to keep pressure fixed to initial hydro Bondi value
    //and copy temperature, adjust rho    
    ldouble temp=calc_PEQ_Tfromurho(uintl,rhol);
    if(temp<TAMB) temp=TAMB;    
    rho = calc_PEQ_rhofromTu(temp,uint);
    #endif

    #ifdef FIX_TEMPERATURE
    //to keep temperature fixed and copy rho
    rho = rhoout;
    uint=calc_PEQ_ufromTrho(TAMB,rho);
    #endif

    //velocities
    ldouble uconl[4]={0.,url,0.,0.};
    ldouble ucon[4]={0.,0.,0.,0.};
    conv_vels(uconl,uconl,VELPRIM,VEL4,geoml.gg,geoml.GG);

    //by default fixed mdot
    ldouble mdot=rhol*uconl[1]*geoml.gdet;
    ucon[1]=mdot/rho/geom.gdet;
    conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);
    
    pp[0]=rho;
    pp[1]=uint;
    pp[2]=ucon[1];
    pp[3]=ucon[2];
    pp[4]=ucon[3];

    //velocity fixed instead of flat mdot
    #ifdef FIX_VELBONDI
    pp[2]=get_u(pproblem1,VX,ix,iy,iz);
    pp[3]=get_u(pproblem1,VY,ix,iy,iz);
    pp[4]=get_u(pproblem1,VZ,ix,iy,iz);
    #endif
    pp[5]=calc_Sfromu(rho,uint);

    #ifdef RADIATION
    //outflow / no inflow
    
    ldouble El,E,Nf;    
    //last but one cell
    url=get_u(p,FX0,NX-1,iy,iz);
    El=get_u(p,EE0,NX-1,iy,iz);
    #ifdef NCOMPTONIZATION
    Nf=get_u(p,NF0,NX-1,iy,iz);
    #endif
    uconl[0]=uconl[2]=uconl[3]=0.;
    uconl[1]=url;
    conv_vels(uconl,uconl,VELPRIM,VEL4,geoml.gg,geoml.GG);

    //flat luminosity
    ldouble Rijl[4][4];
    calc_Rij(&get_u(p,0,NX-1,iy,iz),&geoml,Rijl); 
    ldouble Rtr = Rijl[1][0]*geoml.gdet/geom.gdet;
    E=Rtr/(4./3.*uconl[0]*uconl[1]+1./3.*geom.gg[0][1]);

    //velocity to BL
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
	E=Eout;
	Nf=pp[NF0]=calc_NFfromE(Eout);
      }

    //velocity back to MYCOORDS
    trans2_coco(geomBL.xxvec,ucon,ucon,BLCOORDS,MYCOORDS);
    conv_vels(ucon,ucon,VEL4,VELPRIMRAD,geom.gg,geom.GG);

    pp[EE0]=E*geomBLl.gdet/geomBL.gdet;
    pp[FX0]=ucon[1];
    pp[FY0]=ucon[2];
    pp[FZ0]=ucon[3]; 
    #ifdef NCOMPTONIZATION
    pp[NF0]=Nf*geomBLl.gdet/geomBL.gdet;
    #endif

#endif	 
    
    p2u(pp,uu,&geom);
 
    //printf("%d > %e %e %e %e\n",ix,url,uconl[1],ucon[1],pp[VX]);getch();
    return 0.;
  }
 else if(ix<global_ix1) //outflow near BH or at inner boundaries
   {
     iix=global_ix1;
     iiy=iy;
     iiz=iz;
   
     //copying primitives with gdet taken into account
     for(iv=0;iv<NV;iv++)
       { 
	 if(1 && iv==VX)
	   {
	     //first cell
	     ldouble urr=get_u(p,VX,global_ix1,iy,iz);
	     ldouble rhor=get_u(p,RHO,global_ix1,iy,iz);
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
  
