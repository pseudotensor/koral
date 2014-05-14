//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  


struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geoml;
fill_geometry(NX-1,iy,iz,&geoml);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

struct geometry geomBLl;
fill_geometry_arb(NX-1,iy,iz,&geomBLl,KERRCOORDS);


/**********************/

//radius
if(ix>=NX) //analytical solution at rout only
  {
    ldouble rho,uint,ur,url,rhol;

    rho=get_u(pproblem1,RHO,ix,iy,iz);
    uint=get_u(pproblem1,UU,ix,iy,iz);

    //last but one cell
    url=get_u(p,VX,NX-1,iy,iz);
    rhol=get_u(p,RHO,NX-1,iy,iz);
    ldouble uconl[4]={0.,url,0.,0.};
    conv_vels(uconl,uconl,VELPRIM,VEL4,geoml.gg,geoml.GG);
    trans2_coco(geoml.xxvec,uconl,uconl,MYCOORDS,BLCOORDS);
    ldouble mdot=rhol*uconl[1]*geomBLl.xx*geomBLl.xx;

    //BL velocity in the ghost cell
    ldouble ucon[4]={0.,0.,0.,0.};
    ucon[1]=mdot/rho/geomBL.xx/geomBL.xx;
    conv_vels(ucon,ucon,VEL4,VELPRIM,geomBL.gg,geomBL.GG);

    
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
	E=El;//ERADRES*get_u(p,UU,NX-1,iy,iz);
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
#endif	 
    
    //transforming primitives from BL to MYCOORDS
    trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,geomBL.xxvec,&geomBL,&geom);
    p2u(pp,uu,&geom);
 
    //printf("%d > %e %e %e %e\n",ix,url,uconl[1],ucon[1],pp[VX]);getch();
    return 0.;
  }
 else if(ix<0) //outflow near BH
   {
     ldouble gdet=get_g(g,3,4,ix,iy,iz);  
     ldouble gg[4][5],eup[4][4],elo[4][4],ggsrc[4][5];

     iix=0;
     iiy=iy;
     iiz=iz;
     pick_g(ix,iy,iz,gg);
     pick_T(emuup,ix,iy,iz,eup);
     pick_T(emulo,ix,iy,iz,elo);

     pick_g(iix,iiy,iiz,ggsrc);
     gdet_src=get_g(g,3,4,iix,iiy,iiz);  
     gdet_bc=get_g(g,3,4,ix,iy,iz);        
     ldouble rsrc=get_x(iix,0);
     ldouble rbc=get_x(ix,0);

     ldouble Fx,Fy,Fz,rho,rho0,Tgas0,E,uint,ur,Tgas,Trad,r,prad,pgas,ut,vx;

     //copying primitives with gdet taken into account
     for(iv=0;iv<NV;iv++)
       { 
	 /*
	 if(iv==2 || iv==7)
	   pp[iv]=get_u(p,iv,iix,iiy,iiz)*(1.-(rsrc-rbc)/(.5*(rsrc+rbc)));
	 else if(iv==3 || iv==4 || iv==8 || iv==9)
	   pp[iv]=get_u(p,iv,iix,iiy,iiz)*(1.+(rsrc-rbc)/(.5*(rsrc+rbc)));
	 else 
	   pp[iv]=get_u(p,iv,iix,iiy,iiz);//gdet_src/gdet_bc;

	 //following ~r**-1.5 scaling
	 if(iv==0)
	   pp[iv]=get_u(p,iv,iix,iiy,iiz)*pow(rsrc/rbc,1.5);
	 if(iv==1)
	   pp[iv]=get_u(p,iv,iix,iiy,iiz)*pow(rsrc/rbc,1.5*GAMMA);
	 if(iv==2)
	   pp[iv]=vx;
	 */
	
	 pp[iv]=get_u(p,iv,iix,iiy,iiz);

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
  
