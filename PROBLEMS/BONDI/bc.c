//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  

gdet_bc=get_g(g,3,4,ix,iy,iz);  
gdet_src=get_g(g,3,4,iix,iiy,iiz);
ldouble gg[4][5],GG[4][5],ggsrc[4][5],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
pick_T(emuup,ix,iy,iz,eup);
pick_T(emulo,ix,iy,iz,elo);
ldouble xx=get_x(ix,0);

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);


struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);


/**********************/

//radius
//  if(ix>=NX || ix<0) //analytical solution on both sides
if(ix>=NX) //analytical solution at rout only
  {
    ldouble rho,uint,ur,E;

    rho=get_u(pproblem1,RHO,ix,iy,iz);
    uint=get_u(pproblem1,UU,ix,iy,iz);
    ur=get_u(pproblem1,VX,ix,iy,iz);

    //four-vel in BL
    ldouble ucon[4]={0.,ur,0.,0.};
    conv_vels(ucon,ucon,VEL4,VELPRIM,geomBL.gg,geomBL.GG);

    //rad. four-vel in BL
    ldouble urfcon[4]={0.,0.,0.,0.};
    conv_vels(urfcon,urfcon,VEL4,VELPRIM,geomBL.gg,geomBL.GG);

    pp[0]=rho;
    pp[1]=uint;
    pp[2]=ucon[1];
    pp[3]=ucon[2];
    pp[4]=ucon[3];
    pp[5]=calc_Sfromu(rho,uint);
#ifdef RADIATION
    E=get_u(pproblem1,EE0,ix,iy,iz);
    pp[6]=E;
    pp[7]=urfcon[1];
    pp[8]=urfcon[2];
    pp[9]=urfcon[3]; 
#endif	 

    p2u(pp,uu,&geom);

    //transforming primitives from BL to MYCOORDS
    trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,geomBL.xxvec,&geomBL,&geom);

     
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

	 //unchanged primitives
	 pp[iv]=get_u(p,iv,iix,iiy,iiz);

       }

     /*  
     //free-fall acc.to analytical solution from the ix=0 cell
     ldouble Vinfall=sqrtl(2./rbc)*(1.-2./rbc);
     ldouble Dratio=((rsrc*rsrc*sqrtl(2./rsrc*(1.-2./rsrc)))) /
     ((rbc*rbc*sqrtl(2./rbc*(1.-2./rbc))));
     ldouble Eratio=(pow(rsrc*rsrc*sqrt(2./rsrc),GAMMA)*pow(1.-2./rsrc,(GAMMA+1.)/4.)) /
     (pow(rbc*rbc*sqrt(2./rbc),GAMMA)*pow(1.-2./rbc,(GAMMA+1.)/4.));
     ldouble rhosrc=get_u(p,0,iix,iiy,iiz);
     ldouble uintsrc=get_u(p,1,iix,iiy,iiz);
     ldouble Vsrc=get_u(p,2,iix,iiy,iiz);
     ldouble Wsrc=1./sqrtl(1.-Vsrc*Vsrc*ggsrc[1][1]);
     ldouble Dsrc=rhosrc*Wsrc;
     ldouble Esrc=uintsrc*Wsrc;
     ldouble V=Vsrc*(sqrtl(2./rbc)*(1.-2./rbc))/(sqrtl(2./rsrc)*(1.-2./rsrc));
     ldouble W=1./sqrtl(1.-V*V*gg[1][1]);
     ldouble D=Dsrc*Dratio;
     ldouble E=Esrc*Eratio;
     ldouble rho=D/W;
     ldouble uint=E/W;      
     pp[0]=rho; pp[1]=uint; pp[2]=V; pp[3]=pp[4]=0.;
     */
      
     /*
     //analytical solution
     ldouble r=get_x(ix,0);
     ldouble D=PAR_D/(r*r*sqrtl(2./r*(1.-2./r)));
     ldouble E=PAR_E/(pow(r*r*sqrt(2./r),GAMMA)*pow(1.-2./r,(GAMMA+1.)/4.));
     ldouble V=sqrtl(2./r)*(1.-2./r);
     ldouble W=1./sqrtl(1.-V*V*gg[1][1]);
     ldouble rho=D/W;
     ldouble uint=E/W;
     pp[0]=rho; pp[1]=uint; pp[2]=-V; pp[3]=pp[4]=0.;
     */
     
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
  
