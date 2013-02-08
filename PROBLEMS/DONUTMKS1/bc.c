//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  
ldouble xxvec[4],xxvecBL[4],xx,yy,zz;

get_xx(ix,iy,iz,xxvec);
coco_N(xxvec,xxvecBL,MYCOORDS,BLCOORDS);
xx=xxvecBL[1];
yy=xxvecBL[2];
zz=xxvecBL[3];
gdet_bc=get_g(g,3,4,ix,iy,iz);  
//gdet_src=get_g(g,3,4,iix,iiy,iiz);
ldouble gg[4][5],GG[4][5],ggsrc[4][5],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
pick_T(emuup,ix,iy,iz,eup);
pick_T(emulo,ix,iy,iz,elo);
//working in BL
ldouble ggBL[4][5],GGBL[4][5];
calc_g_arb(xxvecBL,ggBL,KERRCOORDS);
calc_G_arb(xxvecBL,GGBL,KERRCOORDS);
ldouble eupBL[4][4],eloBL[4][4];
ldouble tupBL[4][4],tloBL[4][4];
calc_tetrades(ggBL,tupBL,tloBL,KERRCOORDS);
calc_ZAMOes(ggBL,eupBL,eloBL,KERRCOORDS);
/**********************/


//radius
if(ix>=NX) //analytical solution at rout only
  {
    ldouble podpierd=-(GGBL[0][0]-2.*ELL*GGBL[0][3]+ELL*ELL*GGBL[3][3]);
    ldouble ut=-1./sqrt(podpierd);
    ut/=UTPOT;
    ldouble uint,Vphi,rho,Vr;
    ldouble D,E,W,eps,uT,uphi,uPhi;
    if(ut<-1 || podpierd<0.|| NODONUT)
      {
	//ambient
	set_hdatmosphere(pp,xxvec,gg,GG,0);

	
	//BL free-fall velocity
	ldouble ucon[4];
	ldouble r=xx;
	ucon[0]=0.;
	ucon[1]=-sqrtl(2./r)*(1.-2./r);
	ucon[2]=ucon[3]=0.;
	conv_vels(ucon,ucon,VEL3,VEL4,ggBL,GGBL);
	trans2_coco(xxvecBL,ucon,ucon,BLCOORDS,MYCOORDS);
	conv_vels(ucon,ucon,VEL4,VELPRIM,gg,GG);
	pp[2]=ucon[1];
	pp[3]=ucon[2];
	pp[4]=ucon[3];
	
      }
    else
      {
	ldouble h=-1./ut;
	ldouble eps=(h-1.)/GAMMA;
	rho=powl(eps*(GAMMA-1.)/KKK,1./(GAMMA-1.));
	uint=rho*eps;
	uphi=-ELL*ut;
	uT=GGBL[0][0]*ut+GGBL[0][3]*uphi;
	uPhi=GGBL[3][3]*uphi+GGBL[0][3]*ut;
	Vphi=uPhi/uT;
	Vr=0.;

	//4-velocity in BL
	ldouble ucon[4]={0.,-Vr,0.,Vphi};
	conv_vels(ucon,ucon,VEL3,VEL4,ggBL,GGBL);
	trans2_coco(xxvecBL,ucon,ucon,BLCOORDS,MYCOORDS);
	conv_vels(ucon,ucon,VEL4,VELPRIM,gg,GG);
   
	pp[2]=ucon[1]; 
	pp[3]=ucon[2];
	pp[4]=ucon[3];

	//density etc.
	pp[0]=rho; pp[1]=uint; 
      }     

   
    pp[5]=calc_Sfromu(pp[0],pp[1]);
    
    //testing if primitives make sense
    check_floors_hd(pp,VELPRIM,gg,GG);
    //end of floor section

    p2u(pp,uu,gg,GG);

    return 0.;
  }
 else if(ix<0) //outflow near BH
   {
     iix=0;
     iiy=iy;
     iiz=iz;

     ldouble r=xx; //BL radius
     ldouble xx0[4];
     get_xx(iix,iiy,iiz,xx0);
     coco_N(xx0,xx0,MYCOORDS,BLCOORDS);
     ldouble r0=xx0[1];
     
     pp[0]=get_u(p,0,iix,iiy,iiz)*pow(r/r0,-1.5);
     pp[1]=get_u(p,1,iix,iiy,iiz)*pow(r/r0,-2.5);

     //copying MLCOORDS velocities
     for(iv=2;iv<NV;iv++)
       { 
	 //unchanged primitives
	 pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }

     if(pp[2]>0.) pp[2]=0.;

     //testing if interpolated primitives make sense
     check_floors_hd(pp,VELPRIM,gg,GG);
     //end of floor section

     p2u(pp,uu,gg,GG);
     return 0;
   }

//reflections in theta 
if(iy<0.) //spin axis
  {      
    iiy=-iy-1;
    iiz=iz;
    iix=ix;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  
    for(iv=0;iv<NV;iv++)
      {
	if(iv==3)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,gg,GG);
    //end of floor section

    p2u(pp,uu,gg,GG);
    return 0;
  }
if(iy>=NY) //equatorial plane
  {
    iiy=NY-(iy-NY)-1;
    iiz=iz;
    iix=ix;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  
  	  
    for(iv=0;iv<NV;iv++)
      {
	if(iv==3)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }
    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,gg,GG);
    //end of floor section

    p2u(pp,uu,gg,GG); 
    return 0; 
  }
   
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

//testing if interpolated primitives make sense
check_floors_hd(pp,VELPRIM,gg,GG);
//end of floor section
 
return 0;

