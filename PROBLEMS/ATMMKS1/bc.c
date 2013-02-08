//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

int iix,iiy,iiz,iv;  

ldouble xx,yy,zz,xxvecBL[4],xxvec[4];
get_xx(ix,iy,iz,xxvec);
coco_N(xxvec,xxvecBL,MYCOORDS,BLCOORDS);
xx=xxvecBL[1];
yy=xxvecBL[2];
zz=xxvecBL[3];

ldouble gg[4][5],eup[4][4],elo[4][4],GG[4][5];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);

//KERR metric
ldouble ggBL[4][5],GGBL[4][5];
calc_g_arb(xxvecBL,ggBL,KERRCOORDS);
calc_G_arb(xxvecBL,GGBL,KERRCOORDS);


/**********************/

//radius
if(ix>=NX) 
  {
    //ambient
    set_hdatmosphere(pp,xxvec,gg,GG,0);

    //BL free-fall velocity
    ldouble ucon[4];
    ldouble r=xx;
    ucon[1]=-sqrtl(2./r)*(1.-2./r);
    ucon[2]=ucon[3]=0.;
    conv_vels(ucon,ucon,VEL3,VEL4,ggBL,GGBL);
    trans2_coco(xxvecBL,ucon,ucon,BLCOORDS,MYCOORDS);
    conv_vels(ucon,ucon,VEL4,VELPRIM,gg,GG);
    pp[2]=ucon[1];
    pp[3]=ucon[2];
    pp[4]=ucon[3];

    //calculating entropy
    pp[5]=calc_Sfromu(pp[0],pp[1]);
	      
    //testing if interpolated primitives make sense
    check_floors_hd(pp,VELPRIM,gg,GG);
    //end of floor section

    p2u(pp,uu,gg,GG);

    return 0.;
  }
 else if(ix<0) //outflow at the inner edge
   {      
     iix=0;
     iiy=iy;
     iiz=iz;

     ldouble r=xxvecBL[1];
     ldouble r0,xx0[4];
     get_xx(iix,iiy,iiz,xx0);
     coco_N(xx0,xx0,MYCOORDS,BLCOORDS);
     r0=xx0[1];

     
     pp[0]=get_u(p,0,iix,iiy,iiz)*pow(r/r0,-1.5);
     pp[1]=get_u(p,1,iix,iiy,iiz)*pow(r/r0,-2.5);

     //copying MYCOORDS velocities
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

   
//symmetric in phi:
while(iz<0)    iz+=NZ;
while(iz>=NZ)    iz-=NZ; 
//symmetric in theta:
while(iy<0)    iy+=NY;
while(iy>=NY)    iy-=NY; 

for(iv=0;iv<NV;iv++)
  {
    uu[iv]=get_u(u,iv,ix,iy,iz);
    pp[iv]=get_u(p,iv,ix,iy,iz);      
  }
  
return 0;
  
