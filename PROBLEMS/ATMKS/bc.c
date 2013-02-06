//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

int iix,iiy,iiz,iv;  	  
ldouble gg[4][5],ggsrc[4][5],eup[4][4],elo[4][4],GG[4][5];
ldouble xx,yy,zz,xxvec[4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);

get_xx(ix,iy,iz,xxvec);
//coco_N(xxvec,xxvec,KSCOORDS,BLCOORDS);
xx=xxvec[1];
yy=xxvec[2];
zz=xxvec[3];


/**********************/

//radius
if(ix>=NX) 
  {
    //ambient
    set_hdatmosphere(pp,xxvec,gg,GG,0);

    //BL inflow velocity
    ldouble ucon[4];
    ldouble r=xx;
    ucon[1]=-sqrtl(2./r)*(1.-2./r);
    ucon[2]=ucon[3]=0.;
    trans2_coco(xxvec,ucon,ucon,BLCOORDS,MYCOORDS);
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

     //copying primitives with gdet taken into account
     for(iv=0;iv<NV;iv++)
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
  
