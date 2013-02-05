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
    pp[0]=PAR_D;
    pp[1]=PAR_U;

    //calculating entropy
    pp[5]=calc_Sfromu(pp[0],pp[1]);	      

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
  
