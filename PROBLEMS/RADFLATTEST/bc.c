//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);



/**********************/

//radius
if(ix>=NX) //total boundary, properties of the galaxy
  {    
    for(iv=0;iv<NV;iv++)
      { 
	pp[iv]=get_u(p,iv,NX-1,iiy,iiz);
      } 

    if(pp[VX]<0.) pp[VX]=0.;
    if(pp[FX0]<0.) pp[FX0]=0.;

    p2u(pp,uu,&geom);
 
    return 0.;
  }

 else if(ix<0) //outflow near BH
   {
     for(iv=0;iv<NVMHD;iv++)
       { 
	 pp[iv]=get_u(p,iv,0,iiy,iiz);
       } 

     /*
     if(geom.xx<0.3)
       pp[0]=1000.;
     else
       pp[0]=.001;
     */
     pp[0]=.001;
     pp[1]=.1;
     pp[2]=0.;
     pp[3]=0.;
     pp[4]=0.;
     pp[5]=calc_Sfromu(pp[0],pp[1]);

#ifdef RADIATION
     pp[6]=.1;
     pp[7]=0.9;
     pp[8]=0.;
     pp[9]=0.;
#endif	 
    
     p2u(pp,uu,&geom);
 
     return 0;
   }

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
  
