//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  

gdet_bc=get_g(g,3,4,ix,iy,iz);  
gdet_src=get_g(g,3,4,iix,iiy,iiz);
ldouble gg[4][5],GG[4][5],ggsrc[4][5],tup[4][4],tlo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
pick_T(tmuup,ix,iy,iz,tup);
pick_T(tmulo,ix,iy,iz,tlo);
ldouble xx=get_x(ix,0);

/**********************/

//radius
//  if(ix>=NX || ix<0) //analytical solution on both sides
if(ix>=NX || ix<0) //analytical solution at rout only
  {
    

ldouble rho0=0.1;
ldouble r0=100.;
ldouble u0=0.1;
ldouble r=xx;

 ldouble rho=rho0*r0/r;
 ldouble uint=u0*r0*r0/r/r/(1.-2./r);

 pp[0]=rho;
 pp[1]=uint;

    pp[2]=0.;
    pp[3]=0.;
    pp[4]=0.;
    pp[5]=calc_Sfromu(pp[0],pp[1]);


    p2u(pp,uu,gg,GG);
    
    return 0.;
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
  
