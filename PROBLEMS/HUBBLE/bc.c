//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  

/**********************/
gdet_bc=get_g(g,3,4,ix,iy,iz);  
ldouble gg[4][5],GG[4][5],ggsrc[4][5],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
pick_g(ix,iy,iz,GG);
pick_T(emuup,ix,iy,iz,eup);
pick_T(emulo,ix,iy,iz,elo);
ldouble xx=get_x(ix,0);
ldouble xx0,xx1;

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

if(ix<0)
  {
    xx0=get_x(0,0);
    xx1=get_x(1,0);
    for(iv=0;iv<NV;iv++)
      {
	pp[iv]=get_u(p,iv,0,iy,iz)+(get_u(p,iv,1,iy,iz)-get_u(p,iv,0,iy,iz))*(xx-xx0)/(xx1-xx0);
      }

    //    printf("%d %e %e %e %e %e\n",ix,get_u(p,iv,1,iy,iz),get_u(p,iv,0,iy,iz),xx,xx0,xx1);getchar();
    p2u(pp,uu,&geom);
    return 0;
  }

if(ix>=NX)
  {
    xx0=get_x(NX-1,0);
    xx1=get_x(NX-2,0);
    for(iv=0;iv<NV;iv++)
      {
	pp[iv]=get_u(p,iv,NX-1,iy,iz)+(get_u(p,iv,NX-2,iy,iz)-get_u(p,iv,NX-1,iy,iz))*(xx-xx0)/(xx1-xx0);
      }
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
    pp[iv]=get_u(p,iv,iix,iiy,iiz);      
  }

p2u(pp,uu,&geom);

  
return 0;
  
