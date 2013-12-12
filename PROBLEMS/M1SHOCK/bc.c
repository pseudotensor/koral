//KORAL - bc.c
//returns problem specific BC
//**********************
//included in:
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) 
//from problem.c
//**********************

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  
gdet_bc=get_g(g,3,4,ix,iy,iz);  
//gdet_src=get_g(g,3,4,iix,iiy,iiz);
ldouble gg[4][5],GG[4][5],ggsrc[4][5],tup[4][4],tlo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
pick_T(tmuup,ix,iy,iz,tup);
pick_T(tmulo,ix,iy,iz,tlo);
ldouble xx=get_x(ix,0);
ldouble Fx,Fy,Fz;

/**********************/
/**********************/
/**********************/
if(ix<0 || ix>=NX)
  {
    Fx=Fz=Fy=0.;
    pp[0]=1.;
    pp[1]=calc_PEQ_ufromTrho(TEMPIN,pp[RHO]);
    ldouble ucon[4]={0.,0.,0.,0.};
    conv_vels(ucon,ucon,VEL3,VELPRIM,gg,GG);
    pp[2]=ucon[1]; 
    pp[3]=ucon[2];
    pp[4]=ucon[3]; 

    pp[5]=calc_Sfromu(pp[RHO],pp[UU]);

    pp[EE0]=calc_LTE_EfromT(TEMPOUT);
    pp[FX0]=Fx;
    pp[FY0]=Fy;
    pp[FZ0]=Fz; 

    prad_ff2lab(pp,pp,&geom);
    p2u(pp,uu,&geom);	

    return 0.;
  }

iix=ix;
iiz=iz;
iiy=iy;
//periodic
while(iiy<0)    iiy+=NY;
while(iiy>=NY)    iiy-=NY; 
//periodic
while(iiz<0)    iiz+=NZ;
while(iiz>=NZ)    iiz-=NZ; 
 
for(iv=0;iv<NV;iv++)
  {
    pp[iv]=get_u(p,iv,iix,iiy,iiz);      
  }

//converting to conserved
p2u(pp,uu,&geom);
  
return 0;
  
