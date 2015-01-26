ldouble rho,mx,my,mz,m,E,uint,Fx,Fy,Fz,pLTE;  
ldouble xx,yy,zz;
ldouble uu[NV],xxvec[4],xxvecBL[4];
ldouble pp[NV],ppback[NV],T;

/***********************************************/
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

/***********************************************/
//background
pp[RHO]=RHOZERO;
pp[UU]=calc_PEQ_ufromTrho(TEMPZERO,RHOZERO);
pp[VX]=0.;
pp[VY]=0.;
pp[VZ]=0.;

#ifdef SPHERE
ldouble w1=exp(-((geom.xx-XBLOB1)*(geom.xx-XBLOB1) + (geom.yy-YBLOB1)*(geom.yy-YBLOB1) )/SIZEBLOB1/SIZEBLOB1);
pp[RHO]*=(1. + BLOBMAG1*w1);
#endif

#ifdef SQUARE
if(geom.xx>0.2*(MAXX-MINX) && geom.xx<0.5*(MAXX-MINX) && geom.yy>0.2*(MAXY-MINY) && geom.yy<0.5*(MAXY-MINY))
  pp[RHO]*=BLOBMAG1;
#endif





#ifdef RADIATION
ldouble temp=calc_PEQ_Tfromurho(pp[UU],pp[RHO]);

pp[EE0]=calc_LTE_EfromT(temp);
pp[FX0]=pp[FY0]=pp[FZ0]=0.;


#endif


/***********************************************/
//entropy
pp[5]=calc_Sfromu(pp[0],pp[1]);
p2u(pp,uu,&geom);
/***********************************************/

int iv;

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }


//entropy
update_entropy(ix,iy,iz,0);
set_cflag(0,ix,iy,iz,0);
