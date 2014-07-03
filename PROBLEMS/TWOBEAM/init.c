ldouble rho,mx,my,mz,m,E,uint,Fx,Fy,Fz,pLTE;  
ldouble xx,yy,zz;
ldouble uu[NV],xxvec[4],xxvecBL[4];
ldouble pp[NV],ppback[NV],T;

/***********************************************/
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

/***********************************************/
//background
pp[RHO]=1.0;
pp[UU]=1.0;
pp[VX]=0.;
pp[VY]=0.;
pp[VZ]=0.;



#ifdef RADIATION

pp[EE0]=1.0;
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
