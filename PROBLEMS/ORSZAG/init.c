
ldouble uu[NV];
ldouble pp[NV];

//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

/***********************************************/

ldouble x=geom.xx;
ldouble y=geom.yy;

pp[RHO]=25./(36.*M_PI);
pp[UU]=(5./12/M_PI)/GAMMAM1/CSCALE/CSCALE;
pp[VX]=-sin(2.*M_PI*y)/CSCALE;
pp[VY]=sin(2.*M_PI*x)/CSCALE;
pp[VZ]=0.;

#ifdef MAGNFIELD
//TODO: put A_z here and convert to B's lateron
pp[B1]=0.;
pp[B2]=0.;
pp[B3]=(BZERO/4./M_PI)*cos(4.*M_PI*x)+(BZERO/2./M_PI)*cos(2.*M_PI*y);
pp[B3]/=CSCALE;
#endif
//entropy
pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU]);

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
