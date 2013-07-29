
ldouble uu[NV];
ldouble pp[NV];

//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

/***********************************************/

pp[RHO]=1.;
pp[UU]=0.001;
pp[VX]=pp[VY]=pp[VZ]=0.;

#ifdef MAGNFIELD
pp[B1]=pp[B2]=pp[B3]=0.;
#endif

//entropy
pp[5]=calc_Sfromu(pp[RHO],pp[UU]);

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
