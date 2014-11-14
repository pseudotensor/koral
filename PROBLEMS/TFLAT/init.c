//sets the initial conditions
//called from a loop going over ix,iy,iz

/***********************************************/
//structure of geometry
//holds coordinates and metric, see ko.h
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

/***********************************************/
//vectors of primitives and conserved
ldouble pp[NV],uu[NV];

/***********************************************/
//domain
pp[RHO]=RHOAMB;
pp[UU]=UUAMB;
pp[VZ]=0.;
pp[VY]=0.;
pp[VX]=0.;

/***********************************************/
//calculate entropy from rho & uint
pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU]);

/***********************************************/
//convert primitives to conserved
p2u(pp,uu,&geom);

/***********************************************/
//save to memory
int iv;
for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }
