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
ldouble dist = sqrt((geom.xx - BCX)*(geom.xx - BCX) + 
		    (geom.yy - BCY)*(geom.yy - BCY) + 
		    (geom.zz - BCZ)*(geom.zz - BCZ)) ;
pp[RHO]=RHOAMB + (BRHO - RHOAMB) * exp(-dist*dist / BW / BW);
//printf("%e %e %e\n",dist, exp(-dist*dist / BW / BW),(BRHO - RHOAMB) * exp(-dist*dist / BW / BW)); getch();
pp[UU]=UUAMB;
pp[VZ]=BVZ;
pp[VY]=BVY;
pp[VX]=BVX;

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
