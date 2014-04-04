//sets the initial conditions
//called from a loop going over ix,iy,iz

/***********************************************/
//structure of geometry
//holds coordinates and metric, see ko.h
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
struct geometry geomcart;
fill_geometry_arb(ix,iy,iz,&geomcart,MINKCOORDS);


/***********************************************/
//vectors of primitives and conserved
ldouble pp[NV],uu[NV];

/***********************************************/
//ambient medium
ldouble BCX = RBLOB;
ldouble BCZ = 0.;
ldouble dist = sqrt((geomcart.xx-BCX)*(geomcart.xx-BCX)+(geomcart.zz-BCZ)*(geomcart.zz-BCZ));
//ldouble omega = sqrt(1./(RBLOB*RBLOB*RBLOB));
ldouble omega = sqrt(1./(geom.xx*geom.xx*geom.xx));
pp[RHO]=RHOAMB + (BRHO - RHOAMB) * exp(-dist*dist / BW / BW );
pp[UU]=UUAMB;
ldouble OMEGA1=sqrt(omega*omega/(-geom.gg[0][0])/(1+omega*omega*geom.gg[3][3]/geom.gg[0][0]));
if(pp[RHO]>2.*RHOAMB)
  {
    pp[VZ]=OMEGA1;
    pp[VY]=0.;
    pp[VX]=0.;
  }
 else
   pp[VY]=pp[VX]=pp[VZ]=0.;

//magn field
#ifdef MAGNFIELD
pp[B1]=pp[B2]=pp[B3]=0.; 

pp[B3]=my_max(0.,pp[RHO]-.1*BRHO);
    
#endif

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
