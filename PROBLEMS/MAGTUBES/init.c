
ldouble uu[NV];
ldouble pp[NV];

//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

/***********************************************/
if(geom.xx<0.)
  {
    pp[RHO]=RHOL;
    pp[UU]=PL/GAMMAM1;
    pp[VX]=VXL;
    pp[VY]=VYL;
    pp[VZ]=VZL;

#ifdef MAGNFIELD
    pp[B1]=0.;
    pp[B2]=BYL;
    pp[B3]=BZL;
#endif
  }
 else
  {
    pp[RHO]=RHOR;
    pp[UU]=PR/GAMMAM1;
    pp[VX]=VXR;
    pp[VY]=VYR;
    pp[VZ]=VZR;
    
#ifdef MAGNFIELD
    pp[B1]=0.;
    pp[B2]=BYR;
    pp[B3]=BZR;
#endif
   }
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
