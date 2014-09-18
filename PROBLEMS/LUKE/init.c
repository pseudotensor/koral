ldouble uu[NV];
ldouble pp[NV];

/***********************************************/
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

/***********************************************/
//background
pp[RHO]=RHO_AMB + 100.*RHO_AMB*exp(-(geom.xx*geom.xx+geom.yy*geom.yy+geom.zz*geom.zz)/0.02);
pp[UU]=U_AMB*RHO_AMB/pp[RHO];
pp[VX]=0.;
if(pp[RHO]>1.1*RHO_AMB) pp[VX]=0.1;
pp[VY]=0.;
pp[VZ]=0.;

/***********************************************/
//entropy
pp[5]=calc_Sfromu(pp[0],pp[1]);
p2u(pp,uu,&geom);
/***********************************************/

//print_primitives(pp); getch();

int iv;

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }
