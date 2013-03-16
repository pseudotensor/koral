
int ix,iy,iz;

ldouble pp[NV],uu[NV];

ix=IXDOT;
iy=IYDOT;
iz=IZDOT;

ldouble gg[4][5],GG[4][5];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

pp[0]=1.;
pp[1]=1.;
pp[2]=0.;
pp[3]=0.;
pp[4]=0.;
pp[5]=calc_Sfromu(pp[0],pp[1]);
pp[6]=LTEFACTOR*calc_LTE_Efromurho(pp[0],pp[1]);
pp[7]=0.;
pp[8]=0.;
pp[9]=0.;

if(ix==IXDOT && iy==IYDOT && iz==IZDOT)
  {
    if(NZ==1)
      pp[6]*=100.;
    else
      pp[6]*=10000.;
    pp[8]=FYDOT*pp[6];
  }

prad_ff2lab(pp,pp,&geom);

p2u(pp,uu,gg,GG);

/**************************/

int iv;

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }

//entropy
update_entropy(ix,iy,iz,0);
