
#ifndef FORGETDOTS

int ix,iy,iz;
/**************************/

int iv;
int irf=0;
ldouble pp[NV],uu[NV];

pp[0]=1.;
pp[1]=1.;
pp[2]=0.;
pp[3]=0.;
pp[4]=0.;
pp[5]=calc_Sfromu(pp[0],pp[1]);

ix=IXDOT1;
iy=IYDOT1;
iz=IZDOT1;

ldouble gg[4][5],GG[4][5];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);


pp[6]=LTEFACTOR*calc_LTE_Efromurho(pp[0],pp[1]);
pp[6]*=1000.;
pp[FX(irf)]=FXDOT1*pp[6];
pp[FY(irf)]=FYDOT1*pp[6];
pp[FZ(irf)]=FZDOT1*pp[6];

for(irf=1;irf<NRF;irf++)
  {
    pp[EE(irf)]=EEFLOOR;
    pp[FX(irf)]=0.;
    pp[FY(irf)]=0.;
    pp[FZ(irf)]=0.;
  }


prad_ff2lab(pp,pp,&geom);

p2u(pp,uu,&geom);


#ifdef MULTIRADFLUID

//print_Nvector(uu,NV);
redistribute_radfluids(pp,uu,&geom);
//print_Nvector(uu,NV);getchar();
u2p_rad(uu,pp,&geom,&irf);

#endif

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }

//entropy
update_entropy(ix,iy,iz,0);


ix=IXDOT2;
iy=IYDOT2;
iz=IZDOT2;
irf=0;

pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);

fill_geometry(ix,iy,iz,&geom);


pp[6]=LTEFACTOR*calc_LTE_Efromurho(pp[0],pp[1]);
pp[6]*=1000.;
pp[FX(irf)]=FXDOT2*pp[6];
pp[FY(irf)]=FYDOT2*pp[6];
pp[FZ(irf)]=FZDOT2*pp[6];

for(irf=1;irf<NRF;irf++)
  {
    pp[EE(irf)]=EEFLOOR;
    pp[FX(irf)]=0.;
    pp[FY(irf)]=0.;
    pp[FZ(irf)]=0.;
  }


prad_ff2lab(pp,pp,&geom);

p2u(pp,uu,&geom);


#ifdef MULTIRADFLUID

//print_Nvector(uu,NV);
redistribute_radfluids(pp,uu,&geom);
//print_Nvector(uu,NV);getchar();
u2p_rad(uu,pp,&geom,&irf);

#endif


for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }

//entropy
update_entropy(ix,iy,iz,0);

#endif
