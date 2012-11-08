//KORAL - p2u.c
//primitives to conserved conversion

#include "ko.h"

//**********************************************************************
//**********************************************************************
//**********************************************************************
//primitive to conserved converter
int
p2u(ldouble *p, ldouble *u, ldouble g[][5], ldouble eup[][4], ldouble elo[][4])
{
  ldouble gtt=g[0][0];
  ldouble gtph=g[0][3];
  ldouble grr=g[1][1];
  ldouble gthth=g[2][2];
  ldouble gphph=g[3][3];

  ldouble gdet=g[3][4];

  ldouble rho=p[0];
  ldouble uu=p[1];
  ldouble vr=p[2];
  ldouble vth=p[3];
  ldouble vph=p[4];
  ldouble S=p[5];

  //************************************
  //************************************
  //************************************
  //radiation part
  //************************************
  //************************************
  //************************************

#ifdef RADIATION

  ldouble E=p[6];
  ldouble F[3]={p[7],p[8],p[9]};
  ldouble Rij[4][4];

  calc_Rij(p,Rij);
  boost22_ff2zamo(Rij,Rij,p,g,eup);
  trans22_zamo2lab(Rij,Rij,g,elo);  
  indices_2221(Rij,Rij,g);

  u[6]=Rij[0][0]; //R^t_t
  u[7]=Rij[0][1]; //R^t_i
  u[8]=Rij[0][2];
  u[9]=Rij[0][3];
 
#endif

  //************************************
  //************************************
  //************************************
  //hydro part
  //************************************
  //************************************
  //************************************
 
  ldouble ut2=-1./(gtt + 2.*vph*gtph + vr*vr*grr + vph*vph*gphph + vth*vth*gthth);
  if(ut2<0.)
    {
      my_err("ut2.lt.0 in p2u\n"); ut2=0.;
    }

  ldouble ut=sqrtl(ut2);
  ldouble rhout = rho*ut;
  ldouble Sut;

  if(uu<0. || rho<0.)
    Sut=S*ut;
  else
    Sut=calc_Sfromu(rho,uu)*ut;

  ldouble Tttt=rhout*(1+ut*(gtt+vph*gtph))+GAMMA*uu*ut2*(gtt+vph*gtph)+uu*(GAMMA-1.);  
  ldouble Ttr=(rho+GAMMA*uu)*ut2*vr*grr;
  ldouble Ttth=(rho+GAMMA*uu)*ut2*vth*gthth;
  ldouble Ttph=(rho+GAMMA*uu)*ut2*(gtph+vph*gphph);

  u[0]=rhout;
  u[1]=Tttt;
  u[2]=Ttr;
  u[3]=Ttth;
  u[4]=Ttph;
  u[5]=Sut;

  return 0.;
}

