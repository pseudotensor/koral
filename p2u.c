//KORAL - p2u.c
//primitives to conserved conversion

#include "ko.h"

//**********************************************************************
//**********************************************************************
//**********************************************************************
//primitive to conserved converter
int
p2u(ldouble *p, ldouble *u, ldouble g[][5], ldouble G[][5])
{
  ldouble gtt=g[0][0];
  ldouble gtph=g[0][3];
  ldouble grr=g[1][1];
  ldouble gthth=g[2][2];
  ldouble gphph=g[3][3];

  ldouble gdet=g[3][4];

  ldouble rho=p[0];
  ldouble uu=p[1];
  ldouble vcon[4],ucon[4],ucov[4];
  vcon[1]=p[2];
  vcon[2]=p[3];
  vcon[3]=p[4];
  vcon[0]=0.;
  ldouble S=p[5];

  //converting to 4-velocity
  conv_vels(vcon,ucon,VELR,VEL4,g,G);
  //conv_vels(vcon,ucon,VEL3,VEL4,g,G);
  indices_21(ucon,ucov,g);

  //************************************
  //************************************
  //************************************
  //radiation part
  //************************************
  //************************************
  //************************************

#ifdef RADIATION
  
  p2u_rad(p,u,g,G);
 
#endif

  //************************************
  //************************************
  //************************************
  //hydro part
  //************************************
  //************************************
  //************************************
 

  //TODO: generalize
  ldouble ut=ucon[0];
  ldouble vr=p[2];
  ldouble vth=p[3];
  ldouble vph=p[4];

  ldouble rhout = rho*ut;
  ldouble Sut;

  if(uu<0. || rho<0.)
    Sut=S*ut;
  else
    Sut=calc_Sfromu(rho,uu)*ut;

  ldouble pre=(GAMMA-1.)*uu;
  ldouble w=rho+uu+pre;
  ldouble Tttt=rhout + w*ucon[0]*ucov[0] + pre;
  ldouble Ttr =w*ucon[0]*ucov[1];
  ldouble Ttth =w*ucon[0]*ucov[2];
  ldouble Ttph =w*ucon[0]*ucov[3];
   
  u[0]=rhout;
  u[1]=Tttt;
  u[2]=Ttr;
  u[3]=Ttth;
  u[4]=Ttph;
  u[5]=Sut;

  return 0.;
}

/********************************************************/
/**** converts radiative primitives xs************************/
/********************************************************/
/********************************************************/
int p2u_rad(ldouble *p,ldouble *u,ldouble g[][5],ldouble G[][5])
{
  int i,j;
  ldouble Erf=p[6];

  //relative four-velocity
  ldouble urf[4];
  urf[0]=0.;
  urf[1]=p[7];
  urf[2]=p[8];
  urf[3]=p[9];

  //converting to lab four-velocity
  conv_vels(urf,urf,VELR,VEL4,g,G);
  
  ldouble Rtop[4];
  Rtop[0]=4./3.*Erf*urf[0]*urf[0] + 1./3.*Erf*G[0][0]; //R^t_t
  Rtop[1]=4./3.*Erf*urf[0]*urf[1] + 1./3.*Erf*G[0][1];
  Rtop[2]=4./3.*Erf*urf[0]*urf[2] + 1./3.*Erf*G[0][2];
  Rtop[3]=4./3.*Erf*urf[0]*urf[3] + 1./3.*Erf*G[0][3];

  indices_21(Rtop,Rtop,g); //R^t_mu

  u[6]=Rtop[0]; //R^t_t
  u[7]=Rtop[1]; //R^t_i
  u[8]=Rtop[2];
  u[9]=Rtop[3];

  return 0;
}
