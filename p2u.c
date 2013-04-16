//KORAL - p2u.c
//primitives to conserved conversion

#include "ko.h"

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates conserved in given cell using global array p[]
int
calc_conserved(int ix,int iy,int iz)
{
  int iv;
  ldouble uu[NV],pp[NV];
  ldouble gg[4][5],GG[4][5],tlo[4][4],tup[4][4];
  
  pick_g(ix,iy,iz,gg);
  pick_G(ix,iy,iz,GG);

  for(iv=0;iv<NV;iv++)
    {
      pp[iv]=get_u(p,iv,ix,iy,iz);
    }

  p2u(pp,uu,gg,GG);


  for(iv=0;iv<NV;iv++)
    {
      set_u(u,iv,ix,iy,iz,uu[iv]);
    }

  return 0;
}
 

//**********************************************************************
//**********************************************************************
//**********************************************************************
//primitive to conserved converter
int
p2u(ldouble *p, ldouble *u, ldouble g[][5], ldouble G[][5])
{
  ldouble rho=p[0];
  ldouble uu=p[1];
  ldouble vcon[4],ucon[4],ucov[4];
  vcon[1]=p[2];
  vcon[2]=p[3];
  vcon[3]=p[4];
  vcon[0]=0.;
  ldouble S=p[5];

  //converting to 4-velocity

  conv_vels(vcon,ucon,VELPRIM,VEL4,g,G);
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
 

  ldouble ut=ucon[0];
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

#ifdef EDDINGTON_APR
  int irf,ii;
  int pp10[10];
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],(*tlo)[4],(*tup)[4];
  gg=geom->gg;
  GG=geom->GG;
  tlo=geom->tlo;
  tup=geom->tup;

  ldouble Rij[4][4];
  for(ii=0;ii<NVHD;ii++)
    {
      pp10[ii]=pp[ii];
      uu10[ii]=uu[ii];
    }
  
  for(irf=0;irf<NRF;irf++)
    {
      pp10[EE(0)]=p[EE(irf)];
      pp10[FX(0)]=p[FX(irf)];
      pp10[FY(0)]=p[FY(irf)];
      pp10[FZ(0)]=p[FZ(irf)];

      calc_Rij_ff(pp10,Rij);  
      trans22_on2cc(Rij,Rij,tlo);  
      boost22_ff2lab(Rij,Rij,pp,gg,GG); 
      indices_2221(Rij,Rij,gg);  

      u[EE(irf)]=Rij[0][0];
      u[FX(irf)]=Rij[0][1];
      u[FY(irf)]=Rij[0][2];
      u[FZ(irf)]=Rij[0][3];
    }
  return 0;
#endif

#ifdef LABRADFLUXES
  
  u[6]=p[6]; //R^t_t
  u[7]=p[7]; //R^t_i
  u[8]=p[8];
  u[9]=p[9];
  return 0;
 
#endif
  
  int i,j,irf;

  for(irf=0;irf<NRF;irf++)
    {
      ldouble Erf=p[EE(irf)];

      //relative four-velocity
      ldouble urf[4];
      urf[0]=0.;
      urf[1]=p[FX(irf)];
      urf[2]=p[FY(irf)];
      urf[3]=p[FZ(irf)];

      //converting to lab four-velocity
      conv_vels(urf,urf,VELPRIMRAD,VEL4,g,G);
  
      ldouble Rtop[4];
      Rtop[0]=4./3.*Erf*urf[0]*urf[0] + 1./3.*Erf*G[0][0]; //R^t_t
      Rtop[1]=4./3.*Erf*urf[0]*urf[1] + 1./3.*Erf*G[0][1];
      Rtop[2]=4./3.*Erf*urf[0]*urf[2] + 1./3.*Erf*G[0][2];
      Rtop[3]=4./3.*Erf*urf[0]*urf[3] + 1./3.*Erf*G[0][3];

      indices_21(Rtop,Rtop,g); //R^t_mu

      u[EE(irf)]=Rtop[0]; //R^t_t
      u[FX(irf)]=Rtop[1]; //R^t_i
      u[FY(irf)]=Rtop[2];
      u[FZ(irf)]=Rtop[3];
    }

  return 0;
}
