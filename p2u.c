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
  
  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);

  for(iv=0;iv<NV;iv++)
    {
      pp[iv]=get_u(p,iv,ix,iy,iz);
    }

  p2u(pp,uu,&geom);

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
p2u(ldouble *p, ldouble *u, void *ggg)
{
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  gdet=geom->gdet;
  GG=geom->GG;
  gdetu=gdet;

#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif


  ldouble rho=p[0];
  ldouble uu=p[1];
  ldouble vcon[4],vcov[4],ucon[4],ucov[4];
  ldouble bcon[4]={0.,0.,0.,0.},bcov[4]={0.,0.,0.,0.},bsq=0.;
  vcon[1]=p[2];
  vcon[2]=p[3];
  vcon[3]=p[4];
  vcon[0]=0.;
  ldouble S=p[5];

  //converting to 4-velocity

  conv_vels(vcon,ucon,VELPRIM,VEL4,gg,GG);
  //indices_21(ucon,ucov,gg);
  //print_4vector(ucov);
  conv_velscov(vcon,ucov,VELPRIM,VEL4,gg,GG);
  //print_4vector(ucov);

#ifdef MAGNFIELD
  calc_bcon_4vel(p,ucon,ucov,bcon);
  indices_21(bcon,bcov,gg); 
  bsq = dot(bcon,bcov);
#endif

  //************************************
  //************************************
  //************************************
  //radiation part
  //************************************
  //************************************
  //************************************

#ifdef RADIATION
  
  p2u_rad(p,u,ggg);
 
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

  //S=pp[5] updated appropriately in u2p_hot, u2p_entropy and floors
  Sut=S*ut;

  ldouble pre=(GAMMA-1.)*uu; 
  ldouble w=rho+uu+pre;
  ldouble eta=w+bsq;
  ldouble ptot=pre+0.5*bsq;

  //~T^i_j 
  ldouble Tttt=rhout + eta*ucon[0]*ucov[0] + ptot - bcon[0]*bcov[0];
  ldouble Ttr =eta*ucon[0]*ucov[1] - bcon[0]*bcov[1];
  ldouble Ttth =eta*ucon[0]*ucov[2] - bcon[0]*bcov[2];
  ldouble Ttph =eta*ucon[0]*ucov[3] - bcon[0]*bcov[3];

  /*
  printf("%e %e %e %e -> %.20e\n",eta*ucon[0]*ucov[3] - bcon[0]*bcov[3],eta,ucon[0],ucov[3],dot(ucon,ucov));
  print_4vector(vcon);
  print_4vector(ucon);
  print_metric(gg);
  
  
 
 int i,j,k;
  /*for(i=0;i<4;i++)
    {
      ucov[i]=0.;
      for(k=0;k<4;k++)
	{
	  ucov[i]+=ucon[k]*gg[i][k];
	  if(i==3) printf("+ %e (%d %d)\n",ucon[k]*gg[i][k],i,k);
	}	  
    }
 
 //print_4vector(ucov);

 indices_21(vcon,vcov,gg);
 //print_4vector(vcov);

 // print_4vector(vcon);getch();
 ldouble qsq=0.;
 for(i=1;i<4;i++)
   for(j=1;j<4;j++)
     qsq+=vcon[i]*vcon[j]*gg[i][j];
 ldouble gamma2=1.+qsq;
 ldouble alpha2=-1./GG[0][0];
 ldouble gamma=sqrt(gamma2);
 ldouble alpha=sqrt(alpha2);
 for(i=0;i<4;i++)
   ucov[i]=vcov[i]-alpha*gamma*delta(0,i);

  print_4vector(ucov);
  */

	 

  u[0]=gdetu*rhout;
  u[1]=gdetu*Tttt;
  u[2]=gdetu*Ttr;
  u[3]=gdetu*Ttth;
  u[4]=gdetu*Ttph;
  u[5]=gdetu*Sut;


#ifdef TRACER
  ldouble tracerut=p[TRA]*ut;
  u[TRA]= gdetu*tracerut;
#endif

  //************************************
  //************************************
  //************************************
  //magnetic part
  //************************************
  //************************************
  //************************************
 
#ifdef MAGNFIELD
  u[B1]=gdetu*p[B1];
  u[B2]=gdetu*p[B2];
  u[B3]=gdetu*p[B3];
#endif
  
 

  return 0.;
}

/********************************************************/
/**** converts radiative primitives xs************************/
/********************************************************/
/********************************************************/
int p2u_rad(ldouble *pp,ldouble *uu,void *ggg)
{
  int i,j,irf;

   struct geometry *geom
   = (struct geometry *) ggg;

   ldouble (*gg)[5],(*GG)[5],(*tlo)[4],(*tup)[4],gdet,gdetu;
  gg=geom->gg;
  GG=geom->GG;
  tlo=geom->tlo;
  tup=geom->tup;
  gdet=geom->gdet;gdetu=gdet;

#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif
 
#ifdef LABRADFLUXES
  
  uu[EE0]=gdetu*pp[EE0]; //R^t_t
  uu[FX0]=gdetu*pp[FX0]; //R^t_i
  uu[FY0]=gdetu*pp[FY0];
  uu[FZ0]=gdetu*pp[FZ0];
  return 0;
 
#endif

#ifdef EDDINGTON_APR
  int ii,jj;
  ldouble Rij[4][4],h[4][4];
  ldouble ucov[4],ucon[4]={0,pp[2],pp[3],pp[4]};
  conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);
  indices_21(ucon,ucov,gg);
  ldouble EE=pp[EE0];
  ldouble Fcon[4]={0.,pp[FX0],pp[FY0],pp[FZ0]};
  Fcon[0]=-1./ucov[0]*(Fcon[1]*ucov[1]+Fcon[2]*ucov[2]+Fcon[3]*ucov[3]); //F^0 u_0 = - F^i u_i
  //projection tensor
  for(ii=0;ii<4;ii++)
    for(jj=0;jj<4;jj++)
      h[ii][jj]=GG[ii][jj] + ucon[ii]*ucon[jj];
  //Fragile's formula
  for(ii=0;ii<4;ii++)
    for(jj=0;jj<4;jj++)
      Rij[ii][jj]=EE*ucon[ii]*ucon[jj] + Fcon[ii]*ucon[jj] + Fcon[jj]*ucon[ii] + 1./3.*EE*delta(ii,jj)*h[ii][jj];
  indices_2221(Rij,Rij,gg);

  //  print_Nvector(p,NV);
  //  print_4vector(ucon);
  //  print_4vector(Fcon);
      
  uu[EE0]=gdetu*Rij[0][0];
  uu[FX0]=gdetu*Rij[0][1];
  uu[FY0]=gdetu*Rij[0][2];
  uu[FZ0]=gdetu*Rij[0][3];

  return 0;
#endif
  
  //M1
  for(irf=0;irf<NRF;irf++)
    {
      ldouble Erf=pp[EE(irf)];

      //relative four-velocity
      ldouble urf[4];
      urf[0]=0.;
      urf[1]=pp[FX(irf)];
      urf[2]=pp[FY(irf)];
      urf[3]=pp[FZ(irf)];

      //converting to lab four-velocity
      conv_vels(urf,urf,VELPRIMRAD,VEL4,gg,GG);
  
      ldouble Rtopp[4];
      Rtopp[0]=4./3.*Erf*urf[0]*urf[0] + 1./3.*Erf*GG[0][0]; //R^t_t
      Rtopp[1]=4./3.*Erf*urf[0]*urf[1] + 1./3.*Erf*GG[0][1];
      Rtopp[2]=4./3.*Erf*urf[0]*urf[2] + 1./3.*Erf*GG[0][2];
      Rtopp[3]=4./3.*Erf*urf[0]*urf[3] + 1./3.*Erf*GG[0][3];

      indices_21(Rtopp,Rtopp,gg); //R^t_mu

      uu[EE(irf)]=gdetu*Rtopp[0]; //R^t_t
      uu[FX(irf)]=gdetu*Rtopp[1]; //R^t_i
      uu[FY(irf)]=gdetu*Rtopp[2];
      uu[FZ(irf)]=gdetu*Rtopp[3];
    }

  return 0;
}
