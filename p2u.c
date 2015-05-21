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
  p2u_mhd(p,u,ggg);

#ifdef RADIATION
  p2u_rad(p,u,ggg);
#endif

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//primitive to conserved converter
int
p2u_mhd(ldouble *p, ldouble *u, void *ggg)
{

#ifdef NONRELMHD
  p2u_mhd_nonrel(p,u,ggg);
  return 0;
#endif

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

  conv_vels_both(vcon,ucon,ucov,VELPRIM,VEL4,gg,GG);

#ifdef MAGNFIELD
  calc_bcon_4vel(p,ucon,ucov,bcon);
  indices_21(bcon,bcov,gg); 
  bsq = dot(bcon,bcov);
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

  //S=pp[5] updated appropriately in u2p_hot, u2p_entropy and floors so from outside
  Sut=S*ut;

  ldouble pre=(GAMMA-1.)*uu; 
  ldouble w=rho+uu+pre;
  ldouble eta=w+bsq;
  ldouble etap = uu+pre+bsq; //eta-rho
  ldouble ptot=pre+0.5*bsq;
  
  //this computes utp1=1+u_t
  ldouble utp1;
  utp1=calc_utp1(vcon,ucon,geom);

  ldouble Tttt=etap*ucon[0]*ucov[0] + rho*ucon[0]*utp1 + ptot - bcon[0]*bcov[0];
  ldouble Ttr =eta*ucon[0]*ucov[1] - bcon[0]*bcov[1];
  ldouble Ttth =eta*ucon[0]*ucov[2] - bcon[0]*bcov[2];
  ldouble Ttph =eta*ucon[0]*ucov[3] - bcon[0]*bcov[3];

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

//**********************************************************************
//**********************************************************************
//**********************************************************************
//primitive to conserved converter - non-relativistic!
int
p2u_mhd_nonrel(ldouble *p, ldouble *u, void *ggg)
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

  conv_vels_both(vcon,ucon,ucov,VELPRIM,VEL4,gg,GG);

  ldouble v2=dot3nr(ucon,ucov);

#ifdef MAGNFIELD
  calc_bcon_4vel(p,ucon,ucov,bcon);
  indices_21(bcon,bcov,gg); 
  bsq = dot(bcon,bcov);
#endif


  //************************************
  //************************************
  //************************************
  //hydro part
  //************************************
  //************************************
  //************************************
 
  ldouble Ttt=-(uu + bsq/2. + rho*v2/2.);
  ldouble Tttt=Ttt;
  ldouble Ttr =rho*ucov[1];
  ldouble Ttth =rho*ucov[2];
  ldouble Ttph =rho*ucov[3];

  u[0]=gdetu*rho;
  u[1]=gdetu*Tttt;
  u[2]=gdetu*Ttr;
  u[3]=gdetu*Ttth;
  u[4]=gdetu*Ttph;
  u[5]=gdetu*S;


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
 
  //M1
  //for(irf=0;irf<NRF;irf++)
    {
      ldouble Erf=pp[EE];

      //relative four-velocity
      ldouble urf[4];
      urf[0]=0.;
      urf[1]=pp[FX];
      urf[2]=pp[FY];
      urf[3]=pp[FZ];

      //converting to lab four-velocity
      conv_vels(urf,urf,VELPRIMRAD,VEL4,gg,GG);
  
      ldouble Rtopp[4];
      Rtopp[0]=4./3.*Erf*urf[0]*urf[0] + 1./3.*Erf*GG[0][0]; //R^t_t
      Rtopp[1]=4./3.*Erf*urf[0]*urf[1] + 1./3.*Erf*GG[0][1];
      Rtopp[2]=4./3.*Erf*urf[0]*urf[2] + 1./3.*Erf*GG[0][2];
      Rtopp[3]=4./3.*Erf*urf[0]*urf[3] + 1./3.*Erf*GG[0][3];

      indices_21(Rtopp,Rtopp,gg); //R^t_mu

      uu[EE]=gdetu*Rtopp[0]; //R^t_t
      uu[FX]=gdetu*Rtopp[1]; //R^t_i
      uu[FY]=gdetu*Rtopp[2];
      uu[FZ]=gdetu*Rtopp[3];

      #ifdef NCOMPTONIZATION
      uu[NF]=gdetu*pp[NF]*urf[0];
      #endif
    }

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//takes primitives and calculates quantities that are averaged for the avg file
int
p2avg(int ix,int iy,int iz,ldouble *avg)
{
  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  struct geometry geomout;
  fill_geometry_arb(ix,iy,iz,&geomout,OUTCOORDS);
  struct geometry geoml;
  fill_geometry_face(ix,iy,iz,0,&geoml);
  struct geometry geomoutl;
  fill_geometry_face_arb(ix,iy,iz,0,&geomoutl,OUTCOORDS);

 
  int iv,iv2;ldouble pp[NV],uu[NV];
  for(iv=0;iv<NV;iv++)
    {
      uu[iv]=get_u(u,iv,ix,iy,iz); //conserved 
      pp[iv]=get_u(p,iv,ix,iy,iz); //primitives 

      avg[iv]=pp[iv]; //first NV slots in pavg are regular primitives in MYCOORDS!
   }

  //primitives to OUTCOORDS
#ifdef RADIATION
  trans_prad_coco(pp, pp, MYCOORDS,OUTCOORDS, geom.xxvec,&geom,&geomout);
#endif
  trans_pmhd_coco(pp, pp, MYCOORDS,OUTCOORDS, geom.xxvec,&geom,&geomout);

  ldouble (*gg)[5],(*GG)[5];
  gg=geomout.gg;
  GG=geomout.GG;

  //four-vectors etc
  ldouble rho=pp[0];
  ldouble uint=pp[1];
  ldouble Tgas=(GAMMA-1.)*uint*MU_GAS*M_PROTON/K_BOLTZ/rho; 
  ldouble vcon[4],vcov[4],ucon[4],ucov[4];
  ldouble bcon[4]={0.,0.,0.,0.},bcov[4]={0.,0.,0.,0.},bsq=0.;
  vcon[1]=pp[2];
  vcon[2]=pp[3];
  vcon[3]=pp[4];
  vcon[0]=0.;
  ldouble S=pp[5];
  conv_vels_both(vcon,ucon,ucov,VELPRIM,VEL4,gg,GG); 
#ifdef MAGNFIELD
  calc_bcon_4vel(pp,ucon,ucov,bcon);
  indices_21(bcon,bcov,gg); 
  bsq = dot(bcon,bcov);
#endif
  
  //hydro stress-energy
  ldouble Tij[4][4];
  calc_Tij(pp,&geomout,Tij);
  indices_2221(Tij,Tij,gg);

  //#ifdef BHDISK_PROBLEMTYPE
  //avg already in OUTCOORDS
  avg[AVGBSQ]=bsq;
  avg[AVGTGAS]=Tgas;
  for(iv=0;iv<4;iv++)
    avg[AVGUCON(iv)]=ucon[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGUCOV(iv)]=ucov[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGBCON(iv)]=bcon[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGBCOV(iv)]=bcov[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGRHOUCON(iv)]=rho*ucon[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGRHOUCOV(iv)]=rho*ucov[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGUUUCON(iv)]=uint*ucon[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGUUCOV(iv)]=uint*ucov[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGBSQUCON(iv)]=bsq*ucon[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGBSQUCOV(iv)]=bsq*ucov[iv];
  for(iv=0;iv<4;iv++)
    for(iv2=0;iv2<4;iv2++)
      avg[AVGRHOUCONUCOV(iv,iv2)]=rho*ucon[iv]*ucov[iv2];
  for(iv=0;iv<4;iv++)
    for(iv2=0;iv2<4;iv2++)
      avg[AVGUUUCONUCOV(iv,iv2)]=uint*ucon[iv]*ucov[iv2];
  for(iv=0;iv<4;iv++)
    for(iv2=0;iv2<4;iv2++)
      avg[AVGBSQUCONUCOV(iv,iv2)]=bsq*ucon[iv]*ucov[iv2];
  for(iv=0;iv<4;iv++)
    for(iv2=0;iv2<4;iv2++)
      avg[AVGBCONBCOV(iv,iv2)]=bcon[iv]*bcov[iv2];
  for(iv=0;iv<4;iv++)
    avg[AVGWUCON(iv)]=(rho+uint+bsq/2)*ucon[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGWUCON(iv)]=(rho+uint+bsq/2)*ucon[iv];

  //fluxes at faces, including the diffusive part
  //conserved fluxes at left face in MYCOORDS
  for(iv=0;iv<NV;iv++)
    {
      avg[AVGFLUXXL(iv)]=get_ub(flbx,iv,ix,iy,iz,0);
      avg[AVGFLUXYL(iv)]=get_ub(flby,iv,ix,iy,iz,1);
      avg[AVGFLUXZL(iv)]=get_ub(flbz,iv,ix,iy,iz,2);
    }

  //converting rest-mass flux to BLCOORDS 
  ldouble vector[4];
  //primitives and conserved at left faces - used to fill missing time-component
  ldouble uface[NV],pface[NV],fd_uLl[NV],fd_uRl[NV];
  int i;
  for(i=0;i<NV;i++)
    {
      fd_uLl[i]=get_ub(pbLx,i,ix,iy,iz,0);
      fd_uRl[i]=get_ub(pbRx,i,ix,iy,iz,0);
      pface[i]=.5*(fd_uLl[i]+fd_uRl[i]);
    }
  p2u(pface,uface,&geoml);

  //rest-mass flux in radius
  vector[0]=uface[RHO]; //rho ut gdet
  vector[1]=get_ub(flbx,iv,ix,iy,iz,0); //rho ur gdet
  vector[2]=0.; //unimportant within Kerr-Shield
  vector[3]=0.;
  trans2_coco(geoml.xxvec,vector,vector,MYCOORDS, OUTCOORDS); //this does not work?
  avg[AVGRHOURDIFF]=vector[1];
   
#ifdef RADIATION
  ldouble Rtt,Ehat,ugas[4];
  //calc_ff_Rtt(pp,&Rtt,ugas,&geomout); //this very slow - why?

  /*
  ldouble Rijp[4][4];
  calc_Rij_M1(pp,&geomout,Rijp);
 
  indices_2221(Rijp,Rijp,geomout.gg);
  Rtt=0.;
  int i1,i2;
  for(i1=0;i1<4;i1++)
    for(i2=0;i2<4;i2++)
      Rtt+=-Rijp[i1][i2]*ucon[i2]*ucov[i1];
  
  Ehat=-Rtt;       
  */

  ldouble Rij[4][4];
  calc_Rij(pp,&geomout,Rij);
  indices_2221(Rij,Rij,geomout.gg);

  //Ehat calculation from Rij, gas velocity still in ucon, ucov
  Rtt=0.;
  int i1,i2;
  for(i1=0;i1<4;i1++)
    for(i2=0;i2<4;i2++)
    Rtt+=-Rij[i1][i2]*ucon[i2]*ucov[i1];

  Ehat=-Rtt;       

  vcon[1]=pp[FX];
  vcon[2]=pp[FY];
  vcon[3]=pp[FZ];
  vcon[0]=0.;  
  conv_vels_both(vcon,ucon,ucov,VELPRIM,VEL4,gg,GG); 

  /*
  //fluid frame Ghat
  ldouble Gi[4],Gic[4];
  calc_Gi(pp,&geomout,Gi,0); 

  ldouble kappaes=calc_kappaes(pp,&geomout);
  ldouble uffcon[4]={1.,0.,0.,0.};
  calc_Compt_Gi(pp,&geomout,Gic,Ehat,Tgas,kappaes,uffcon);
  */

  //four fource
  ldouble Gi[4],Gic[4],Giff[4],Gicff[4];
  //calc_Gi(pp,&geomout,Gi,1); 
  //boost2_lab2ff(Gi,Giff,pp,geomout.gg,geomout.GG);
  calc_Gi(pp,&geomout,Giff,0); 

  
#if defined(COMPTONIZATION) || defined(NCOMPTONIZATION)
  //uwaga! boost sprawia, ze znaki fluid frame Compt i abs part rozne! spojrzec dlaczego!

  //test - directly in ff
  ucon[1]=ucon[2]=ucon[3]=0.;  ucon[0]=1.;
  ldouble kappaes=calc_kappaes(pp,&geomout);
  calc_Compt_Gi(pp,&geomout,Gicff,Ehat,Tgas,kappaes,ucon);
  
  //in lab frame + boost:
  //   calc_Compt_Gi(pp,&geomout,Gic,Ehat,Tgas,kappaes,ucon);
  //boost2_lab2ff(Gic,Gicff,pp,geomout.gg,geomout.GG);

  


#endif 

  //radiation temperature
  ldouble Thatrad;
  #ifdef NCOMPTONIZATION //number of photons conserved
  Thatrad = calc_ncompt_Thatrad(pp,&geomout,Ehat);
  #else //thermal comptonization
  Thatrad = calc_LTE_TfromE(Ehat);
  #endif
  
  for(iv=0;iv<4;iv++)
    avg[AVGURFCON(iv)]=ucon[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGURFCOV(iv)]=ucov[iv];
   
  avg[AVGEHAT]=Ehat;
  avg[AVGTRAD]=Thatrad;

  for(iv=0;iv<4;iv++)
    for(iv2=0;iv2<4;iv2++)
      avg[AVGRIJ(iv,iv2)]=Rij[iv][iv2];

  for(iv=0;iv<4;iv++)
    avg[AVGEHATUCON(iv)]=Ehat*ucon[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGEHATUCOV(iv)]=Ehat*ucov[iv];

  for(iv=0;iv<4;iv++)
    avg[AVGGHAT(iv)]=Giff[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGGHATCOMPT(iv)]=Gicff[iv];

#endif

  //#endif

  return 0.;
}

//this computes utp1=1+u_t , which for nonrelativistic cases is ~0 .  If computed as 1+u_t, then residual will be large error if small residual.
ldouble
calc_utp1(ldouble *vcon, ldouble *ucon, void *ggg)
{
   struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  gdet=geom->gdet;
  GG=geom->GG;
  gdetu=gdet;

  ldouble utp1;
  if(VELPRIM==VELR) //based on vELR
    {
      int i,j;
      ldouble qsq=0.;
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  qsq+=vcon[i]*vcon[j]*gg[i][j];
      ldouble gamma2=(1.+qsq);
      ldouble alphasq=(-1./GG[0][0]);
      ldouble alpgam=sqrt(alphasq*gamma2);
      //\beta^i \beta_i / \alpha^2 = g^{ti} g_{ti}
      ldouble betasqoalphasq=gg[0][1]*GG[0][1] + gg[0][2]*GG[0][2] + gg[0][3]*GG[0][3]; 

      ldouble ud0tilde = 0.0;
      SLOOPA(j) ud0tilde += vcon[j]*gg[0][j]; // \tilde{u}_t = \tilde{u}^i g_{ti} since \tilde{u}^t=0
      utp1= ud0tilde + (geom->gttpert - alphasq*(betasqoalphasq + qsq))/(1.0+alpgam);
    }
  else //based on ucon[]
    {
      int i,j,k;
      // 3-velocity in coordinate basis
      ldouble vconp[4];
      SLOOPA(j) vconp[j]=ucon[j]/ucon[0];

      ldouble plus1gv00=geom->gttpert;
      ldouble vsq=geom->gttpert;
      SLOOPA(j) vsq+=2.0*geom->gg[0][j]*vconp[j];
      SLOOP(j,k) vsq+=geom->gg[j][k]*vconp[j]*vconp[k];

      ldouble gvtt=geom->gg[0][0];
      
      ldouble alpha=0.0;
      SLOOPA(j) alpha+=geom->gg[j][0]*ucon[j];

      ldouble uu0 = ucon[0];

      utp1 = alpha + ((1.0-gvtt)*plus1gv00 - uu0*uu0*vsq*gvtt*gvtt)/(1.0-gvtt*uu0);
    }
  return utp1;
}
