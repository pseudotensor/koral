//KORAL - physics.c
//some problem independent physics

#include "ko.h"
//*************************************************
//calculates left and right wave speeds at cell center
//*************************************************
int
calc_wavespeeds_lr_pure(ldouble *pp,void *ggg,ldouble *aaa)
{
  struct geometry *geom
    = (struct geometry *) ggg;
  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;

  gg=geom->gg;
  GG=geom->GG;

  int iv;
  
  ldouble axhdl,axhdr,ayhdl,ayhdr,azhdl,azhdr;
  ldouble axl,axr,ayl,ayr,azl,azr;
  axl=axr=ayl=ayr=azl=azr=1.;
  
  ldouble ucon[4],ucov[4],cst1,cst2,cst3,cst4;
  

  //**********************************************************************
  //***** hydro: speed of sound ******************************************
  //**********************************************************************
	      
  ldouble rho=pp[RHO];
  ldouble uu=pp[UU];
 
  ldouble pre=(GAMMA-1.)*uu;
  ldouble cs2=GAMMA*pre/(rho+uu+pre);

  //test
  //cs2*=4.;
  
  if(cs2<0.) cs2=0.;

  //**********************************************************************
  //***** other stuff ****************************************************
  //**********************************************************************

  for(iv=1;iv<4;iv++)
    ucon[iv]=pp[1+iv];
  ucon[0]=0.;
  conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);
  indices_21(ucon,ucov,gg);

  //**********************************************************************
  //algorithm from HARM to transform the fluid frame wavespeed into lab frame
  //TODO IMMEDIATELY: generalize! ???
  //**********************************************************************

  ldouble Acov[4],Acon[4],Bcov[4],Bcon[4],Asq,Bsq,Au,Bu,AB,Au2,Bu2,AuBu,A,B,discr,wspeed2;

  //**********************************************************************
  //**********************************************************************
  //x
  Acov[0]=0.;
  Acov[1]=1.;
  Acov[2]=0.;
  Acov[3]=0.;
  indices_12(Acov,Acon,GG);
   
  Bcov[0]=1.;
  Bcov[1]=0.;
  Bcov[2]=0.;
  Bcov[3]=0.;
  indices_12(Bcov,Bcon,GG);

  Asq = dot(Acon,Acov);
  Bsq = dot(Bcon,Bcov);
  Au = dot(Acov, ucon);
  Bu = dot(Bcov, ucon);
  AB = dot(Acon, Bcov);
  Au2 = Au * Au;
  Bu2 = Bu * Bu;
  AuBu = Au * Bu;

  //hydro
  wspeed2=cs2;
  B = 2. * (AuBu * (1.0-wspeed2)  - AB*wspeed2);
  A = Bu2 * (1.0 - wspeed2) - Bsq * wspeed2;
  discr = 4.0 * wspeed2 * ((AB * AB - Asq * Bsq) * wspeed2 + (2.0 * AB * Au * Bu - Asq * Bu2 - Bsq * Au2) * (wspeed2 - 1.0));

  if(discr<0.) {printf("x discr in wavespeeds lt 0\n"); discr=0.;}
  discr = sqrt(discr);
  cst1 = -(-B + discr) / (2. * A);
  cst2 = -(-B - discr) / (2. * A);  
  if(cst2>cst1)
    {
      axhdl=cst1;  axhdr=cst2;
    }
  else
    {
      axhdr=cst1;  axhdl=cst2;
    }

  //**********************************************************************
  //**********************************************************************
  //y
  Acov[0]=0.;
  Acov[1]=0.;
  Acov[2]=1.;
  Acov[3]=0.;
  indices_12(Acov,Acon,GG);
  
  Asq = dot(Acon,Acov);
  Bsq = dot(Bcon,Bcov);
  Au = dot(Acov, ucon);
  Bu = dot(Bcov, ucon);
  AB = dot(Acon, Bcov);
  Au2 = Au * Au;
  Bu2 = Bu * Bu;
  AuBu = Au * Bu;

  //hydro
  wspeed2=cs2;
  B = 2. * (AuBu * (1.0-wspeed2)  - AB*wspeed2);
  A = Bu2 * (1.0 - wspeed2) - Bsq * wspeed2;
  discr = 4.0 * wspeed2 * ((AB * AB - Asq * Bsq) * wspeed2 + (2.0 * AB * Au * Bu - Asq * Bu2 - Bsq * Au2) * (wspeed2 - 1.0));
  if(discr<0.) {printf("x discr in wavespeeds lt 0\n"); discr=0.;}
  discr = sqrt(discr);
  cst1 = -(-B + discr) / (2. * A);
  cst2 = -(-B - discr) / (2. * A);  
  if(cst2>cst1)
    {
      ayhdl=cst1;  ayhdr=cst2;
    }
  else
    {
      ayhdr=cst1;  ayhdl=cst2;
    }
 
       
  //**********************************************************************
  //**********************************************************************
  //z
  Acov[0]=0.;
  Acov[1]=0.;
  Acov[2]=0.;
  Acov[3]=1.;
  indices_12(Acov,Acon,GG);
   
  Asq = dot(Acon,Acov);
  Bsq = dot(Bcon,Bcov);
  Au = dot(Acov, ucon);
  Bu = dot(Bcov, ucon);
  AB = dot(Acon, Bcov);
  Au2 = Au * Au;
  Bu2 = Bu * Bu;
  AuBu = Au * Bu;

  //hydro
  wspeed2=cs2;
  B = 2. * (AuBu * (1.0-wspeed2)  - AB*wspeed2);
  A = Bu2 * (1.0 - wspeed2) - Bsq * wspeed2;
  discr = 4.0 * wspeed2 * ((AB * AB - Asq * Bsq) * wspeed2 + (2.0 * AB * Au * Bu - Asq * Bu2 - Bsq * Au2) * (wspeed2 - 1.0));
  if(discr<0.) {printf("y discr in wavespeeds lt 0\n"); discr=0.;}
  discr = sqrt(discr);
  cst1 = -(-B + discr) / (2. * A);
  cst2 = -(-B - discr) / (2. * A);  
  if(cst2>cst1)
    {
      azhdl=cst1;  azhdr=cst2;
    }
  else
    {
      azhdr=cst1;  azhdl=cst2;
    }

#ifdef RADIATION
  //**********************************************************************
  //***** radiation: characteristic wave speed ***************************
  //**********************************************************************

  ldouble aval[6];
  int verbose=0;

  //physical size of the cell
  ldouble dx[3];
  ldouble xx[4]={0.,geom->xx,geom->yy,geom->zz};
  
  //ix,iy,iz has the indices of the face, so the depth taken from left/right
  dx[0]=my_max(get_size_x(geom->ix,0)*sqrt(gg[1][1]),get_size_x(geom->ix+1,0)*sqrt(gg[1][1]));
  dx[1]=my_max(get_size_x(geom->iy,1)*sqrt(gg[2][2]),get_size_x(geom->iy+1,1)*sqrt(gg[2][2]));
  dx[2]=my_max(get_size_x(geom->iz,2)*sqrt(gg[3][3]),get_size_x(geom->iz+1,2)*sqrt(gg[3][3]));
  ldouble tautot[3];
  calc_tautot(pp,xx,dx,tautot);

#ifndef MULTIRADFLUID
  calc_rad_wavespeeds(pp,geom,tautot,aval,verbose);

#ifdef FULLRADWAVESPEEDS
  aval[0]=aval[1]=1./sqrt(gg[1][1]);
  aval[2]=aval[3]=1./sqrt(gg[2][2]);
  aval[4]=aval[5]=1./sqrt(gg[3][3]);
#endif

#else //multifluid
  calc_rad_wavespeeds_mf_total(pp,gg,GG,tautot,aval);
#endif

  axl=aval[0];
  axr=aval[1];
  ayl=aval[2];
  ayr=aval[3];
  azl=aval[4];
  azr=aval[5];

#endif

  //**********************************************************************
  //**********************************************************************
    
  //combining rad and hydro velocities when SIMPLEVISCOSITY with total pressure
#ifdef ENFORCERADWAVESPEEDS
#ifdef SIMPLEHDVISCOSITY

 //based on local opacities
  ldouble xxvecBL[4];
  coco_N(xx,xxvecBL,MYCOORDS,BLCOORDS);
  //damping in radius
  ldouble fdampr=step_function(xxvecBL[1]-RMINVISC,RMINVISC/5.);

  dx[0]=1.;dx[1]=1.;dx[2]=1.;//opt depth per Rg
  calc_tautot(pp,xx,dx,tautot);
  ldouble Rij[4][4];
  calc_Rij(pp,ggg,Rij);
  // skip boosts what reasonable when gas v<<1
  // boost22_lab2ff(Rij,Rij,pp,gg,GG);
  // trans22_cc2on(Rij,Rij,tup);
  tautot[0]+=1.e-50;
  ldouble PARAM=TAUSUPPRESSPARAM;
  //to damp less than add pressure
  PARAM/=10.;
  ldouble fdamptau=exp(-PARAM/tautot[0]/tautot[0]);
  ldouble fdamp=fdamptau*fdampr;

  axhdl -= axl*fdamp;
  axhdr += axr*fdamp;
  ayhdl -= ayl*fdamp;
  ayhdr += ayr*fdamp;
  azhdl -= azl*fdamp;
  azhdr += azr*fdamp;
  
#endif
#endif

  //zeroing 'co-going' velocities
  if(axhdl>0.) axhdl=0.;
  if(axhdr<0.) axhdr=0.;
  if(ayhdl>0.) ayhdl=0.;
  if(ayhdr<0.) ayhdr=0.;
  if(azhdl>0.) azhdl=0.;
  if(azhdr<0.) azhdr=0.;
  if(axl>0.) axl=0.;
  if(axr<0.) axr=0.;
  if(ayl>0.) ayl=0.;
  if(ayr<0.) ayr=0.;
  if(azl>0.) azl=0.;
  if(azr<0.) azr=0.;

  //saving and passing up
  aaa[0]=axhdl;
  aaa[1]=axhdr;
  aaa[2]=ayhdl;
  aaa[3]=ayhdr;
  aaa[4]=azhdl;
  aaa[5]=azhdr;
  aaa[6]=axl;
  aaa[7]=axr;
  aaa[8]=ayl;
  aaa[9]=ayr;
  aaa[10]=azl;
  aaa[11]=azr;

  return 0;
}

//*************************************************
//calculates left and right wave speeds at cell center
//*************************************************
int
calc_wavespeeds_lr(int ix, int iy, int iz,ldouble *aaa)
{
  ldouble (*gg)[5],(*GG)[5];

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  
  //temporary using local arrays
  gg=geom.gg;
  GG=geom.GG;

  ldouble pp[NV];
  int iv;

  //picking up primitives 
  for(iv=0;iv<NV;iv++)
    pp[iv]=get_u(p,iv,ix,iy,iz);

  calc_wavespeeds_lr_pure(pp,&geom,aaa);

  /*
  conv_velsinprims(pp,VELPRIM,VEL3,gg,GG);
  if(fabs(pp[2])>.2 || (ix==NX-1 && iy==NY-1))
    {
      printf("%d %d %d\n %e %e %e\n",ix,iy,iz,pp[2],pp[3],pp[4]);
      print_Nvector(aaa,6);
      getchar();
    }
  */
  return 0;
}

//***************************************
//returns othersource terms for all conserved quantities
//***************************************
int f_other_source_term_arb(ldouble *pp,void *ggg,ldouble *ss)
{
  int i;

  struct geometry *geom
    = (struct geometry *) ggg;
  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;

  int ix,iy,iz;
  ix=geom->ix;
  iy=geom->iy;
  iz=geom->iz;
 
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet;gdetu=gdet;
#if (GDETIN==0) //no metric determinant inside derivatives
  gdetu=1.;
#endif

  ldouble rho=pp[RHO];
  ldouble u=pp[1];  
  ldouble E=pp[6];  
  ldouble pr=(GAMMA-1.)*(u);
  ldouble T=pr*MU_GAS*M_PROTON/K_BOLTZ/rho;
  
  ldouble Gi[4];
  calc_Gi(pp,ggg,Gi);

  ldouble ucon[4]={0.,pp[2],pp[3],pp[4]},ucov[4];
  conv_vels(ucon,ucon,VELPRIM,VEL4,geom->gg,geom->GG);
   indices_21(ucon,ucov,geom->gg);

  ldouble entropy_source= -1./T*dot(ucov,Gi);

  for(i=0;i<NV;i++)
    ss[i]=0.;
  
  //  ss[5]=entropy_source;


  return 0;
}


//***************************************
//returns geometrical source terms for all conserved quantities
//***************************************
int f_other_source_term(int ix, int iy, int iz,ldouble *ss)
{
  int i;
  ldouble pp[NV];  

  for(i=0;i<NV;i++)
    pp[i]=get_u(p,i,ix,iy,iz);  
  
  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);

  f_other_source_term_arb(pp,&geom,ss);

  return 0;
}

//***************************************
//returns geometrical source terms for all conserved quantities
//***************************************
int f_metric_source_term_arb(ldouble *pp,void *ggg,ldouble *ss)
{
  int i;

  struct geometry *geom
    = (struct geometry *) ggg;
  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;

  int ix,iy,iz;
  ix=geom->ix;
  iy=geom->iy;
  iz=geom->iz;
 
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet;gdetu=gdet;
#if (GDETIN==0) //no metric determinant inside derivatives
  gdetu=1.;
#endif

  ldouble dlgdet[3];
  dlgdet[0]=gg[0][4]; //D[gdet,x1]/gdet
  dlgdet[1]=gg[1][4]; //D[gdet,x2]/gdet
  dlgdet[2]=gg[2][4]; //D[gdet,x3]/gdet
  
  ldouble ut;
  ldouble T[4][4];

  //calculating stress energy tensor components
  calc_Tij(pp,geom,T);
  indices_2221(T,T,gg);

  int ii, jj;
  for(ii=0;ii<4;ii++)
    for(jj=0;jj<4;jj++)
      {
	if(isnan(T[ii][jj])) 
	  {
	    printf("%d %d %e\n",ii,jj,T[ii][jj]);
	    my_err("nan in metric_source_terms\n");
	  }
      }
 
  ldouble rho=pp[RHO];
  ldouble u=pp[UU];
  ldouble vcon[4],ucon[4];
  vcon[1]=pp[2];
  vcon[2]=pp[3];
  vcon[3]=pp[4];
  ldouble S=pp[5];

  //converting to 4-velocity
  conv_vels(vcon,ucon,VELPRIM,VEL4,gg,GG);
  
  int k,l,iv;
  for(iv=0;iv<NV;iv++)
    ss[iv]=0.;

  /***************************************************/
#ifdef RADIATION
  /***************************************************/

#ifndef MULTIRADFLUID
  ldouble Rij[4][4];
  calc_Rij(pp,geom,Rij); //R^ij
  indices_2221(Rij,Rij,gg); //R^i_j

  //terms with Christoffels
  for(k=0;k<4;k++)
    for(l=0;l<4;l++)
      {
	ss[1]+=gdetu*T[k][l]*get_gKr(l,0,k,ix,iy,iz);
	ss[2]+=gdetu*T[k][l]*get_gKr(l,1,k,ix,iy,iz);
	ss[3]+=gdetu*T[k][l]*get_gKr(l,2,k,ix,iy,iz);
	ss[4]+=gdetu*T[k][l]*get_gKr(l,3,k,ix,iy,iz);
	ss[6]+=gdetu*Rij[k][l]*get_gKr(l,0,k,ix,iy,iz);
	ss[7]+=gdetu*Rij[k][l]*get_gKr(l,1,k,ix,iy,iz);
	ss[8]+=gdetu*Rij[k][l]*get_gKr(l,2,k,ix,iy,iz);
	ss[9]+=gdetu*Rij[k][l]*get_gKr(l,3,k,ix,iy,iz);
      }

  //terms with dloggdet
#if (GDETIN==0)
  for(l=1;l<4;l++)
    {
      ss[0]+=-dlgdet[l-1]*rho*ucon[l];
      ss[1]+=-dlgdet[l-1]*(T[l][0]+rho*ucon[l]);
      ss[2]+=-dlgdet[l-1]*(T[l][1]);
      ss[3]+=-dlgdet[l-1]*(T[l][2]);
      ss[4]+=-dlgdet[l-1]*(T[l][3]);
      ss[5]+=-dlgdet[l-1]*S*ucon[l];
      ss[6]+=-dlgdet[l-1]*(Rij[l][0]);
      ss[7]+=-dlgdet[l-1]*(Rij[l][1]);
      ss[8]+=-dlgdet[l-1]*(Rij[l][2]);
      ss[9]+=-dlgdet[l-1]*(Rij[l][3]);
    }
#endif

#else
  int irf;
  ldouble Rij[NRF][4][4];
  calc_Rij_mf(pp,gg,GG,Rij); //R^ij
  for(ii=0;ii<NRF;ii++)
    indices_2221(Rij[ii],Rij[ii],gg); //R^i_jl

  //terms with Christoffels
  //hydro first
  for(k=0;k<4;k++)
    for(l=0;l<4;l++)
      {
	ss[1]+=gdetu*T[k][l]*get_gKr(l,0,k,ix,iy,iz);
	ss[2]+=gdetu*T[k][l]*get_gKr(l,1,k,ix,iy,iz);
	ss[3]+=gdetu*T[k][l]*get_gKr(l,2,k,ix,iy,iz);
	ss[4]+=gdetu*T[k][l]*get_gKr(l,3,k,ix,iy,iz);	 
      }
  //now radiation
  for(irf=0;irf<NRF;irf++)
    for(k=0;k<4;k++)
      for(l=0;l<4;l++)
	{
	  ss[EE(irf)]+=gdetu*Rij[irf][k][l]*get_gKr(l,0,k,ix,iy,iz);
	  ss[FX(irf)]+=gdetu*Rij[irf][k][l]*get_gKr(l,1,k,ix,iy,iz);
	  ss[FY(irf)]+=gdetu*Rij[irf][k][l]*get_gKr(l,2,k,ix,iy,iz);
	  ss[FZ(irf)]+=gdetu*Rij[irf][k][l]*get_gKr(l,3,k,ix,iy,iz);
	}

  //terms with dloggdet
#if (GDETIN==0)
  //hydro first
  for(l=1;l<4;l++)
    {
      ss[0]+=-dlgdet[l-1]*rho*ucon[l];
      ss[1]+=-dlgdet[l-1]*(T[l][0]+rho*ucon[l]);
      ss[2]+=-dlgdet[l-1]*(T[l][1]);
      ss[3]+=-dlgdet[l-1]*(T[l][2]);
      ss[4]+=-dlgdet[l-1]*(T[l][3]);
      ss[5]+=-dlgdet[l-1]*S*ucon[l];
    }

  //rad now
  for(irf=0;irf<NRF;irf++)
    for(l=1;l<4;l++)
      {
	ss[EE(irf)]+=-dlgdet[l-1]*(Rij[irf][l][0]);
	ss[FX(irf)]+=-dlgdet[l-1]*(Rij[irf][l][1]);
	ss[FY(irf)]+=-dlgdet[l-1]*(Rij[irf][l][2]);
	ss[FZ(irf)]+=-dlgdet[l-1]*(Rij[irf][l][3]);
      }
#endif

#endif

  /***************************************************/
#else //pure hydro
  /***************************************************/

  //terms with Christoffels
  for(k=0;k<4;k++)
    for(l=0;l<4;l++)
      {
	ss[1]+=gdetu*T[k][l]*get_gKr(l,0,k,ix,iy,iz);
	ss[2]+=gdetu*T[k][l]*get_gKr(l,1,k,ix,iy,iz);
	ss[3]+=gdetu*T[k][l]*get_gKr(l,2,k,ix,iy,iz);
	ss[4]+=gdetu*T[k][l]*get_gKr(l,3,k,ix,iy,iz);
      }

  //terms with dloggdet  
#if (GDETIN==0)
  for(l=1;l<4;l++)
    {
      ss[0]+=-dlgdet[l-1]*rho*ucon[l];
      ss[1]+=-dlgdet[l-1]*(T[l][0]+rho*ucon[l]);
      ss[2]+=-dlgdet[l-1]*(T[l][1]);
      ss[3]+=-dlgdet[l-1]*(T[l][2]);
      ss[4]+=-dlgdet[l-1]*(T[l][3]);
      ss[5]+=-dlgdet[l-1]*S*ucon[l];
    }   
#endif

  /***************************************************/
#endif
  /***************************************************/

  return 0;
}

//***************************************
//returns geometrical source terms for all conserved quantities
//***************************************
int f_metric_source_term(int ix, int iy, int iz,ldouble *ss)
{
  int i;
  ldouble pp[NV];  

  for(i=0;i<NV;i++)
    pp[i]=get_u(p,i,ix,iy,iz);  
  
  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);

  f_metric_source_term_arb(pp,&geom,ss);

  return 0;
}

//***************************************
// calculates fluxes at faces
//***************************************
ldouble f_flux_prime( ldouble *pp, int idim, int ix, int iy, int iz,ldouble *ff)
{  
  int iv;
  for(iv=0;iv<NV;iv++)
    ff[iv]=0.;

  //picking up metric from a cell face  
  struct geometry geom;
  fill_geometry_face(ix,iy,iz,idim,&geom);
 
  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom.gg;
  GG=geom.GG;
  gdet=geom.gdet;gdetu=gdet;
#if (GDETIN==0) //no metric determinant inside derivatives
  gdetu=1.;
#endif

  //calculating Tij
  ldouble T[4][4];

  calc_Tij(pp,&geom,T);
  indices_2221(T,T,gg);

  //primitives
  ldouble rho=pp[RHO];
  ldouble u=pp[UU];
  ldouble vcon[4],ucon[4];
  vcon[1]=pp[2];
  vcon[2]=pp[3];
  vcon[3]=pp[4];
  ldouble S=pp[5];

  //converting to 4-velocity
  //TODO: introduce structure of state
  conv_vels(vcon,ucon,VELPRIM,VEL4,gg,GG);

  //4-velocity
  ldouble u1=ucon[1];
  ldouble u2=ucon[2];
  ldouble u3=ucon[3];
  ldouble gam=GAMMA;

  int ii, jj, irf;
  for(ii=0;ii<4;ii++)
    for(jj=0;jj<4;jj++)
      {
	if(isnan(T[ii][jj])) 
	  {
	    printf("nan tmunu: %d %d %e\n",ii,jj,T[ii][jj]);
	    print_Nvector(pp,NV);
	    my_err("nan in flux_prime\n");
	  }
      }
 
#ifdef RADIATION

#ifndef MULTIRADFLUID
  ldouble Rij[4][4];
  calc_Rij(pp,&geom,Rij); //R^ij
  indices_2221(Rij,Rij,gg); //R^i_j
#else
  ldouble Rij[NRF][4][4];
  calc_Rij_mf(pp,gg,GG,Rij); //R^ij

  for(ii=0;ii<NRF;ii++)
    {
      indices_2221(Rij[ii],Rij[ii],gg); //R^i_j
    }
 #endif

  //to move gdet in/out derivative:
  //here, up in metric source terms, in u2p and p2u, as well as in finite.c with del4[], and below in calc_entropy

   if(idim==0) //x
    {
      ff[0]= gdetu*rho*u1;

      ff[1]= gdetu*(T[1][0]+rho*u1);

      ff[2]= gdetu*(T[1][1]);

      ff[3]= gdetu*(T[1][2]);
 
      ff[4]= gdetu*(T[1][3]);

      ff[5]= gdetu*S*u1;

#ifndef MULTIRADFLUID
      ff[6]= gdetu*Rij[1][0];
      
      ff[7]= gdetu*Rij[1][1];
      
      ff[8]= gdetu*Rij[1][2];
      
      ff[9]= gdetu*Rij[1][3];
#else
      for(irf=0;irf<NRF;irf++)
	{
	  ff[EE(irf)]=gdetu*Rij[irf][1][0];
	  ff[FX(irf)]=gdetu*Rij[irf][1][1];
	  ff[FY(irf)]=gdetu*Rij[irf][1][2];
	  ff[FZ(irf)]=gdetu*Rij[irf][1][3];
	}
#endif
    }  
  if(idim==1) //y
    {
      ff[0]= gdetu*rho*u2;

      ff[1]= gdetu*(T[2][0]+rho*u2);

      ff[2]= gdetu*(T[2][1]);

      ff[3]= gdetu*(T[2][2]);

      ff[4]= gdetu*(T[2][3]);

      ff[5]= gdetu*S*u2;
 
#ifndef MULTIRADFLUID
      ff[6]= gdetu*Rij[2][0];
      
      ff[7]= gdetu*Rij[2][1];
      
      ff[8]= gdetu*Rij[2][2];
      
      ff[9]= gdetu*Rij[2][3];
#else
      for(irf=0;irf<NRF;irf++)
	{
	  ff[EE(irf)]=gdetu*Rij[irf][2][0];
	  ff[FX(irf)]=gdetu*Rij[irf][2][1];
	  ff[FY(irf)]=gdetu*Rij[irf][2][2];
	  ff[FZ(irf)]=gdetu*Rij[irf][2][3];
	}
#endif
    }  
  if(idim==2) //z
    {
      ff[0]= gdetu*rho*u3;

      ff[1]= gdetu*(T[3][0]+rho*u3);
 
      ff[2]= gdetu*(T[3][1]);

      ff[3]= gdetu*(T[3][2]);
 
      ff[4]= gdetu*(T[3][3]);

      ff[5]= gdetu*S*u3;
 
#ifndef MULTIRADFLUID
      ff[6]= gdetu*Rij[3][0];
      
      ff[7]= gdetu*Rij[3][1];
      
      ff[8]= gdetu*Rij[3][2];
       
      ff[9]= gdetu*Rij[3][3];
#else
      for(irf=0;irf<NRF;irf++)
	{
	  ff[EE(irf)]=gdetu*Rij[irf][3][0];
	  ff[FX(irf)]=gdetu*Rij[irf][3][1];
	  ff[FY(irf)]=gdetu*Rij[irf][3][2];
	  ff[FZ(irf)]=gdetu*Rij[irf][3][3];
	}
#endif
    } 

#else //pure hydro

  if(idim==0) //x
    {
      ff[0]= gdetu*rho*u1;

      ff[1]= gdetu*(T[1][0]+rho*u1);

      ff[2]= gdetu*(T[1][1]);

      ff[3]= gdetu*(T[1][2]);
 
      ff[4]= gdetu*(T[1][3]);

      ff[5]= gdetu*S*u1;
    }  
  if(idim==1) //y
    {
      ff[0]= gdetu*rho*u2;

      ff[1]= gdetu*(T[2][0]+rho*u2);

      ff[2]= gdetu*(T[2][1]);

      ff[3]= gdetu*(T[2][2]);

      ff[4]= gdetu*(T[2][3]);

      ff[5]= gdetu*S*u2;
    }  
  if(idim==2) //z
    {
      ff[0]= gdetu*rho*u3;

      ff[1]= gdetu*(T[3][0]+rho*u3);
 
      ff[2]= gdetu*(T[3][1]);

      ff[3]= gdetu*(T[3][2]);
 
      ff[4]= gdetu*(T[3][3]);

      ff[5]= gdetu*S*u3;
    } 

#endif
  /*

  if(isinf(ff[5]))
    {
      printf("%e %e %e \n",S,u1,u2); getchar();
    }
  */
  return 0.;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates energy-momentum tensor components basing on vector of primitivies p and given metric g
//returns T^munu
int
calc_Tij(ldouble *pp, void* ggg, ldouble T[][4])
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  int iv,i,j;
  ldouble rho=pp[RHO];
  ldouble uu=pp[UU];
  ldouble ucon[4],ucov[4];
  
  //converts to 4-velocity
  for(iv=1;iv<4;iv++)
    ucon[iv]=pp[1+iv];
  ucon[0]=0.;
  conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);
  //indices_21(ucon,ucov,gg);

  ldouble w=rho+GAMMA*uu;
  ldouble p=(GAMMA-1.)*uu;

  /*
  //this old mixed formulation worked even at the polar axis but provides wrong Tij there
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      T[i][j]=w*ucon[i]*ucov[j]+p*delta(i,j);
  */

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      T[i][j]=w*ucon[i]*ucon[j]+p*GG[i][j];

#if (HDVISCOSITY!=NOVISCOSITY)
  ldouble Tvisc[4][4];
  calc_visc_Tij(pp,ggg,Tvisc);

  //test
  /*
  //limitting maximal allowed viscous change
  ldouble maxviscchange=MAXVISCCHANGE,change=-1.,viscdamp;

  ldouble maxspatial=-1.;
  for(i=1;i<4;i++)
    for(j=1;j<4;j++)
      {
	if(fabs(T[i][j])>maxspatial) maxspatial=fabs(T[i][j]);
      }

  for(i=1;i<4;i++)
    for(j=1;j<4;j++)
      {
	ldouble d1=fabs(Tvisc[i][j])/maxspatial;
	if(d1>change) 
	  change=d1;
      }

  if(change>maxviscchange) 
    viscdamp=maxviscchange/change;
  else
    viscdamp=1.;
  */

  ldouble viscdamp=1.0;
	
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {	
	T[i][j]+=viscdamp*Tvisc[i][j];
      }

#endif  

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates viscouss energy-momentum tensor components
//returns T^munu
int
calc_visc_Tij(ldouble *pp, void* ggg, ldouble T[][4])
{
  int i,j;
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],(*tlo)[4],(*tup)[4];
  gg=geom->gg;
  GG=geom->GG;
  tlo=geom->tlo;
  tup=geom->tup;

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      T[i][j]=0.;

#if (HDVISCOSITY==SHEARVISCOSITY)
  //calculating shear
  ldouble shear[4][4],shearon[4][4];
  calc_shear_comoving(geom->ix,geom->iy,geom->iz,shear,0);
  //calc_shear_lab(geom->ix,geom->iy,geom->iz,shear,0);
  indices_1122(shear,shear,geom->GG);
  
  //to ortonormal
  trans22_cc2on(shear,shearon,geom->tup);

 
  //calculating the viscosity coefficient 
  ldouble xxvec[4]={0.,geom->xx,geom->yy,geom->zz};
  ldouble xxvecBL[4];
  coco_N(xxvec,xxvecBL,MYCOORDS,BLCOORDS);
  
  ldouble fdampr=step_function(xxvecBL[1]-RMINVISC,RMINVISC/10.);
  ldouble pgas=(GAMMA-1.)*pp[UU];

  ldouble eta;
  ldouble Omk = 1./sqrt(xxvecBL[1]*xxvecBL[1]*xxvecBL[1]);

  //kind of alpha p
  eta = ALPHAHDVISC * pgas / Omk;
  
  //damping in radius
  eta*=fdampr;
  
  /*
  if(PROBLEM==30 || PROBLEM==42) //RADNT & RVDONUTIC to overcome huge gradients near rout
    if(geom->ix>=NX-2)
      eta = 0.;  
  */

  //limiting
  ldouble maxspatial=-1.;
  for(i=1;i<4;i++)
    for(j=1;j<4;j++)
      {
	if(fabs(shearon[i][j])>maxspatial)
	  maxspatial=fabs(shearon[i][j]);
      }

  ldouble param=1./3.;
  if(2.*eta*maxspatial > param)
    {
      printf("limiting hd eta: %e->%e at (%d %d %d)\n",2.*eta*maxspatial,param,geom->ix,geom->iy,geom->iz);
      eta = param/2./maxspatial;
    }

  //test  
#ifdef SHEARVISCOSITYONLYRPHI
  //zeroing not rphi components
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	if((i==1 && j==3) || (i==3 && j==1)) 
	  ;
	else
	  shear[i][j]=0.;
      }
#endif

  //to lab frame - only if comoving shear
  boost22_ff2lab(shear,shear,pp,geom->gg,geom->GG);

  //multiply by viscosity to get viscous tensor
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	T[i][j]= -2. * eta * shear[i][j];
      }

  /*
  //comparing with SIMPLEVISC
  print_tensor(T);

  ldouble Ttemp[4][4];

  //test
  for(i=0;i<4;i++)
  for(j=0;j<4;j++)
  {
  Ttemp[i][j]=0.;
  }
  
  Ttemp[1][3]=ALPHAHDVISC*pgas*3.;
  Ttemp[3][1]=ALPHAHDVISC*pgas*3.;

  trans22_on2cc(Ttemp,Ttemp,tlo);
  boost22_ff2lab(Ttemp,Ttemp,pp,gg,GG); 
 
  print_tensor(Ttemp);
  if((xxvecBL[1]>20. && xxvecBL[1]<22.) && geom->iy==NY-1)getchar();
  */

#endif 

#if (HDVISCOSITY==SIMPLEVISCOSITY)  
  ldouble xxvec[4]={0.,geom->xx,geom->yy,geom->zz};
  ldouble xxvecBL[4];
  coco_N(xxvec,xxvecBL,MYCOORDS,BLCOORDS);
  
  ldouble fdampr=step_function(xxvecBL[1]-RMINVISC,RMINVISC/5.);
  //  if(xxvecBL[1]<RMINVISC)) fdampr=0.;
  ldouble pgas=(GAMMA-1.)*pp[UU];

#ifdef RADIATION
#ifdef ALPHATOTALPRESSURE


  /*
  //suppressing radiation pressure by looking at diagonal fluid-frame rad.pressure components  
  ldouble Rij[4][4];
  calc_Rij(pp,ggg,Rij);
  boost22_lab2ff(Rij,Rij,pp,gg,GG);
  trans22_cc2on(Rij,Rij,tup);
  ldouble diag[3]={Rij[1][1],Rij[2][2],Rij[3][3]};
  ldouble maxdiag=my_max(diag[0],my_max(diag[1],diag[2]))/Rij[0][0];
  ldouble THINSUPPPARAM = 0.00001; //prad=0 for maxdiag=0.4
  ldouble prad=1./3.*Rij[0][0]*exp(-(maxdiag-1./3.)*(maxdiag-1./3.)/THINSUPPPARAM);  
  if(geom->ix==NX-1)
    {
      printf("ix: %d %d\n",geom->ix,geom->iy);
      print_Nvector(pp,NV);
      print_tensor(Rij);
      printf("%f\n",maxdiag);
      printf("prad: %e -> %e\n",1./3.*Rij[0][0],
	     1./3.*Rij[0][0]*exp(-(maxdiag-1./3.)*(maxdiag-1./3.)/THINSUPPPARAM));
      getchar();
    }
  */

  //based on local opacities
  ldouble xx[4]={0.,geom->xx,geom->yy,geom->zz};
  ldouble dx[3]={1.,1.,1.}; //opt depth per Rg
  ldouble tautot[3];
  calc_tautot(pp,xx,dx,tautot);

  ldouble Rij[4][4];
  calc_Rij(pp,ggg,Rij);
  // skip boosts what reasonable when gas v<<1
  boost22_lab2ff(Rij,Rij,pp,gg,GG);
  trans22_cc2on(Rij,Rij,tup);
  
  tautot[0]+=1.e-50;

  ldouble PARAM=TAUSUPPRESSPARAM;
  ldouble fdamptau=exp(-PARAM/tautot[0]/tautot[0]);
  ldouble prad=1./3.*Rij[0][0]*fdamptau;

  /*
 if(geom->iy==NY-1)
    {
      printf("ix: %d %d r: %f\n",geom->ix,geom->iy,xxvecBL[1]);
      print_tensor(Rij);
      printf("(%f) fdamps: %e %e\n",tautot[0],fdamptau,fdampr);
      getchar();
    }
  */
  

  ldouble p=pgas+prad;
#else
  ldouble p=pgas;
#endif //ALPHATOTALPRESSURE
#else //RADIATION

  ldouble p=pgas;
  
#endif //RADIATION

  //damping radially
  p*=fdampr;

  T[1][3]=ALPHAHDVISC*p;
  T[3][1]=ALPHAHDVISC*p;

  trans22_on2cc(T,T,tlo);
  boost22_ff2lab(T,T,pp,gg,GG); 


#endif //SIMPLEVISCOSITY

  return 0;

}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//entropy-related routines

ldouble
calc_ufromS(ldouble S,ldouble rho)
{  
  return pow((pow(rho,1./(GAMMAM1)+1.)*exp(S/rho)),GAMMAM1)/(GAMMA-1.);
}

ldouble
calc_Sfromu(ldouble rho,ldouble u)
{
  ldouble indexn=1.0/GAMMAM1;
  return rho*log(pow(GAMMAM1*u,indexn)/pow(rho,indexn+1.0));

  //return rho*log(pow((GAMMAM1*u/rho),1./GAMMAM1)/rho);
  
  //HARM - gives the same result
  //ldouble indexn=1.0/GAMMAM1;
  //printf("entr: %e %e\n",ret,rho*log(pow(GAMMAM1*u,indexn)/pow(rho,indexn+1.0)));getchar();

}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//updates entropy (p[5]) basing on new primitives or stays with the old one if entropy u2p solver was involved
int
update_entropy(int ix,int iy,int iz,int u2pflag)
{
  ldouble gg[4][5],GG[4][5];
  pick_G(ix,iy,iz,GG);
  pick_g(ix,iy,iz,gg);
  ldouble gdet=gg[3][4],gdetu=gdet;;
#if (GDETIN==0) //no metric determinant inside derivatives
  gdetu=1.;
#endif

  ldouble ucon[4],ut,S,Sut,rho,uu;
  int iv;

  //density and energy density
  rho=get_u(p,0,ix,iy,iz);
  uu=get_u(p,1,ix,iy,iz);

  //converts to 4-velocity
  for(iv=1;iv<4;iv++)
    ucon[iv]=get_u(p,iv+1,ix,iy,iz);
  ucon[0]=0.;
  conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);
  ut=ucon[0];

  //u2p_hot worked
  if(u2pflag==0 && uu>0. && rho>0.)
    {
      S=calc_Sfromu(rho,uu);      
      set_u(p,5,ix,iy,iz,S);
      set_u(u,5,ix,iy,iz,S*ut*gdetu); 
    }

  /*
  //u2p_hot didn't work
  //keeping Sut, updating pp[5]
  else
    {
      Sut=get_u(u,5,ix,iy,iz);
      S=Sut/ut;

      //ldouble uint=calc_ufromS(S,rho);


      //      set_u(p,5,ix,iy,iz,S);
    }
  */
  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates shear tensor sigma_ij in the comoving frame of gas or radiation
//whichvel == 0 -> using gas velocity
//whichvel == 1 -> using radiative velocity
int
calc_shear_comoving(int ix,int iy,int iz,ldouble S[][4],int hdorrad)
{
  int i,j,k,iv;

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);

  ldouble (*gg)[5],(*GG)[5],(*tlo)[4],(*tup)[4];
  gg=geom.gg;
  GG=geom.GG;
  tlo=geom.tlo;
  tup=geom.tup;

  //gas or radiation
  int istart,whichvel;
  if(hdorrad==0)
    {
      whichvel=VELPRIM;
      istart=VX;
    }
  else if(hdorrad==1)
    {
      whichvel=VELPRIMRAD;
      istart=FX(0);
    }
 
  //let's start with derivatives
  ldouble du[4][4]; //du_i,j
  ldouble du2[4][4]; //du^i,j

  //time derivatives are zero in the comoving frame!
  for(i=0;i<4;i++)
    {
      du[i][0] = 0.;
      du2[i][0] = 0.;
    }

  //spatial derivatives

  //not to go out of bounds - ghost cell should not use this anyway
  while(ix<=-NG) ix++;
  while(iy<=-NG) iy++;
  while(iz<=-NG) iz++;
  while(ix>=NX+NG-1) ix--;
  while(iy>=NY+NG-1) iy--;
  while(iz>=NZ+NG-1) iz--; 

  ldouble ppm1[NV],ppp1[NV],pp[NV];
  ldouble ggm1[4][5],GGm1[4][5];
  ldouble ggp1[4][5],GGp1[4][5];
  ldouble xxvecm1[4],xxvec[4],xxvecp1[4];
  ldouble uconm1[4],uconp1[4],ucon[4];
  ldouble ucovm1[4],ucovp1[4],ucov[4];
  int idim;

  //four-velocity at cell
  get_xx(ix,iy,iz,xxvec);
  for(iv=0;iv<NV;iv++)
    {
      pp[iv]=get_u(p,iv,ix,iy,iz);
    }
  ucon[1]=pp[istart];  ucon[2]=pp[istart+1];  ucon[3]=pp[istart+2];
  conv_vels(ucon,ucon,whichvel,VEL4,gg,GG);  
  boost2_lab2rf(ucon,ucon,pp,gg,GG); //should give (-1,0,0,0)
  indices_21(ucon,ucov,gg);
   
  //derivatives
  for(idim=1;idim<4;idim++)
    {
      
      if(idim==1)
	{
	  get_xx(ix-1,iy,iz,xxvecm1);
	  get_xx(ix+1,iy,iz,xxvecp1);

	  for(iv=0;iv<NV;iv++)
	    {
	      ppm1[iv]=get_u(p,iv,ix-1,iy,iz);
	      ppp1[iv]=get_u(p,iv,ix+1,iy,iz);
	    }

	  pick_g(ix-1,iy,iz,ggm1);  pick_G(ix-1,iy,iz,GGm1);
	  pick_g(ix+1,iy,iz,ggp1);  pick_G(ix+1,iy,iz,GGp1);
	}

      if(idim==2)
	{
	  get_xx(ix,iy-1,iz,xxvecm1);
	  get_xx(ix,iy+1,iz,xxvecp1);

	  for(iv=0;iv<NV;iv++)
	    {
	      ppm1[iv]=get_u(p,iv,ix,iy-1,iz);
	      ppp1[iv]=get_u(p,iv,ix,iy+1,iz);
	    }

	  pick_g(ix,iy-1,iz,ggm1);  pick_G(ix,iy-1,iz,GGm1);
	  pick_g(ix,iy+1,iz,ggp1);  pick_G(ix,iy+1,iz,GGp1);
	}

     if(idim==3)
	{
	  get_xx(ix,iy,iz-1,xxvecm1);
	  get_xx(ix,iy,iz+1,xxvecp1);

	  for(iv=0;iv<NV;iv++)
	    {
	      ppm1[iv]=get_u(p,iv,ix,iy,iz-1);
	      ppp1[iv]=get_u(p,iv,ix,iy,iz+1);
	    }

	  pick_g(ix,iy,iz-1,ggm1);  pick_G(ix,iy,iz-1,GGm1);
	  pick_g(ix,iy,iz+1,ggp1);  pick_G(ix,iy,iz+1,GGp1);
	}

     //calculating four velocity

     uconm1[1]=ppm1[istart];  uconm1[2]=ppm1[istart+1];  uconm1[3]=ppm1[istart+2];
     uconp1[1]=ppp1[istart];  uconp1[2]=ppp1[istart+1];  uconp1[3]=ppp1[istart+2];

     conv_vels(uconm1,uconm1,whichvel,VEL4,ggm1,GGm1);
     conv_vels(uconp1,uconp1,whichvel,VEL4,ggp1,GGp1);

     //boosting with respect to local cell
     boost2_lab2rf(uconm1,uconm1,pp,ggm1,GGm1); 
     boost2_lab2rf(uconp1,uconp1,pp,ggp1,GGp1);    

     indices_21(uconm1,ucovm1,ggm1);
     indices_21(uconp1,ucovp1,ggp1);
  
     for(i=0;i<4;i++) //second order (ucon=0)
       {
	 du[i][idim]=(ucovp1[i]-ucovm1[i]) / (xxvecp1[idim] - xxvecm1[idim]);
	 du2[i][idim]=(uconp1[i]-uconm1[i]) / (xxvecp1[idim] - xxvecm1[idim]);
       }
    }

  //covariant derivative tensor du_i;j
  ldouble dcu[4][4];
  //covariant derivative tensor du^i;j - only for expansion
  ldouble dcu2[4][4];
  ldouble Krsum;

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Krsum=0.;
	  for(k=0;k<4;k++)
	    Krsum+=get_gKr(k,i,j,ix,iy,iz)*ucov[k];

	  dcu[i][j] = du[i][j] - Krsum;
	}
    
      //only diagonal terms for expansion
      Krsum=0.;
      for(k=0;k<4;k++)
	Krsum+=get_gKr(i,i,k,ix,iy,iz)*ucon[k];
    
      dcu2[i][i] = du2[i][i] + Krsum; 
    }

  //expansion
  ldouble theta=0.;
  for(i=0;i<4;i++)
    theta+=dcu2[i][i];

  //projection tensors P11=P_ij, P21=P^i_j
  ldouble P11[4][4],P21[4][4];
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	P11[i][j] = gg[i][j] + ucov[i]*ucov[j];
	P21[i][j] = delta(i,j) + ucon[i]*ucov[j];
      }

  //the shear tensor sigma_ij
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	ldouble sum1,sum2;
	sum1=sum2=0.;
	for(k=0;k<4;k++)
	  {
	    sum1+=dcu[i][k]*P21[k][j];
	    sum2+=dcu[j][k]*P21[k][i];
	  }
	S[i][j] = 0.5*(sum1+sum2) - 1./3.*theta*P11[i][j];
      }

  return 0;
  //returns shear in sigma_ij in the radiation rest frame!
}
    
   

  

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates shear tensor sigma_ij in the lab frame
//whichvel == 0 -> using gas velocity
//whichvel == 1 -> using radiative velocity
int
calc_shear_lab(int ix,int iy,int iz,ldouble S[][4],int hdorrad)
{
  int i,j,k,iv;

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);

  ldouble (*gg)[5],(*GG)[5],(*tlo)[4],(*tup)[4];
  gg=geom.gg;
  GG=geom.GG;
  tlo=geom.tlo;
  tup=geom.tup;

  //let's start with derivatives
  ldouble du[4][4]; //du_i,j
  ldouble du2[4][4]; //du^i,j

  //time derivatives
  ldouble ucontm1[4],ucovtm1[4],ucontm2[4],ucovtm2[4];

  int istart,whichvel;
  if(hdorrad==0)
    {
      whichvel=VELPRIM;
      istart=VX;
    }
  else if(hdorrad==1)
    {
      whichvel=VELPRIMRAD;
      istart=FX(0);
    }

  ucontm1[0]=ucontm2[0]=0.; //time component will be calculated
  ucontm1[1]=get_u(ptm1,istart,ix,iy,iz);
  ucontm1[2]=get_u(ptm1,istart+1,ix,iy,iz);
  ucontm1[3]=get_u(ptm1,istart+2,ix,iy,iz);
  ucontm2[1]=get_u(ptm2,istart,ix,iy,iz);
  ucontm2[2]=get_u(ptm2,istart+1,ix,iy,iz);
  ucontm2[3]=get_u(ptm2,istart+2,ix,iy,iz);

  conv_vels(ucontm1,ucontm1,whichvel,VEL4,gg,GG);
  conv_vels(ucontm2,ucontm2,whichvel,VEL4,gg,GG);

  indices_21(ucontm1,ucovtm1,gg);
  indices_21(ucontm2,ucovtm2,gg);

  for(i=0;i<4;i++)
    {
#ifndef ZEROTIMEINSHEAR
      if(fabs(ttm1-ttm2) < SMALL)
	{
	  du[i][0] = 0.;
	  du2[i][0] = 0.;
	}
      else
	{
	  du[i][0]=(ucovtm1[i]-ucovtm2[i])/(ttm1-ttm2);
	  du2[i][0]=(ucontm1[i]-ucontm2[i])/(ttm1-ttm2);
	}
#else //force d/dt = 0 in shear
      du[i][0] = 0.;
      du2[i][0] = 0.;
#endif
    }

  //  if(ix==NX-5) {print_4vector(ucontm1);print_4vector(ucontm2);printf("du/dt: %e %e %e %e t12: %e %e\n",du2[0][0],du2[1][0],du2[2][0],du2[3][0],ttm2,ttm1);getchar();}

  //spatial derivatives

  //not to go out of bounds - ghost cell should not use this anyway
  while(ix<=-NG) ix++;
  while(iy<=-NG) iy++;
  while(iz<=-NG) iz++;
  while(ix>=NX+NG-1) ix--;
  while(iy>=NY+NG-1) iy--;
  while(iz>=NZ+NG-1) iz--; 

  ldouble ppm1[NV],ppp1[NV],pp[NV];
  ldouble ggm1[4][5],GGm1[4][5];
  ldouble ggp1[4][5],GGp1[4][5];
  ldouble xxvecm1[4],xxvec[4],xxvecp1[4];
  ldouble uconm1[4],uconp1[4],ucon[4];
  ldouble ucovm1[4],ucovp1[4],ucov[4];
  int idim;

  //four-velocity at cell
  get_xx(ix,iy,iz,xxvec);
  for(iv=0;iv<NV;iv++)
    {
      pp[iv]=get_u(p,iv,ix,iy,iz);
    }
  ucon[1]=pp[istart];  ucon[2]=pp[istart+1];  ucon[3]=pp[istart+2];
  conv_vels(ucon,ucon,whichvel,VEL4,gg,GG);  
  indices_21(ucon,ucov,gg);
   
  //derivatives
  for(idim=1;idim<4;idim++)
    {
      
      if(idim==1)
	{
	  get_xx(ix-1,iy,iz,xxvecm1);
	  get_xx(ix+1,iy,iz,xxvecp1);

	  for(iv=0;iv<NV;iv++)
	    {
	      ppm1[iv]=get_u(p,iv,ix-1,iy,iz);
	      ppp1[iv]=get_u(p,iv,ix+1,iy,iz);
	    }

	  pick_g(ix-1,iy,iz,ggm1);  pick_G(ix-1,iy,iz,GGm1);
	  pick_g(ix+1,iy,iz,ggp1);  pick_G(ix+1,iy,iz,GGp1);
	}

      if(idim==2)
	{
	  get_xx(ix,iy-1,iz,xxvecm1);
	  get_xx(ix,iy+1,iz,xxvecp1);

	  for(iv=0;iv<NV;iv++)
	    {
	      ppm1[iv]=get_u(p,iv,ix,iy-1,iz);
	      ppp1[iv]=get_u(p,iv,ix,iy+1,iz);
	    }

	  pick_g(ix,iy-1,iz,ggm1);  pick_G(ix,iy-1,iz,GGm1);
	  pick_g(ix,iy+1,iz,ggp1);  pick_G(ix,iy+1,iz,GGp1);
	}

     if(idim==3)
	{
	  get_xx(ix,iy,iz-1,xxvecm1);
	  get_xx(ix,iy,iz+1,xxvecp1);

	  for(iv=0;iv<NV;iv++)
	    {
	      ppm1[iv]=get_u(p,iv,ix,iy,iz-1);
	      ppp1[iv]=get_u(p,iv,ix,iy,iz+1);
	    }

	  pick_g(ix,iy,iz-1,ggm1);  pick_G(ix,iy,iz-1,GGm1);
	  pick_g(ix,iy,iz+1,ggp1);  pick_G(ix,iy,iz+1,GGp1);
	}

     //calculating four velocity

     uconm1[1]=ppm1[istart];  uconm1[2]=ppm1[istart+1];  uconm1[3]=ppm1[istart+2];
     uconp1[1]=ppp1[istart];  uconp1[2]=ppp1[istart+1];  uconp1[3]=ppp1[istart+2];

     conv_vels(uconm1,uconm1,whichvel,VEL4,ggm1,GGm1);
     conv_vels(uconp1,uconp1,whichvel,VEL4,ggp1,GGp1);

     indices_21(uconm1,ucovm1,ggm1);
     indices_21(uconp1,ucovp1,ggp1);
  
     for(i=0;i<4;i++)
       {
	 du[i][idim]=(ucovp1[i]-ucovm1[i]) / (xxvecp1[idim] - xxvecm1[idim]);
	 du2[i][idim]=(uconp1[i]-uconm1[i]) / (xxvecp1[idim] - xxvecm1[idim]);
       }
    }

  //covariant derivative tensor du_i;j
  ldouble dcu[4][4];
  //covariant derivative tensor du^i;j - only for expansion
  ldouble dcu2[4][4];
  ldouble Krsum;

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Krsum=0.;
	  for(k=0;k<4;k++)
	    Krsum+=get_gKr(k,i,j,ix,iy,iz)*ucov[k];

	  dcu[i][j] = du[i][j] - Krsum;
	}
    
      //only diagonal terms for expansion
      Krsum=0.;
      for(k=0;k<4;k++)
	Krsum+=get_gKr(i,i,k,ix,iy,iz)*ucon[k];
    
      dcu2[i][i] = du2[i][i] + Krsum; 
    }

  //expansion
  ldouble theta=0.;
  for(i=0;i<4;i++)
    theta+=dcu2[i][i];

  //projection tensors P11=P_ij, P21=P^i_j
  ldouble P11[4][4],P21[4][4];
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	P11[i][j] = gg[i][j] + ucov[i]*ucov[j];
	P21[i][j] = delta(i,j) + ucon[i]*ucov[j];
      }

  //the shear tensor sigma_ij
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	ldouble sum1,sum2;
	sum1=sum2=0.;
	for(k=0;k<4;k++)
	  {
	    sum1+=dcu[i][k]*P21[k][j];
	    sum2+=dcu[j][k]*P21[k][i];
	  }
	S[i][j] = 0.5*(sum1+sum2) - 1./3.*theta*P11[i][j];
      }

  return 0;
}
    
   

  

