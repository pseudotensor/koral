//KORAL - physics.c
//some problem independent physics

#include "ko.h"

//*************************************************
//calculates left and right wave speeds at cell center
//*************************************************
int
calc_wavespeeds_lr_core(ldouble *ucon,ldouble GG[][5],ldouble *aret,ldouble wspeed2,int idim)
{
  ldouble Acov[4],Acon[4],Bcov[4],Bcon[4],Asq,Bsq,Au,Bu,AB,Au2,Bu2,AuBu,A,B,discr,cst1,cst2;
  
  Acov[0]=0.;
  Acov[1]=0.;
  Acov[2]=0.;
  Acov[3]=0.;

  Acov[idim+1]=1.;

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

  //wspeed2=cs2;
  B = 2. * (AuBu * (1.0-wspeed2)  - AB*wspeed2);
  A = Bu2 * (1.0 - wspeed2) - Bsq * wspeed2;
  discr = 4.0 * wspeed2 * ((AB * AB - Asq * Bsq) * wspeed2 + (2.0 * AB * Au * Bu - Asq * Bu2 - Bsq * Au2) * (wspeed2 - 1.0));

  if(discr<0.) {printf("discr in wavespeeds lt 0\n"); return -1;}
  discr = sqrt(discr);
  cst1 = -(-B + discr) / (2. * A);
  cst2 = -(-B - discr) / (2. * A);  
  if(cst2>cst1)
    {
      aret[0]=cst1;  aret[1]=cst2;
    }
  else
    {
      aret[0]=cst1;  aret[1]=cst2;
    }
  return 0;
}

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
  
  ldouble utcon[4],ucon[4],ucov[4],cst1,cst2,cst3,cst4;
  ldouble bcon[4],bcov[4],bsq;
  ldouble cs2,va2,EF,EEloc; 
  ldouble rho,uu,pre;

  //**********************************************************************
  //***** four velocity **************************************************
  //**********************************************************************

  for(iv=1;iv<4;iv++)
    utcon[iv]=pp[1+iv];
  utcon[0]=0.;
  conv_vels_both(utcon,ucon,ucov,VELPRIM,VEL4,gg,GG);
 
  //**********************************************************************
  //***** hydro: speed of sound ******************************************
  //**********************************************************************
	      
  rho=pp[RHO];
  uu=pp[UU];
 
  pre=(GAMMA-1.)*uu;
  cs2=GAMMA*pre/(rho+uu+pre);
  if(cs2<0.) cs2=0.;

  //**********************************************************************
  //***** magn: alvenic speed ****** *************************************
  //**********************************************************************
	      
  va2=0.;

#ifdef MAGNFIELD
  calc_bcon_4vel(pp,ucon,ucov,bcon);
  indices_21(bcon,bcov,gg); 
  bsq = dot(bcon,bcov);
  EF = rho + GAMMA*uu;
  EEloc = bsq + EF ;
  va2 = bsq/EEloc ;
  if(va2<0.) va2=0.;
#endif

  //**********************************************************************
  //***** mhd: fast magnetosonic speed ***********************************
  //**********************************************************************

  ldouble vtot2; //total characteristic velocity
  vtot2=cs2 + va2 - cs2*va2;

  //**********************************************************************
  //algorithm from HARM to transform the fluid frame wavespeed into lab frame
  //**********************************************************************

  ldouble aret[2];
  int ret;
  ret=calc_wavespeeds_lr_core(ucon,GG,aret,vtot2,0);
  if(ret<0) {printf("error wsx at %d | %d | %d\n",geom->ix,geom->iy,geom->iz);}
  axhdl=aret[0];
  axhdr=aret[1];
  
  ret=calc_wavespeeds_lr_core(ucon,GG,aret,vtot2,1);
  if(ret<0) {printf("error wsy at %d | %d | %d\n",geom->ix,geom->iy,geom->iz);}
  ayhdl=aret[0];
  ayhdr=aret[1];
  
  ret=calc_wavespeeds_lr_core(ucon,GG,aret,vtot2,2);
  if(ret<0) {printf("error wsz at %d | %d | %d\n",geom->ix,geom->iy,geom->iz);}
  azhdl=aret[0];
  azhdr=aret[1];
 
#ifdef RADIATION
  //**********************************************************************
  //***** radiation: characteristic wave speed ***************************
  //**********************************************************************

  ldouble aval[6];
  int verbose=0;

  //physical size of the cell
  ldouble dx[3];
  ldouble xx[4]={0.,geom->xx,geom->yy,geom->zz};
  
  //ix,iy,iz could be the indices of a face, so the depth taken from left/right
  dx[0]=my_max(get_size_x(geom->ix,0)*sqrt(gg[1][1]),get_size_x(geom->ix+1,0)*sqrt(gg[1][1]));
  dx[1]=my_max(get_size_x(geom->iy,1)*sqrt(gg[2][2]),get_size_x(geom->iy+1,1)*sqrt(gg[2][2]));
  dx[2]=my_max(get_size_x(geom->iz,2)*sqrt(gg[3][3]),get_size_x(geom->iz+1,2)*sqrt(gg[3][3]));
  ldouble tautot[3];
  calc_tautot(pp,geom,dx,tautot);

  //M1
  calc_rad_wavespeeds(pp,geom,tautot,aval,verbose);

  axl=aval[0];
  axr=aval[1];
  ayl=aval[2];
  ayr=aval[3];
  azl=aval[4];
  azr=aval[5];

#endif

#ifdef OVERWRITERADWAVESPEEDSWITHHD
  axl=axhdl;
  axr=axhdr;
  ayl=ayhdl;
  ayr=ayhdr;
  azl=azhdl;
  azr=azhdr;
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
  //hd:
  aaa[0]=axhdl;
  aaa[1]=axhdr;
  aaa[2]=ayhdl;
  aaa[3]=ayhdr;
  aaa[4]=azhdl;
  aaa[5]=azhdr;
  //rad:
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

  for(i=0;i<NV;i++)
    ss[i]=0.;

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
    {
      pp[i]=get_u(p,i,ix,iy,iz);  
    }
  
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

  ldouble Rij[4][4];
  calc_Rij(pp,geom,Rij); //R^ij
  indices_2221(Rij,Rij,gg); //R^i_j
  ldouble urfcon[4];
  urfcon[0]=0.;
  urfcon[1]=pp[FX0];
  urfcon[2]=pp[FY0];
  urfcon[3]=pp[FZ0];
  conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,gg,GG);

  //terms with Christoffels
  for(k=0;k<4;k++)
    for(l=0;l<4;l++)
      {
	ss[1]+=gdetu*T[k][l]*get_gKr(l,0,k,ix,iy,iz);
	ss[2]+=gdetu*T[k][l]*get_gKr(l,1,k,ix,iy,iz);
	ss[3]+=gdetu*T[k][l]*get_gKr(l,2,k,ix,iy,iz);
	ss[4]+=gdetu*T[k][l]*get_gKr(l,3,k,ix,iy,iz);
	ss[EE0]+=gdetu*Rij[k][l]*get_gKr(l,0,k,ix,iy,iz);
	ss[FX0]+=gdetu*Rij[k][l]*get_gKr(l,1,k,ix,iy,iz);
	ss[FY0]+=gdetu*Rij[k][l]*get_gKr(l,2,k,ix,iy,iz);
	ss[FZ0]+=gdetu*Rij[k][l]*get_gKr(l,3,k,ix,iy,iz);
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
      ss[EE0]+=-dlgdet[l-1]*(Rij[l][0]);
      ss[FX0]+=-dlgdet[l-1]*(Rij[l][1]);
      ss[FY0]+=-dlgdet[l-1]*(Rij[l][2]);
      ss[FZ0]+=-dlgdet[l-1]*(Rij[l][3]);
      #ifdef NCOMPTONIZATION
      ss[NF0]+=-dlgdet[l-1]*pp[NF0]*urfcon[l];
      #endif
    }


#endif //GDETIN

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
int f_general_source_term_arb(ldouble *pp,void *ggg,ldouble *ss)
{
  int i;

  struct geometry *geom
    = (struct geometry *) ggg;

  int ix,iy,iz,iv;
  ix=geom->ix;
  iy=geom->iy;
  iz=geom->iz;
 
  PLOOP(iv) ss[iv]=0.;

  /***************************************************/

  return 0;
}


//***************************************
//returns geometrical source terms for all conserved quantities
//***************************************
int f_metric_source_term(int ix, int iy, int iz,ldouble *ss)
{
  int i;

  //ldouble pp[NV];  
  //  for(i=0;i<NV;i++)
  //pp[i]=get_u(p,i,ix,iy,iz);  
  
  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);

  f_metric_source_term_arb(&get_u(p,0,ix,iy,iz),&geom,ss);

  return 0;
}

//***************************************
//returns general source terms for all conserved quantities
//***************************************
int f_general_source_term(int ix, int iy, int iz,ldouble *ss)
{
  int i;

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);

  f_general_source_term_arb(&get_u(p,0,ix,iy,iz),&geom,ss);

  return 0;
}

//***************************************
// calculates fluxes at faces
//***************************************
int f_flux_prime( ldouble *pp, int idim, int ix, int iy, int iz,ldouble *ff,int lr)
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
#ifdef TRACER
  ldouble tracer=pp[TRA];
#endif
  ldouble rho=pp[RHO];
  ldouble u=pp[UU];

  ldouble vcon[4],ucon[4],ucov[4],bcon[4],bcov[4],bsq=0.;


  vcon[1]=pp[2];
  vcon[2]=pp[3];
  vcon[3]=pp[4];
  ldouble S=pp[5];

  //converting to 4-velocity
  conv_vels_both(vcon,ucon,ucov,VELPRIM,VEL4,gg,GG);

#ifdef MAGNFIELD
  calc_bcon_4vel(pp,ucon,ucov,bcon);
  indices_21(bcon,bcov,gg); 
  bsq = dot(bcon,bcov);
#endif


  ldouble pre=(GAMMA-1.)*u; 
  ldouble w=rho+u+pre;
  ldouble eta=w+bsq;
  ldouble etap = u+pre+bsq; //eta-rho


  int ii, jj, irf;
  for(ii=0;ii<4;ii++)
    for(jj=0;jj<4;jj++)
      {
	if(isnan(T[ii][jj])) 
	  {
	    printf("nan tmunu: %d %d %e at %d %d %d\n",ii,jj,T[ii][jj],ix,iy,iz);
	    print_Nvector(pp,NV);
	    my_err("nan in flux_prime\n");
	  }
      }


  ldouble utp1=calc_utp1(vcon,ucon,&geom);

  //fluxes

  //hydro
  ff[0]= gdetu*rho*ucon[idim+1];
  //ff[1]= gdetu*(T[idim+1][0]+rho*ucon[idim+1]);
  //to avoid slow cancellation:
  ff[1]= gdetu*(etap*ucon[idim+1]*ucov[0] + rho*ucon[idim+1]*utp1);
#ifdef MAGNFIELD
  ff[1]+= - gdetu*bcon[idim+1]*bcov[0];
#endif

  ff[2]= gdetu*(T[idim+1][1]);
  ff[3]= gdetu*(T[idim+1][2]); 
  ff[4]= gdetu*(T[idim+1][3]);
  ff[5]= gdetu*S*ucon[idim+1];

#ifdef TRACER
  ff[TRA]= gdetu*tracer*ucon[idim+1]; 
#endif

#ifdef MAGNFIELD
  //magnetic
  ff[B1]=gdetu*(bcon[1]*ucon[idim+1] - bcon[idim+1]*ucon[1]);
  ff[B2]=gdetu*(bcon[2]*ucon[idim+1] - bcon[idim+1]*ucon[2]);
  ff[B3]=gdetu*(bcon[3]*ucon[idim+1] - bcon[idim+1]*ucon[3]);

  //test
  //ff[B1]=ff[B2]=ff[B3]=0.;
#endif

  //radiation
#ifdef RADIATION

  if(RADCLOSURE==VETCLOSURE)
    {
      //to pick up intensities from the cell center correspondingly to left/right biased fluxes
   
      if(lr==0)//left biased
	{
	  geom.par=lr;
	  //if(idim==0) geom.ix--;
	  //if(idim==1) geom.iy--;
	  //if(idim==2) geom.iz--;
	}
	
      //fill_geometry(ix,iy,iz,&geom);
      
    }
 
  f_flux_prime_rad(pp,idim,&geom,ff);
#endif

  return 0;
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
  ldouble utcon[4],ucon[4],ucov[4];  
  ldouble bcon[4],bcov[4],bsq=0.;
  
  //converts to 4-velocity
  for(iv=1;iv<4;iv++)
    utcon[iv]=pp[1+iv];
  utcon[0]=0.;
  conv_vels_both(utcon,ucon,ucov,VELPRIM,VEL4,gg,GG);

#ifdef MAGNFIELD
  calc_bcon_4vel(pp,ucon,ucov,bcon);
  indices_21(bcon,bcov,gg); 
  bsq = dot(bcon,bcov);
#else
  bcon[0]=bcon[1]=bcon[2]=bcon[3]=0.;
  bsq=0.;
#endif
  
  ldouble p=(GAMMA-1.)*uu; 
  ldouble w=rho+uu+p;
  ldouble eta=w+bsq;
  ldouble ptot=p+0.5*bsq;
  
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      T[i][j]=eta*ucon[i]*ucon[j] + ptot*GG[i][j] - bcon[i]*bcon[j];

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

  return 0;
}

/*******************************************/

ldouble calc_PEQ_rhofromTu(ldouble T,ldouble u)
{
  ldouble p=GAMMAM1*u;
  ldouble rho = p/K_BOLTZ/T*MU_GAS*M_PROTON;

  return rho;
}

ldouble calc_PEQ_csfromT(ldouble T)
{
  //ldouble p=K_BOLTZ*rho*T/MU_GAS/M_PROTON;	      
  ldouble cs = sqrt(GAMMA*K_BOLTZ*T/MU_GAS/M_PROTON);
  return cs;
}

ldouble calc_PEQ_ufromTrho(ldouble T,ldouble rho)
{
  ldouble p=K_BOLTZ*rho*T/MU_GAS/M_PROTON;	      
  ldouble u=p/(GAMMA-1.);
  return u;
}

ldouble calc_PEQ_Tfromurho(ldouble u,ldouble rho)
{
  ldouble p=u*(GAMMA-1.);
  ldouble T=p/(K_BOLTZ*rho/MU_GAS/M_PROTON);
  return T;
}
