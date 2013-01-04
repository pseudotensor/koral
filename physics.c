//KORAL - physics.c
//some problem independent physics

#include "ko.h"

//*************************************************
//calculates left and right wave speeds at cell center
//*************************************************
int
calc_wavespeeds_lr(int ix, int iy, int iz,ldouble *aaa)
{
  //coordinates
  ldouble xx[4];
  xx[0]=0.;
  xx[1]=get_x(ix,0);
  xx[2]=get_x(iy,1);
  xx[3]=get_x(iz,2);

  //metric
  ldouble g00=get_g(g,0,0,ix,iy,iz);
  ldouble g03=get_g(g,0,3,ix,iy,iz);
  ldouble g30=g03;
  ldouble g11=get_g(g,1,1,ix,iy,iz);
  ldouble g22=get_g(g,2,2,ix,iy,iz);  
  ldouble g33=get_g(g,3,3,ix,iy,iz);

  //inversed metric
  ldouble G00=get_g(G,0,0,ix,iy,iz);
  ldouble G03=get_g(G,0,3,ix,iy,iz);
  ldouble G11=get_g(G,1,1,ix,iy,iz);
  ldouble G22=get_g(G,2,2,ix,iy,iz);  
  ldouble G33=get_g(G,3,3,ix,iy,iz);
  ldouble G30=G03;

  //extent of the cell
  //TODO
  ldouble dx[3];
  dx[0]=get_size_x(ix,0)*sqrt(g11);
  dx[1]=get_size_x(iy,1)*sqrt(g22);
  dx[2]=get_size_x(iz,2)*sqrt(g33);   

  ldouble pp[NV],uuu[NV];
  int iv;
  
  ldouble axhdl,axhdr,ayhdl,ayhdr,azhdl,azhdr;
  ldouble axl,axr,ayl,ayr,azl,azr;
  axl=axr=ayl=ayr=azl=azr=1.;
  
  ldouble ucon[4],ucov[4],cst1,cst2,cst3,cst4;

  //picking up primitives 
  for(iv=0;iv<NV;iv++)
    pp[iv]=get_u(p,iv,ix,iy,iz);

  //**********************************************************************
  //***** hydro: speed of sound ******************************************
  //**********************************************************************
	      
  ldouble rho=pp[0];
  ldouble uu=pp[1];
  ldouble vr=pp[2];
  ldouble vth=pp[3];
  ldouble vph=pp[4];

  ldouble pre=(GAMMA-1.)*uu;
  ldouble cs2=GAMMA*pre/(rho+uu+pre);

  if(cs2<0.) cs2=0.;

#ifdef RADIATION
  //**********************************************************************
  //***** radiation: characteristic wave speed ***************************
  //**********************************************************************

  ldouble tautot[3];
  calc_tautot(pp,xx,dx,tautot);

  ldouble E=get_u(p,6,ix,iy,iz);
  ldouble nx=get_u(p,7,ix,iy,iz)/E ;
  ldouble ny=get_u(p,8,ix,iy,iz)/E;
  ldouble nz=get_u(p,9,ix,iy,iz)/E;
  ldouble nl=sqrtl(nx*nx+ny*ny+nz*nz);

  //solving for radiation wave speeds
  max_eigen_lr_radx(nx,ny,nz,0,&axl,&axr);  
  max_eigen_lr_radx(nx,ny,nz,1,&ayl,&ayr);
  max_eigen_lr_radx(nx,ny,nz,2,&azl,&azr);
  
  axr=my_min(axr, 4./3./tautot[0]);
  axl=my_max(axl, -4./3./tautot[0]);
  ayr=my_min(ayr, 4./3./tautot[1]);
  ayl=my_max(ayl, -4./3./tautot[1]);
  azr=my_min(azr, 4./3./tautot[2]);
  azl=my_max(azl, -4./3./tautot[2]);
#endif

  //**********************************************************************
  //***** other stuff ****************************************************
  //**********************************************************************
  ldouble ut2=-1./(g00 + 2.*vph*g03 + vr*vr*g11 + vth*vth*g22 + vph*vph*g33 );
  if(ut2<0.) my_err("ut2.lt.0 in wavespeeds\n");
  
  ldouble ut=sqrtl(ut2);  
  ldouble ur=vr*ut;
  ldouble uth=vth*ut;
  ldouble uph=vph*ut;

  ucon[0]=ut;
  ucon[1]=ur;
  ucon[2]=uth;
  ucon[3]=uph;
  
  ucov[0]=g00*ucon[0]+g03*ucon[3];
  ucov[1]=g11*ucon[1];
  ucov[2]=g22*ucon[2];
  ucov[3]=g33*ucon[3]+g03*ucon[0];


  //**********************************************************************
  //algorithm from HARM to transform the fluid frame wavespeed into lab frame
  //**********************************************************************

  ldouble Acov[4],Acon[4],Bcov[4],Bcon[4],Asq,Bsq,Au,Bu,AB,Au2,Bu2,AuBu,A,B,discr,wspeed2;

  //**********************************************************************
  //**********************************************************************
  //x
  Acov[0]=0.;
  Acov[1]=1.;
  Acov[2]=0.;
  Acov[3]=0.;
  Acon[0]=G00*Acov[0]+G03*Acov[3];
  Acon[1]=G11*Acov[1];
  Acon[2]=G22*Acov[2];
  Acon[3]=G33*Acov[3]+G03*Acov[0];
  
  Bcov[0]=1.;
  Bcov[1]=0.;
  Bcov[2]=0.;
  Bcov[3]=0.;
  Bcon[0]=G00*Bcov[0]+G03*Bcov[3];
  Bcon[1]=G11*Bcov[1];
  Bcon[2]=G22*Bcov[2];
  Bcon[3]=G33*Bcov[3]+G03*Bcov[0];


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
  discr = sqrtl(discr);
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

#ifdef RADIATION
  wspeed2=axl*axl;
  B = 2. * (AuBu * (1.0-wspeed2)  - AB*wspeed2);
  A = Bu2 * (1.0 - wspeed2) - Bsq * wspeed2;
  discr = 4.0 * wspeed2 * ((AB * AB - Asq * Bsq) * wspeed2 + (2.0 * AB * Au * Bu - Asq * Bu2 - Bsq * Au2) * (wspeed2 - 1.0));
  if(discr<0.) {printf("x1discr in ravespeeds lt 0\n"); discr=0.;}
  discr = sqrtl(discr);
  cst1 = -(-B + discr) / (2. * A);
  cst2 = -(-B - discr) / (2. * A);  

  wspeed2=axr*axr;
  B = 2. * (AuBu * (1.0-wspeed2)  - AB*wspeed2);
  A = Bu2 * (1.0 - wspeed2) - Bsq * wspeed2;
  discr = 4.0 * wspeed2 * ((AB * AB - Asq * Bsq) * wspeed2 + (2.0 * AB * Au * Bu - Asq * Bu2 - Bsq * Au2) * (wspeed2 - 1.0));
  if(discr<0.) {printf("x2discr in ravespeeds lt 0\n"); discr=0.;}
  discr = sqrtl(discr);
  cst3 = -(-B + discr) / (2. * A);
  cst4 = -(-B - discr) / (2. * A);  

  axl = my_min(my_min(cst1,cst2),my_min(cst3,cst4));
  axr = my_max(my_max(cst1,cst2),my_max(cst3,cst4));
#endif
     
  //**********************************************************************
  //**********************************************************************
  //y
  Acov[0]=0.;
  Acov[1]=0.;
  Acov[2]=1.;
  Acov[3]=0.;
  Acon[0]=G00*Acov[0]+G03*Acov[3];
  Acon[1]=G11*Acov[1];
  Acon[2]=G22*Acov[2];
  Acon[3]=G33*Acov[3]+G03*Acov[0];
  
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
  discr = sqrtl(discr);
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
 
#ifdef RADIATION
  wspeed2=ayl*ayl;
  B = 2. * (AuBu * (1.0-wspeed2)  - AB*wspeed2);
  A = Bu2 * (1.0 - wspeed2) - Bsq * wspeed2;
  discr = 4.0 * wspeed2 * ((AB * AB - Asq * Bsq) * wspeed2 + (2.0 * AB * Au * Bu - Asq * Bu2 - Bsq * Au2) * (wspeed2 - 1.0));
  if(discr<0.) {printf("y1discr in ravespeeds lt 0\n"); discr=0.;}
  discr = sqrtl(discr);
  cst1 = -(-B + discr) / (2. * A);
  cst2 = -(-B - discr) / (2. * A);  

  wspeed2=ayr*ayr;
  B = 2. * (AuBu * (1.0-wspeed2)  - AB*wspeed2);
  A = Bu2 * (1.0 - wspeed2) - Bsq * wspeed2;
  discr = 4.0 * wspeed2 * ((AB * AB - Asq * Bsq) * wspeed2 + (2.0 * AB * Au * Bu - Asq * Bu2 - Bsq * Au2) * (wspeed2 - 1.0));
  if(discr<0.) {printf("y1discr in ravespeeds lt 0\n"); discr=0.;}
  discr = sqrtl(discr);
  cst3 = -(-B + discr) / (2. * A);
  cst4 = -(-B - discr) / (2. * A);  

  ayl = my_min(my_min(cst1,cst2),my_min(cst3,cst4));
  ayr = my_max(my_max(cst1,cst2),my_max(cst3,cst4));
#endif
     
   
  //**********************************************************************
  //**********************************************************************
  //z
  Acov[0]=0.;
  Acov[1]=0.;
  Acov[2]=0.;
  Acov[3]=1.;
  Acon[0]=G00*Acov[0]+G03*Acov[3];
  Acon[1]=G11*Acov[1];
  Acon[2]=G22*Acov[2];
  Acon[3]=G33*Acov[3]+G03*Acov[0];
  
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
  discr = sqrtl(discr);
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
  wspeed2=azl*azl;
  B = 2. * (AuBu * (1.0-wspeed2)  - AB*wspeed2);
  A = Bu2 * (1.0 - wspeed2) - Bsq * wspeed2;
  discr = 4.0 * wspeed2 * ((AB * AB - Asq * Bsq) * wspeed2 + (2.0 * AB * Au * Bu - Asq * Bu2 - Bsq * Au2) * (wspeed2 - 1.0));
  if(discr<0.) {printf("z1discr in ravespeeds lt 0\n"); discr=0.;}
  discr = sqrtl(discr);
  cst1 = -(-B + discr) / (2. * A);
  cst2 = -(-B - discr) / (2. * A);  

  wspeed2=azr*azr;
  B = 2. * (AuBu * (1.0-wspeed2)  - AB*wspeed2);
  A = Bu2 * (1.0 - wspeed2) - Bsq * wspeed2;
  discr = 4.0 * wspeed2 * ((AB * AB - Asq * Bsq) * wspeed2 + (2.0 * AB * Au * Bu - Asq * Bu2 - Bsq * Au2) * (wspeed2 - 1.0));
  if(discr<0.) {printf("z2discr in ravespeeds lt 0\n"); discr=0.;}
  discr = sqrtl(discr);
  cst3 = -(-B + discr) / (2. * A);
  cst4 = -(-B - discr) / (2. * A);  

  azl = my_min(my_min(cst1,cst2),my_min(cst3,cst4));
  azr = my_max(my_max(cst1,cst2),my_max(cst3,cst4));
#endif
  //**********************************************************************
  //**********************************************************************
    
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

//***************************************
//returns external source terms
//***************************************
int f_source_term(int ix, int iy, int iz,ldouble *ss)
{
  int iv;
  for(iv=0;iv<NV;iv++)
    ss[iv]=0.;

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

  ldouble gg[4][5],ggxl[4][5],ggxr[4][5];

  pick_g(ix,iy,iz,gg);
  ldouble gdet=gg[3][4];
  pick_gb(ix,iy,iz,0,ggxl);
  pick_gb(ix+1,iy,iz,0,ggxr);

  ldouble g00=gg[0][0];
  ldouble g03=gg[0][3];
  ldouble g11=gg[1][1];
  ldouble g22=gg[2][2];
  ldouble g33=gg[3][3];

  ldouble dlgdet[3];
  dlgdet[0]=gg[0][4]; //D[gdet,x1]/gdet
  dlgdet[1]=gg[1][4]; //D[gdet,x2]/gdet
  dlgdet[2]=gg[2][4]; //D[gdet,x3]/gdet
  
  ldouble gam=GAMMA;
  ldouble ut;
  ldouble T[4][4];

  //calculating stress energy tensor components
  calc_Tmunu(pp,gg,T,&ut);

  int ii, jj;
  for(ii=0;ii<4;ii++)
    for(jj=0;jj<4;jj++)
      {
	if(isnan(T[ii][jj])) 
	  {
	    printf("%d %d %Le\n",ii,jj,T[ii][jj]);
	    my_err("nan in metric_source_terms\n");
	  }
      }
 
  ldouble rho=pp[0];
  ldouble u=pp[1];
  ldouble v1=pp[2];
  ldouble v2=pp[3];
  ldouble v3=pp[4];
  ldouble S=pp[5];

  ldouble ucon[4];
  ldouble ut2=-1./(g00 + 2.*v3*g03 + v1*v1*g11 + v3*v3*g33 + v2*v2*g22);
  if(ut2<0.) 
    my_err("ut2.lt.0 in metric_source_terms\n");
  ut=sqrtl(ut2);

  ucon[0]=ut;
  ucon[1]=v1*ut;
  ucon[2]=v2*ut;
  ucon[3]=v3*ut;
 
  int k,l,iv;

  for(iv=0;iv<NV;iv++)
    ss[iv]=0.;

  /***************************************************/
#ifdef RADIATION
  /***************************************************/

  ldouble tup[4][4],tlo[4][4];
  pick_T(tmuup,ix,iy,iz,tup);
  pick_T(tmulo,ix,iy,iz,tlo);
  ldouble Rij[4][4];
  calc_Rij_ff(pp,Rij);
  //boost22_ff2zamo(Rij,Rij,p,g,eup);
  //trans22_zamo2lab(Rij,Rij,g,elo);  
  trans22_on2cc(Rij,Rij,gg,tlo);
  boost22_ff2lab(Rij,Rij,pp,gg);
  indices_2221(Rij,Rij,gg);

  //terms with Christoffels
  for(k=0;k<4;k++)
    for(l=0;l<4;l++)
      {
	ss[1]+=T[k][l]*get_gKr(l,0,k,ix,iy,iz);
	ss[2]+=T[k][l]*get_gKr(l,1,k,ix,iy,iz);
	ss[3]+=T[k][l]*get_gKr(l,2,k,ix,iy,iz);
	ss[4]+=T[k][l]*get_gKr(l,3,k,ix,iy,iz);
	ss[6]+=Rij[k][l]*get_gKr(l,0,k,ix,iy,iz);
	ss[7]+=Rij[k][l]*get_gKr(l,1,k,ix,iy,iz);
	ss[8]+=Rij[k][l]*get_gKr(l,2,k,ix,iy,iz);
	ss[9]+=Rij[k][l]*get_gKr(l,3,k,ix,iy,iz);
      }

  //terms with dloggdet
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

  /***************************************************/
#else
  /***************************************************/

  //terms with Christoffels
  for(k=0;k<4;k++)
    for(l=0;l<4;l++)
      {
	ss[1]+=T[k][l]*get_gKr(l,0,k,ix,iy,iz);
	ss[2]+=T[k][l]*get_gKr(l,1,k,ix,iy,iz);
	ss[3]+=T[k][l]*get_gKr(l,2,k,ix,iy,iz);
	ss[4]+=T[k][l]*get_gKr(l,3,k,ix,iy,iz);
      }

  //terms with dloggdet  
  for(l=1;l<4;l++)
    {
      ss[0]+=-dlgdet[l-1]*rho*ucon[l];
      ss[1]+=-dlgdet[l-1]*(T[l][0]+rho*ucon[l]);
      ss[2]+=-dlgdet[l-1]*(T[l][1]);
      ss[3]+=-dlgdet[l-1]*(T[l][2]);
      ss[4]+=-dlgdet[l-1]*(T[l][3]);
      ss[5]+=-dlgdet[l-1]*S*ucon[l];
    }   

  /***************************************************/
#endif
  /***************************************************/

  return 0;
}

//***************************************
// calculates fluxes
//***************************************
ldouble f_flux_prime( ldouble *pp, int idim, int ix, int iy, int iz,ldouble *ff)
{  
  int iv;
  for(iv=0;iv<NV;iv++)
    ff[iv]=0.;

  //picking up metric from a cell face  
  ldouble gg[4][5];
  pick_gb(ix,iy,iz,idim,gg);
  ldouble gdet=gg[3][4];
  
  //calculating Tmunu
  ldouble ut;
  ldouble T[4][4];
  calc_Tmunu(pp,gg,T,&ut);

  //primitives
  ldouble rho=pp[0];
  ldouble u=pp[1];
  ldouble v1=pp[2];
  ldouble v2=pp[3];
  ldouble v3=pp[4];
  ldouble S=pp[5];

  if(isnan(S)) 
    {
      printf("nan S: %d %d %d S: %Le\n",ix,iy,iz,S);
      my_err("nan S in flux prime\n");
    }

  ldouble u1=v1*ut;
  ldouble u2=v2*ut;
  ldouble u3=v3*ut;
  ldouble gam=GAMMA;

  int ii, jj;
  for(ii=0;ii<4;ii++)
    for(jj=0;jj<4;jj++)
      {
	if(isnan(T[ii][jj])) 
	  {
	    printf("nan tmunu: %d %d %Le\n",ii,jj,T[ii][jj]);
	    print_Nvector(pp,NV);
	    my_err("nan in flux_prime\n");
	  }
      }
 
#ifdef RADIATION
  ldouble tup[4][4],tlo[4][4];
  pick_T(tmuup,ix,iy,iz,tup);
  pick_T(tmulo,ix,iy,iz,tlo);
  ldouble Rij[4][4];
  calc_Rij_ff(pp,Rij);
  //boost22_ff2zamo(Rij,Rij,p,g,eup);
  //trans22_zamo2lab(Rij,Rij,g,elo);  
  trans22_on2cc(Rij,Rij,gg,tlo);
  boost22_ff2lab(Rij,Rij,pp,gg);
  indices_2221(Rij,Rij,gg);

  //to move gdet in/out derivative:
  //here, up in metric source terms, in u2p and p2u, as well as in finite.c with del4[]

   if(idim==0) //x
    {
      ff[0]= rho*u1;

      ff[1]= (T[1][0]+rho*u1);

      ff[2]= (T[1][1]);

      ff[3]= (T[1][2]);
 
      ff[4]= (T[1][3]);

      ff[5]= S*u1;

      ff[6]= Rij[1][0];
      
      ff[7]= Rij[1][1];
      
      ff[8]= Rij[1][2];
      
      ff[9]= Rij[1][3];
    }  
  if(idim==1) //y
    {
      ff[0]= rho*u2;

      ff[1]= (T[2][0]+rho*u2);

      ff[2]= (T[2][1]);

      ff[3]= (T[2][2]);

      ff[4]= (T[2][3]);

      ff[5]= S*u2;
 
      ff[6]= Rij[2][0];
      
      ff[7]= Rij[2][1];
      
      ff[8]= Rij[2][2];
      
      ff[9]= Rij[2][3];
    }  
  if(idim==2) //z
    {
      ff[0]= rho*u3;

      ff[1]= (T[3][0]+rho*u3);
 
      ff[2]= (T[3][1]);

      ff[3]= (T[3][2]);
 
      ff[4]= (T[3][3]);

      ff[5]= S*u3;
 
      ff[6]= Rij[3][0];
      
      ff[7]= Rij[3][1];
      
      ff[8]= Rij[3][2];
       
      ff[9]= Rij[3][3];
    } 

#else //pure hydro

  if(idim==0) //x
    {
      ff[0]= rho*u1;

      ff[1]= (T[1][0]+rho*u1);

      ff[2]= (T[1][1]);

      ff[3]= (T[1][2]);
 
      ff[4]= (T[1][3]);

      ff[5]= S*u1;
    }  
  if(idim==1) //y
    {
      ff[0]= rho*u2;

      ff[1]= (T[2][0]+rho*u2);

      ff[2]= (T[2][1]);

      ff[3]= (T[2][2]);

      ff[4]= (T[2][3]);

      ff[5]= S*u2;
    }  
  if(idim==2) //z
    {
      ff[0]= rho*u3;

      ff[1]= (T[3][0]+rho*u3);
 
      ff[2]= (T[3][1]);

      ff[3]= (T[3][2]);
 
      ff[4]= (T[3][3]);

      ff[5]= S*u3;
    } 

#endif

  return 0.;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//entropy-related routines

ldouble
calc_ufromS(ldouble S,ldouble rho)
{  
  return powl((powl(rho,1./(GAMMAM1)+1.)*exp(S/rho)),GAMMAM1)/(GAMMA-1.);
}

ldouble
calc_Sfromu(ldouble rho,ldouble u)
{
  ldouble ret = rho*log(pow((GAMMAM1*u/rho),1./GAMMAM1)/rho);
  return ret;
}

