
int SPHboundary(ldouble *pp, void *ggg, void *gggBL)
{

  struct geometry *geom
     = (struct geometry *) ggg;
   
   struct geometry *geomBL
     = (struct geometry *) gggBL;
 
  ldouble th=geomBL->xxvec[2];
  ldouble r=geomBL->xxvec[1];

  //ambient
  set_hdatmosphere(pp,geom->xxvec,geom->gg,geom->GG,0);
#ifdef RADIATION
  set_radatmosphere(pp,geom->xxvec,geom->gg,geom->GG,0);
#endif

  int ix,iy,iz;
  ix=geom->ix;
  iy=geom->iy;
  iz=geom->iz;

  int iys,izs;  //shifted to fit the projection array

  //temporary, to avoid problems in corners
  //  if(iy<0) iy=0;
  //if(iy>=NY) iy=NY-1;
  //if(iz<0) iz=0;
  //if(iz>=NZ) iz=NZ-1;

  iys=iy+NG;
  izs=iz+NG;

  //printf("%d %d %d > %d %d\n",ix,iy,iz,PROCID,TI);

  double rhoamb,rho0,rho1,vr0,vr1,vth0,vth1,vph0,vph1,temp0,temp1;
  double rho,temp,vr,vth,vph;
  double ucon[4];

  rhoamb=pp[RHO];
  rho0=rhoCGS2GU(SPHdata0[SPHprojection[iys][izs][0]][SPHprojection[iys][izs][1]][0]);
  rho1=rhoCGS2GU(SPHdata1[SPHprojection[iys][izs][0]][SPHprojection[iys][izs][1]][0]);

  //verify if we have non-zero SPH input with density exceeding the ambient density
  if(rho0<rhoamb && rho1<rhoamb)
    return -1; //outside disk;
   
  temp0=tempCGS2GU(SPHdata0[SPHprojection[iys][izs][0]][SPHprojection[iys][izs][1]][1]);
  temp1=tempCGS2GU(SPHdata1[SPHprojection[iys][izs][0]][SPHprojection[iys][izs][1]][1]);

  vr0=SPHdata0[SPHprojection[iys][izs][0]][SPHprojection[iys][izs][1]][2]/CCC0;
  vr1=SPHdata1[SPHprojection[iys][izs][0]][SPHprojection[iys][izs][1]][2]/CCC0;

  #ifdef SCALEVR
  vr0*=sqrt(SPHRADIUS/geomBL->xx);
  vr1*=sqrt(SPHRADIUS/geomBL->xx);
  #endif

  ldouble rcgs=geomBL->xx * GMC2;

  //angular velocities given in rad/s
  vth0=SPHdata0[SPHprojection[iys][izs][0]][SPHprojection[iys][izs][1]][3]*GMC3;
  vth1=SPHdata1[SPHprojection[iys][izs][0]][SPHprojection[iys][izs][1]][3]*GMC3;

  vph0=SPHdata0[SPHprojection[iys][izs][0]][SPHprojection[iys][izs][1]][4]*GMC3;
  vph1=SPHdata1[SPHprojection[iys][izs][0]][SPHprojection[iys][izs][1]][4]*GMC3;

  //interpolation
  rho = rho1 - (rho1 - rho0) * (SPHtime1 - global_time) / (SPHtime1 - SPHtime0); 
  temp = temp1 - (temp1 - temp0) * (SPHtime1 - global_time) / (SPHtime1 - SPHtime0); 
  vr = vr1 - (vr1 - vr0) * (SPHtime1 - global_time) / (SPHtime1 - SPHtime0); 
  vth = vth1 - (vth1 - vth0) * (SPHtime1 - global_time) / (SPHtime1 - SPHtime0); 
  vph = vph1 - (vph1 - vph0) * (SPHtime1 - global_time) / (SPHtime1 - SPHtime0); 

  //verify again
  if(rho<rhoamb)
    return -1; //outside disk;
   
  
  ucon[1]=vr;
  ucon[2]=vth;
  ucon[3]=vph;
    
  conv_vels(ucon,ucon,VEL3,VELPRIM,geomBL->gg,geomBL->GG);

  //if(geom->iy==TNY/2 && geom->ix==TNX) printf("%d > %e %e %e %e %e\n",geom->iz,rho,temp,vr,vth,vph);
 
  pp[0]=rho;
  pp[1]=calc_PEQ_ufromTrho(temp,rho);
  ldouble uint=pp[1];
  pp[2]=ucon[1]; 
  pp[3]=ucon[2];
  pp[4]=ucon[3];

  //print_primitives(pp);getchar();

#ifdef MAGNFIELD//setting them zero not to break the following coordinate transformation
  pp[B1]=pp[B2]=pp[B3]=0.; 
#endif

#ifdef RADIATION
  //split pressure into rad and gas assuming LTE
  ldouble Eatm=pp[EE0];
  ldouble E,Fx, Fy,Fz;
  //distributing pressure
  ldouble P,aaa,bbb;
  P=GAMMAM1*uint;
  //solving for T satisfying P=pgas+prad=bbb T + aaa T^4
  aaa=4.*SIGMA_RAD/3.;
  bbb=K_BOLTZ*rho/MU_GAS/M_PROTON;
  ldouble naw1=cbrt(9*aaa*Power(bbb,2) - Sqrt(3)*Sqrt(27*Power(aaa,2)*Power(bbb,4) + 256*Power(aaa,3)*Power(P,3)));
  ldouble T4=-Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))/2. + Sqrt((4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 - naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa) + (2*bbb)/(aaa*Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))))/2.;
    
  E=calc_LTE_EfromT(T4);
  Fx=Fy=Fz=0.;
  uint=calc_PEQ_ufromTrho(T4,rho);

  pp[UU]=uint;//my_max(uint,uintatm);
  pp[EE0]=E;//my_max(E,Eatm);

  pp[FX0]=Fx;
  pp[FY0]=Fy;
  pp[FZ0]=Fz;

  //temp
  // if(geom->iy==NY/2) printf("%f %e\n",E/3. / (2./3.*uint),T4);
 
  //transforming from BL lab radiative primitives to code non-ortonormal primitives
  prad_ff2lab(pp,pp,geomBL);
    
#endif

  //transforming primitives from BL to MYCOORDS
  trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,geomBL->xxvec,gggBL,ggg);
    
 
    
#ifdef MAGNFIELD 
  //artificially impose poloidal magnetic field
  ldouble Pgas,Prad;
  Pgas=GAMMAM1*uint;
#ifdef RADIATION
  Prad=E/3.;
#else
  Prad=0.;
#endif
  
  //radial
  //pp[B1]=sqrt(2.*(Pgas+Prad)*MAGBETA)/sqrt(geom->gg[1][1])*sin(theq2/thmax*M_PI*2.)*cos(MAGNOMEGA*global_time);
  //if(global_time < VERTBTIME)
  //pp[B2]=sqrt(2.*(Pgas+Prad)*MAGBETA/2.)/sqrt(geom->gg[2][2])*sin(theq2/thmax*M_PI)*cos(MAGNOMEGA*global_time);
 
  pp[B1]=sqrt(2.*(Pgas+Prad)*MAGBETA)/sqrt(geom->gg[1][1])*cos(MAGNOMEGA*global_time);
  if(global_time < VERTBTIME)
    pp[B2]=sqrt(2.*(Pgas+Prad)*MAGBETA)/sqrt(geom->gg[2][2])*cos(MAGNOMEGA*global_time);
#endif

  
#ifdef PERTMAGN //perturb to break axisymmetry
pp[UU]*=1.+PERTMAGN*sin(10.*2.*M_PI*(MAXZ-geomBL->zz)/(MAXZ-MINZ));
#endif
 
//if(TI==NTX-1 && TJ==NTY/2 && TK==0 && iy==NY-1) printf("%d > %d > %e %e\n",TK,iz,pp[RHO],pp[B1]);
// if(TI==NTX-1 && TJ==NTY/2 && TK==NTZ-1 && iy==NY-1) printf("%d > %d > %e %e\n",TK,iz,pp[RHO],pp[B1]);

 pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU]);

  return 0;
}

