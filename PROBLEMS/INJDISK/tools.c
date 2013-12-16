int diskatboundary_ktorus(ldouble *pp, void *ggg, void *gggBL)
{
  struct geometry *geom
     = (struct geometry *) ggg;
   
   struct geometry *geomBL
     = (struct geometry *) gggBL;
 
   //torus as in Kato+04
   ldouble Rsph=ROUT;
   ldouble Rcyl=ROUT*sin(geomBL->yy);
   ldouble ell0=sqrt(DISKRCIR*DISKRCIR*DISKRCIR)/(DISKRCIR-2.);
   ldouble ell = ell0*pow(Rcyl/DISKRCIR,ELLA);
   ldouble Psi0 = -1./(DISKRCIR-2.);
   ldouble PsiT0 = Psi0 + 1./(2.*(1.-ELLA))*pow(ell0/DISKRCIR,2.);
   ldouble Psi = -1./(Rsph-2.);
   ldouble PsiT = Psi + 1./(2.*(1.-ELLA))*pow(ell/Rcyl,2.);
   ldouble podpierd = 1. - GAMMA/(VSZERO*VSZERO)*(PsiT-PsiT0)/(NPOLI+1.);

   if(podpierd<0. || Rcyl<6.) 
     return -1; //outside disk;
   
   ldouble rho,mx,my,mz,m,E,uint,pgas,Fx,Fy,Fz,pLTE;  
   ldouble uu[NV], ppback[NV],T;
   ldouble Vphi,Vr;
   ldouble D,W,eps,uT,uphi,uPhi;
   
   rho = DISKRHO * pow(podpierd,NPOLI);

   pgas = DISKRHO * VSZERO * VSZERO / GAMMA * pow(rho/DISKRHO,1.+1./NPOLI);
   uint = pgas / GAMMAM1;
   Vphi=ell/Rcyl/Rcyl;

   //3-velocity in BL transformed to MYCOORDS
   ldouble ucon[4]={0.,DISKVR,0.,Vphi};
   conv_vels(ucon,ucon,VEL3,VELPRIM,geomBL->gg,geomBL->GG);
   
   ldouble rhoatm=pp[RHO];
   ldouble uintatm =pp[UU];
   pp[0]=my_max(rho,rhoatm); 
   pp[1]=my_max(uint,uintatm);
   pp[2]=ucon[1]; 
   pp[3]=ucon[2];
   pp[4]=ucon[3];

   //   print_primitives(pp);getchar();

#ifdef MAGNFIELD//setting them zero not to break the following coordinate transformation
    pp[B1]=pp[B2]=pp[B3]=0.; 
#endif

#ifdef RADIATION
    ldouble Eatm=pp[EE0];
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

    pp[UU]=my_max(uint,uintatm);
    pp[EE0]=my_max(E,Eatm);

    pp[FX0]=Fx;
    pp[FY0]=Fy;
    pp[FZ0]=Fz;

 
    //transforming from BL lab radiative primitives to code non-ortonormal primitives
    prad_ff2lab(pp,pp,geomBL);
    
#endif

    //transforming primitives from BL to MYCOORDS
    trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,geomBL->xxvec,gggBL,ggg);
    
#ifdef MAGNFIELD 
    //artificiall impose poloidal magnetic field
    ldouble Pgas,Prad;
    Pgas=GAMMAM1*uint;
#ifdef RADIATION
    Prad=E/3.;
#else
    Prad=0.;
#endif
    pp[B2]=sqrt((Pgas+Prad)*MAGBETA)/sqrt(geom->gg[2][2])*cos(MAGNOMEGA*global_time);   

#endif
  
    return 0;
}


int diskatboundary(ldouble *pp, void *ggg, void *gggBL)
{
  struct geometry *geom
     = (struct geometry *) ggg;
   
   struct geometry *geomBL
     = (struct geometry *) gggBL;
 
  ldouble th=geomBL->xxvec[2];
  ldouble r=geomBL->xxvec[1];
  ldouble theq=fabs(M_PI/2.-th);
  ldouble thmax=DISKHR*M_PI/2.;

  if(theq>thmax)
    return -1; //outside disk;
   
  ldouble rho,uint,uphi,Om,ucon[4],temp;

  rho=DISKRHO*pow(1.-pow(theq/thmax,2.),3.);
  temp=DISKTEMP*(1.-pow(theq/thmax,2.));
		   
  uint=calc_PEQ_ufromTrho(temp,rho);
  Om=sqrt(DISKRCIR)/r/r;

  ldouble rhoatm=pp[RHO];
  ldouble uintatm =pp[UU];

  ucon[1]=DISKVR;
  ucon[2]=0.;
  ucon[3]=Om;
    
  conv_vels_ut(ucon,ucon,VEL4,VELPRIM,geomBL->gg,geomBL->GG);
   
  pp[0]=my_max(rho,rhoatm);
  pp[1]=my_max(uint,uintatm);
  pp[2]=ucon[1]; 
  pp[3]=ucon[2];
  pp[4]=ucon[3];

  //print_primitives(pp);getchar();

#ifdef MAGNFIELD//setting them zero not to break the following coordinate transformation
  pp[B1]=pp[B2]=pp[B3]=0.; 
#endif

#ifdef RADIATION
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

  pp[UU]=my_max(uint,uintatm);
  pp[EE0]=my_max(E,Eatm);

  pp[FX0]=Fx;
  pp[FY0]=Fy;
  pp[FZ0]=Fz;

 
  //transforming from BL lab radiative primitives to code non-ortonormal primitives
  prad_ff2lab(pp,pp,geomBL);
    
#endif

  //transforming primitives from BL to MYCOORDS
  trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,geomBL->xxvec,gggBL,ggg);
    
 
    
#ifdef MAGNFIELD 
  //artificiall impose poloidal magnetic field
  ldouble Pgas,Prad;
  Pgas=GAMMAM1*uint;
#ifdef RADIATION
  Prad=E/3.;
#else
  Prad=0.;
#endif
  pp[B2]=sqrt((Pgas+Prad)*MAGBETA)/sqrt(geom->gg[2][2])*cos(MAGNOMEGA*global_time);
  
  /*
  //MYCOORDS vector potential to calculate B's
  ldouble Acov[4];
  Acov[0]=Acov[1]=Acov[2]=0.;
  Acov[3]=my_max(pow(pp[RHO]*geomBL.xx*geomBL.xx/4.e-20,2.)-0.02,0.)*sqrt(1.e-23)*pow(sin(fabs(geomBL.yy)),4.);

  pp[B1]=0.;
  pp[B2]=0.;
  pp[B3]=Acov[3];
  */
#endif
  
  return 0;
}
