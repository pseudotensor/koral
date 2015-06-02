//definitions
ldouble rho,mx,my,mz,m,E,uint,pgas,Fx,Fy,Fz,pLTE,ell;  
ldouble uu[NV], pp[NV],ppback[NV],T;
ldouble Vphi,Vr;
ldouble D,W,eps,uT,uphi,uPhi;

//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

//coordinates
ldouble r=geomBL.xx;
ldouble th=geomBL.yy;

//ambient
set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
#ifdef RADIATION
set_radatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
#endif

if(get_u(pproblem1,0,ix,iy,iz)>pp[RHO]) //sph input if density higher than photosphere
  {
    ldouble rho,uint,vr,vth,vph,temp,ucon[4];

    rho=rhoCGS2GU(get_u(pproblem1,0,ix,iy,iz));
    temp=tempCGS2GU(get_u(pproblem1,1,ix,iy,iz));
    uint=calc_PEQ_ufromTrho(temp,rho);
    vr=tempCGS2GU(get_u(pproblem1,2,ix,iy,iz))/CCC0;
    vth=tempCGS2GU(get_u(pproblem1,3,ix,iy,iz))*GMC3;
    vph=tempCGS2GU(get_u(pproblem1,4,ix,iy,iz))*GMC3;
    //if(iy==TNY/2) printf("%d %d > %e %e %e\n",ix,iz,vr,vth,vph);  
    ucon[1]=vr;
    ucon[2]=vth;
    ucon[3]=vph;
    
    conv_vels(ucon,ucon,VEL3,VELPRIM,geomBL.gg,geomBL.GG);

    pp[0]=rho;
    pp[1]=uint;
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

  //temperature
  //if(geom.iy==TNY/2) printf("%d %d > %e %e %e %e %e\n",geom.iy,geom.iz,E, Eatm, E/3. / (2./3.*uint),temp,T4);
 
  //transforming from BL lab radiative primitives to code non-ortonormal primitives
  prad_ff2lab(pp,pp,geomBL);
    
#endif

  //transforming primitives from BL to MYCOORDS
  trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,geomBL.xxvec,&geomBL,&geom);
    
 
    
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
  //pp[B1]=sqrt(2.*(Pgas+Prad)*MAGBETA)/sqrt(geom.gg[1][1])*sin(theq2/thmax*M_PI*2.)*cos(MAGNOMEGA*global_time);
  //if(global_time < VERTBTIME)
  //pp[B2]=sqrt(2.*(Pgas+Prad)*MAGBETA/2.)/sqrt(geom.gg[2][2])*sin(theq2/thmax*M_PI)*cos(MAGNOMEGA*global_time);
 
  pp[B1]=sqrt(2.*(Pgas+Prad)*MAGBETA)/sqrt(geom.gg[1][1])*cos(MAGNOMEGA*global_time);
  if(global_time < VERTBTIME)
    pp[B2]=sqrt(2.*(Pgas+Prad)*MAGBETA)/sqrt(geom.gg[2][2])*cos(MAGNOMEGA*global_time);
#endif

  
#ifdef PERTMAGN //perturb to break axisymmetry
pp[UU]*=1.+PERTMAGN*sin(10.*2.*M_PI*(MAXZ-geomBL.zz)/(MAXZ-MINZ));
#endif
 
//if(TI==NTX-1 && TJ==NTY/2 && TK==0 && iy==NY-1) printf("%d > %d > %e %e\n",TK,iz,pp[RHO],pp[B1]);
// if(TI==NTX-1 && TJ==NTY/2 && TK==NTZ-1 && iy==NY-1) printf("%d > %d > %e %e\n",TK,iz,pp[RHO],pp[B1]);

 pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU]);


  }
 else //leave the atmosphere
   {
   
   }

//entropy
pp[5]=calc_Sfromu(pp[0],pp[1]);
//to conserved
p2u(pp,uu,&geom);

/***********************************************/

int iv;

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }

//entropy
update_entropy(ix,iy,iz,0);
set_cflag(0,ix,iy,iz,0);
