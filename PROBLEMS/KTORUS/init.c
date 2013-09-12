
ldouble rho,mx,my,mz,m,E,uint,pgas,Fx,Fy,Fz,pLTE;  
ldouble uu[NV], pp[NV],ppback[NV],T;
ldouble Vphi,Vr;
ldouble D,W,eps,uT,uphi,uPhi;

//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

//torus as in Kato+04
ldouble Rsph=geomBL.xx;
ldouble Rcyl=geomBL.xx*sin(geomBL.yy);
ldouble ell0=sqrt(RZERO*RZERO*RZERO)/(RZERO-2.);
ldouble ell = ell0*pow(Rcyl/RZERO,ELLA);
ldouble Psi0 = -1./(RZERO-2.);
ldouble PsiT0 = Psi0 + 1./(2.*(1.-ELLA))*pow(ell0/RZERO,2.);
ldouble Psi = -1./(Rsph-2.);
ldouble PsiT = Psi + 1./(2.*(1.-ELLA))*pow(ell/Rcyl,2.);
ldouble podpierd = 1. - GAMMA/(VSZERO*VSZERO)*(PsiT-PsiT0)/(NPOLI+1.);

//if(geom.iy==NY/2){printf("%e %e %e\n",geomBL.xx,geomBL.yy,(PsiT-PsiT0));getch();}

if(podpierd<0. || Rcyl<10. || 1)// outside donut
  {
    //ambient
    set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
    pp[VZ]=0.01/geom.xx/sqrt(geom.xx);
#ifdef RADIATION
    set_radatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
#endif
  }
 else //inside donut
   {
    //ambient
    set_hdatmosphere(ppback,geom.xxvec,geom.gg,geom.GG,0);
#ifdef RADIATION
    set_radatmosphere(ppback,geom.xxvec,geom.gg,geom.GG,0);
#endif

    rho = RHOZERO * pow(podpierd,NPOLI);

    pgas = RHOZERO * VSZERO * VSZERO / GAMMA * pow(rho/RHOZERO,1.+1./NPOLI);
    uint = pgas / GAMMAM1;
    Vphi=ell/Rcyl/Rcyl;

    //3-velocity in BL transformed to rel-velocity
    ldouble ucon[4]={0.,0.,0.,Vphi};
    conv_vels(ucon,ucon,VEL3,VELPRIM,geomBL.gg,geomBL.GG);
   
    pp[0]=my_max(rho,ppback[0]); 
    pp[1]=my_max(uint,ppback[1]);
    pp[2]=ucon[1]; 
    pp[3]=ucon[2];
    pp[4]=ucon[3];

#ifdef MAGNFIELD//setting them zero not to break the following coordinate transformation
    pp[B1]=pp[B2]=pp[B3]=0.; 
#endif

#ifdef RADIATION
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

    pp[UU]=my_max(uint,ppback[1]);
    pp[EE0]=my_max(E,ppback[EE0]);

    pp[FX0]=Fx;
    pp[FY0]=Fy;
    pp[FZ0]=Fz;

    //transforming from BL lab radiative primitives to code non-ortonormal primitives
    prad_ff2lab(pp,pp,&geomBL);

#endif

    //transforming primitives from BL to MYCOORDS
    trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,geomBL.xxvec,&geomBL,&geom);
    
#ifdef MAGNFIELD 
    //MYCOORDS vector potential to calculate B's
    ldouble Acov[4];
    Acov[0]=Acov[1]=Acov[2]=0.;
    Acov[3]=my_max(pow(pp[RHO]*geomBL.xx*geomBL.xx/4.e-20,2.)-0.02,0.)*sqrt(1.e-23)*pow(sin(fabs(geomBL.yy)),4.);

    pp[B1]=0.;
    pp[B2]=0.;
    pp[B3]=Acov[3];
#endif

  }

/*
#ifdef RADIATION
    //setting up radiative energy density to residual value - 
    //should immetiadelly increase sucking out the gas internal energy in opt.thick regions
    pp[EE0]=pp[UU]/1.e6;
    pp[FX0]=pp[FY0]=pp[FZ0]=0.;    
#endif
*/

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
