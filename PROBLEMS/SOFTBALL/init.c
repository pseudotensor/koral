//sets the initial conditions
//called from a loop going over ix,iy,iz

/***********************************************/
//structure of geometry
//holds coordinates and metric, see ko.h
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
struct geometry geomcart;
fill_geometry_arb(ix,iy,iz,&geomcart,MINKCOORDS);
struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);


/***********************************************/
//vectors of primitives and conserved
ldouble pp[NV],uu[NV];

/***********************************************/
//angular momentum of the torus:
ldouble L=3.8;

ldouble W = (1./2.)*log(-(geom.gg[0][0]*geom.gg[3][3])/(geom.gg[3][3]+L*L*geom.gg[0][0]));
ldouble rin=4.57599;
ldouble Win = -0.0416192; 
ldouble w=exp(-(W-Win));

//OS: I didn't know how to manipulate the size of the torus in your formulae, so I used my old ones:
//ldouble podpierd=-(geom.GG[0][0]-2.*L*geom.GG[0][3]+L*L*geom.GG[3][3]);
//ldouble ut=-1./sqrt(podpierd);
//ut/=0.985; //determines the torus size, the closer to 1, the bigger the torus
//ldouble w=-1./ut;
ldouble epsilon = (w-1.)/GAMMA;  //OS: dot after 1.
ldouble vmichel = get_u(pproblem1,0,ix,iy,iz);

//when constructing rad-pressure supported torus we want to have pressure like for a gamma=4/3 gas, because radiation pressure has effective gamma = 4/3 
ldouble effgamma=GAMMA;
#ifdef RADIATION
effgamma=4./3.;
#endif

if(epsilon>0. && geomBL.xx>rin) //OS: interior of the torus
  {



    ldouble kappa = TORUSENTR; //OS: entropy constant, 0.01 gave temperature < 1e5 what was a bit too low
    //density
    ldouble rho0=powl((effgamma-1)*epsilon/kappa,1/(effgamma-1)); //OS: without the if(w>1.) condition rho0 could be NaN
    pp[RHO]=rho0;

    //~pressure
    ldouble uu0 = kappa * pow(rho0, effgamma) / (effgamma - 1.); //OS: you forgot to define pressure which must be consistent with the torus model
    pp[UU]=uu0;

    //angular velocity
    ldouble omega = -L*(geom.gg[0][0]/geom.gg[3][3]);
    ldouble OMEGA1=sqrt(-omega*omega/(geom.gg[0][0]+omega*omega*geom.gg[3][3]));

    pp[VZ]=OMEGA1;
    pp[VY]=0.;
    pp[VX]=0.;

    
    //superimposing vmichel
    #ifdef UNPERTURBED
    vmichel=0.;
    #endif

    ldouble vmichel4vel[4]={0.,0.06*vmichel,0.,0.};
    trans2_coco(geom.xxvec,vmichel4vel,vmichel4vel,BLCOORDS,MYCOORDS);
    pp[VX]=vmichel4vel[1];
    pp[VY]=vmichel4vel[2];
    
    //just in case VELPRIM!=VEL4
    conv_velsinprims(pp,VEL4,VELPRIM,geom.gg,geom.GG);

    
#ifdef RADIATION
    //taking original torus pressure and distributing it between gas and radiation to satisfy local thermal equilibrium:
    //solving for T satisfying P=pgas+prad=bbb T + aaa T^4

    ldouble uint=pp[UU];
    ldouble rho=pp[RHO];

    ldouble P,aaa,bbb;
    P=(effgamma-1.)*uint;
    aaa=4.*SIGMA_RAD/3.;
    bbb=K_BOLTZ*rho/MU_GAS/M_PROTON;
    ldouble naw1=cbrt(9*aaa*Power(bbb,2) - Sqrt(3)*Sqrt(27*Power(aaa,2)*Power(bbb,4) + 256*Power(aaa,3)*Power(P,3)));

    //effective temperature of radiation
    ldouble T4=-Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))/2. + Sqrt((4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 - naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa) + (2*bbb)/(aaa*Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))))/2.;

    //radiation energy density
    ldouble E=calc_LTE_EfromT(T4); 
    
    //fluid frame fluxes
    ldouble Fx,Fy,Fz;
    Fx=Fy=Fz=0.;

    //corrected gas pressure / internal energy density
    pp[UU]=calc_PEQ_ufromTrho(T4,rho);

    //writing to primitives
    pp[EE0]=E;
    pp[FX0]=Fx;
    pp[FY0]=Fy;
    pp[FZ0]=Fz;

    //transforming from BL lab radiative primitives to code non-ortonormal primitives
    prad_ff2lab(pp,pp,&geomBL);
#endif

  }
 else //OS: atmosphere outside the torus
   {
     pp[RHO]=RHOAMB;
     pp[UU]=UUAMB;
     pp[VY]=pp[VX]=pp[VZ]=0.;

#ifdef RADIATION
    set_radatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
#endif
   }

/***********************************************/
//calculate entropy from rho & uint
pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU]);

/***********************************************/
//convert primitives to conserved
p2u(pp,uu,&geom);

/***********************************************/
//save to memory
int iv;
for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }
