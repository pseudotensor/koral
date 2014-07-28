//sets the initial conditions
//called from a loop going over ix,iy,iz

/***********************************************/
//structure of geometry
//holds coordinates and metric, see ko.h
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
struct geometry geomcart;
fill_geometry_arb(ix,iy,iz,&geomcart,MINKCOORDS);


/***********************************************/
//vectors of primitives and conserved
ldouble pp[NV],uu[NV];

/***********************************************/
//angular momentum of the torus:
ldouble L=3.8;

ldouble W = (1./2.)*log(-(geom.gg[0][0]*geom.gg[3][3])/(geom.gg[3][3]+L*L*geom.gg[0][0]));
//ldouble rin=4;
ldouble Win = -0.0416192; 
ldouble w=exp(-(W-Win));

//OS: I didn't know how to manipulate the size of the torus in your formulae, so I used my old ones:
//ldouble podpierd=-(geom.GG[0][0]-2.*L*geom.GG[0][3]+L*L*geom.GG[3][3]);
//ldouble ut=-1./sqrt(podpierd);
//ut/=0.985; //determines the torus size, the closer to 1, the bigger the torus
//ldouble w=-1./ut;
ldouble epsilon = (w-1.)/GAMMA;  //OS: dot after 1.
ldouble vmichel = get_u(pproblem1,0,ix,iy,iz);
if(epsilon>0.) //OS: interior of the torus
  {
    ldouble kappa = 1.; //OS: entropy constant, 0.01 gave temperature < 1e5 what was a bit too low
    //density
    ldouble rho0=powl((GAMMA-1)*epsilon/kappa,1/(GAMMA-1)); //OS: without the if(w>1.) condition rho0 could be NaN
    pp[RHO]=rho0;

    //~pressure
    ldouble uu0 = kappa * pow(rho0, GAMMA) / (GAMMA - 1.); //OS: you forgot to define pressure which must be consistent with the torus model
    pp[UU]=uu0;

    //angular velocity
    ldouble omega = -L*(geom.gg[0][0]/geom.gg[3][3]);
    ldouble OMEGA1=sqrt(-omega*omega/(geom.gg[0][0]+omega*omega*geom.gg[3][3]));

    pp[VZ]=OMEGA1;
    pp[VY]=0.;
    pp[VX]=0.;

    //superimposing vmichel
    ldouble vmichel4vel[4]={0.,vmichel,0.,0.};
    trans2_coco(geom.xxvec,vmichel4vel,vmichel4vel,BLCOORDS,MYCOORDS);
    pp[VX]=vmichel4vel[1];
    pp[VY]=vmichel4vel[2];

    //just in case VELPRIM!=VEL4
    conv_velsinprims(pp,VEL4,VELPRIM,geom.gg,geom.GG);
  }
 else //OS: atmosphere outside the torus
   {
     pp[RHO]=RHOAMB;
     pp[UU]=UUAMB;
     pp[VY]=pp[VX]=pp[VZ]=0.;
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
