/*
//Torus Definition
ldouble L=3.8;
ldouble W = (1./2.)*log(-(geom.gg[0][0]*geom.gg[3][3])/(geom.gg[3][3]+L*L*geom.gg[0][0]));
ldouble Win = -0.0416192; 
ldouble w=exp(-(W-Win));
ldouble epsilon = (w-1.)/GAMMA;
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

    //just in case VELPRIM!=VEL4
    conv_velsinprims(pp,VEL4,VELPRIM,geom.gg,geom.GG);
  }
 else //OS: atmosphere outside the torus
   {
     pp[RHO]=RHOAMB;
     pp[UU]=UUAMB;
     pp[VY]=pp[VX]=pp[VZ]=0.;
   }
*/

ldouble Gamma_mo = GAMMA-1.;
//calculate quantities needed for michel solution
//We are fixing critical density and polytropic constant
ldouble vc2      = kappa*GAMMA*powl(rhocri,(Gamma_mo))/(1.+kappa*GAMMA/(Gamma_mo)*powl(rhocri,(Gamma_mo)));
ldouble uxcri    = sqrt(vc2/(1.+3.*vc2));
ldouble rcri     = 0.5/powl(uxcri,2);
ldouble c1       = rho0*uxcri*powl(rcri,2);
ldouble c3       = powl((1.+kp*GAMMA*powl(rhocri,(Gamma_mo))/(Gamma_mo)),2.)*(1.-2./rcri+powl(uxcri,2.));

//root finder
for(ii=0,ii<Nloop_0,ii++)
  {        
      int iv;
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2]; 

      struct geometry geom;
      fill_geometry(ix,iy,iz,&geom);
      struct geometry geomBL;
      fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

      double radius = geomBL.xx;

      if (radius < rcri)
      pp[VX] = uxcri + 1.
      else
      pp[VX] = powl(10,-5)
        
        ldouble hru = kappa*GAMMA/(Gamma_mo)*powl((c1/(pp[VX]*powl(radius,2))),(Gamma_mo));
        ldouble fu  = powl((1+hru),2)*(1.-2.*1./radius+powl(pp[VX],2.)) - c3;
        ldouble dfdu = 2.*(1.+hru)*((1-GAMMA)*hru*(1.-2./radius+powl(pp[VX],2))+powl(pp[VX],2)*(1+hru))/pp[VX];
        vmichel == fmax(powl(10,-5),pp[VX]-fu/dfdu);
    }
}
set_u(pproblem1,0,ix,iy,iz,vmichel);
