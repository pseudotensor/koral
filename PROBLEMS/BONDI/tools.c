

/*
//precalculates hydro Bondi solution, saves it to 
int ix,iy,iz;


int
func (double r, const double y[], double f[],
      void *params)
{
  ldouble v=y[0];

  ldouble Be=(5.-3.*GAMMA)/(4.*(GAMMA-1.))/RBONDI;
  
  ldouble cs2=GAMMAM1*(Be+1./r-1./2.*v*v);

  ldouble dvdr=v/(v*v - cs2)*(2./r*cs2-1./r/r);

  //printf("%e %e %e %e\n",r,v,dvdr,cs2);getch();
  
  f[0] = dvdr;

  return GSL_SUCCESS;
}



#pragma omp parallel for private(ix,iy,iz) schedule (dynamic)
for(iz=0;iz<NZ;iz++)
  {
    for(iy=0;iy<NY;iy++)
      {
	for(ix=-NGCX;ix<NX+NGCX;ix++)
	  {
	    ldouble rho,E,uint,Tgas,Trad,r,prad,pgas,ut,vx,Be,Kappa,urs,cs2;

	    struct geometry geom;
	    fill_geometry(ix,iy,iz,&geom);

	    struct geometry geomBL;
	    fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

	    //initial point
	    //at RBONDI
	    r=RMAX;
	    Be=(5.-3.*GAMMA)/(4.*(GAMMA-1.))/RBONDI;
	    urs=-sqrt(1./2./RBONDI);

	    gsl_odeiv2_system sys = {func, NULL, 1, NULL};

	     int i;
	    double r0=RBONDI, r1 = geomBL.xx, ur;
	    if(r1>RBONDI)
	      ur=(1.-1.e-6)*urs;
	    else
	      ur=(1.+1.e-6)*urs;

	    gsl_odeiv2_driver * d = 
	      gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,
					     my_sign(r1-r0)*1e-6, 1e-6, 0.0);
	    int status = gsl_odeiv2_driver_apply (d, &r0, r1, &ur);

	    if (status != GSL_SUCCESS)
	      {
		printf ("error, return value=%d\n", status);
		break;
	      }
	  
	    gsl_odeiv2_driver_free (d);

	    //at given cell
	    r=geomBL.xx;
	    rho=rhoCGS2GU(-MDOT*calc_mdotEdd()/(4.*Pi*lenGU2CGS(r)*lenGU2CGS(r)*velGU2CGS(ur)));
	    cs2=GAMMAM1*(Be+1./r-1./2.*ur*ur);
	    pgas=rho/GAMMA*cs2;
	    uint=pgas/GAMMAM1;
	    prad=PRADGAS*pgas;
	    E=prad*3.;

	    set_u(pproblem1,RHO,ix,iy,iz,rho);
	    set_u(pproblem1,UU,ix,iy,iz,uint);
	    set_u(pproblem1,VX,ix,iy,iz,ur);
	    set_u(pproblem1,VY,ix,iy,iz,0.);
	    set_u(pproblem1,VZ,ix,iy,iz,0.);
#ifdef RADIATION
	    set_u(pproblem1,EE0,ix,iy,iz,E);
	    set_u(pproblem1,FX0,ix,iy,iz,0.);
	    set_u(pproblem1,FY0,ix,iy,iz,0.);
	    set_u(pproblem1,FZ0,ix,iy,iz,0.);
#endif


	  }
      }
  }

*/
