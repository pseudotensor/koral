//precalculate the initial dipolar magnetic field to use in the ghost cells - stored in pproblem1

#ifdef MAGNFIELD


{
  int ii
#pragma omp parallel for
 for(ii=0;ii<Nloop_02;ii++) //domain and ghost cells
	{
	  int ix, iy, iz;
	  ix=loop_02[ii][0];
	  iy=loop_02[ii][1];
	  iz=loop_02[ii][2]; 
	    //geomBLetries
	    struct geometry geomBL;
	    fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

	    ldouble r=geomBL.xx;
	    ldouble th=geomBL.yy;
	    ldouble ph=geomBL.zz;

	    //vector potential

	    set_u(pproblem1,B1,ix,iy,iz,0.);
	    set_u(pproblem1,B2,ix,iy,iz,0.);
#ifdef DIPOLE
	    set_u(pproblem1,B3,ix,iy,iz,sin(th)*sin(th)/r);
#endif
#ifdef MONOPOLE
	    set_u(pproblem1,B3,ix,iy,iz,1.-cos(th)); 
#endif
	  }
}


      if(PROCID==0) {printf("Calculating magn. field in prepinit.c ... (1) ");fflush(stdout);}
#pragma omp parallel
      calc_BfromA(pproblem1,0);

#pragma omp parallel
      {
	int ii;
 for(ii=0;ii<Nloop_02;ii++) //domain and ghost cells
	{
	  int ix, iy, iz;
	  ix=loop_02[ii][0];
	  iy=loop_02[ii][1];
	  iz=loop_02[ii][2]; 

	  ldouble fac=1.;
	  #ifdef BETANORMFACTOR
	  fac=BETANORMFACTOR;
#endif

	    int iv;
	    set_u(pproblem1,B1,ix,iy,iz,fac*get_u(pvecpot,1,ix,iy,iz));
	    set_u(pproblem1,B2,ix,iy,iz,fac*get_u(pvecpot,2,ix,iy,iz));
	    set_u(pproblem1,B3,ix,iy,iz,fac*get_u(pvecpot,3,ix,iy,iz));

	  }
      }
if(PROCID==0) {printf("Done.\n");fflush(stdout);}
#endif
