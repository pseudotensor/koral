//scales magnetic pressure so MAXBETA = pmag/(pgas+prad)
int ix,iy,iz;
ldouble fac;

#ifdef MAGNFIELD
if(PROCID==0) {printf("Renormalizing magnetic field... ");fflush(stdout);}


//manual normalization - verify beta!
#ifdef BETANORMFACTOR
fac=BETANORMFACTOR;
#endif

#pragma omp parallel for private(ix,iy,iz) schedule (dynamic)
for(iz=0;iz<NZ;iz++)
  {
    for(iy=0;iy<NY;iy++)
      {
	for(ix=0;ix<NX;ix++)
	  {


	    /***********************************************/
	    ldouble pp[NV],uu[NV];
	    struct geometry geom;
	    fill_geometry(ix,iy,iz,&geom);
	    struct geometry geomBL;
	    fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);
	    int iv;

	    //if(ix==0) {printf("%d %d > %e %e\n",ix,iy,get_u(p,B1,ix,iy,iz),2*cos(geomBL.yy)/geomBL.xx/geomBL.xx/geomBL.xx);getch();}

	    PLOOP(iv)
	      pp[iv]=get_u(p,iv,ix,iy,iz);
	    
	    pp[B1]*=fac;
	    pp[B2]*=fac;
	    pp[B3]*=fac;

	    p2u(pp,uu,&geom);

	    PLOOP(iv)
	    {
	      set_u(p,iv,ix,iy,iz,pp[iv]);
	      set_u(u,iv,ix,iy,iz,uu[iv]);
	    }
	  }
      }
  }


if(PROCID==0) {printf("done!\n");fflush(stdout);}
#endif
	    
