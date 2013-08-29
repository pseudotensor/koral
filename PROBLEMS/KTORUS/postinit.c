//scales magnetic pressure so MAXBETA = pmag/(pgas+prad)
int ix,iy,iz;

#ifdef MAGNFIELD
ldouble maxbeta=0.;
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
	    fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
 
	    int iv;
	    PLOOP(iv)
	      pp[iv]=get_u(p,iv,ix,iy,iz);

	    ldouble bcon[4],bcov[4],bsq,pmag;
	    calc_bcon_prim(pp,bcon,&geom);
	    indices_21(bcon,bcov,geom.gg); 
	    bsq = dot(bcon,bcov);
	    pmag = bsq/2.;
	    
	    ldouble pgas=GAMMAM1*pp[UU];
	    ldouble ptot=pgas;

#ifdef RADIATION

		ldouble Rtt,Ehat,ucon[4],prad;
		calc_ff_Rtt(pp,&Rtt,ucon,&geom);
		Ehat=-Rtt; 
		prad=Ehat/3.;
		ptot+=prad;
		
#endif

	    if(geom.iy==NY/2)
	      {
#pragma omp critical
		if(pmag/ptot>maxbeta) maxbeta=pmag/ptot;
	      }

	  }
      }
  }

ldouble fac=sqrt(MAXBETA/maxbeta);
printf("rescaling magn.fields by %e (%e)\n",fac,maxbeta);

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
	    int iv;

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

/*
printf("adjusting gas pressure...\n",fac);

#pragma omp parallel for private(ix,iy,iz) schedule (dynamic)
for(iz=0;iz<NZ;iz++)
  {
    for(iy=0;iy<NY;iy++)
      {
	for(ix=0;ix<NX;ix++)
	  {

	    ldouble pp[NV],uu[NV];
	    struct geometry geom;
	    fill_geometry(ix,iy,iz,&geom);
	    int iv;

	    PLOOP(iv)
	      pp[iv]=get_u(p,iv,ix,iy,iz);
	    
	    ldouble pgas=GAMMAM1*pp[UU];

	    ldouble bcon[4],bcov[4],bsq,pmag;
	    calc_bcon_prim(pp,bcon,&geom);
	    indices_21(bcon,bcov,geom.gg); 
	    bsq = dot(bcon,bcov);
	    pmag = bsq/2.;

	    pgas-=pmag;

	    pp[UU]=pgas;

	    p2u(pp,uu,&geom);

	    PLOOP(iv)
	    {
	      set_u(p,iv,ix,iy,iz,pp[iv]);
	      set_u(u,iv,ix,iy,iz,uu[iv]);
	    }

	  }
      }
  }
*/
#endif
	    
