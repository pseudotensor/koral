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

	    //donut formulae
	    ldouble podpierd=-(geomBL.GG[0][0]-2.*ELL*geomBL.GG[0][3]+ELL*ELL*geomBL.GG[3][3]);
	    ldouble ut=-1./sqrt(podpierd);
	    ut/=UTPOT; //rescales rin
	    ut*=1.001; //not to account for suface

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
	    if(ut<-1 || podpierd<0.) //outside donut
	      {
		;
	      }
	    else
	      {
		ldouble Rtt,Ehat,ucon[4],prad;
		calc_ff_Rtt(pp,&Rtt,ucon,&geom);
		Ehat=-Rtt; 
		prad=Ehat/3.;
		ptot+=prad;
		//if(geom.ix==NX/2 && geom.iy==NY/2) printf("%e %e %e\n",pgas,prad,ptot);

#pragma omp critical
		if(pmag/ptot>maxbeta) maxbeta=pmag/ptot;
	      }
	    #endif
	  }
      }
  }

ldouble fac=sqrt(MAXBETA/maxbeta);

printf("rescaling magn.fields by %e\n",fac);

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
	    
