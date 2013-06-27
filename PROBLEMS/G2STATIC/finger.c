
//imposes initial condition where cloud is not

#ifdef TRACER
//skipping finger because I am not evolving cells (SKIP_CELLS)

int ix,iy,iz,ii,iv;
ldouble trace,rhoorg;
ldouble pp[NV],uu[NV];
/**************************/

#pragma omp parallel for private(ix,iy,iz,iv) schedule (dynamic)
for(ii=0;ii<Nloop_0;ii++) //domain only
  {
    ix=loop_0[ii][0];
    iy=loop_0[ii][1];
    iz=loop_0[ii][2]; 
 
    trace=get_u(p,TRA,ix,iy,iz);
    rhoorg=get_u(p,0,ix,iy,iz);

    if(trace<MINTRACE)
      {
	struct geometry geom;
	fill_geometry(ix,iy,iz,&geom);
	for(iv=0;iv<NV;iv++)
	  {
	    pp[iv]=get_u(pinit,iv,ix,iy,iz);	
	  }
	//hd floors
	check_floors_hd(pp,VELPRIM,&geom);
	p2u(pp,uu,&geom);
	for(iv=0;iv<NV;iv++)
	  {
	    set_u(u,iv,ix,iy,iz,uu[iv]);
	    set_u(p,iv,ix,iy,iz,pp[iv]);
	  }
     }
  }
#endif
