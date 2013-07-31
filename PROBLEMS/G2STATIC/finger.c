
//imposes initial condition where cloud is not

#ifdef TRACER

int ix,iy,iz,ii,iv;
/**************************/

#pragma omp parallel for private(ix,iy,iz,iv) schedule (dynamic)
for(ii=0;ii<Nloop_0;ii++) //domain only
  {
    ldouble trace,f;
    ldouble pp[NV],ppdisk[NV],uu[NV];
    struct geometry geom;
    ix=loop_0[ii][0];
    iy=loop_0[ii][1];
    iz=loop_0[ii][2]; 
  
    trace=get_u(p,TRA,ix,iy,iz);
  
    if(trace<MINTRACE)
      {
	fill_geometry(ix,iy,iz,&geom);
	f = trace / MINTRACE;
	for(iv=0;iv<NV;iv++)
	  {
	    ppdisk[iv]=get_u(pproblem2,iv,ix,iy,iz);	
	    pp[iv]=get_u(p,iv,ix,iy,iz);

	    pp[iv]=ppdisk[iv]*(1.-f) + pp[iv]*f;
	  }

	//hd floors
	check_floors_hd(pp,VELPRIM,&geom);
	p2u(pp,uu,&geom);
	for(iv=0;iv<NV;iv++)
	  {
	    if(iv==TRA) continue; //do not overwrite the tracer field
	    set_u(u,iv,ix,iy,iz,uu[iv]);
	    set_u(p,iv,ix,iy,iz,pp[iv]);
	  }
     }
  }
#endif
