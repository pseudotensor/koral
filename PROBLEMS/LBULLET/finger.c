


int ix,iy,iz,ii,iv;
/**************************/

//#pragma omp parallel for private(ix,iy,iz,iv) schedule (dynamic)
for(ii=0;ii<Nloop_0;ii++) //domain only
  {
    ldouble pp[NV],uu[NV];
    struct geometry geom;
    ix=loop_0[ii][0];
    iy=loop_0[ii][1];
    iz=loop_0[ii][2]; 
 
    ldouble vu = get_u(pproblem1,VZ,ix,iy,iz); 
  
    if(fabs(vu)<1.e-50 && abs((iy+TOJ)-TNY/2)>TNY/4)
      {
    	fill_geometry(ix,iy,iz,&geom);
    	for(iv=0;iv<NV;iv++)
    	  {
    	    pp[iv]=get_u(pproblem1,iv,ix,iy,iz);	
    	  }

    	//hd floors
    	//check_floors_hd(pp,VELPRIM,&geom);
    	p2u(pp,uu,&geom);
     
        for(iv=0;iv<NV;iv++){
         set_u(p,iv,ix,iy,iz,pp[iv]);
         set_u(u,iv,ix,iy,iz,uu[iv]);
      }
  }
  }
