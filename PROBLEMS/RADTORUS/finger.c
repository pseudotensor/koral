
//imposes initial condition where cloud is not

#ifdef FIXEDALLBUTTEMP

int ix,iy,iz,ii,iv;
/**************************/

#pragma omp parallel for private(ix,iy,iz,iv) schedule (dynamic)
for(ii=0;ii<Nloop_0;ii++) //domain only
  {
    struct geometry geom;
    ix=loop_0[ii][0];
    iy=loop_0[ii][1];
    iz=loop_0[ii][2]; 
 
    fill_geometry(ix,iy,iz,&geom);

    ldouble temp = calc_PEQ_Tfromurho(get_u(p,UU,ix,iy,iz),get_u(p,RHO,ix,iy,iz));
	
    //resetting all gas quantities but for temperature
    for(iv=0;iv<NVMHD;iv++)
      set_u(p,iv,ix,iy,iz,get_u(pinit,iv,ix,iy,iz));
    
    set_u(p,UU,ix,iy,iz,calc_PEQ_ufromTrho(temp,get_u(p,RHO,ix,iy,iz)));

    p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
  }
#endif
