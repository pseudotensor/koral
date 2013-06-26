
//imposes initial condition where cloud is not

#ifdef TRACER
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
	

	//#include "init.c"

	/*
	//modified init.c

	ldouble rho,mx,my,mz,m,E,uint,Fx,Fy,Fz,pLTE;  
	ldouble xx,yy,zz;
	ldouble uu[NV],xxvec[4],xxvecBL[4];
	ldouble pp[NV],ppback[NV],T;

	struct geometry geom;
	fill_geometry(ix,iy,iz,&geom);
	struct geometry geomBL;
	fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

	//Sgr A* atmosphere
	set_sgradisk(pp,geom.xxvec,&geom,&geomBL);

	//TODO!!!!

	//increase rho
	ldouble mag=CLMAG;
	ldouble factor=(1.+mag*exp(-dist*dist/CLWIDTH/CLWIDTH));
	ldouble atmrho = pp[0];
	ldouble clrho = (factor-1.)*atmrho;
	pp[0] =atmrho+clrho;

#ifdef TRACER
	//tracer : fraction of cloud gas 
	pp[TRA]=clrho/pp[0];
#endif

	//velocity
	ldouble OmKep = 1./sqrt(geomBL.xx*geomBL.xx*geomBL.xx);

	//ldouble ucon[4]={0.,0,+OmKep,0.};
	ldouble ucon[4]={0.,0,0,-OmKep};

	conv_vels(ucon,ucon,VEL3,VEL4,geomBL.gg,geomBL.GG);
	trans2_coco(geomBL.xxvec,ucon,ucon,KERRCOORDS,MYCOORDS);
	conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);

	pp[2]=(pp[2]*atmrho + ucon[1]*clrho ) / (atmrho + clrho);
	pp[3]=(pp[3]*atmrho + ucon[2]*clrho ) / (atmrho + clrho);
	pp[4]=(pp[4]*atmrho + ucon[3]*clrho ) / (atmrho + clrho);    


	//entropy
	pp[5]=calc_Sfromu(pp[0],pp[1]);
	//hd floors
	check_floors_hd(pp,VELPRIM,&geom);
	//to conserved
	p2u(pp,uu,&geom);

	int iv;

	for(iv=0;iv<NV;iv++)
	  {
	    set_u(u,iv,ix,iy,iz,uu[iv]);
	    set_u(p,iv,ix,iy,iz,pp[iv]);
	  }

	//entropy
	update_entropy(ix,iy,iz,0);
	set_cflag(0,ix,iy,iz,0);
*/
      }
  }
#endif
