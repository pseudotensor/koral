//reads cloud from the file / creates it and asigns densities / velocities to pproblem[]

/***********************************************/
#pragma omp parallel for private(ix,iy,iz) schedule (dynamic)
for(iz=0;iz<NZ;iz++)
  {
    for(iy=0;iy<NY;iy++)
      {
	for(ix=0;ix<NX;ix++)
	  {
	    /***********************************************/
	    ldouble pp[NV];
	    struct geometry geom;
	    fill_geometry(ix,iy,iz,&geom);
	    struct geometry geomBL;
	    fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

	    //faked cloud
	    ldouble clix,cliy,cliz;
	    clix=NX*0.8;
	    cliy=NY/2;
	    cliz=NZ*.25;

	    ldouble clxx[4];
	    get_xx(clix,cliy,cliz,clxx);

	    //to cartesian
	    ldouble xxmink[4],clxxmink[4];
	    coco_N(geom.xxvec,xxmink,MYCOORDS,MINKCOORDS);
	    coco_N(clxx,clxxmink,MYCOORDS,MINKCOORDS);

	    //distance from the center
	    ldouble dist = sqrt((xxmink[1]-clxxmink[1])*(xxmink[1]-clxxmink[1])+
				(xxmink[2]-clxxmink[2])*(xxmink[2]-clxxmink[2])+
				(xxmink[3]-clxxmink[3])*(xxmink[3]-clxxmink[3]));


	    //increase rho
	    /*
	    ldouble mag=CLRHO;
	    ldouble factor=(1.+mag*exp(-dist*dist/CLWIDTH/CLWIDTH));
	    ldouble atmrho = pp[0];
	    ldouble clrho = (factor-1.)*atmrho;
	    */
	    pp[0] =CLRHO*exp(-dist*dist/CLWIDTH/CLWIDTH);

	    //velocity
	    ldouble OmKep = 1./sqrt(geomBL.xx*geomBL.xx*geomBL.xx);

	    //ldouble ucon[4]={0.,0,+OmKep,0.};
	    ldouble ucon[4]={0.,0,0,-OmKep};

	    conv_vels(ucon,ucon,VEL3,VEL4,geomBL.gg,geomBL.GG);
	    trans2_coco(geomBL.xxvec,ucon,ucon,KERRCOORDS,MYCOORDS);
	    conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);

	    pp[2]=ucon[1];
	    pp[3]=ucon[2];
	    pp[4]=ucon[3];

	    
	    int iv;

	    for(iv=0;iv<NV;iv++)
	      {
		set_u(pproblem,iv,ix,iy,iz,pp[iv]);
	      }

	    //uint & tr & entropy not used
	  }
      }
  }

