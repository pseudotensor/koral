//reads cloud from the file / creates it and asigns densities / velocities to pproblem[]

/***********************************************/
/***********************************************/
/***********************************************/
/***********************************************/
/***********************************************/
/***********************************************/
/** cloud from file smoothed with Gaussian *****/
/** kernel *************************************/
/***********************************************/
#ifdef G2CLOUD
/*
//uses p[0] for temporary storage of the 
//p[2..4] for storage of cartesian velocities
int ii;
#pragma omp parallel for private(ix,iy,iz) schedule (static)
  for(ii=0;ii<Nloop_0;ii++) //domain only
    {
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2]; 
      
      set_u(p,0,ix,iy,iz,0.);
      set_u(p,2,ix,iy,iz,0.);
      set_u(p,3,ix,iy,iz,0.);
      set_u(p,4,ix,iy,iz,0.);
    }
*/

//holds cloud particles in memory
ldouble **clpar = (ldouble**)malloc(sizeof(ldouble*));
clpar[0]=(ldouble*)malloc(6*sizeof(ldouble));
int npar=0;

//reads from file
FILE* cloud=fopen("cloud.dat","r");
ldouble x[4],xsph[4],v[4];
while(fscanf(cloud,"%lf %lf %lf %lf %lf %lf\n",&x[1],&x[2],&x[3],&v[1],&v[2],&v[3])!=EOF)
  {    
    //if particle too close - to get rid of alone particleas ahead of the cloud
    if(sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3])<RMINFORPART)
      continue;

    npar++;    

    if(npar>0)
      {
	clpar=(ldouble**)realloc(clpar,(npar)*sizeof(ldouble*));
	clpar[npar-1]=(ldouble*)malloc(6*sizeof(ldouble));
      }
    
   /*
    printf("%e %e %e\n",x[1],x[2],x[3]);
    coco_N(x,xsph,MINKCOORDS,SPHCOORDS);
    coco_N(x,x,MINKCOORDS,MYCOORDS);
    
    //searching for cell that encompasses this one
    if(x[1]<MINX || x[1]>MAXX) continue;
    if(x[2]<MINY || x[2]>MAXY) continue;
    if(x[3]<MINZ || x[3]>MAXZ) continue;
    
    ix=iy=iz=0;
    do
      ix++;
    while(get_xb(ix,0)<x[1] && iz<=NX);
    ix--;
    do
      iy++;
    while(get_xb(iy,1)<x[2] && iy<=NY);
    iy--;
    do
      iz++;
    while(get_xb(iz,2)<x[3] && iz<=NZ);
    iz--;
    
    printf("%e %e %e -> %d %d %d\n",x[1],x[2],x[3],ix,iy,iz);
    printf("[%e %e]\n",get_xb(ix,0),get_xb(ix+1,0));
    printf("[%e %e]\n",get_xb(iy,1),get_xb(iy+1,1));
    printf("[%e %e]\n",get_xb(iz,2),get_xb(iz+1,2));
    getchar();

    //converting velocitities to spherical
    ldouble ucon[4];
    ldouble r=xsph[1];
    ldouble th=xsph[2];
    ldouble ph=xsph[3];
    
    ucon[1]=cos(ph)*sin(th)*v[1]+sin(ph)*sin(th)*v[2]+cos(th)*v[3];
    ucon[2]=1./r*cos(ph)*cos(th)*v[1]+1./r*sin(ph)*cos(th)*v[2]-1./r*sin(th)*v[3];
    ucon[3]=-1./r*sin(ph)/sin(th)*v[1]+1./r*cos(ph)/sin(th)*v[2];

    //velocity to MYCOORDS velocity
    ldouble gg[4][5],GG[4][5],ggBL[4][5],GGBL[4][5];
    calc_g_arb(x,gg,MYCOORDS);
    calc_G_arb(x,GG,MYCOORDS);
    calc_g_arb(xsph,ggBL,BLCOORDS);
    calc_G_arb(xsph,GGBL,BLCOORDS);

    conv_vels(ucon,ucon,VEL3,VEL4,ggBL,GGBL);
    trans2(xsph,ucon,ucon,BLCOORDS,MYCOORDS);
    conv_vels(ucon,ucon,VEL4,VELPRIM,gg,GG);
    */

     
    //saving to cloud particle temporary table
    clpar[npar-1][0]=x[1];
    clpar[npar-1][1]=x[2];
    clpar[npar-1][2]=x[3];
    clpar[npar-1][3]=v[1];
    clpar[npar-1][4]=v[2];
    clpar[npar-1][5]=v[3];
  }

//cloud read, now going through cells and applying smoothing kernel

int ii;
#pragma omp parallel for private(ix,iy,iz) schedule (static)
  for(ii=0;ii<Nloop_0;ii++) //domain only
    {
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2]; 

      struct geometry geom;
      fill_geometry(ix,iy,iz,&geom);
      struct geometry geomBL;
      fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
      struct geometry geomXYZ;
      fill_geometry_arb(ix,iy,iz,&geomXYZ,MINKCOORDS);
      ldouble x,y,z;
      x=geomXYZ.xx;
      y=geomXYZ.yy;
      z=geomXYZ.zz;

      set_u(pproblem,0,ix,iy,iz,0.);
      set_u(pproblem,2,ix,iy,iz,0.);
      set_u(pproblem,3,ix,iy,iz,0.);
      set_u(pproblem,4,ix,iy,iz,0.);

      int icl;
      //over cloud particles
      for(icl=0;icl<npar;icl++)
	{
	  ldouble dist=
	    sqrt(pow(x-clpar[icl][0],2.)+pow(y-clpar[icl][1],2.)+pow(z-clpar[icl][2],2.));

	  ldouble kern=
	    exp(-dist*dist / KERNELWIDTH / KERNELWIDTH);

	  set_u(pproblem,0,ix,iy,iz,
		get_u(pproblem,0,ix,iy,iz)+kern); //assumes equal weights of particles

	  set_u(pproblem,2,ix,iy,iz,
		get_u(pproblem,2,ix,iy,iz)+kern*clpar[icl][3]); //x-momentum
	  
	  set_u(pproblem,3,ix,iy,iz,
		get_u(pproblem,3,ix,iy,iz)+kern*clpar[icl][4]); //y-momentum
	  
	  set_u(pproblem,4,ix,iy,iz,
		get_u(pproblem,4,ix,iy,iz)+kern*clpar[icl][5]); //z-momentum

	}
      //if(ix==13 && iy==10 && iz==21) getchar();
      //momentum -> velocity

      if(get_u(pproblem,0,ix,iy,iz)>SMALL)
	{
	  set_u(pproblem,2,ix,iy,iz,
		get_u(pproblem,2,ix,iy,iz)/get_u(pproblem,0,ix,iy,iz));
	  set_u(pproblem,3,ix,iy,iz,
		get_u(pproblem,3,ix,iy,iz)/get_u(pproblem,0,ix,iy,iz));
	  set_u(pproblem,4,ix,iy,iz,
		get_u(pproblem,4,ix,iy,iz)/get_u(pproblem,0,ix,iy,iz));
	}
      else
	{
	  set_u(pproblem,2,ix,iy,iz,0.);
	  set_u(pproblem,3,ix,iy,iz,0.);
	  set_u(pproblem,4,ix,iy,iz,0.);
	}
	  
      
      /*
      if(iy==NY/2 && ix==NX/2)
	{
	  printf("1 %d %d %d-> %e | %e %e %e\n",ix,iy,iz
		 ,get_u(pproblem,0,ix,iy,iz)
		 ,get_u(pproblem,2,ix,iy,iz)
		 ,get_u(pproblem,3,ix,iy,iz)
		 ,get_u(pproblem,4,ix,iy,iz));
	}
      */
       
      //converting velocitities to spherical
      ldouble v[4]={0.,get_u(pproblem,2,ix,iy,iz),get_u(pproblem,3,ix,iy,iz),get_u(pproblem,4,ix,iy,iz)};
      ldouble ucon[4];
      ldouble r=geomBL.xx;
      ldouble th=geomBL.yy;
      ldouble ph=geomBL.zz;

   
      ucon[1]=cos(ph)*sin(th)*v[1]+sin(ph)*sin(th)*v[2]+cos(th)*v[3];
      ucon[2]=1./r*cos(ph)*cos(th)*v[1]+1./r*sin(ph)*cos(th)*v[2]-1./r*sin(th)*v[3];
      ucon[3]=-1./r*sin(ph)/sin(th)*v[1]+1./r*cos(ph)/sin(th)*v[2];

      //velocity to MYCOORDS velocity
      conv_vels(ucon,ucon,VEL3,VEL4,geomBL.gg,geomBL.GG);
      trans2_coco(geomBL.xxvec,ucon,ucon,BLCOORDS,MYCOORDS);
      conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);

      //saving to pproblem
      set_u(pproblem,2,ix,iy,iz,ucon[1]);
      set_u(pproblem,3,ix,iy,iz,ucon[2]);
      set_u(pproblem,4,ix,iy,iz,ucon[3]);

      /*
      if(iy==NY/2 && ix==NX/2)
	{
	  printf("2 %d %d %d-> %e | %e %e %e\n",ix,iy,iz
		 ,get_u(pproblem,0,ix,iy,iz)
		 ,get_u(pproblem,2,ix,iy,iz)
		 ,get_u(pproblem,3,ix,iy,iz)
		 ,get_u(pproblem,4,ix,iy,iz));
	  getchar();
	}
      */
      
      //density normalized later
    }

  //calculating total mass in pproblem
  ldouble xx[4],dx[3],mass,rho,gdet;
  mass=0.;
  for(iz=0;iz<NZ;iz++)
    {
      for(iy=0;iy<NY;iy++)
	{
	  for(ix=0;ix<NX;ix++)
	    {
	      get_xx(ix,iy,iz,xx);
	      dx[0]=get_size_x(ix,0);
	      dx[1]=get_size_x(iy,1);
	      dx[2]=get_size_x(iz,2);
	      gdet=calc_gdet(xx);
	      rho=get_u(pproblem,0,ix,iy,iz);
	      //mass in cgs:
	      mass+=rhoGU2CGS(rho)*lenGU2CGS(lenGU2CGS(lenGU2CGS(dx[0]*dx[1]*dx[2]*gdet)));
	    }
	}
    }

//normalizing cloud
//printf("MASSCLIUD: %e %e\n",MASSCLOUD,mass); getch();
ldouble scale=MASSCLOUD/mass;
for(iz=0;iz<NZ;iz++)
  for(iy=0;iy<NY;iy++)
    for(ix=0;ix<NX;ix++)
      set_u(pproblem,0,ix,iy,iz,get_u(pproblem,0,ix,iy,iz)*scale);
	    
#endif

/***********************************************/
/***********************************************/
/***********************************************/
/***********************************************/
/***********************************************/
#ifdef SPHCLOUD
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
#endif
