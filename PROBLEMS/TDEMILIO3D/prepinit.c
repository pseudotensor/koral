int gix,giy,giz;
int ix,iy,iz;
int kk;
//SPH data handling required everywhere
//reading in the coordinates and loking for the projection SPH -> KORAL
    //reading in the SPH coordinates in theta,phi
    FILE* coords=fopen("coordBL.dat","r");
    int cellcount=0,ret;
    ldouble sphth,sphph,th,ph;
    while(fscanf(coords,"%*d %d %d %*lf %lf %lf\n",&iy,&iz,&sphth,&sphph)!=EOF)
      {
	cellcount++;
	SPHcoordsY[iy]=sphth;
	SPHcoordsZ[iz]=sphph;
	//printf("%d %d %f %f \n",iy,iz,sphth,sphph);
      }
    fclose(coords);

    if(cellcount!=SPHTNY*SPHTNZ)
      my_err("coordBL.dat incompatible with SPH input\n");

    //now calculating the radial coordinates
     for(kk=0;kk<SPHTNX;kk++)
       {
	 ldouble r=exp(log(SPHRMIN) + (ldouble)kk/(ldouble)SPHTNX* (log(SPHRMAX) - log(SPHRMIN)));
	 //Exp[Log[RMIN] + (ir - 1)/NR (Log[RMAX] - Log[RMIN])];
	 SPHcoordsX[kk]=r;
       }

    struct geometry geomBL;
	  
    //searchnig for projection
    //for each SPH grid point find the closest point on the global KORAL grid
	    
    //in radius
    ldouble dist=BIG,distloc;
    int kbest=-1;
    for (kk=0;kk<SPHTNX;kk++)
      {
	ldouble rsph = SPHcoordsX[kk];
	dist=BIG;
	for(gix=0;gix<TNX;gix++)
	  {
	    int ix=gix-TOI;
 
	    ldouble xx[4]={0.,0.5*(calc_xb(ix,0)+calc_xb(ix+1,0)),0.,0.};
	    coco_N(xx,xx,MYCOORDS,BLCOORDS);
	    ldouble rkoral=xx[1];
	    distloc=fabs(rsph-rkoral);
	    if(distloc<dist) {kbest=gix;dist=distloc;}
	  }
	SPHprojX[kk]=kbest;
      }


    //in theta
    dist=BIG;
    kbest=-1;
    for (j=0;j<SPHTNY;j++)
      {
	ldouble thsph = SPHcoordsY[j];
	dist=BIG;
	for(giy=0;giy<TNY;giy++)
	  {
	     int iy=giy-TOJ;
	   
	    ldouble xx[4]={0.,0.,0.5*(calc_xb(iy,1)+calc_xb(iy+1,1)),0.};
	    coco_N(xx,xx,MYCOORDS,BLCOORDS);
	    ldouble thkoral=xx[2];
	    distloc=fabs(thsph-thkoral);
	    if(distloc<dist) {kbest=giy;dist=distloc;}
	  }
	SPHprojY[j]=kbest;
      }


    //in phi
    dist=BIG;
    kbest=-1;
    for (j=0;j<SPHTNZ;j++)
      {
	ldouble phsph = SPHcoordsZ[j];
	dist=BIG;
	for(giz=0;giz<TNZ;giz++)
	  {
	    int iz=giz-TOK;
	    ldouble xx[4]={0.,0.,0.,0.5*(calc_xb(iz,2)+calc_xb(iz+1,2))};
	    coco_N(xx,xx,MYCOORDS,BLCOORDS);
	    ldouble phkoral=xx[3];
	    distloc=fabs(phsph-phkoral);
	    if(distloc<dist) {kbest=giz;dist=distloc;}
	  }

	SPHprojZ[j]=kbest;
      }

if(PROCID==10)
  {
	for(kk=0;kk<SPHTNZ;kk++)
	  printf("%d> %d %e %d\n",PROCID,kk,SPHcoordsZ[kk],SPHprojZ[kk]);
  }


    //reading in the file
    char fname[100];
    sprintf(fname,"reddwarf3D");
    fhandle_problem1=fopen(fname,"rb");
    if(fhandle_problem1==NULL) {printf("no SPH2KORAL file %s!\n",fname); exit(1);}

    int it,ic;
    ldouble datain[5];

    if(PROCID==0) printf("reading in  SPH data from %s\n",fname);

    int count[NX][NY][NZ];
    
    for(i=0;i<NX;i++)
      for(j=0;j<NY;j++)
	for(kk=0;kk<NZ;kk++)
	  count[i][j][kk]=0;

    for(iz=0;iz<SPHTNZ;iz++)
      for(iy=0;iy<SPHTNY;iy++)
	for(ix=0;ix<SPHTNX;ix++)
	  {
	    i=SPHprojX[ix];
	    j=SPHprojY[iy];
	    kk=SPHprojZ[iz];

	    ret=fread(datain,sizeof(double),5,fhandle_problem1);
	  
	    if(if_indomain(i,j,kk))
	      {
		int sx,sy,sz;
		ldouble rho=datain[0];
		for(sx=-SPHSMEARX;sx<=SPHSMEARX;sx++)
		  for(sy=-SPHSMEARY;sy<=SPHSMEARY;sy++)
		    for(sz=-SPHSMEARZ;sz<=SPHSMEARZ;sz++)
		      {
			int iix,iiy,iiz;
			iix=i+sx;
			iiy=j+sy;
			iiz=kk+sz;
			if(if_indomain(iix,iiy,iiz))
			  {
			    set_u(pproblem1,0,i+sx,j+sy,kk+sz,rho+get_u(pproblem1,0,i+sx,j+sy,kk+sz));
			    int ivv;
			    for(ivv=1;ivv<5;ivv++)
			      set_u(pproblem1,ivv,i+sx,j+sy,kk+sz,rho*datain[ivv]+get_u(pproblem1,ivv,i+sx,j+sy,kk+sz));
			    count[i+sx][j+sy][kk+sz]++;
			  }
		      }

		//		printf("%d %d %d\n",ix,iy,iz);
	      }
	  }

    for(i=0;i<NX;i++)
      for(j=0;j<NY;j++)
	for(kk=0;kk<NZ;kk++)
	  {
	    if(count[i][j][kk]>0)
	      {
		ldouble rhotot=get_u(pproblem1,0,i,j,kk);
		set_u(pproblem1,0,i,j,kk,get_u(pproblem1,0,i,j,kk)/count[i][j][kk]);
		for(ix=1;ix<5;ix++)
		  set_u(pproblem1,ix,i,j,kk,get_u(pproblem1,ix,i,j,kk)/rhotot);
	      }

	    if(count[i][j][kk]<0 ||  get_u(pproblem1,RHO,i,j,kk)<SPHRHOCUT)
	      set_u(pproblem1,RHO,i,j,kk,-1.); //negative density if no sph data at this cell
	  }

    
    //printing at the equatorial plane
/*
    iy=NY/2;
    ix=50;
    for(iz=0;iz<NZ;iz++)
      {
	fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);
	  
	printf("%d %d %d > %f %f %f > %e\n",
	       ix, iy,iz,
	       geomBL.xx,
	       geomBL.yy,
	       geomBL.zz,
	       get_u(pproblem1,RHO,ix,iy,iz));
      }

exit(1);
*/

