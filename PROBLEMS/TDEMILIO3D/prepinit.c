int gix,giy,giz;
int ix,iy,iz;
int k;

    //reading in the SPH file
    char fname[100];
    sprintf(fname,"reddwarf3D");
    fhandle_problem1=fopen(fname,"rb");
    if(fhandle_problem1==NULL) {printf("no SPH2KORAL file %s!\n",fname); exit(1);}

    int it,ic;
    ldouble datain[8];

    if(PROCID==0) printf("reading in  SPH data from %s\n",fname);

    int count[NX][NY][NZ];
    
    for(i=0;i<NX;i++)
      for(j=0;j<NY;j++)
	for(k=0;k<NZ;k++)
	  count[i][j][k]=0;

while(fread(datain,sizeof(double),8,fhandle_problem1)>0)
  {
    ldouble r=datain[0]/GMC2;
    ldouble th=datain[1];
    ldouble ph=datain[2];
    int gi[3];

    find_globalindex(r,th,ph,gi);
    
    i=gi[0];
    j=gi[1];
    k=gi[2];
	    
    //printf("%e %e %e > %d %d %d\n",r,th,ph,i,j,k);
    //getch();
    //continue;

    if(if_indomain(i,j,k))
      {
	int sx,sy,sz;
	ldouble rho=datain[3];
	for(sx=-SPHSMEARX;sx<=SPHSMEARX;sx++)
	  for(sy=-SPHSMEARY;sy<=SPHSMEARY;sy++)
	    for(sz=-SPHSMEARZ;sz<=SPHSMEARZ;sz++)
	      {
		int iix,iiy,iiz;
		iix=i+sx;
		iiy=j+sy;
		iiz=k+sz;
		if(if_indomain(iix,iiy,iiz))
		  {
		    set_u(pproblem1,0,i+sx,j+sy,k+sz,rho+get_u(pproblem1,0,i+sx,j+sy,k+sz));
		    int ivv;
		    for(ivv=4;ivv<8;ivv++)
		      set_u(pproblem1,ivv-3,i+sx,j+sy,k+sz,rho*datain[ivv]+get_u(pproblem1,ivv-3,i+sx,j+sy,k+sz));
		    count[i+sx][j+sy][k+sz]++;
		  }
	      }
      }
  };
    

for(i=0;i<NX;i++)
  for(j=0;j<NY;j++)
    for(k=0;k<NZ;k++)
      {
	if(count[i][j][k]>0)
	  {
	    ldouble rhotot=get_u(pproblem1,0,i,j,k);
	    set_u(pproblem1,0,i,j,k,get_u(pproblem1,0,i,j,k)/count[i][j][k]);
	    for(ix=1;ix<5;ix++)
	      set_u(pproblem1,ix,i,j,k,get_u(pproblem1,ix,i,j,k)/rhotot);
	  }

	if(count[i][j][k]<0 ||  get_u(pproblem1,RHO,i,j,k)<SPHRHOCUT)
	  set_u(pproblem1,RHO,i,j,k,-1.); //negative density if no sph data at this cell
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

