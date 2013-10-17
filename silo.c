//KORAL - silo.c
//routines for writing a silo file with quadratic mesh
//used both on the go and separately

#include "ko.h"
#include <silo.h>
#include <string.h>

/*********************************************/
/* writes silo file in dumps
/*********************************************/
int fprint_silofile(ldouble time, int num, char* folder, char* prefix)
{
  char bufor[50];
  sprintf(bufor,"%s/%s%04d.silo",folder,prefix,num);
 
  DBfile *file = NULL;/* The Silo file pointer */
  char *coordnames[3];/* Names of the coordinates */
  ldouble *nodex;/* The coordinate arrays */
  ldouble *nodey;
  ldouble *nodez;
  ldouble *coordinates[3];/* The array of coordinatearrays */
  int dimensions[3];/* The number of nodes */
   
  /* Create the Silo file */
  file = DBCreate(bufor, DB_CLOBBER, DB_LOCAL, NULL,DB_PDB);

  /* Name the coordinate axes ‘X’ and ‘Y’ */
  coordnames[0] = strdup("X");
  coordnames[1] = strdup("Y");
  coordnames[2] = strdup("Z");
  
  /* Give the cartesian coordinates of the mesh */
  int ix,iy,iz,iv;
  ldouble pp[NV],xxvec[4],xxveccar[4],xxvecsph[4];

  int nx=NX;
  int ny=NY;
  int nz=NZ;

#ifdef FULLPHI //printing one more cell in phi to close the sphere
  nz++;
#endif

  nodex=(ldouble *)malloc(nx*ny*nz*sizeof(ldouble));
  nodey=(ldouble *)malloc(nx*ny*nz*sizeof(ldouble));
  nodez=(ldouble *)malloc(nx*ny*nz*sizeof(ldouble));

  /* get the primitives */
  ldouble *rho = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *uint = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *temp = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  #ifdef TRACER
  ldouble *tracer = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  #endif 

  ldouble *vx = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *vy = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *vz = (ldouble*)malloc(nx*ny*nz*sizeof(double));

  #ifdef MAGNFIELD
  ldouble *bsq = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Bx = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *By = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Bz = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  #endif

  #ifdef RADIATION
  ldouble *Erad = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Fx = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Fy = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Fz = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  #endif

  for(iz=0;iz<nz;iz++)
    {
      for(iy=0;iy<ny;iy++)
	{
	  for(ix=0;ix<nx;ix++)
	    {
	      int iix,iiy,iiz;
	      iix=ix;
	      iiy=iy;
	      iiz=iz;
	      if(iiz>=NZ) iiz-=NZ;

	      struct geometry geom;
	      fill_geometry(iix,iiy,iiz,&geom);
	      struct geometry geomout;
	      fill_geometry_arb(iix,iiy,iiz,&geomout,OUTCOORDS);
	      
	      get_xx(iix,iiy,iiz,xxvec);
	      coco_N(xxvec,xxvecsph,MYCOORDS,SPHCOORDS);
	      coco_N(xxvec,xxveccar,MYCOORDS,MINKCOORDS);


	      ldouble r=xxvecsph[1];
	      ldouble th=xxvecsph[2];
	      ldouble ph=xxvecsph[3];

	      //if(OUTCOORDS==BLCOORDS && geomout.xx<r_horizon_BL(BHSPIN))
	      //continue;

	      int nodalindex=iz*(ny*nx) + iy*nx + ix;
	      for(iv=0;iv<NV;iv++)
		{
		  pp[iv]=get_u(p,iv,iix,iiy,iiz);
		}


	      //coordinates
	      nodex[nodalindex]=xxveccar[1];
	      nodey[nodalindex]=xxveccar[2];
	      nodez[nodalindex]=xxveccar[3];

	      //primitives to OUTCOORDS
              #ifdef RADIATION
	      trans_prad_coco(pp, pp, MYCOORDS,OUTCOORDS, xxvec,&geom,&geomout);
              #endif
	      trans_pmhd_coco(pp, pp, MYCOORDS,OUTCOORDS, xxvec,&geom,&geomout);

	      
	      //velocities
	      ldouble vel[4];

	      if(doingavg==0) //using snapshot date
		{
		  rho[nodalindex]=pp[RHO];
		  uint[nodalindex]=pp[UU];
		  vel[1]=pp[VX];
		  vel[2]=pp[VY];
		  vel[3]=pp[VZ];
		  conv_vels(vel,vel,VELPRIM,VEL4,geomout.gg,geomout.GG);						  
		}
	      else //using averaged data
		{
		  rho[nodalindex]=get_uavg(pavg,RHO,ix,iy,iz);
		  uint[nodalindex]=get_uavg(pavg,UU,ix,iy,iz);
		  vel[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  vel[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  vel[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);		  
		}

	      #ifdef CGSOUTPUT
	      rho[nodalindex]=rhoGU2CGS(rho[nodalindex]);
	      uint[nodalindex]=endenGU2CGS(uint[nodalindex]);
	      #endif

	      temp[nodalindex]=calc_PEQ_Tfromurho(uint[nodalindex],rho[nodalindex]);
	      
	      #ifdef TRACER
	      tracer[nodalindex]=pp[TRA];
	      #endif

	      //outvel - non-ortonormal VEL4
	      vx[nodalindex]=vel[1];
	      vy[nodalindex]=vel[2];
	      vz[nodalindex]=vel[3];

	      //transform to cartesian
	      if (MYCOORDS==SCHWCOORDS || MYCOORDS==KSCOORDS || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS || MYCOORDS==MKS1COORDS || MYCOORDS==MKS2COORDS)
		{
		  vel[2]*=r;
		  vel[3]*=r*sin(th);
		  
		  vx[nodalindex] = sin(th)*cos(ph)*vel[1] 
		    + cos(th)*cos(ph)*vel[2]
		    - sin(ph)*vel[3];

		  vy[nodalindex] = sin(th)*sin(ph)*vel[1] 
		    + cos(th)*sin(ph)*vel[2]
		    + cos(ph)*vel[3];

		  vz[nodalindex] = cos(th)*vel[1] 
		    - sin(th)*vel[2];
		}
		

	      #ifdef MAGNFIELD
	      //magnetic field
	      
	      ldouble bcon[4],bcov[4];

	      if(doingavg==0)
		{
		  calc_bcon_prim(pp,bcon,&geomout);
		  indices_21(bcon,bcov,geomout.gg); 
		  bsq[nodalindex] = dot(bcon,bcov);
		}
	      else
		{
		  bcon[1]=get_uavg(pavg,AVGBCON(1),ix,iy,iz);
		  bcon[2]=get_uavg(pavg,AVGBCON(2),ix,iy,iz);
		  bcon[3]=get_uavg(pavg,AVGBCON(3),ix,iy,iz);
		  bsq[nodalindex]=get_uavg(pavg,AVGBSQ,ix,iy,iz);
		}
	      	      
	      //to ortonormal	      
	      //trans2_cc2on(bcon,bcon,geomout.tup);

	      Bx[nodalindex]=bcon[1];
	      By[nodalindex]=bcon[2];
	      Bz[nodalindex]=bcon[3];

	      //transform to cartesian
	      if (MYCOORDS==SCHWCOORDS || MYCOORDS==KSCOORDS || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS || MYCOORDS==MKS1COORDS || MYCOORDS==MKS2COORDS)
		{
		  bcon[2]*=r;
		  bcon[3]*=r*sin(th);

		  Bx[nodalindex] = sin(th)*cos(ph)*bcon[1] 
		    + cos(th)*cos(ph)*bcon[2]
		    - sin(ph)*bcon[3];

		  By[nodalindex] = sin(th)*sin(ph)*bcon[1] 
		    + cos(th)*sin(ph)*bcon[2]
		    + cos(ph)*bcon[3];

		  Bz[nodalindex] = cos(th)*bcon[1] 
		    - sin(th)*bcon[2];
		}
	      #endif

	      #ifdef RADIATION

	      ldouble Rtt,Ehat,ugas[4],rvel[4],Rij[4][4];

	      if(doingavg==0)
		{
		  calc_ff_Rtt(pp,&Rtt,ugas,&geomout);
		  Ehat=-Rtt;  	      							  
		  //prad_lab2on(pp,pp,&geomout);
		  //rvel[1]=pp[FX0];
		  //rvel[2]=pp[FY0];
		  //rvel[3]=pp[FZ0];
		  //conv_vels(rvel,rvel,VELPRIM,VEL4,geomout.gg,geomout.GG);
		  calc_Rij(pp,&geomout,Rij); //calculates R^munu in OUTCOORDS
		}
	      else
		{
		  Ehat=get_uavg(pavg,AVGEHAT,ix,iy,iz);
		  Rij[0][0]=get_uavg(pavg,AVGRIJ(0,0),ix,iy,iz);
		  Rij[0][1]=get_uavg(pavg,AVGRIJ(0,1),ix,iy,iz);
		  Rij[0][2]=get_uavg(pavg,AVGRIJ(0,2),ix,iy,iz);
		  Rij[0][3]=get_uavg(pavg,AVGRIJ(0,3),ix,iy,iz);		  
		}
	      	
	      Erad[nodalindex]=Ehat;
	      Fx[nodalindex]=Rij[0][1];
	      Fy[nodalindex]=Rij[0][2];
	      Fz[nodalindex]=Rij[0][3];

	       //transform to cartesian
	      if (MYCOORDS==SCHWCOORDS || MYCOORDS==KSCOORDS || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS || MYCOORDS==MKS1COORDS || MYCOORDS==MKS2COORDS)
		{
		  Rij[0][2]*=r;
		  Rij[0][3]*=r*sin(th);

		  Fx[nodalindex] = sin(th)*cos(ph)*Rij[0][1] 
		    + cos(th)*cos(ph)*Rij[0][2]
		    - sin(ph)*Rij[0][3];

		  Fy[nodalindex] = sin(th)*sin(ph)*Rij[0][1] 
		    + cos(th)*sin(ph)*Rij[0][2]
		    + cos(ph)*Rij[0][3];

		  Fz[nodalindex] = cos(th)*Rij[0][1] 
		    - sin(th)*Rij[0][2];
		}
	      #endif
	  
	      /*
	      printf("%d %d %d | %e %e %e | %e %e %e | %e %e | %e %e %e\n",ix,iy,iz,
		     xxvec[1],xxvec[2],xxvec[3],
		     xxveccar[1],xxveccar[2],xxveccar[3],
		     rho[nodalindex],temp[nodalindex],
		     vx[nodalindex],vy[nodalindex],vz[nodalindex]);getchar();
	      */
	    }
	}
    }
  
  /* assign grid */
  int ndim;
  if(ny==1 && nz==1) //1d
    {
      ndim=1;
      dimensions[0]=nx;
      coordinates[0]=nodex;
    }
  else if(nz==1) //2d
    {
      ndim=2;

      dimensions[0] = nx;
      dimensions[1] = ny;

      coordinates[0] = nodex;
      coordinates[1] = nodey;
#ifdef SILO2D_XZPLANE
      coordinates[1] = nodez;
#endif
    }
  else if(ny==1) //2d, switch order
    {
      ndim=2;
      
      /* How many nodes in each direction? */
      dimensions[0] = nx;
      dimensions[1] = nz;

      /* Assign coordinates to coordinates array */
      coordinates[0] = nodex;
      coordinates[1] = nodez; 
    }
  else //3d
    {
      ndim=3;
      
      /* How many nodes in each direction? */
      dimensions[0] = nx;
      dimensions[1] = ny;
      dimensions[2] = nz;

      /* Assign coordinates to coordinates array */
      coordinates[0] = nodex;
      coordinates[1] = nodey;
      coordinates[2] = nodez;
    }      
     
  /* Write out the mesh to the file */
  DBPutQuadmesh(file, "mesh1", coordnames, coordinates,
  		dimensions, ndim, DB_DOUBLE, DB_NONCOLLINEAR, NULL);

  /* Write scalars */
  DBoptlist *optList = DBMakeOptlist(1);
  DBAddOption(optList, DBOPT_DTIME, (void*)&time);
  DBPutQuadvar1(file, "rho","mesh1", rho,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);

  DBPutQuadvar1(file, "uint","mesh1", uint,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);

  DBPutQuadvar1(file, "temp","mesh1", temp,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);

  #ifdef TRACER
  DBPutQuadvar1(file, "tracer","mesh1", tracer,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);
  #endif

  #ifdef RADIATION
  DBPutQuadvar1(file, "erad","mesh1", Erad,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);
  #endif

  #ifdef MAGNFIELD
  DBPutQuadvar1(file, "bsq","mesh1", bsq,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);
  #endif


  /* Write vectors */
  optList = DBMakeOptlist(1);
  DBAddOption(optList, DBOPT_DTIME, (void*)&time);
  char *names[3];  
  ldouble *vels[3];

  //velocity
  vels[0]=vx;
  vels[1]=vy;
  vels[2]=vz;
#ifdef SILO2D_XZPLANE
  vels[1]=vz;
#endif

  names[0] = strdup("vel1");
  names[1] = strdup("vel2");
  names[2] = strdup("vel3");
  DBPutQuadvar(file, "velocity","mesh1", 3, names, vels, 
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);

  #ifdef MAGNFIELD
  //magn field
  vels[0]=Bx;
  vels[1]=By;
  vels[2]=Bz;
#ifdef SILO2D_XZPLANE
  vels[1]=Bz;
#endif
  names[0] = strdup("B1");
  names[1] = strdup("B2");
  names[2] = strdup("B3");
  DBPutQuadvar(file, "magn_field","mesh1", 3, names, vels, 
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);
  #endif

  #ifdef RADIATION 
  //radiative flux
  vels[0]=Fx;
  vels[1]=Fy;
  vels[2]=Fz;
#ifdef SILO2D_XZPLANE
  vels[1]=Fz;
#endif
  names[0] = strdup("F1");
  names[1] = strdup("F2");
  names[2] = strdup("F3");
  DBPutQuadvar(file, "rad_flux","mesh1", 3, names, vels, 
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);
  #endif
 
  /* Close the Silo file */
  DBClose(file);

  free(nodex);
  free(nodey);
  free(nodez);

  free(rho);
  free(uint);
  free(temp);
  #ifdef TRACER
  free(tracer);
  #endif

  #ifdef RADIATION
  free(Erad);
  free(Fx);
  free(Fy);
  free(Fz);
  #endif

  #ifdef MAGNFIELD
  free(bsq);
  free(Bx);
  free(By);
  free(Bz);
  #endif
 
  
  free(vx);
  free(vy);
  free(vz);

  return (0);
}
