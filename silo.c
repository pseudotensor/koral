//KORAL - silo.c
//routines for writing a silo file with quadratic mesh
//used both on the go and separately

#include "ko.h"
#include <silo.h>
#include <string.h>

/*********************************************/
/* writes silo file in dumps
/*********************************************/
int fprint_silofile(ldouble time, int num, char* folder)
{
  char bufor[50];
  sprintf(bufor,"%s/sil%04d.silo",folder,num);
 
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
  ldouble pp[NV],uu[NV],xxvec[4],xxveccar[4],xxvecsph[4];

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
  ldouble *temp = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  #ifdef TRACER
  ldouble *tracer = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  #endif 

  ldouble *vx = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *vy = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *vz = (ldouble*)malloc(nx*ny*nz*sizeof(double));

  #ifdef MAGNFIELD
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

	      int nodalindex=iz*(ny*nx) + iy*nx + ix;
	      for(iv=0;iv<NV;iv++)
		{
		  uu[iv]=get_u(u,iv,iix,iiy,iiz);
		  pp[iv]=get_u(p,iv,iix,iiy,iiz);
		}

	      get_xx(iix,iiy,iiz,xxvec);
	      coco_N(xxvec,xxveccar,MYCOORDS,MINKCOORDS);

	      //coordinates
	      nodex[nodalindex]=xxveccar[1];
	      nodey[nodalindex]=xxveccar[2];
	      nodez[nodalindex]=xxveccar[3];

	      //primitives to OUTCOORDS
              #ifdef RADIATION
	      trans_prad_coco(pp, pp, MYCOORDS,OUTCOORDS, xxvec,&geom,&geomout);
              #endif
	      trans_pmhd_coco(pp, pp, MYCOORDS,OUTCOORDS, xxvec,&geom,&geomout);

	      //scalars
	      #ifndef CGSOUTPUT
	      rho[nodalindex]=pp[RHO];
	      #else
	      rho[nodalindex]=rhoGU2CGS(pp[RHO]);
	      #endif
	      temp[nodalindex]=calc_PEQ_Tfromurho(pp[UU],pp[RHO]);
	      
	      #ifdef TRACER
	      tracer[nodalindex]=pp[TRA];
	      #endif

	      //velocities
	      ldouble vel[4]={0,pp[VX],pp[VY],pp[VZ]};	

	      conv_vels(vel,vel,VELPRIM,VEL4,geomout.gg,geomout.GG);						  
	      trans2_cc2on(vel,vel,geomout.tup);

	      //outvel - ortonormal VEL4
	      vx[nodalindex]=vel[1];
	      vy[nodalindex]=vel[2];
	      vz[nodalindex]=vel[3];

	      //transform to cartesian
	      if (MYCOORDS==SCHWCOORDS || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS || MYCOORDS==MKS1COORDS)
		{
		  coco_N(xxvec,xxvecsph,MYCOORDS,SPHCOORDS);
		  ldouble r=xxvecsph[1];
		  ldouble th=xxvecsph[2];
		  ldouble ph=xxvecsph[3];

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
	      ldouble bcon[4];
	      calc_bcon_prim(pp,bcon,&geomout);
	      
	      //to ortonormal	      
	      trans2_cc2on(bcon,bcon,geomout.tup);

	      Bx[nodalindex]=bcon[1];
	      By[nodalindex]=bcon[2];
	      Bz[nodalindex]=bcon[3];

	      //transform to cartesian
	      if (MYCOORDS==SCHWCOORDS || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS || MYCOORDS==MKS1COORDS)
		{
		  coco_N(xxvec,xxvecsph,MYCOORDS,SPHCOORDS);
		  ldouble r=xxvecsph[1];
		  ldouble th=xxvecsph[2];
		  ldouble ph=xxvecsph[3];

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
	      prad_lab2on(pp,pp,&geomout);

	      Erad[nodalindex]=pp[EE0];

	      ldouble rvel[4]={0,pp[FX0],pp[FY0],pp[FZ0]};	

	      conv_vels(rvel,rvel,VELPRIM,VEL4,geomout.gg,geomout.GG);						  
	      trans2_cc2on(rvel,rvel,geomout.tup);

	      //outvel - ortonormal VEL4
	      Fx[nodalindex]=rvel[1];
	      Fy[nodalindex]=rvel[2];
	      Fz[nodalindex]=rvel[3];

	       //transform to cartesian
	      if (MYCOORDS==SCHWCOORDS || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS || MYCOORDS==MKS1COORDS)
		{
		  coco_N(xxvec,xxvecsph,MYCOORDS,SPHCOORDS);
		  ldouble r=xxvecsph[1];
		  ldouble th=xxvecsph[2];
		  ldouble ph=xxvecsph[3];

		  Fx[nodalindex] = sin(th)*cos(ph)*rvel[1] 
		    + cos(th)*cos(ph)*rvel[2]
		    - sin(ph)*rvel[3];

		  Fy[nodalindex] = sin(th)*sin(ph)*rvel[1] 
		    + cos(th)*sin(ph)*rvel[2]
		    + cos(ph)*rvel[3];

		  Fz[nodalindex] = cos(th)*rvel[1] 
		    - sin(th)*rvel[2];
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

  DBPutQuadvar1(file, "temp","mesh1", temp,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);
  #ifdef TRACER
  DBPutQuadvar1(file, "tracer","mesh1", tracer,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_NODECENT, optList);
  #endif
  #ifdef RADIATION
  DBPutQuadvar1(file, "Erad","mesh1", Erad,
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
  free(temp);
  #ifdef TRACER
  free(tracer);
  #endif
  
  free(vx);
  free(vy);
  free(vz);

  return (0);
}
