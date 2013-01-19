//KORAL - misc.c
//misceleanous routines

#include "ko.h"

//**********************************************************************
//**********************************************************************
//**********************************************************************
//integrates mass in the domain, assumes flat spacetime and uniform mesh
ldouble
calc_totalmass()
{
  int ix,iy,iz;
  ldouble dx,dy,dz,mass,rho;
  dx=(MAXX-MINX)/NX;
  dy=(MAXY-MINY)/NY;
  dz=(MAXZ-MINZ)/NZ;
  mass=0.;
  for(iz=0;iz<NZ;iz++)
    {
      for(iy=0;iy<NY;iy++)
	{
	  for(ix=0;ix<NX;ix++)
	    {
	      rho=get_u(u,0,ix,iy,iz);
	      mass+=rho*dx*dy*dz;

	    }
	}
    }
  return mass;
}
	  
//**********************************************************************
//**********************************************************************
//**********************************************************************
//calls gnuplot for 2d
int
convert_out2gif_2d(char *fname,char *fname2,int niter,ldouble t)
{
  FILE *fgnu=fopen("plot.gp","w");
  char bufor[50];

  //PROBLEMS/XXX/out2gid_2d.c
  #include PR_OUT2GIF_2D

  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  int i=system("gnuplot plot.gp ");
  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calls gnuplot for 1d
int
convert_out2gif_1d(char *fname,char *fname2,int niter,ldouble t)
{
  FILE *fgnu=fopen("plot.gp","w");
  char bufor[50];

  //PROBLEMS/XXX/out2gid_1d.c
  #include PR_OUT2GIF_1D

  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  int i=system("gnuplot plot.gp ");
  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
ldouble my_min(ldouble a, ldouble b)
{
  if(a<b) return a;
  else return b;
}

ldouble my_sign(ldouble x)
{
  if(x>0.) return 1.;
  if(x<0.) return -1.;
  if(x==0.) return 0.;
  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//initializes arrays
int
initialize_arrays()
{
  int i,j,k;

  
  //grid 
  x=(ldouble*)malloc((NX+NY+NZ+6*NG)*sizeof(ldouble));
  xb=(ldouble*)malloc((NX+1+NY+1+NZ+1+6*NG)*sizeof(ldouble));

  //primitives at cell centers
  p=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  //temporary
  pt0=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  //primitives at cell centers after reconstruction
  px=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  py=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  pz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));

  //conserved averages
  u=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  //source terms at cell centers
  s=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));

  //conserved at cell centers
  cellflag=(int*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NFLAGS*sizeof(int));
 
  //metric at cell centers
  g=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*gSIZE*sizeof(ldouble));
  //metric at cell centers
  G=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*gSIZE*sizeof(ldouble));
  //Kristofels at cell centers
  gKr=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*64*sizeof(ldouble));

  //LNRF basis one-forms
  emuup=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //LNRF basis vectors
  emulo=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //LNRF basis one-forms bx
  emuupbx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //LNRF basis vectors bx
  emulobx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //LNRF basis one-forms by
  emuupby=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*16*sizeof(ldouble));
  //LNRF basis vectors by
  emuloby=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*16*sizeof(ldouble));
  //LNRF basis one-forms bz
  emuupbz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*16*sizeof(ldouble));
  //LNRF basis vectors bz
  emulobz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*16*sizeof(ldouble));

  //ortonormal tetrad one-forms
  tmuup=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //ortonormal tetrad vectors
  tmulo=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //ortonormal tetrad one-forms bx
  tmuupbx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //ortonormal tetrad vectors bx
  tmulobx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG)*(NZ+2*NG)*16*sizeof(ldouble));
  //ortonormal tetrad one-forms by
  tmuupby=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*16*sizeof(ldouble));
  //ortonormal tetrad vectors by
  tmuloby=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*16*sizeof(ldouble));
  //ortonormal tetrad one-forms bz
  tmuupbz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*16*sizeof(ldouble));
  //ortonormal tetrad vectors bz
  tmulobz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*16*sizeof(ldouble));

  //left-interpolated primitives at cell x-faces
  pbLx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG )*(NZ+2*NG)*NV*sizeof(ldouble));
  //right-interpolated primitives at cell x-faces
  pbRx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  //left-interpolated primitives at cell y-faces
  pbLy=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*NV*sizeof(ldouble));
  //right-interpolated primitives at cell y-faces
  pbRy=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*NV*sizeof(ldouble));
  //left-interpolated primitives at cell z-faces
  pbLz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*NV*sizeof(ldouble));
  //right-interpolated primitives at cell z-faces
  pbRz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*NV*sizeof(ldouble));

   //left-interpolated primitives at cell x-faces
  sbLx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG )*(NZ+2*NG)*NV*sizeof(ldouble));
  //right-interpolated primitives at cell x-faces
  sbRx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  //left-interpolated primitives at cell y-faces
  sbLy=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*NV*sizeof(ldouble));
  //right-interpolated primitives at cell y-faces
  sbRy=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*NV*sizeof(ldouble));
  //left-interpolated primitives at cell z-faces
  sbLz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*NV*sizeof(ldouble));
  //right-interpolated primitives at cell z-faces
  sbRz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*NV*sizeof(ldouble));

 //left-interpolated conserved at cell x-faces - unused
  ubLx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG )*(NZ+2*NG)*NV*sizeof(ldouble));
  //right-interpolated conserved at cell x-faces
  ubRx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  //left-interpolated conserved at cell y-faces
  ubLy=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*NV*sizeof(ldouble));
  //right-interpolated conserved at cell y-faces
  ubRy=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*NV*sizeof(ldouble));
  //left-interpolated conserved at cell z-faces
  ubLz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*NV*sizeof(ldouble));
  //right-interpolated conserved at cell z-faces
  ubRz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*NV*sizeof(ldouble));

  //corrected flux at x faces
  flbx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG )*(NZ+2*NG)*NV*sizeof(ldouble));
  //corrected flux at x faces
  flby=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1 )*(NZ+2*NG)*NV*sizeof(ldouble));
  //corrected flux at x faces
  flbz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG )*(NZ+2*NG+1)*NV*sizeof(ldouble));
  //flux based on left-interpolated conserved at cell x-faces
  flLx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG )*(NZ+2*NG)*NV*sizeof(ldouble));
  //flux based on right-interpolated conserved at cell x-faces
  flRx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  //flux based on left-interpolated conserved at cell y-faces
  flLy=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*NV*sizeof(ldouble));
  //flux based on right-interpolated conserved at cell y-faces
  flRy=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*NV*sizeof(ldouble));
  //flux based on left-interpolated conserved at cell z-faces
  flLz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*NV*sizeof(ldouble));
  //flux based on right-interpolated conserved at cell z-faces
  flRz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*NV*sizeof(ldouble));


  //Krzysie at cell x-faces
  gKrbx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG )*(NZ+2*NG)*64*sizeof(ldouble));
   //Krzysie at cell x-faces
  gKrby=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1 )*(NZ+2*NG)*64*sizeof(ldouble));
 //Krzysie at cell x-faces
  gKrbz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG )*(NZ+2*NG+1)*64*sizeof(ldouble));

  //metric at cell x-faces
  gbx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG )*(NZ+2*NG)*gSIZE*sizeof(ldouble));
  //metric at cell y-faces
  gby=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*gSIZE*sizeof(ldouble));
  //metric at cell z-faces
  gbz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*gSIZE*sizeof(ldouble));
  //metric at cell x-faces
  Gbx=(ldouble*)malloc((NX+2*NG+1)*(NY+2*NG )*(NZ+2*NG)*gSIZE*sizeof(ldouble));
  //metric at cell y-faces
  Gby=(ldouble*)malloc((NX+2*NG)*(NY+2*NG+1)*(NZ+2*NG)*gSIZE*sizeof(ldouble));
  //metric at cell z-faces
  Gbz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG+1)*gSIZE*sizeof(ldouble));


  //indices of the ghost cells
  gcidx=(int**)malloc(3*sizeof(int*));
  gcidx[0]=(int*)malloc(2*NG*sizeof(int));
  gcidx[1]=(int*)malloc(2*NG*sizeof(int));
  gcidx[2]=(int*)malloc(2*NG*sizeof(int));

  //auxiliary primitive arrays
  ut0=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  ut1=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  ut2=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  ut3=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  ut4=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  du=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  u_bak=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  u_step1=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  u_step2=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  u_step3=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));
  u_step4=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*NV*sizeof(ldouble));

  //wavespeeds hd and rad - max(al,ar)
  ahdx=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  ahdy=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  ahdz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  aradx=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  arady=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  aradz=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));

  //wavespeeds hd and rad - leftgoing
  ahdxl=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  ahdyl=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  ahdzl=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  aradxl=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  aradyl=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  aradzl=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));

  //wavespeeds hd and rad - lrightgoing
  ahdxr=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  ahdyr=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  ahdzr=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  aradxr=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  aradyr=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));
  aradzr=(ldouble*)malloc((NX+2*NG)*(NY+2*NG)*(NZ+2*NG)*sizeof(ldouble));

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
int
free_arrays()
{
  int i1,i,j;
  free(cellflag);
  free(x);
  free(xb);
  free(p);
  free(pt0);
  free(px);
  free(py);
  free(pz);
  free(u);
  free(g);
  free(G);
  free(gKr);
  free(gKrbx);
  free(gKrby);
  free(gKrbz);
  free(emuup);
  free(emulo);
  free(emuupbx);
  free(emulobx);
  free(emuupby);
  free(emuloby);
  free(emuupbz);
  free(emulobz);
  free(tmuup);
  free(tmulo);
  free(tmuupbx);
  free(tmulobx);
  free(tmuupby);
  free(tmuloby);
  free(tmuupbz);
  free(tmulobz);
  free(pbLx);
  free(pbRx);
  free(pbLy);
  free(pbRy);
  free(pbLz);
  free(pbRz);
  free(sbLx);
  free(sbRx);
  free(sbLy);
  free(sbRy);
  free(sbLz);
  free(sbRz);
  free(ubLx);
  free(ubRx);
  free(ubLy);
  free(ubRy);
  free(ubLz);
  free(ubRz);
  free(flbx);
  free(flby);
  free(flbz);
  free(flLx);
  free(flRx);
  free(flLy);
  free(flRy);
  free(flLz);
  free(flRz);
  free(gbx);
  free(gby);
  free(gbz);
  free(Gbx);
  free(Gby);
  free(Gbz);
 
  free(gcidx[0]);
  free(gcidx[1]);
  free(gcidx[2]);
  free(gcidx);

  free(ut0);
  free(ut1);
  free(ut2);
  free(ut3);
  free(ut4);
  free(du);
  free(u_bak);
  free(u_step1);
  free(u_step2);
  free(u_step3);
  free(u_step4);
  free(aradx);
  free(arady);
  free(aradz);
  free(ahdx);
  free(ahdy);
  free(ahdz);
  free(aradxl);
  free(aradyl);
  free(aradzl);
  free(ahdxl);
  free(ahdyl);
  free(ahdzl);
  free(aradxr);
  free(aradyr);
  free(aradzr);
  free(ahdxr);
  free(ahdyr);
  free(ahdzr);
 
  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//inverse 4by4 matrix
int
inverse_44matrix(ldouble a[][4], ldouble ia[][4])
{

  ldouble mat[16],dst[16];
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      mat[i*4+j]=a[i][j];

  ldouble	tmp[12]; ldouble	src[16]; ldouble det;
  /* transpose matrix */
  for (i = 0; i <4; i++)
    {
      src[i]=mat[i*4];
      src[i+4]=mat[i*4+1];
      src[i+8]=mat[i*4+2];
      src[i+12]=mat[i*4+3];
    }
  /* calculate pairs for first 8 elements (cofactors) */
  tmp[0] = src[10] * src[15];
  tmp[1] = src[11] * src[14];
  tmp[2] = src[9] * src[15];
  tmp[3] = src[11] * src[13]; 
  tmp[4] = src[9] * src[14]; 
  tmp[5] = src[10] * src[13];
  tmp[6] = src[8] * src[15];
  tmp[7] = src[11] * src[12];
  tmp[8] = src[8] * src[14];
  tmp[9] = src[10] * src[12];
  tmp[10] = src[8] * src[13];
  tmp[11] = src[9] * src[12];
  /* calculate first 8 elements (cofactors) */
  dst[0] = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7]; 
  dst[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];
  dst[1] = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7]; 
  dst[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7]; 
  dst[2] = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];
  dst[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7]; 
  dst[3] = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6]; 
  dst[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6]; 
  dst[4] = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3]; 
  dst[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3]; 
  dst[5] = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3]; 
  dst[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];
  dst[6] = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3]; 
  dst[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];
  dst[7] = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];
  dst[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];
  /* calculate pairs for second 8 elements (cofactors) */
  tmp[0] = src[2]*src[7]; 
  tmp[1] = src[3]*src[6];
  tmp[2] = src[1]*src[7];
  tmp[3] = src[3]*src[5]; 
  tmp[4] = src[1]*src[6];
  tmp[5] = src[2]*src[5];
  tmp[6] = src[0]*src[7];
  tmp[7] = src[3]*src[4];
  tmp[8] = src[0]*src[6];
  tmp[9] = src[2]*src[4];
  tmp[10] = src[0]*src[5];
  tmp[11] = src[1]*src[4];
  /* calculate second 8 elements (cofactors) */
  dst[8] = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15]; 
  dst[8] -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];
  dst[9] = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15]; 
  dst[9] -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15]; 
  dst[10] = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];
  dst[10]-= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15]; 
  dst[11] = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];
  dst[11]-= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14]; 
  dst[12] = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];
  dst[12]-= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10]; 
  dst[13] = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10]; 
  dst[13]-= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8]; 
  dst[14] = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8]; 
  dst[14]-= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9]; 
  dst[15] = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9]; 
  dst[15]-= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8];
  /* calculate determinant */
  det=src[0]*dst[0]+src[1]*dst[1]+src[2]*dst[2]+src[3]*dst[3];

 
  /* calculate matrix inverse */
  det = 1/det; 

  if(isnan(det))
    my_err("det in inverse 4x4 zero\n");

  for (j = 0; j < 16; j++)
    dst[j] *= det;

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      ia[i][j]= dst[i*4+j];

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//prints error message and gets chars
int
my_err(char *message)
{
  char bufor[200];
  sprintf(bufor,"|err| : %s\n",message);
  printf("%s",bufor);
  getchar();
  return 0;  
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
int
getch()
{
  getchar();
  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//returns location of the horizone in BL
ldouble
r_horizon_BL(ldouble a)
{
  return 1.+sqrtl(1-a*a);
}

//returns location of the co-rotating marginally bound orbit in BL
ldouble
r_mbound_BL(ldouble a)
{
  return 2.*(1.-a/2.+sqrtl(1.-a));
}

//returns location of the photon orbit in BL
ldouble
r_photon_BL(ldouble a)
{
  return 2.*(1.-cosl(2./3.*acosl(-a)));
}
