//KORAL - misc.c
//misceleanous routines

#include "ko.h"
	  
//**********************************************************************
//**********************************************************************
//**********************************************************************
//calls gnuplot for 2d
int
convert_out2gif_2d(char *fname,char *fname2,int niter,ldouble t)
{
  char buf1[50],buf2[50];  
  sprintf(buf1,"plot.gp.%d",PROCID);
  FILE *fgnu=fopen(buf1,"w");


  //PROBLEMS/XXX/out2gid_2d.c
  #include PR_OUT2GIF_2D

  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  sprintf(buf2,"gnuplot %s",buf1);
  int i=system(buf2);
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

//**********************************************************************
//**********************************************************************
//**********************************************************************
ldouble my_min_N(ldouble *v,int N)
{
  ldouble min=v[0];
  int i;
  for(i=1;i<N;i++)
    if(v[i]<min) min=v[i];
  return min;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
ldouble my_max_N(ldouble *v,int N)
{
  ldouble max=v[0];
  int i;
  for(i=1;i<N;i++)
    if(v[i]>max) max=v[i];
  return max;
}

ldouble my_sign(ldouble x)
{
  if(x>0.) return 1.;
  if(x<0.) return -1.;
  if(x==0.) return 1.;
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
  
  /********* basic arrays, used for postprocessing ***********/

  //grid 
  x=(ldouble*)malloc((NX+NY+NZ+6*NG)*sizeof(ldouble));
  xb=(ldouble*)malloc((NX+1+NY+1+NZ+1+6*NG)*sizeof(ldouble));

  //primitives at cell centers
  p=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));

  //quantities to average in time
  pavg=(ldouble*)malloc((SX)*(SY)*(SZ)*(NV+NAVGVARS)*sizeof(ldouble));
  
  //conserved averages
  u=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));

  //flags at cell centers
  cellflag=(int*)malloc((SX)*(SY)*(SZ)*NFLAGS*sizeof(int));
 

  //metric at cell centers
  g=(ldouble*)malloc((SX)*(SY)*(SZMET)*gSIZE*sizeof(ldouble));
  //metric at cell centers
  G=(ldouble*)malloc((SX)*(SY)*(SZMET)*gSIZE*sizeof(ldouble));
  //Kristofels at cell centers
  gKr=(ldouble*)malloc((SX)*(SY)*(SZMET)*64*sizeof(ldouble));

  /********* extra arrays, used only for time evolution ***********/

  if(doingpostproc==0)
    {
     
      //primitives at cell centers at initial state - used for fixed boundary conditions
      pinit=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
 
      //primitives at cell centers for other uses 
      pproblem1=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
      pproblem2=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));

      //arrays for temporary use (e.g., vector potential, mimic_dynamo)
      ptemp1=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
      pvecpot=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));

      //primitives at cell centers in previous time steps
      //ptm1=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
      //ptm2=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));

      //buffer for sending/receiving messages
#ifdef MPI
      msgbufs=(ldouble**)malloc(MPIMSGBUFSIZE*sizeof(ldouble*));
      for(i=0;i<MPIMSGBUFSIZE;i++)
	msgbufs[i]=(ldouble*)malloc(my_max3(NX*NY*NV*NG,NY*NZ*NV*NG,NZ*NX*NV*NG)*sizeof(ldouble));
#endif

#ifdef MAGNFIELD
      //electromotive force at corners
      emf=(ldouble*)malloc((NX+1)*(NY+1)*(NZ+1)*3*sizeof(ldouble));
#endif

      //metric at cell x-faces
      gbx=(ldouble*)malloc((SX+1)*(SY)*(SZMET)*gSIZE*sizeof(ldouble));
      //metric at cell y-faces
      gby=(ldouble*)malloc((SX)*(SY+1)*(SZMET)*gSIZE*sizeof(ldouble));
      //metric at cell z-faces
      gbz=(ldouble*)malloc((SX)*(SY)*(SZMET+1)*gSIZE*sizeof(ldouble));
      //metric at cell x-faces
      Gbx=(ldouble*)malloc((SX+1)*(SY)*(SZMET)*gSIZE*sizeof(ldouble));
      //metric at cell y-faces
      Gby=(ldouble*)malloc((SX)*(SY+1)*(SZMET)*gSIZE*sizeof(ldouble));
      //metric at cell z-faces
      Gbz=(ldouble*)malloc((SX)*(SY)*(SZMET+1)*gSIZE*sizeof(ldouble));

      
      //LNRF basis one-forms
      emuup=(ldouble*)malloc((SX)*(SY)*(SZMET)*16*sizeof(ldouble));
      //LNRF basis vectors
      emulo=(ldouble*)malloc((SX)*(SY)*(SZMET)*16*sizeof(ldouble));

      //ortonormal tetrad one-forms
      tmuup=(ldouble*)malloc((SX)*(SY)*(SZMET)*16*sizeof(ldouble));
      //ortonormal tetrad vectors
      tmulo=(ldouble*)malloc((SX)*(SY)*(SZMET)*16*sizeof(ldouble));

      //left-interpolated primitives at cell x-faces
      pbLx=(ldouble*)malloc((SX+1)*(SY )*(SZ)*NV*sizeof(ldouble));
      //right-interpolated primitives at cell x-faces
      pbRx=(ldouble*)malloc((SX+1)*(SY)*(SZ)*NV*sizeof(ldouble));
      //left-interpolated primitives at cell y-faces
      pbLy=(ldouble*)malloc((SX)*(SY+1)*(SZ)*NV*sizeof(ldouble));
      //right-interpolated primitives at cell y-faces
      pbRy=(ldouble*)malloc((SX)*(SY+1)*(SZ)*NV*sizeof(ldouble));
      //left-interpolated primitives at cell z-faces
      pbLz=(ldouble*)malloc((SX)*(SY)*(SZ+1)*NV*sizeof(ldouble));
      //right-interpolated primitives at cell z-faces
      pbRz=(ldouble*)malloc((SX)*(SY)*(SZ+1)*NV*sizeof(ldouble));


      //corrected flux at x faces
      flbx=(ldouble*)malloc((SX+1)*(SY)*(SZ)*NV*sizeof(ldouble));
      //corrected flux at x faces
      flby=(ldouble*)malloc((SX)*(SY+1)*(SZ)*NV*sizeof(ldouble));
      //corrected flux at x faces
      flbz=(ldouble*)malloc((SX)*(SY)*(SZ+1)*NV*sizeof(ldouble));

      //flux based on left-interpolated conserved at cell x-faces
      flLx=(ldouble*)malloc((SX+1)*(SY )*(SZ)*NV*sizeof(ldouble));
      //flux based on right-interpolated conserved at cell x-faces
      flRx=(ldouble*)malloc((SX+1)*(SY)*(SZ)*NV*sizeof(ldouble));
      //flux based on left-interpolated conserved at cell y-faces
      flLy=(ldouble*)malloc((SX)*(SY+1)*(SZ)*NV*sizeof(ldouble));
      //flux based on right-interpolated conserved at cell y-faces
      flRy=(ldouble*)malloc((SX)*(SY+1)*(SZ)*NV*sizeof(ldouble));
      //flux based on left-interpolated conserved at cell z-faces
      flLz=(ldouble*)malloc((SX)*(SY)*(SZ+1)*NV*sizeof(ldouble));
      //flux based on right-interpolated conserved at cell z-faces
      flRz=(ldouble*)malloc((SX)*(SY)*(SZ+1)*NV*sizeof(ldouble));




      //auxiliary primitive arrays
      ut0=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
      ut1=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
      ut2=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
      ut3=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
      dut0=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
      dut1=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
      dut2=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
      dut3=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
      drt0=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
      drt1=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
      drt2=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
      drt3=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
      uforget=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
      u_bak_fixup=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
      p_bak_fixup=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
      u_bak_subzone=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));
      p_bak_subzone=(ldouble*)malloc((SX)*(SY)*(SZ)*NV*sizeof(ldouble));

      //wavespeeds hd and rad - max(al,ar)
      ahdx=(ldouble*)malloc((SX)*(SY)*(SZ)*sizeof(ldouble));
      ahdy=(ldouble*)malloc((SX)*(SY)*(SZ)*sizeof(ldouble));
      ahdz=(ldouble*)malloc((SX)*(SY)*(SZ)*sizeof(ldouble));
      aradx=(ldouble*)malloc((SX)*(SY)*(SZ)*sizeof(ldouble));
      arady=(ldouble*)malloc((SX)*(SY)*(SZ)*sizeof(ldouble));
      aradz=(ldouble*)malloc((SX)*(SY)*(SZ)*sizeof(ldouble));

      //wavespeeds hd and rad - leftgoing
      ahdxl=(ldouble*)malloc((SX)*(SY)*(SZ)*sizeof(ldouble));
      ahdyl=(ldouble*)malloc((SX)*(SY)*(SZ)*sizeof(ldouble));
      ahdzl=(ldouble*)malloc((SX)*(SY)*(SZ)*sizeof(ldouble));
      aradxl=(ldouble*)malloc((SX)*(SY)*(SZ)*sizeof(ldouble));
      aradyl=(ldouble*)malloc((SX)*(SY)*(SZ)*sizeof(ldouble));
      aradzl=(ldouble*)malloc((SX)*(SY)*(SZ)*sizeof(ldouble));

      //wavespeeds hd and rad - lrightgoing
      ahdxr=(ldouble*)malloc((SX)*(SY)*(SZ)*sizeof(ldouble));
      ahdyr=(ldouble*)malloc((SX)*(SY)*(SZ)*sizeof(ldouble));
      ahdzr=(ldouble*)malloc((SX)*(SY)*(SZ)*sizeof(ldouble));
      aradxr=(ldouble*)malloc((SX)*(SY)*(SZ)*sizeof(ldouble));
      aradyr=(ldouble*)malloc((SX)*(SY)*(SZ)*sizeof(ldouble));
      aradzr=(ldouble*)malloc((SX)*(SY)*(SZ)*sizeof(ldouble));

      //timesteps required by each cell
      cell_tstepden=(ldouble*)malloc((SX)*(SY)*(SZ)*sizeof(ldouble));
      cell_dt=(ldouble*)malloc((SX)*(SY)*(SZ)*sizeof(ldouble));
    }
  
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
  free(pinit);
  free(pproblem1);
  free(pproblem2);
  free(ptemp1);
  free(pvecpot);

#ifdef MAGNFIELD
  free(emf);
#endif

  free(ptm1);
  free(ptm2);
  for(i=0;i<12;i++)
  free(msgbufs[i]);
  free(msgbufs);
  free(px);
  free(py);
  free(pz);
  free(u);
  free(g);
  free(G);
  free(gKr);
  //free(gKrbx);
  //free(gKrby);
  //free(gKrbz);
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
 
  //free(ubLx);
  //free(ubRx);
  //free(ubLy);
  //free(ubRy);
  //free(ubLz);
  //free(ubRz);
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
 

  free(ut0);
  free(ut1);
  free(ut2);
  free(ut3);
  free(dut0);
  free(dut1);
  free(dut2);
  free(dut3);
  free(drt0);
  free(drt1);
  free(drt2);
  free(drt3);
  free(uforget);
  free(u_bak_fixup);
  free(p_bak_fixup);
  free(u_bak_subzone);
  free(p_bak_subzone);

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
  free(cell_dt);
  free(cell_tstepden);

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//inverse general matrix
int
inverse_matrix(ldouble *a, ldouble *ia, int N)
{
  /*
  gsl_matrix_view m 
    = gsl_matrix_view_array (a, N, N);
  gsl_matrix_view im 
    = gsl_matrix_view_array (ia, N, N);
  
  gsl_permutation * p = gsl_permutation_alloc (N);

  int s;

  gsl_linalg_LU_decomp (&m.matrix, p, &s);

  gsl_linalg_LU_invert (&m.matrix, p, &im.matrix);
 
  gsl_permutation_free(p);
  */

  //exit(0);
  //printf("inversing... ");

  gsl_matrix *m
    = gsl_matrix_alloc (N, N);
  gsl_matrix *im
    = gsl_matrix_alloc (N, N);
  int i,j;
  for(i=0;i<N;i++)
      for(j=0;j<N;j++)
	  gsl_matrix_set(m,i,j,a[i*N+j]);

  gsl_permutation * p = gsl_permutation_alloc (N);

  int s;

  gsl_linalg_LU_decomp (m, p, &s);

  gsl_linalg_LU_invert (m, p, im);
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      ia[i*N+j]=gsl_matrix_get(im,i,j);

  gsl_matrix_free(m);
  gsl_matrix_free(im);
  gsl_permutation_free(p);
  

//printf("done! \n");
  /*
  gsl_matrix *im
    = gsl_matrix_alloc (N, N);
  int i,j;
  for(i=0;i<N;i++)
    {
      for(j=0;j<N;j++)
	{
	  gsl_matrix_set(im,i,j,a[i*N+j]);
	  //printf("%e ",a[i*N+j]);
	}
      //printf("\n");
    }

  gsl_linalg_cholesky_decomp (im);

  gsl_linalg_cholesky_invert (im);

  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      ia[i*N+j]=gsl_matrix_get(im,i,j);
   
  gsl_matrix_free(im);
  */
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
    //    my_err("det in inverse 4x4 zero\n");
    return -1;

  for (j = 0; j < 16; j++)
    dst[j] *= det;

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	ia[i][j]= dst[i*4+j];
	if(isnan(ia[i][j])) return -1;
      }

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//determinant of 4by4 matrix
ldouble
determinant_44matrix(ldouble a[][4])
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

  return det;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//Heavyside step function around 0
//1       /-------
//       | 
//0 ____/0
//x9 determines the sharpness and says where step function equals 0.95
ldouble
step_function(ldouble x,ldouble x9)
{
  ldouble k=1.47222/x9;
  return 1./(1.+exp(-2.*k*x));
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//prints error message and gets chars
int
my_err(char *message)
{
  if(PROCID==0)
    {
      char bufor[200];
      sprintf(bufor,"|err| : %s\n",message);
      printf("%s",bufor);
      getchar();
    }
  return 0;  
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//prints warning message and returns
int
my_warning(char *message)
{
  if(PROCID==0)
    {
      char bufor[200];
      sprintf(bufor,"|warning| : %s\n",message);
      printf("%s",bufor);
    }
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
//* atan2 in [0,2pi] ***************************************************
//**********************************************************************
ldouble
my_atan2(ldouble y, ldouble x)
{
  ldouble res=atan2(y,x);
  if(res<0.) res+=2.*M_PI;
   
  return res;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
/* Arrange the N elements of ARRAY in random order.
   Only effective if N is much smaller than RAND_MAX;
   if this may not be the case, use a better random
   number generator. 

http://benpfaff.org/writings/clc/shuffle.html

*/

void shuffle_loop(int **array, size_t n)
{
  if (n > 1) {
    size_t i;
    for (i = 0; i < n - 1; i++) {
      size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
      int t[3] = {array[j][0],array[j][1],array[j][2]};
      array[j][0] = array[i][0];
      array[j][1] = array[i][1];
      array[j][2] = array[i][2];
      array[i][0] = t[0];
      array[i][1] = t[1];
      array[i][2] = t[2];
    }
  }
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates eigen values of a symmetric 4x4 matrix
ldouble
calc_eigen_3x3symm(ldouble g[][3], ldouble *ev)
{
  int verbose=0;

  double matrix[]={g[0][0],g[0][1],g[0][2],
		   g[1][0],g[1][1],g[1][2],
		   g[2][0],g[2][1],g[2][2]};		       

  gsl_matrix_view m = gsl_matrix_view_array (matrix, 3, 3);     
  gsl_vector *eval = gsl_vector_alloc (3);
  gsl_eigen_symm_workspace * w = gsl_eigen_symm_alloc (3);       
  gsl_eigen_symm (&m.matrix, eval, w);     
  gsl_eigen_symm_free (w);
       
  int i,j;
     
  for (i = 0; i < 3; i++)
    {
      double eval_i 
	= gsl_vector_get (eval, i);

      ev[i]=eval_i;
    }

  gsl_vector_free (eval);

  return my_max(my_max(fabs(ev[0]),fabs(ev[1])),fabs(ev[2]));
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//verifies the eigenvalues are positive, resets negative to zeros if not
int
make_matrixposdef_3x3symm(ldouble g[][3])
{
  int verbose=0;

  double matrix[]={g[0][0],g[0][1],g[0][2],
		   g[1][0],g[1][1],g[1][2],
		   g[2][0],g[2][1],g[2][2]};		       

  gsl_matrix_view m = gsl_matrix_view_array (matrix, 3, 3);     
  gsl_vector *evals = gsl_vector_alloc (3);
  gsl_matrix *evecs = gsl_matrix_alloc(3,3); //matrix of eigenvectors
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);       
  gsl_eigen_symmv (&m.matrix, evals, evecs, w);     
  gsl_eigen_symmv_free (w);

  int i,j;

  if(gsl_vector_get(evals,0)>=0. && gsl_vector_get(evals,1)>=0. && gsl_vector_get(evals,2)>=0.) return 0; //all is well

  ldouble maxev=my_max(my_max(fabs(gsl_vector_get(evals,0)),fabs(gsl_vector_get(evals,1))),fabs(gsl_vector_get(evals,2)));
  if(maxev<0.) printf("no positive e-v\n");

  //zeroing negative eigenvalues
  ldouble zero=0.;
  if(gsl_vector_get(evals,0)<0.)
    gsl_vector_set(evals,0,zero*maxev);
  if(gsl_vector_get(evals,1)<0.)
    gsl_vector_set(evals,1,zero*maxev);
  if(gsl_vector_get(evals,2)<0.)
    gsl_vector_set(evals,2,zero*maxev);

  ldouble aP[3*3],aiP[3*3];

  gsl_matrix *D = gsl_matrix_alloc(3,3);
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      {
	if(i==j)
	  gsl_matrix_set(D, i,j,gsl_vector_get(evals,i));
	else
	  gsl_matrix_set(D, i,j,0.);

	aP[i*3+j]=gsl_matrix_get(evecs, i, j);
     }
	  
  inverse_matrix(aP,aiP,3);

  gsl_matrix_view P = gsl_matrix_view_array (aP, 3, 3);
  gsl_matrix_view iP = gsl_matrix_view_array (aiP, 3, 3);
  gsl_matrix *PD = gsl_matrix_alloc(3,3);

  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, &P.matrix, D,
                  0.0, PD);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, PD, &iP.matrix,
                  0.0, &P.matrix);

  //returning modified matrix
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      g[i][j]=gsl_matrix_get(&P.matrix, i, j);

  //verbose
  /*
  for(i = 0; i < 3; i++)
    printf("%e ",gsl_vector_get(evals,i));
  printf("\n");
  printf("\n");
	 
  for(i = 0; i < 3; i++)
    {
      for(j = 0; j < 3; j++)
	printf("%e ",g[i][j]);

      printf("\n");
    }
  printf("\n");
  for(i = 0; i < 3; i++)
    {
      for(j = 0; j < 3; j++)
	printf("%e ",gsl_matrix_get(&P.matrix,i,j));

      printf("\n");
    }

  getch();
  */

  gsl_vector_free (evals);
  gsl_matrix_free (evecs);
  gsl_matrix_free (D);
  gsl_matrix_free (PD);
 
  return 1;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates eigen values of a symmetric 4x4 matrix
ldouble
calc_eigen_4x4symm(ldouble g[][4], ldouble *ev)
{
  int verbose=0;
  double matrix[]={g[0][0],g[0][1],g[0][2],g[0][3],
		   g[1][0],g[1][1],g[1][2],g[1][3],
		   g[2][0],g[2][1],g[2][2],g[2][3],
		   g[3][0],g[3][1],g[3][2],g[3][3]};		       

  gsl_matrix_view m = gsl_matrix_view_array (matrix, 4, 4);     
  gsl_vector *eval = gsl_vector_alloc (4);
  gsl_eigen_symm_workspace * w = gsl_eigen_symm_alloc (4);       
  gsl_eigen_symm (&m.matrix, eval, w);     
  gsl_eigen_symm_free (w);
       
  int i,j;
     
  for (i = 0; i < 4; i++)
    {
      double eval_i 
	= gsl_vector_get (eval, i);

      ev[i]=eval_i;
    }

  gsl_vector_free (eval);

  return my_max(my_max(fabs(ev[0]),fabs(ev[1])),my_max(fabs(ev[2]),fabs(ev[3])));
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates eigen values of a non-symmetric 4x4 matrix
ldouble
calc_eigen_4x4(ldouble g[][4], ldouble *ev)
{
  int verbose=0;
  double matrix[]={g[0][0],g[0][1],g[0][2],g[0][3],
		   g[1][0],g[1][1],g[1][2],g[1][3],
		   g[2][0],g[2][1],g[2][2],g[2][3],
		   g[3][0],g[3][1],g[3][2],g[3][3]};		       

  gsl_matrix_view m = gsl_matrix_view_array (matrix, 4, 4);     
  gsl_vector_complex *eval = gsl_vector_complex_alloc (4);
  gsl_eigen_nonsymm_workspace * w = gsl_eigen_nonsymm_alloc (4);       
  gsl_eigen_nonsymm (&m.matrix, eval, w);     
  gsl_eigen_nonsymm_free (w);
       
  int i,j;
     
  for (i = 0; i < 4; i++)
    {
      gsl_complex eval_i 
           = gsl_vector_complex_get (eval, i);
      ev[i]=GSL_REAL(eval_i);
    }

  gsl_vector_complex_free (eval);

  return my_max(my_max(fabs(ev[0]),fabs(ev[1])),my_max(fabs(ev[2]),fabs(ev[3])));
}

/*****************************************************************/
/* prints tensor to screen *****************************************/
/*****************************************************************/
int
print_tensor(ldouble T[][4])
{
  int i;
  printf("============\n");
  for(i=0;i<4;i++)
    printf("%10.16e %10.16e %10.16e %10.16e\n",T[0][i],T[1][i],T[2][i],T[3][i]);
  printf("============\n");
  return 0;  
}

/*****************************************************************/
/* prints metric to screen *****************************************/
/*****************************************************************/
int
print_metric(ldouble T[][5])
{
  int i;
  printf("============\n");
  for(i=0;i<4;i++)
    printf("%10.3e %10.3e %10.3e %10.3e\n",T[i][0],T[i][1],T[i][2],T[i][3]);
  printf("============\n");
  return 0;  
}

/*****************************************************************/
/* prints 4vector to screen ****************************************/
/*****************************************************************/
int
print_4vector(ldouble v[4])
{
  int i;
  printf("============\n");
  printf("%10.8e %10.8e %10.8e %10.8e\n",v[0],v[1],v[2],v[3]);
  printf("============\n");
  return 0;  
}

/*****************************************************************/
/* prints Nvvector to screen *************************************/
/*****************************************************************/
int
print_NVvector(ldouble *v)
{
  print_Nvector(v,NV);
  return 0;
}

/*****************************************************************/
/* prints primitives to screen *************************************/
/*****************************************************************/
int
print_primitives(ldouble *p)
{
  printf("\n");
  printf("rho = %.15e\n",p[RHO]);
  printf("ugas = %.15e\n",p[UU]);
  printf("u^1 = %.15e\n",p[VX]);
  printf("u^2 = %.15e\n",p[VY]);
  printf("u^3 = %.15e\n",p[VZ]);
  printf("S = %.15e\n",p[ENTR]);
#ifdef MAGNFIELD
  printf("B^1 = %.15e\n",p[B1]);
  printf("B^2 = %.15e\n",p[B2]);
  printf("B^3 = %.15e\n",p[B3]);
#endif
#ifdef RADIATION
  printf("Erf = %.15e\n",p[EE0]);
  printf("ur^1 = %.15e\n",p[FX0]);
  printf("ur^2 = %.15e\n",p[FY0]);
  printf("ur^3 = %.15e\n",p[FZ0]);
#ifdef NCOMPTONIZATION
  printf("nph = %.15e\n",p[NF0]);
#endif
#endif
   printf("\n");

  return 0;
}

/*****************************************************************/
/* prints primitives to screen *************************************/
/*****************************************************************/
int
print_conserved(ldouble *u)
{
  printf("\n");
  printf("rho u^t = %.15e\n",u[RHO]);
  printf("T^t_t + rho u^t = %.15e\n",u[UU]);
  printf("T^t_1 = %.15e\n",u[VX]);
  printf("T^t_2 = %.15e\n",u[VY]);
  printf("T^t_3 = %.15e\n",u[VZ]);
  printf("S u^t = %.15e\n",u[ENTR]);
#ifdef MAGNFIELD
  printf("B^1 = %.15e\n",u[B1]);
  printf("B^2 = %.15e\n",u[B2]);
  printf("B^3 = %.15e\n",u[B3]);
#endif
#ifdef RADIATION
  printf("R^t_t = %.15e\n",u[EE0]);
  printf("R^t_1 = %.15e\n",u[FX0]);
  printf("R^t_2 = %.15e\n",u[FY0]);
  printf("R^t_3 = %.15e\n",u[FZ0]);
#ifdef NCOMPTONIZATION
  printf("nph u^t = %.15e\n",u[NF0]);
#endif
#endif
  printf("\n");
 
  return 0;
}

/*****************************************************************/
/* prints Nvector to screen ****************************************/
/*****************************************************************/
int
print_Nvector(ldouble v[4],int N)
{
  int i;
  printf("============\n");
  for(i=0;i<N;i++)
  printf("%10.16e ",v[i]);
  printf("\n============\n");
  return 0;  
}


int
my_clock_gettime(void* tsptr)
{
struct timespec *ts
    = (struct timespec *) tsptr;

#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  ts->tv_sec = mts.tv_sec;
  ts->tv_nsec = mts.tv_nsec;

#else
  clock_gettime(CLOCK_REALTIME, ts);
#endif
  return 0;
}


//frequency integrated opacity after appendix of Bell & Lin (ApJ, 427, 987, 1994)
ldouble
opacity_BellLin(ldouble rhoc, ldouble Tc)
{
  ldouble X=0.75;
  ldouble Z=0.02;
      
  //electron scattering
  ldouble opaces=0.2*(1.0+X);

  //Bound-free and free-free - Kramers?
  //opacKr=4.d25*(1.d0+X)*(Z+1.d-3)*rhoc/Tc**3.5
  ldouble opacKr=1.5e20*rhoc/Tc/Tc/sqrt(Tc);

  //H-scattering
  //opacHm=1.1d-25*sqrt(Z*rhoc)*Tc**7.7
  //opacHm=1.d-13*rhoc**0.333333*Tc**4
  ldouble opacHm=1.e-36*cbrt(rhoc)*Tc*Tc*Tc*Tc*Tc*Tc*Tc*Tc*Tc*Tc;

  //molecules
  //opacmol=0.1d0*Z
  ldouble opacmol=1.e-8*cbrt(rhoc*rhoc)*Tc*Tc*Tc;

  //evaporation of metal grains
  ldouble Tc24=Tc*Tc; Tc24*=Tc24; Tc24*=Tc24; Tc24*=Tc24*Tc24;
  ldouble opacmetevap=2.e81*rhoc/Tc24;

  //metal grains
  ldouble opacmet=0.1*sqrt(Tc);

  //evaporation of ice grains
  ldouble opaciceevap=2.e16/Tc/Tc/Tc/Tc/Tc/Tc/Tc;
  
  //ice grains
  ldouble opacice=2.e-4*Tc*Tc;

  #ifdef SIMPLEBELLLIN
  ldouble opacity=1.e0/((1.e0/(opacHm+opacmol))+(1.e0/(0.*opaces+opacKr)));    
  #else
  ldouble opacity=1.e0/((1.e0/(opacHm+opacmol))+(1.e0/(0.*opaces+opacKr)))
    +1.e0/(1.e0/opacice+1.e0/opaciceevap)
    +1.e0/(1.e0/opacmet+1.e0/opacmetevap);
  #endif
  
  return opacity;
}


//decomposes velocities into cartesian
int
decompose_vels(ldouble *pp,int velidx, ldouble v[4],void *ggg,  void *gggBL)
{
  struct geometry *geom
    = (struct geometry *) ggg;
  struct geometry *geomBL
    = (struct geometry *) gggBL;
  int iv;
  int velprim=VELPRIM; if(velidx==FX) velprim=VELPRIMRAD;
  ldouble r,th,ph,ucon[4];
  DLOOPA(iv) ucon[iv]=pp[velidx-1+iv];
  ucon[0]=0.;
  //conv_vels(ucon,ucon,velprim,VEL4,geom->gg,geom->GG);
  //trans2_coco(geom->xxvec,ucon,ucon,MYCOORDS,BLCOORDS);
  //conv_vels(ucon,ucon,VEL4,velprim,geomBL->gg,geomBL->GG);
  r=geomBL->xx;	     th=geomBL->yy;	     ph=geomBL->zz;
  //to ortonormal cartesian
  //ucon[2]*=r;
  //  ucon[3]*=r*sin(th);
  ucon[1]*=sqrt(geom->gg[1][1]);
  ucon[2]*=sqrt(geom->gg[2][2]);
  ucon[3]*=sqrt(geom->gg[3][3]);

  v[1] = sin(th)*cos(ph)*ucon[1] + cos(th)*cos(ph)*ucon[2] - sin(ph)*ucon[3];
  v[2] = sin(th)*sin(ph)*ucon[1] + cos(th)*sin(ph)*ucon[2] + cos(ph)*ucon[3];
  v[3] = cos(th)*ucon[1] - sin(th)*ucon[2];
  return 0;
}

//returns problem specific time limiter
ldouble
get_tsteplimiter()
{
  if(PROBLEM==7) //BONDI
    {
      if(global_time<1.e8)  return 0.01;
      else if(global_time<1.e9)  return 0.03;
      else return 0.05;
    }

    return 0.5;
}
