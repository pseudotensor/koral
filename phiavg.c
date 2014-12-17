#include <stdio.h>
#include <stdlib.h>

#include "ko.h"

int
main
(int argc, char **argv)
{  
  //which files to read
  int no1,no2,nostep,ifavg;
  
  int ifile,itot=0;
  int i,j,k,iv;
  int nx,ny,nz,nv,ret,problem;
  double time;

  if(argc!=5)
    {
      printf("Not enough input arguments. Asks for ./phiavg no1 no2 nostep ifavg\n");
      return -1;
    }
  else
    {      
      no1=atof(argv[1]);
      no2=atof(argv[2]);
      nostep=atof(argv[3]);
      ifavg=atof(argv[4]);
    }

  char folder[100],folderout[100],bufer[100],base[10],fname[100],fnamehead[100],fnameout[100],fnameheadout[100],bufor[500];
  sprintf(folder,"%s","dumps");
  sprintf(folderout,"%s","dumps_phiavg");


  nv=NV;

  if(ifavg)
    {
      sprintf(base,"%s","avg");
      nv+=NAVGVARS;
    }
  else
    sprintf(base,"%s","res");


  printf("phiavg: %d %d %d %d\n",no1,no2,nostep,ifavg);


  //allocate memory
  nx=NX;ny=NY;
  double ***prims,*pp;
  pp=(double*)malloc(nv*sizeof(double));
  prims=(double***)malloc(nx*sizeof(double**));
  for(i=0;i<nx;i++)
    {
      prims[i]=(double**)malloc(ny*sizeof(double*));
      for(j=0;j<ny;j++)
	prims[i][j]=(double*)malloc(nv*sizeof(double));
    }

  prims[0][0][0]=0.;

  for(ifile=no1;ifile<=no2;ifile+=nostep)
    {
      itot++;

      //reading the files
      sprintf(fname,"%s/%s%04d.dat",folder,base,ifile);
      sprintf(fnamehead,"%s/%s%04d.head",folder,base,ifile);
      sprintf(fnameheadout,"%s/%s%04d.head",folderout,base,ifile);
      sprintf(fnameout,"%s/%s%04d.dat",folderout,base,ifile);

      FILE *fdump,*fout;

      

      /***********/
      //header file
      fdump=fopen(fnamehead,"r");

      //reading damp file parameters
      int intpar[6];

      if(ifavg)
	{
	  sprintf(bufor,"cp %s %s\n",fnamehead,fnameheadout);
	  system(bufor);
	  nx=TNX;
	  ny=TNY;
	  nz=TNZ;
	  problem=PROBLEM;
	}
      else
	{
	  ret=fscanf(fdump,"## %d %d %lf %d %d %d %d\n",&intpar[0],&intpar[1],&time,&intpar[2],&intpar[3],&intpar[4],&intpar[5]);
	  problem=intpar[2];
	  nx=intpar[3];
	  ny=intpar[4];
	  nz=intpar[5];
	  fclose(fdump);
      
	  fout=fopen(fnameheadout,"w");
	  sprintf(bufor,"## %5d %5d %10.6e %5d %5d %5d %5d\n",intpar[0],intpar[1],time,intpar[2],intpar[3],intpar[4],intpar[5]);
	  fprintf(fout,"%s",bufor);
	  fprintf(fout,"phi-averaged to NZ=1\n");
	  fclose(fout);
	}
      

      /***********/
      fout=fopen(fnameout,"w");

     
      
      if(nx!=TNX || ny!=TNY || problem!=PROBLEM)
	{
	  printf("PROBLEM files inconsistent with the dump file\n");
	  exit (-1);
	}

      printf("restart file (%s) read no. %d at time: %f of PROBLEM: %d with NXYZ: %d %d %d\n",
	     fname,intpar[0],time,intpar[2],intpar[3],intpar[4],intpar[5]); 


      printf("time: %e resolution: %d x %d x %d NV: %d\n",time,nx,ny,nz,nv);


      //body file
      fdump=fopen(fname,"rb");
 
      for(j=0;j<ny;j++)
	for(i=0;i<nx;i++)
	  for(iv=0;iv<nv;iv++)
	    prims[i][j][iv]=0.;

      for(k=0;k<nz;k++)
	{
	  for(j=0;j<ny;j++)
	    {
	      for(i=0;i<nx;i++)
		{
	    	  int gix,giy,giz;
		  
		  ret=fread(&gix,sizeof(int),1,fdump);
		  ret=fread(&giy,sizeof(int),1,fdump);
		  ret=fread(&giz,sizeof(int),1,fdump);
		  ret=fread(pp,sizeof(double),nv,fdump);


		  if(gix<0 || gix>=nx ||giy<0 || giy>=ny)
		    printf("blont: %d %d %d vs %d %d %d | %d %d %d\n",i,j,k,gix,giy,giz,nx,ny,nz);
		  else
		    for(iv=0;iv<nv;iv++)
		      prims[gix][giy][iv]+=pp[iv];
		}
	    }
	}

      for(j=0;j<ny;j++)
	for(i=0;i<nx;i++)
	  for(iv=0;iv<nv;iv++)
	    prims[i][j][iv]/=(double)nz;

      //print to a file
      k=0;
      for(j=0;j<ny;j++)
	for(i=0;i<nx;i++)
	  {
	    fwrite(&i,sizeof(int),1,fout);
	    fwrite(&j,sizeof(int),1,fout);
	    fwrite(&k,sizeof(int),1,fout);
	    fwrite(&prims[i][j][0],sizeof(ldouble),nv,fout);	    
	  }

	
      fclose(fdump);
      fclose(fout);
    }
  return 0;
}
