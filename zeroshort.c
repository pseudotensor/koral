#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "ko.h"

//temporary
#define GRID_SPACING get_size_x(0,0)      //distance between adjacent grid cells (grid is assumed cubic with equal spacing in all directions)
#define LIGHT_C 1.		//Speed of light in code units
#define STEFAN_BOLTZMANN SIGMA_RAD //stefan boltzmann constant in code units
#define PI 3.14159265358979

#define BINAVERAGE 0


//------------  Subroutines to setup BSP tree data structure  -------------------------

struct solverarg
{
  double F_Start_forward;
  double F_Start_backward;
  double *Intensities;			
  double *F_final_norm;
  double fFinal;
  double *cosang;
};

/**********************/
double
f_stretchFactor (double k, void *argsin)
{
  struct solverarg *args
    = (struct solverarg *) argsin;

  double sumplus,summinus;
  sumplus=summinus=0.;

  int p;
  for(p=0;p<NUMANGLES;p++)
    {
      if(args->cosang[p]>0.)
	sumplus+=args->Intensities[p]*sqrt(1.+(k*k-1.)*args->cosang[p]*args->cosang[p]);
      else
	summinus+=args->Intensities[p]*sqrt(1.+(1./k/k-1.)*args->cosang[p]*args->cosang[p]);
    }
 
  return args->fFinal - 
    (k*args->F_Start_forward + 1./k*args->F_Start_backward)/
    (sumplus + summinus);
}


double
calc_stretchFactor(void *argsin)
{

  int status;
  int iter = 0, max_iter = 100;
  double k = 1., kprev;				
  double f,df,dk;
  
  /*
  for(k=-1.e5;k<-1.e-5;k/=1.1)
    {
      printf("%e %e\n",k,f_stretchFactor(k, argsin));
    }

  exit(1);
  */
  do
    {
      iter++;
      if(iter>100)
	{
	  struct solverarg *args
	    = (struct solverarg *) argsin;
	  printf("calc_stretchFactor() failed to find solution with F/E=%f\n",args->fFinal);
	  
	  /*
	  int p;
	  for(p=0;p<NUMANGLES;p++)
	    printf("%d %e\n",p,args->Intensities[p]);

	  double k;
	  for(k=.5;k<1.;k*=1.01)
	    {
	    printf("%e %e\n",k,f_stretchFactor (k, argsin));

	  
	    int p;double sumplus,summinus;
	    sumplus=summinus=0.;
	    for(p=0;p<NUMANGLES;p++)
	      {
		if(args->cosang[p]>0.)
		  sumplus+=args->Intensities[p]*sqrt(1.+(k*k-1.)*args->cosang[p]*args->cosang[p]);
		else
		  summinus+=args->Intensities[p]*sqrt(1.+(1./k/k-1.)*args->cosang[p]*args->cosang[p]);
	      }
  
	    printf("%e %e %e %e %e\n",(k*args->F_Start_forward + 1./k*args->F_Start_backward),
		   (sumplus + summinus),sumplus,summinus,
		   args->fFinal - 
		   (k*args->F_Start_forward + 1./k*args->F_Start_backward)/
		   (sumplus + summinus));
	    }
	  */


	  exit(-1);
	}

      f=f_stretchFactor(k, argsin);
      df=(f_stretchFactor(k+1.e-6*k, argsin) - f)/(1.e-6*k);

      kprev=k;

      dk=f/df;

      
      if(fabs(dk)>.5*k)
	if(dk>0.) dk=.5*k;
	else dk=-.5*k;
      

      k=k-dk;      

      //printf ("%5d %.7f %.7f\n", iter, k, k-kprev);
      
    }
  while (fabs((k-kprev)/k)>1.e-3);
  //while(iter<10);
  
  return k;
}

double
calc_stretchFactor_gsl(void *argsin)
{

  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *solv;
  double s = 1.;
  double s_lo = 1.e-3, s_hi = 1.e3;
  gsl_function F;

  /*
  for(s=1.e-5;s<1.e5;s*=1.1)
    {
      printf("%e %e\n",s,f_stretchFactor(s, argsin));
    }

  exit(1);
  */

  F.function = &f_stretchFactor;
  F.params = argsin;

  T = gsl_root_fsolver_brent;
  solv = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (solv, &F, s_lo, s_hi);
  
  /*
  printf ("using %s method\n", 
          gsl_root_fsolver_name (solv));

  printf ("%5s [%9s, %9s] %9s %9s\n",
          "iter", "lower", "upper", "root", 
           "err(est)");
  */

  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (solv);
      s = gsl_root_fsolver_root (solv);
      s_lo = gsl_root_fsolver_x_lower (solv);
      s_hi = gsl_root_fsolver_x_upper (solv);
      status = gsl_root_test_interval (s_lo, s_hi,
                                       0, 0.0001);
      
      
      if (status == GSL_SUCCESS)
        printf ("Converged:\n");

      printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
              iter, s_lo, s_hi,
              s,   
              s_hi - s_lo);
      
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free (solv);

  return s;
}



struct bsptree* create_node(int angVal, int iterVal)
{
	struct bsptree *node = (struct bsptree*) malloc( sizeof(struct bsptree) );

	node->angIndex = angVal;
	node->iter = iterVal;

	node->lower = 0;
	node->upper = 0;

	return node;
}


void destroy_tree(struct bsptree *leaf)
{
  if( leaf != 0 )
  {
      destroy_tree(leaf->lower);
      destroy_tree(leaf->upper);
      free( leaf );
  }
}

void insertLower(int angVal, int iterVal, struct bsptree *leaf)
{
	if (leaf->upper == 0)
	{
		leaf->lower = (struct bsptree*) malloc( sizeof(struct bsptree) );

		leaf->lower->lower = 0;
		leaf->lower->upper = 0;
	}

	leaf->lower->angIndex = angVal;
	leaf->lower->iter = iterVal;
}

void insertUpper(int angVal, int iterVal, struct bsptree *leaf)
{
	if (leaf->upper == 0)
	{
		leaf->upper = (struct bsptree*) malloc( sizeof(struct bsptree) );
		leaf->lower->lower = 0;
		leaf->lower->upper = 0;
	}

	leaf->upper->angIndex = angVal;
	leaf->upper->iter = iterVal;

}

// --------------------------- End subroutines for setting up/modifying BSP tree ---------------------

//Subroutine to read in data files
//Returns 1 if success, -1 if failure

int readAngleFiles(double angGridCoords[NUMANGLES][3], double angDualGridCoords[NUMDUALANGLES][3], int dualAdjacency[NUMDUALANGLES][3])
{

  int count = 0;
  double dtot = 0;
  double xval, yval, zval;

  FILE *angFile, *angDualFile, *angDualAdjFile;



  // Open several data files
  angFile = fopen("best-xyz.dat", "r");
  angDualFile = fopen("best-dualtri-xyz.dat", "r");
  angDualAdjFile = fopen("best-dualtri-adjIndex.dat", "r");


  if (!angFile)
    {
      fprintf(stderr, "ERROR:  Unable to open Angle Grid File!\n");
      return -1;
    }
  if (!angDualFile)
    {
      fprintf(stderr, "ERROR:  Unable to open Angle Dual Grid File!\n");
      return -1;
    }
  if (!angDualAdjFile)
    {
      fprintf(stderr, "ERROR:  Unable to open Angle Dual Adjacency Grid File!\n");
      return -1;
    }



  //read in angle grid from file


  while(fscanf(angFile,"%lf %lf %lf",&xval,&yval,&zval) > 0)
    {
      if (count >= NUMANGLES)
	{
	  fprintf(stderr, "ERROR:  Angle Grid File has too many entries! (Expected %d lines)\n", NUMANGLES);
	  return -1;
	}

      angGridCoords[count][0] = xval;
      angGridCoords[count][1] = yval;
      angGridCoords[count][2] = zval;

      dtot = sqrt(xval*xval + yval*yval + zval*zval);

      angGridCoords[count][0] = xval/dtot;
      angGridCoords[count][1] = yval/dtot;
      angGridCoords[count][2] = zval/dtot;


      //		printf("READ = %lf %lf %lf\n",angGridX[count],angGridY[count],angGridZ[count]);  
      count++;
    }  
  fclose(angFile);

  if (count < NUMANGLES)
    {
      fprintf(stderr, "ERROR:  Angle Grid File has too few entries! (Expected %d lines, got %d)\n", NUMANGLES, count);
      return -1;
    }




  //read in angle dual grid from file

  count = 0;
  dtot = 0;
  while(fscanf(angDualFile,"%lf %lf %lf",&xval,&yval,&zval) > 0)
    {
      if (count >= NUMDUALANGLES)
	{
	  fprintf(stderr, "ERROR:  Angle Dual Grid File has too many entries! (Expected %d lines)\n", NUMDUALANGLES);
	  return -1;
	}

      angDualGridCoords[count][0] = xval;
      angDualGridCoords[count][1] = yval;
      angDualGridCoords[count][2] = zval;

      dtot = sqrt(xval*xval + yval*yval + zval*zval);

      angDualGridCoords[count][0] = xval/dtot;
      angDualGridCoords[count][1] = yval/dtot;
      angDualGridCoords[count][2] = zval/dtot;

      //		printf("READ = %lf %lf %lf\n",xval,yval,zval);  
      //		printf("READ = %e %e %e\n",angDualGridX[count],angDualGridY[count],angDualGridZ[count]);  

      count++;
    }  
  fclose(angDualFile);

  if (count < NUMDUALANGLES)
    {
      fprintf(stderr, "ERROR:  Angle Dual Grid File has too few entries! (Expected %d lines, got %d)\n", NUMDUALANGLES, count);
      return -1;
    }




  //read in angle dual adjacency grid from file

  count = 0;
  dtot = 0;

  int index0, index1, index2;

  while(fscanf(angDualAdjFile,"%d %d %d",&index0,&index1,&index2) > 0)
    {
      if (count >= NUMDUALANGLES)
	{
	  fprintf(stderr, "ERROR:  Angle Dual Grid File has too many entries! (Expected %d lines)\n", NUMDUALANGLES);
	  return -1;
	}

      dualAdjacency[count][0] = index0;
      dualAdjacency[count][1] = index1;
      dualAdjacency[count][2] = index2;

      count++;
    }  
  fclose(angDualAdjFile);

  if (count < NUMDUALANGLES)
    {
      fprintf(stderr, "ERROR:  Angle Dual Grid File has too few entries! (Expected %d lines, got %d)\n", NUMDUALANGLES, count);
      return -1;
    }





  return 1;  //successfully read in all files!
}

// Subroutine to calculate the interpolation weights for a spherical 3x3x3 cuboid grid

// We find the 4 bounding corners for each ray hitting the boundary, and their appropriate weights to use in interpolation later on

// Assume existence of getCoord(n1, n2, n3) that returns the r, theta, phi coordinates at grid cell labeled n1, n2, n3



void setupInterpWeights_sph3D(int ix, int iy, int iz, double angGridCoords[NUMANGLES][3], int intersectGridIndices[SX][SY][SZ][NUMANGLES][3][4], double intersectGridWeights[SX][SY][SZ][NUMANGLES][4], double intersectDistances[SX][SY][SZ][NUMANGLES])
{
  int n1_central=1,n2_central=1,n3_central=1; //use as central index, will loop over later
  int delta_n1, delta_n2, delta_n3;

  double coord_values[3];
  double xxvec[4];
  double coord_limits[3][3]; //store the coordinate boundaries in r,theta,phi
  //indices are [coordinate r,th,phi][min,mid,max]


  double posX0, posY0, posZ0; //Coordinates of central cell
  double posR0, posTh0, posPh0;

  int l,m;


  //Read in coordinates
  for (l=0; l < 3; l++)
    {
      for (m=-1; m <= 1; m++)
	{
	  if (l==0)
	    {
	      delta_n1=m;
	      delta_n2=0;
	      delta_n3=0;
	    }
	  if (l==1)
	    {
	      delta_n1=0;
	      delta_n2=m;
	      delta_n3=0;
	    }
	  if (l==2)
	    {
	      delta_n1=0;
	      delta_n2=0;
	      delta_n3=m;
	    }
	  //getCoord(n1_central+delta_n1, n2_central+delta_n2, n3_central+delta_n3, coord_values);
	  get_xx_arb(ix+delta_n1, iy+delta_n2, iz+delta_n3, xxvec, SPHCOORDS);
	  coord_values[0]=xxvec[1];
	  coord_values[1]=xxvec[2];
	  coord_values[2]=xxvec[3];

	  coord_limits[l][m+1]=coord_values[l];
	}
    }

  posR0 = coord_limits[0][1];
  posTh0 = coord_limits[1][1];
  posPh0 = coord_limits[2][1];

  posX0 = posR0*sin(posTh0)*cos(posPh0);
  posY0 = posR0*sin(posTh0)*sin(posPh0);
  posZ0 = posR0*cos(posTh0);

  /*
  if(ix==NX/2 && iy==NY/2 && iz==0)
    { 
      for (l=0; l < 3; l++)
	{
	  for (m=0; m <3; m++)
	    printf("%d %d %e\n", l,m,coord_limits[l][m]);
	}
      exit(1);

    }
  */

  int probeAng;
  for (probeAng = 0; probeAng < NUMANGLES; probeAng++)
    {

      double posX, posY, posZ;
      double posR, posTh, posPh;

      double minL = 0., maxL = 2.0*sqrt(coord_limits[0][2]*coord_limits[0][2] - coord_limits[0][0]*coord_limits[0][0]);


      //Bisect on ray to numerically find intersection location
      int iter;
      for (iter=0; iter < 50; iter++)
	{


	  posX = posX0 + (minL+maxL)/2.0 * (-angGridCoords[probeAng][0]);
	  posY = posY0 + (minL+maxL)/2.0 * (-angGridCoords[probeAng][1]);
	  posZ = posZ0 + (minL+maxL)/2.0 * (-angGridCoords[probeAng][2]);


	  posR = sqrt(posX*posX + posY*posY + posZ*posZ);
	  posTh = my_atan2(sqrt(posX*posX + posY*posY),posZ);
	  posPh = my_atan2(posY,posX);


	  //bisect on path length to find nearest boundary intersection
	  if ((posR > coord_limits[0][2]) || (posR < coord_limits[0][0]) || (posTh > coord_limits[1][2]) || (posTh < coord_limits[1][0]) || (posPh > coord_limits[2][2]) || (posPh < coord_limits[2][0]))
	    {
	      maxL = (minL+maxL)/2.0;
	    }
	  else
	    {
	      minL = (minL+maxL)/2.0;
	    }


	  // printf("%e %e (%e %e %e) (%e %e %e)\n", minL, maxL, posR0, posTh0, posPh0, posR, posTh, posPh);
	  // printf("%e %e (%e %e %e)\n", minL, maxL, posX, posY, posZ);


	}


      posX = posX0 + maxL * (-angGridCoords[probeAng][0]);
      posY = posY0 + maxL * (-angGridCoords[probeAng][1]);
      posZ = posZ0 + maxL * (-angGridCoords[probeAng][2]);

      posR = sqrt(posX*posX + posY*posY + posZ*posZ);
      posTh = my_atan2(sqrt(posX*posX + posY*posY),posZ);
      posPh = my_atan2(posY,posX);

      double lowerR, lowerTh, lowerPh;
      double lowerRIndex, lowerThIndex, lowerPhIndex;

      double upperR, upperTh, upperPh;
      double upperRIndex, upperThIndex, upperPhIndex;


      //R Indices
      if ((posR < coord_limits[0][0]) || (posR > coord_limits[0][2]))
	{
	  lowerRIndex = n1_central;
	  upperRIndex = n1_central;

	  lowerR=coord_limits[0][1];
	  upperR=coord_limits[0][1];
	}
      else if ((posR <= coord_limits[0][1]) && (posR >= coord_limits[0][0]))
	{
	  lowerRIndex = n1_central-1;
	  upperRIndex = n1_central;

	  lowerR=coord_limits[0][0];
	  upperR=coord_limits[0][1];
	}
      else 
	{
	  lowerRIndex = n1_central;
	  upperRIndex = n1_central+1;

	  lowerR=coord_limits[0][1];
	  upperR=coord_limits[0][2];
	}


      //Theta Indices
      if ((posTh < coord_limits[1][0]) || (posTh > coord_limits[1][2]))
	{
	  lowerThIndex = n2_central;
	  upperThIndex = n2_central;

	  lowerTh=coord_limits[1][1];
	  upperTh=coord_limits[1][1];
	}
      else if ((posTh <= coord_limits[1][1]) && (posTh >= coord_limits[1][0]))
	{
	  lowerThIndex = n2_central-1;
	  upperThIndex = n2_central;

	  lowerTh=coord_limits[1][0];
	  upperTh=coord_limits[1][1];
	}
      else
	{
	  lowerThIndex = n2_central;
	  upperThIndex = n2_central+1;

	  lowerTh=coord_limits[1][1];
	  upperTh=coord_limits[1][2];
	}


      //Phi Indices
      if ((posPh < coord_limits[2][0]) || (posPh > coord_limits[2][2]))
	{
	  lowerPhIndex = n3_central;
	  upperPhIndex = n3_central;

	  lowerPh=coord_limits[2][1];
	  upperPh=coord_limits[2][1];
	}
      else if ((posPh <= coord_limits[2][1]) && (posPh >= coord_limits[2][0]))
	{
	  lowerPhIndex = n3_central-1;
	  upperPhIndex = n3_central;

	  lowerPh=coord_limits[2][0];
	  upperPh=coord_limits[2][1];

	}
      else
	{
	  lowerPhIndex = n3_central;
	  upperPhIndex = n3_central+1;

	  lowerPh=coord_limits[2][1];
	  upperPh=coord_limits[2][2];
	}






      //Interpolation weight components
      double w0_low, w0_high;
      double w1_low, w1_high;



      //Get interpolation weights for 4 bounding corners of cube intersection
      if ((posR > coord_limits[0][2]) || (posR < coord_limits[0][0]))
	{

	  int realRIndex;


	  if (posR > coord_limits[0][2])
	    {
	      realRIndex = n1_central+1;
	    }
	  if (posR < coord_limits[0][0])
	    {
	      realRIndex = n1_central-1;
	    }

	  int p;
	  for (p=0; p < 4; p++)
	    {
	      intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][p] = realRIndex;
	    }
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][0] = lowerThIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][1] = lowerThIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][2] = upperThIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][3] = upperThIndex;

	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][0] = lowerPhIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][1] = upperPhIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][2] = lowerPhIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][3] = upperPhIndex;

	  w0_high = (posTh - lowerTh)/(upperTh - lowerTh);
	  w0_low = 1.0 - w0_high;
	  w1_high = (posPh - lowerPh)/(upperPh - lowerPh);
	  w1_low = 1.0 - w1_high;

	}


      if ((posTh > coord_limits[1][2]) || (posTh < coord_limits[1][0]))
	{
	  int realThIndex;

	  if (posTh > coord_limits[1][2])
	    {
	      realThIndex = n2_central+1;
	    }
	  if (posTh < coord_limits[1][0])
	    {
	      realThIndex = n2_central-1;
	    }

	  int p;
	  for (p=0; p < 4; p++)
	    {
	      intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][p] = realThIndex;
	    }
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][0] = lowerRIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][1] = lowerRIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][2] = upperRIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][3] = upperRIndex;

	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][0] = lowerPhIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][1] = upperPhIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][2] = lowerPhIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][3] = upperPhIndex;

	  w0_high = (posR - lowerR)/(upperR - lowerR);
	  w0_low = 1.0 - w0_high;
	  w1_high = (posPh - lowerPh)/(upperPh - lowerPh);
	  w1_low = 1.0 - w1_high;

	}


      if ((posPh > coord_limits[2][2]) || (posPh < coord_limits[2][0]))
	{
	  int realPhIndex;

	  if (posPh > coord_limits[2][2])
	    {
	      realPhIndex = n3_central+1;
	    }
	  if (posPh < coord_limits[2][0])
	    {
	      realPhIndex = n3_central-1;
	    }

	  int p;
	  for (p=0; p < 4; p++)
	    {
	      intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][p] = realPhIndex;
	    }

	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][0] = lowerRIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][1] = lowerRIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][2] = upperRIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][3] = upperRIndex;

	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][0] = lowerThIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][1] = upperThIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][2] = lowerThIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][3] = upperThIndex;

	  w0_high = (posR - lowerR)/(upperR - lowerR);
	  w0_low = 1.0 - w0_high;
	  w1_high = (posTh - lowerTh)/(upperTh - lowerTh);
	  w1_low = 1.0 - w1_high;

	}


      intersectGridWeights[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0] = w0_low*w1_low;
      intersectGridWeights[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1] = w0_low*w1_high;
      intersectGridWeights[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2] = w0_high*w1_low;
      intersectGridWeights[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][3] = w0_high*w1_high;

      intersectDistances[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng] = maxL;



      // printf("!! %e %e %e %e !!\n", w0_low, w0_high, w1_low, w1_high);
      // printf("-- %e %e %e --\n", posX, posY, posZ);



    } //end loop over angles


}

// Subroutine to calculate the interpolation weights for a spherical 3x3x1 cuboid grid
// Assumes symmetry in phi

// We find the 4 bounding corners for each ray hitting the boundary, and their appropriate weights to use in interpolation later on

// Assume existence of getCoord(n1, n2, n3) that returns the r, theta, phi coordinates at grid cell labeled n1, n2, n3



void setupInterpWeights_sph2D(int ix, int iy, int iz, double angGridCoords[NUMANGLES][3], int intersectGridIndices[SX][SY][SZ][NUMANGLES][3][4], double intersectGridWeights[SX][SY][SZ][NUMANGLES][4], double intersectDistances[SX][SY][SZ][NUMANGLES], double intersectGridPhi[SX][SY][SZ][NUMANGLES])
{
  int n1_central=1,n2_central=1,n3_central=1; //use as central index, will loop over later
  int delta_n1, delta_n2, delta_n3;

  double coord_values[3];
  double xxvec[4];
  double coord_limits[3][3]; //store the coordinate boundaries in r,theta,phi
  //indices are [coordinate r,th,phi][min,mid,max]


  double posX0, posY0, posZ0; //Coordinates of central cell
  double posR0, posTh0, posPh0;

  int l,m;


  //Read in coordinates
  for (l=0; l < 3; l++)
    {
      for (m=-1; m <= 1; m++)
	{
	  if (l==0)
	    {
	      delta_n1=m;
	      delta_n2=0;
	      delta_n3=0;
	    }
	  if (l==1)
	    {
	      delta_n1=0;
	      delta_n2=m;
	      delta_n3=0;
	    }
	  if (l==2)
	    {
	      delta_n1=0;
	      delta_n2=0;
	      delta_n3=m;
	    }
	  //getCoord(n1_central+delta_n1, n2_central+delta_n2, n3_central+delta_n3, coord_values);
	  get_xx_arb(ix+delta_n1, iy+delta_n2, iz+delta_n3, xxvec, SPHCOORDS);
	  coord_values[0]=xxvec[1];
	  coord_values[1]=xxvec[2];
	  coord_values[2]=xxvec[3];

	  coord_limits[l][m+1]=coord_values[l];
	}
    }

  posR0 = coord_limits[0][1];
  posTh0 = coord_limits[1][1];
  posPh0 = coord_limits[2][1];

  coord_limits[2][0]=-1.0e10;
  coord_limits[2][2]=1.0e10;


  posX0 = posR0*sin(posTh0)*cos(posPh0);
  posY0 = posR0*sin(posTh0)*sin(posPh0);
  posZ0 = posR0*cos(posTh0);

  /*
  if(ix==NX/2 && iy==NY/2 && iz==0)
    { 
      for (l=0; l < 3; l++)
	{
	  for (m=0; m <3; m++)
	    printf("%d %d %e\n", l,m,coord_limits[l][m]);
	}
      exit(1);

    }
  */

  int probeAng;
  for (probeAng = 0; probeAng < NUMANGLES; probeAng++)
    {

      double posX, posY, posZ;
      double posR, posTh, posPh;

      double minL = 0., maxL = 2.0*sqrt(coord_limits[0][2]*coord_limits[0][2] - coord_limits[0][0]*coord_limits[0][0]);


      //Bisect on ray to numerically find intersection location
      int iter;
      for (iter=0; iter < 50; iter++)
	{


	  posX = posX0 + (minL+maxL)/2.0 * (-angGridCoords[probeAng][0]);
	  posY = posY0 + (minL+maxL)/2.0 * (-angGridCoords[probeAng][1]);
	  posZ = posZ0 + (minL+maxL)/2.0 * (-angGridCoords[probeAng][2]);


	  posR = sqrt(posX*posX + posY*posY + posZ*posZ);
	  posTh = my_atan2(sqrt(posX*posX + posY*posY),posZ);
	  posPh = my_atan2(posY,posX);


	  //bisect on path length to find nearest boundary intersection
	  if ((posR > coord_limits[0][2]) || (posR < coord_limits[0][0]) || (posTh > coord_limits[1][2]) || (posTh < coord_limits[1][0]) )
	    {
	      maxL = (minL+maxL)/2.0;
	    }
	  else
	    {
	      minL = (minL+maxL)/2.0;
	    }


	  // printf("%e %e (%e %e %e) (%e %e %e)\n", minL, maxL, posR0, posTh0, posPh0, posR, posTh, posPh);
	  // printf("%e %e (%e %e %e)\n", minL, maxL, posX, posY, posZ);


	}


      posX = posX0 + maxL * (-angGridCoords[probeAng][0]);
      posY = posY0 + maxL * (-angGridCoords[probeAng][1]);
      posZ = posZ0 + maxL * (-angGridCoords[probeAng][2]);

      posR = sqrt(posX*posX + posY*posY + posZ*posZ);
      posTh = my_atan2(sqrt(posX*posX + posY*posY),posZ);
      posPh = my_atan2(posY,posX);
      intersectGridPhi[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng] = posPh - posPh0;


      double lowerR, lowerTh, lowerPh;
      double lowerRIndex, lowerThIndex, lowerPhIndex;

      double upperR, upperTh, upperPh;
      double upperRIndex, upperThIndex, upperPhIndex;


      //R Indices
      if ((posR < coord_limits[0][0]) || (posR > coord_limits[0][2]))
	{
	  lowerRIndex = n1_central;
	  upperRIndex = n1_central;

	  lowerR=coord_limits[0][1];
	  upperR=coord_limits[0][1];
	}
      else if ((posR <= coord_limits[0][1]) && (posR >= coord_limits[0][0]))
	{
	  lowerRIndex = n1_central-1;
	  upperRIndex = n1_central;

	  lowerR=coord_limits[0][0];
	  upperR=coord_limits[0][1];
	}
      else 
	{
	  lowerRIndex = n1_central;
	  upperRIndex = n1_central+1;

	  lowerR=coord_limits[0][1];
	  upperR=coord_limits[0][2];
	}


      //Theta Indices
      if ((posTh < coord_limits[1][0]) || (posTh > coord_limits[1][2]))
	{
	  lowerThIndex = n2_central;
	  upperThIndex = n2_central;

	  lowerTh=coord_limits[1][1];
	  upperTh=coord_limits[1][1];
	}
      else if ((posTh <= coord_limits[1][1]) && (posTh >= coord_limits[1][0]))
	{
	  lowerThIndex = n2_central-1;
	  upperThIndex = n2_central;

	  lowerTh=coord_limits[1][0];
	  upperTh=coord_limits[1][1];
	}
      else
	{
	  lowerThIndex = n2_central;
	  upperThIndex = n2_central+1;

	  lowerTh=coord_limits[1][1];
	  upperTh=coord_limits[1][2];
	}


      //Phi Indices

	lowerPhIndex = n3_central;
	upperPhIndex = n3_central;








      //Interpolation weight components
      double w0_low, w0_high;
      double w1_low, w1_high;



      //Get interpolation weights for 4 bounding corners of cube intersection
      if ((posR > coord_limits[0][2]) || (posR < coord_limits[0][0]))
	{

	  int realRIndex;


	  if (posR > coord_limits[0][2])
	    {
	      realRIndex = n1_central+1;
	    }
	  if (posR < coord_limits[0][0])
	    {
	      realRIndex = n1_central-1;
	    }

	  int p;
	  for (p=0; p < 4; p++)
	    {
	      intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][p] = realRIndex;
	    }
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][0] = lowerThIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][1] = lowerThIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][2] = upperThIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][3] = upperThIndex;

	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][0] = lowerPhIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][1] = upperPhIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][2] = lowerPhIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][3] = upperPhIndex;

	  w0_high = (posTh - lowerTh)/(upperTh - lowerTh);
	  w0_low = 1.0 - w0_high;
	  w1_high = 1.0;
	  w1_low = 0.0;

	}


      if ((posTh > coord_limits[1][2]) || (posTh < coord_limits[1][0]))
	{
	  int realThIndex;

	  if (posTh > coord_limits[1][2])
	    {
	      realThIndex = n2_central+1;
	    }
	  if (posTh < coord_limits[1][0])
	    {
	      realThIndex = n2_central-1;
	    }

	  int p;
	  for (p=0; p < 4; p++)
	    {
	      intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][p] = realThIndex;
	    }
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][0] = lowerRIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][1] = lowerRIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][2] = upperRIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][3] = upperRIndex;

	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][0] = lowerPhIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][1] = upperPhIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][2] = lowerPhIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][3] = upperPhIndex;

	  w0_high = (posR - lowerR)/(upperR - lowerR);
	  w0_low = 1.0 - w0_high;
	  w1_high = 1.0;
	  w1_low = 0.0;

	}


      intersectGridWeights[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0] = w0_low*w1_low;
      intersectGridWeights[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1] = w0_low*w1_high;
      intersectGridWeights[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2] = w0_high*w1_low;
      intersectGridWeights[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][3] = w0_high*w1_high;

      intersectDistances[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng] = maxL;



      // printf("!! %e %e %e %e !!\n", w0_low, w0_high, w1_low, w1_high);
      // printf("-- %e %e %e --\n", posX, posY, posZ);



    } //end loop over angles


}


// Subroutine to calculate the interpolation weights for a square 3x3x1 cubic grid
// (Assuming symmetry in z)

// We find the 2 bounding edges for each ray hitting the boundary, and their appropriate weights to use in interpolation later on

void setupInterpWeights_cart2D(int ix, int iy, int iz, double angGridCoords[NUMANGLES][3], int intersectGridIndices[SX][SY][SZ][NUMANGLES][3][4], double intersectGridWeights[SX][SY][SZ][NUMANGLES][4], double intersectDistances[SX][SY][SZ][NUMANGLES])
{
  int probeAng;
  for (probeAng = 0; probeAng < NUMANGLES; probeAng++)
    {

      double posX, posY, posZ;
      double minL = 0., maxL = 100.0;


      //Bisect on ray to numerically find intersection location
      int iter;
      for (iter=0; iter < 50; iter++)
	{


	  posX = 1.0 + (minL+maxL)/2.0 * (-angGridCoords[probeAng][0]);
	  posY = 1.0 + (minL+maxL)/2.0 * (-angGridCoords[probeAng][1]);
	  posZ = 1.0 + (minL+maxL)/2.0 * (-angGridCoords[probeAng][2]);



	  //bisect on path length to find nearest boundary intersection
	  if ((posX > 2.0) || (posX < 0.0) || (posY > 2.0) || (posY < 0.0))
	    {
	      maxL = (minL+maxL)/2.0;
	    }
	  else
	    {
	      minL = (minL+maxL)/2.0;
	    }


	  // printf("%e %e (%e %e %e) (%e %e %e)\n", minL, maxL, posR0, posTh0, posPh0, posR, posTh, posPh);
	  // printf("%e %e (%e %e %e)\n", minL, maxL, posX, posY, posZ);


	}


      posX = 1.0 + maxL * (-angGridCoords[probeAng][0]);
      posY = 1.0 + maxL * (-angGridCoords[probeAng][1]);
      posZ = 1.0 + maxL * (-angGridCoords[probeAng][2]);

      int lowerXIndex = floor(posX), upperXIndex = lowerXIndex + 1;
      int lowerYIndex = floor(posY), upperYIndex = lowerYIndex + 1;
      int lowerZIndex = floor(posZ), upperZIndex = lowerZIndex + 1;

      //Interpolation weight components
      double w0_low, w0_high;
      double w1_low, w1_high;



      //Get interpolation weights for 4 bounding corners of cube intersection
      if ((posX > 2.0) || (posX < 0.0))
	{

	  int realXIndex;


	  if (posX > 2.0)
	    {
	      realXIndex = lowerXIndex;
	    }
	  if (posX < 0.0)
	    {
	      realXIndex = upperXIndex;
	    }

	  int p;
	  for (p=0; p < 4; p++)
	    {
	      intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][p] = realXIndex;
	      intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][p] = 1;  //Z index is always 1 (since we assume symmetry in Z)
	    }
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][0] = lowerYIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][1] = lowerYIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][2] = upperYIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][3] = upperYIndex;

	  w0_high = posY - floor(posY);
	  w0_low = 1.0 - w0_high;
	  w1_high = 0.0;
	  w1_low = 1.0;

	}
      else if ((posY > 2.0) || (posY < 0.0))
	{
	  int realYIndex;

	  if (posY > 2.0)
	    {
	      realYIndex = lowerYIndex;
	    }
	  if (posY < 0.0)
	    {
	      realYIndex = upperYIndex;
	    }

	  int p;
	  for (p=0; p < 4; p++)
	    {
	      intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][p] = realYIndex;
	      intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][p] = 1;  //Z index is always 1 (since we assume symmetry in Z)
	    }
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][0] = lowerXIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][1] = lowerXIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][2] = upperXIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][3] = upperXIndex;

	  w0_high = posX - floor(posX);
	  w0_low = 1.0 - w0_high;
	  w1_high = 0.0;
	  w1_low = 1.0;

	}
      else
	{
	  int p;
	  for (p=0; p < 4; p++)
	    {
	      intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][p] = 1.0;  //Z index is always 1 (since we assume symmetry in Z)
	    }
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][0] = lowerXIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][1] = lowerXIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][2] = upperXIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][3] = upperXIndex;

	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][0] = lowerYIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][1] = upperYIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][2] = lowerYIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][3] = upperYIndex;

	  w0_high = posX - floor(posX);
	  w0_low = 1.0 - w0_high;
	  w1_high = posY - floor(posY);
	  w1_low = 1.0 - w1_high;

	}


      intersectGridWeights[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0] = w0_low*w1_low;
      intersectGridWeights[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1] = w0_low*w1_high;
      intersectGridWeights[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2] = w0_high*w1_low;
      intersectGridWeights[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][3] = w0_high*w1_high;

      intersectDistances[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng] = maxL*GRID_SPACING;



      // printf("!! %e %e %e %e !!\n", w0_low, w0_high, w1_low, w1_high);
      // printf("-- %e %e %e --\n", posX, posY, posZ);



    } //end loop over angles


}

void setupInterpWeights_cart3D(int ix,int iy,int iz,double angGridCoords[NUMANGLES][3], int intersectGridIndices[SX][SY][SZ][NUMANGLES][3][4], double intersectGridWeights[SX][SY][SZ][NUMANGLES][4], double intersectDistances[SX][SY][SZ][NUMANGLES])
{
  int probeAng;
  for (probeAng = 0; probeAng < NUMANGLES; probeAng++)
    {

      double posX, posY, posZ;
      double minL = 0., maxL = 2.0;


      //Bisect on ray to numerically find intersection location
      int iter;
      for (iter=0; iter < 50; iter++)
	{


	  posX = 1.0 + (minL+maxL)/2.0 * (-angGridCoords[probeAng][0]);
	  posY = 1.0 + (minL+maxL)/2.0 * (-angGridCoords[probeAng][1]);
	  posZ = 1.0 + (minL+maxL)/2.0 * (-angGridCoords[probeAng][2]);



	  //bisect on path length to find nearest boundary intersection
	  if ((posX > 2.0) || (posX < 0.0) || (posY > 2.0) || (posY < 0.0) || (posZ > 2.0) || (posZ < 0.0))
	    {
	      maxL = (minL+maxL)/2.0;
	    }
	  else
	    {
	      minL = (minL+maxL)/2.0;
	    }


	  //		printf("%e %e (%e %e %e) (%e %e %e)\n", minL, maxL, posR0, posTh0, posPh0, posR, posTh, posPh);
	  //		printf("%e %e (%e %e %e)\n", minL, maxL, posX, posY, posZ);


	}


      posX = 1.0 + maxL * (-angGridCoords[probeAng][0]);
      posY = 1.0 + maxL * (-angGridCoords[probeAng][1]);
      posZ = 1.0 + maxL * (-angGridCoords[probeAng][2]);

      int lowerXIndex = floor(posX), upperXIndex = lowerXIndex + 1;
      int lowerYIndex = floor(posY), upperYIndex = lowerYIndex + 1;
      int lowerZIndex = floor(posZ), upperZIndex = lowerZIndex + 1;

      //Interpolation weight components
      double w0_low, w0_high;
      double w1_low, w1_high;


      if ((posX > 2.0) || (posX < 0.0))
	{

	  int realXIndex;


	  if (posX > 2.0)
	    {
	      realXIndex = lowerXIndex;
	    }
	  if (posX < 0.0)
	    {
	      realXIndex = upperXIndex;
	    }

	  int p;
	  for (p=0; p < 4; p++)
	    {
	      intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][p] = realXIndex;
	    }
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][0] = lowerYIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][1] = lowerYIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][2] = upperYIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][3] = upperYIndex;

	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][0] = lowerZIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][1] = upperZIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][2] = lowerZIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][3] = upperZIndex;

	  w0_high = posY - floor(posY);
	  w0_low = 1.0 - w0_high;
	  w1_high = posZ - floor(posZ);
	  w1_low = 1.0 - w1_high;
		
	}


      if ((posY > 2.0) || (posY < 0.0))
	{
	  int realYIndex;

	  if (posY > 2.0)
	    {
	      realYIndex = lowerYIndex;
	    }
	  if (posY < 0.0)
	    {
	      realYIndex = upperYIndex;
	    }

	  int p;
	  for (p=0; p < 4; p++)
	    {
	      intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][p] = realYIndex;
	    }
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][0] = lowerXIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][1] = lowerXIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][2] = upperXIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][3] = upperXIndex;

	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][0] = lowerZIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][1] = upperZIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][2] = lowerZIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][3] = upperZIndex;

	  w0_high = posX - floor(posX);
	  w0_low = 1.0 - w0_high;
	  w1_high = posZ - floor(posZ);
	  w1_low = 1.0 - w1_high;
		
	}


      if ((posZ > 2.0) || (posZ < 0.0))
	{
	  int realZIndex;

	  if (posZ > 2.0)
	    {
	      realZIndex = lowerZIndex;
	    }
	  if (posZ < 0.0)
	    {
	      realZIndex = upperZIndex;
	    }

	  int p;
	  for (p=0; p < 4; p++)
	    {
	      intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][p] = realZIndex;
	    }
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][0] = lowerXIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][1] = lowerXIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][2] = upperXIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][3] = upperXIndex;

	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][0] = lowerYIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][1] = upperYIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][2] = lowerYIndex;
	  intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][3] = upperYIndex;

	  w0_high = posX - floor(posX);
	  w0_low = 1.0 - w0_high;
	  w1_high = posY - floor(posY);
	  w1_low = 1.0 - w1_high;
		
	}


      intersectGridWeights[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0] = w0_low*w1_low;
      intersectGridWeights[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1] = w0_low*w1_high;
      intersectGridWeights[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2] = w0_high*w1_low;
      intersectGridWeights[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][3] = w0_high*w1_high;

      intersectDistances[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng] = maxL*GRID_SPACING;



      //	printf("!! %e %e %e %e !!\n", w0_low, w0_high, w1_low, w1_high);
      //	printf("-- %e %e %e --\n", posX, posY, posZ);



    } //end loop over angles



}


double I_Solve(double S0, double S1, double I1, double dtau)
{
  //Notation is: point 0 = local point, point 1 = far point, I1 = initial intensity at far point
  if (dtau > 1.0e-3)
    {
      double exptau=exp(-dtau);
      return ((exptau*(S0 - S1*(1.0+dtau)) + S0*(dtau - 1.0) + S1)/dtau + I1*exptau);
    }
  else if (dtau > 0)
    {
      return (S1*dtau + I1*(1.0-dtau));
    }
  else
    {
      return I1;
    }
}


//Subroutine for angle rotation, making use of quaternion rotation matrix

// CAVEAT for case of 180 degrees rotation: there are multiple allowed rotation planes, so we just pick an arbitrary allowed one

void calc_rot_M(double vstart[3], double vfinal[3], double rotM[3][3])
{
	double vstart_norm = sqrt(vstart[0]*vstart[0] + vstart[1]*vstart[1] + vstart[2]*vstart[2]);
	double vfinal_norm = sqrt(vfinal[0]*vfinal[0] + vfinal[1]*vfinal[1] + vfinal[2]*vfinal[2]);

	double nx, ny, nz;  //quaternion imaginary components


	vstart[0] = vstart[0]/vstart_norm;
	vstart[1] = vstart[1]/vstart_norm;
	vstart[2] = vstart[2]/vstart_norm;
	vfinal[0] = vfinal[0]/vfinal_norm;
	vfinal[1] = vfinal[1]/vfinal_norm;
	vfinal[2] = vfinal[2]/vfinal_norm;

	double cosangle = vstart[0]*vfinal[0] + vstart[1]*vfinal[1] + vstart[2]*vfinal[2];
	double sinangle = sqrt(1-cosangle*cosangle);


	//printf("%e %e %e\n", vfinal[0], vfinal[1], vfinal[2]);
//	printf("%e %e\n", cosangle, sinangle);



/*
	vec[0]=vstart[0];
	vec[1]=vstart[1];
	vec[2]=vstart[2];
*/

//	vec[0]=1., vec[1]=1.0, vec[2]=1.0;


	
	double nnorm = sinangle;

	if (sinangle!=0.)
	{
		//Cross product to find normal vector for rotation
		nx = vstart[1]*vfinal[2] - vfinal[1]*vstart[2];
		ny = vstart[2]*vfinal[0] - vfinal[2]*vstart[0];
		nz = vstart[0]*vfinal[1] - vfinal[0]*vstart[1];

		nx = nx/nnorm;
		ny = ny/nnorm;
		nz = nz/nnorm;
	}
	else
	{

		if ((vstart[0] != 0) || (vstart[1] != 0))
		{
			nx = -vstart[1];
			ny = vstart[0];
			nz = 0.0;
		}
		else
		{
			nx = 0.0;
			ny = -vstart[2];
			nz = vstart[1];
		}

		nnorm = sqrt(nx*nx + ny*ny + nz*nz);

		nx = nx/nnorm;
		ny = ny/nnorm;
		nz = nz/nnorm;
	}
	

	//printf("%e %e %e\n", nx, ny, nz);	



	//Quaternion derived rotation matrix
	rotM[0][0] = cosangle + nx*nx*(1-cosangle);
	rotM[0][1] = nx*ny*(1-cosangle) - nz*sinangle;
	rotM[0][2] = nx*nz*(1-cosangle) + ny*sinangle;

	rotM[1][0] = ny*nx*(1-cosangle) + nz*sinangle;
	rotM[1][1] = cosangle + ny*ny*(1-cosangle);
	rotM[1][2] = ny*nz*(1-cosangle) - nx*sinangle;

	rotM[2][0] = nz*nx*(1-cosangle) - ny*sinangle;
	rotM[2][1] = nz*ny*(1-cosangle) + nx*sinangle;
	rotM[2][2] = cosangle + nz*nz*(1-cosangle);



}





// CREATE our BSP angle lookup decision tree
// -------------------------------------------
// Idea is to split angle grid into two equal pieces using a plane aligned with one of the coordinate axes
//
// Decision tree is created recursively, and we cycle the orientation of the splitting plane
//
// The pivot point is taken as median value of either x,y,or z coordinates
// (iter specifies which coordinate we take bisection of, 0=x, 1=y, 2=z, cyclic etc...)
void splitAngGrid(int numAvailAnglesInit, int angGridIndexInit[NUMANGLES][3], int iter, double angGridCoords[NUMANGLES][3], struct bsptree **node)
{
		short angGridBitmaskLow[NUMANGLES];
		short angGridBitmaskHigh[NUMANGLES];

		int numAvailAnglesLow = NUMANGLES;
		int angGridIndexAvailLow[NUMANGLES][3];
		int numAvailAnglesHigh = NUMANGLES;
		int angGridIndexAvailHigh[NUMANGLES][3];

		int i,p,q;
		int sortIndex = iter%3; //cycle through x,y,z coordinates



		if (numAvailAnglesInit <= 0)
		{
			return;
		}	



	
	//Loop through and bisect available angle bins



		numAvailAnglesLow = numAvailAnglesInit/2 + 1;
		numAvailAnglesHigh = numAvailAnglesInit - numAvailAnglesLow;

		
		//Initialize bitmask
		for (i=0; i < NUMANGLES; i++)
		{
			angGridBitmaskLow[i] = 0;
			angGridBitmaskHigh[i] = 0;
		}		


		//Turn on bitmasks
		for (i=0; i < numAvailAnglesInit; i++)
		{

//			printf("WOW %d ", i);

			if (i < numAvailAnglesLow)
			{
				angGridBitmaskLow[angGridIndexInit[i][sortIndex]] = 1;
				angGridIndexAvailLow[i][sortIndex] = angGridIndexInit[i][sortIndex];
			}
			else
			{
				angGridBitmaskHigh[angGridIndexInit[i][sortIndex]] = 1;
				angGridIndexAvailHigh[i - numAvailAnglesLow][sortIndex] = angGridIndexInit[i][sortIndex];
			}

		}




		int count2Low = 0, count2High = 0;

		for (p=0; p < 3; p++)
		{
			if (p == sortIndex)
			{
				continue;
			}


			count2Low = 0, count2High = 0;

			for (i=0; i < numAvailAnglesInit; i++)
			{
				if ((angGridBitmaskLow[angGridIndexInit[i][p]] == 1) && (angGridIndexInit[i][p] != angGridIndexAvailLow[numAvailAnglesLow-1][sortIndex]))
				{
					angGridIndexAvailLow[count2Low][p] = angGridIndexInit[i][p];
					count2Low++;
				}


				if (angGridBitmaskHigh[angGridIndexInit[i][p]] == 1)
				{
					angGridIndexAvailHigh[count2High][p] = angGridIndexInit[i][p];
					count2High++;
				}
			}


		}




		//Take pivot point from last entry of lower list

		(*node) = create_node(angGridIndexAvailLow[numAvailAnglesLow-1][sortIndex], sortIndex);
		splitAngGrid(numAvailAnglesLow-1, angGridIndexAvailLow, iter+1, angGridCoords, &((*node)->lower));
		splitAngGrid(numAvailAnglesHigh, angGridIndexAvailHigh, iter+1, angGridCoords, &((*node)->upper));

}







//Same as above, but now acting on the dual angle grid

void splitDualAngGrid(int numAvailAnglesInit, int angGridIndexInit[NUMDUALANGLES][3], int iter, double angGridCoords[NUMDUALANGLES][3], struct bsptree **node)
{
		short angGridBitmaskLow[NUMDUALANGLES];
		short angGridBitmaskHigh[NUMDUALANGLES];

		int numAvailAnglesLow = NUMDUALANGLES;
		int angGridIndexAvailLow[NUMDUALANGLES][3];
		int numAvailAnglesHigh = NUMDUALANGLES;
		int angGridIndexAvailHigh[NUMDUALANGLES][3];

		int i,p,q;
		int sortIndex = iter%3; //cycle through x,y,z coordinates





		if (numAvailAnglesInit <= 0)
		{
			return;
		}	



	
	//Loop through and bisect available angle bins



		numAvailAnglesLow = numAvailAnglesInit/2 + 1;
		numAvailAnglesHigh = numAvailAnglesInit - numAvailAnglesLow;

		
		//Initialize bitmask
		for (i=0; i < NUMDUALANGLES; i++)
		{
			angGridBitmaskLow[i] = 0;
			angGridBitmaskHigh[i] = 0;
		}		


		//Turn on bitmasks
		for (i=0; i < numAvailAnglesInit; i++)
		{

//			printf("WOW %d ", i);

			if (i < numAvailAnglesLow)
			{
				angGridBitmaskLow[angGridIndexInit[i][sortIndex]] = 1;
				angGridIndexAvailLow[i][sortIndex] = angGridIndexInit[i][sortIndex];
			}
			else
			{
				angGridBitmaskHigh[angGridIndexInit[i][sortIndex]] = 1;
				angGridIndexAvailHigh[i - numAvailAnglesLow][sortIndex] = angGridIndexInit[i][sortIndex];
			}

		}






		int count2Low = 0, count2High = 0;

		for (p=0; p < 3; p++)
		{
			if (p == sortIndex)
			{
				continue;
			}


			count2Low = 0, count2High = 0;

			for (i=0; i < numAvailAnglesInit; i++)
			{
				if ((angGridBitmaskLow[angGridIndexInit[i][p]] == 1) && (angGridIndexInit[i][p] != angGridIndexAvailLow[numAvailAnglesLow-1][sortIndex]))
				{
					angGridIndexAvailLow[count2Low][p] = angGridIndexInit[i][p];
					count2Low++;
				}


				if (angGridBitmaskHigh[angGridIndexInit[i][p]] == 1)
				{
					angGridIndexAvailHigh[count2High][p] = angGridIndexInit[i][p];
					count2High++;
				}
			}


		}





		//Take pivot point from last entry of lower list

		(*node) = create_node(angGridIndexAvailLow[numAvailAnglesLow-1][sortIndex], sortIndex);
		splitDualAngGrid(numAvailAnglesLow-1, angGridIndexAvailLow, iter+1, angGridCoords, &((*node)->lower));
		splitDualAngGrid(numAvailAnglesHigh, angGridIndexAvailHigh, iter+1, angGridCoords, &((*node)->upper));

}







//Traverse our BSP tree to find nearest angle
//
//Tree traversal occurs recursively
//NOTE: Refer to k-d tree, nearest neighbor algorithm as described on wikipedia

void bspGetNearestNeighbor(double targetAng[3], double angGridCoords[NUMANGLES][3], struct bsptree *bspCurrentLoc, double *bestDistance, int *bestIndex)
{
		int currentAngIndex = bspCurrentLoc->angIndex;
		int sortIndex = bspCurrentLoc->iter;
		short doublecheck = 0, firstDecision = 0;


		double distance2 = (targetAng[0] - angGridCoords[currentAngIndex][0])*(targetAng[0] - angGridCoords[currentAngIndex][0])
				+ (targetAng[1] - angGridCoords[currentAngIndex][1])*(targetAng[1] - angGridCoords[currentAngIndex][1])
				+ (targetAng[2] - angGridCoords[currentAngIndex][2])*(targetAng[2] - angGridCoords[currentAngIndex][2]);

		if (distance2 < (*bestDistance))
		{
			*bestDistance = distance2;
			*bestIndex = currentAngIndex;
		}


//		printf("%d | %d | %e - %e\n", currentAngIndex, sortIndex, angGridCoords[currentAngIndex][sortIndex], targetAng[sortIndex]);
		



		//go lower
		if (targetAng[sortIndex] < angGridCoords[currentAngIndex][sortIndex])
//		if (1)
		{
//			printf("smaller\n");
			firstDecision = 0;

			if (bspCurrentLoc->lower != 0) 
			{
				bspGetNearestNeighbor(targetAng, angGridCoords, bspCurrentLoc->lower, bestDistance, bestIndex);

			}


		}
		else //go upper
//		if (1)
		{
//			printf("larger\n");
			firstDecision = 1;

			if (bspCurrentLoc->upper != 0) 
			{
				bspGetNearestNeighbor(targetAng, angGridCoords, bspCurrentLoc->upper, bestDistance, bestIndex);

			}


		}




		//Check if both sides of partition need to be computed
		if ((*bestDistance) > (targetAng[sortIndex] - angGridCoords[currentAngIndex][sortIndex]) * (targetAng[sortIndex] - angGridCoords[currentAngIndex][sortIndex]))
//		if (1)
		{

			if (firstDecision == 0)  //If we first looked at lower, now switch and evaluate upper
			{
				if (bspCurrentLoc->upper != 0) 
				{
					bspGetNearestNeighbor(targetAng, angGridCoords, bspCurrentLoc->upper, bestDistance, bestIndex);
				}
			
			}
			else  //Vice-versa
			{
				if (bspCurrentLoc->lower != 0) 
				{
					bspGetNearestNeighbor(targetAng, angGridCoords, bspCurrentLoc->lower, bestDistance, bestIndex);
				}
			}
		}


	
	return;

}




//Same as above, but acting on dual angle grid

void bspGetNearestDualNeighbor(double targetAng[3], double angGridCoords[NUMDUALANGLES][3], struct bsptree *bspCurrentLoc, double *bestDistance, int *bestIndex)
{
		int currentAngIndex = bspCurrentLoc->angIndex;
		int sortIndex = bspCurrentLoc->iter;
		short doublecheck = 0, firstDecision = 0;


		double distance2 = (targetAng[0] - angGridCoords[currentAngIndex][0])*(targetAng[0] - angGridCoords[currentAngIndex][0])
				+ (targetAng[1] - angGridCoords[currentAngIndex][1])*(targetAng[1] - angGridCoords[currentAngIndex][1])
				+ (targetAng[2] - angGridCoords[currentAngIndex][2])*(targetAng[2] - angGridCoords[currentAngIndex][2]);

		if (distance2 < (*bestDistance))
		{
			*bestDistance = distance2;
			*bestIndex = currentAngIndex;
		}


//		printf("%d | %d | %e - %e\n", currentAngIndex, sortIndex, angGridCoords[currentAngIndex][sortIndex], targetAng[sortIndex]);
		



		//go lower
		if (targetAng[sortIndex] < angGridCoords[currentAngIndex][sortIndex])
//		if (1)
		{
//			printf("smaller\n");
			firstDecision = 0;

			if (bspCurrentLoc->lower != 0) 
			{
				bspGetNearestDualNeighbor(targetAng, angGridCoords, bspCurrentLoc->lower, bestDistance, bestIndex);

			}


		}
		else //go upper
//		if (1)
		{
//			printf("larger\n");
			firstDecision = 1;

			if (bspCurrentLoc->upper != 0) 
			{
				bspGetNearestDualNeighbor(targetAng, angGridCoords, bspCurrentLoc->upper, bestDistance, bestIndex);

			}


		}




		//Check if both sides of partition need to be computed
		if ((*bestDistance) > (targetAng[sortIndex] - angGridCoords[currentAngIndex][sortIndex]) * (targetAng[sortIndex] - angGridCoords[currentAngIndex][sortIndex]))
//		if (1)
		{

			if (firstDecision == 0)  //If we first looked at lower, now switch and evaluate upper
			{
				if (bspCurrentLoc->upper != 0) 
				{
					bspGetNearestDualNeighbor(targetAng, angGridCoords, bspCurrentLoc->upper, bestDistance, bestIndex);
				}
			
			}
			else  //Vice-versa
			{
				if (bspCurrentLoc->lower != 0) 
				{
					bspGetNearestDualNeighbor(targetAng, angGridCoords, bspCurrentLoc->lower, bestDistance, bestIndex);
				}
			}
		}


	
	return;

}









//Do a linear search to find nearest angle
int get_angIndex(double targetAng[3], double angGridCoords[NUMANGLES][3])
{
	double dx, dy, dz, dtot2;
	//double dmin = 1.0e10;
	double dmax=-1.0e10;
	int imin = 0;
	int i;

	//Find closest grid cell
	for (i = 0; i < NUMANGLES; i++)
	{

		dtot2=targetAng[0]*angGridCoords[i][0] + targetAng[1]*angGridCoords[i][1] + targetAng[2]*angGridCoords[i][2];

		if (dtot2 > dmax)
		{
			dmax = dtot2;
			imin = i;
		}
	}

	return imin;
}


//Do a linear search to find nearest dual angle
int get_angDualIndex(double targetAng[3], double angGridCoords[NUMDUALANGLES][3])
{
	double dx, dy, dz, dtot2;
	//double dmin = 1.0e10;
	double dmax=-1.0e10;
	int imin = 0;
	int i;

	//Find closest grid cell
	for (i = 0; i < NUMDUALANGLES; i++)
	{

		dtot2=targetAng[0]*angGridCoords[i][0] + targetAng[1]*angGridCoords[i][1] + targetAng[2]*angGridCoords[i][2];

		if (dtot2 > dmax)
		{
			dmax = dtot2;
			imin = i;
		}
	}

	return imin;
}






//Express vector targetAng as linear combination of 3 input angles, as specified by indices angIndex[0]-[2]
void linComb(double targetAng[3], double angGridCoords[NUMANGLES][3], int angIndex[3], double interp_coeff[3])
{
	double M[3][3];
	double inv_M[3][3];
	double det = 0.;

	double lc1, lc2, lc3;



	//specify initial matrix of 3 independent basis vectors
	M[0][0] = angGridCoords[angIndex[0]][0];
	M[0][1] = angGridCoords[angIndex[1]][0];
	M[0][2] = angGridCoords[angIndex[2]][0];

	M[1][0] = angGridCoords[angIndex[0]][1];
	M[1][1] = angGridCoords[angIndex[1]][1];
	M[1][2] = angGridCoords[angIndex[2]][1];

	M[2][0] = angGridCoords[angIndex[0]][2];
	M[2][1] = angGridCoords[angIndex[1]][2];
	M[2][2] = angGridCoords[angIndex[2]][2];






	//invert Matrix
	det = (M[0][0]*M[1][1]*M[2][2] + M[0][1]*M[1][2]*M[2][0] + M[0][2]*M[1][0]*M[2][1] - M[2][0]*M[1][1]*M[0][2] - M[2][1]*M[1][2]*M[0][0] - M[2][2]*M[1][0]*M[0][1]); 

	//printf("det=%e\n",det);

	if (det==0)
	{
		printf("ERROR: colinear basis vectors used! Det = 0!\n");
		return;
	}

	inv_M[0][0] = (M[1][1]*M[2][2] - M[2][1]*M[1][2]);
	inv_M[0][1] = (M[0][2]*M[2][1] - M[2][2]*M[0][1]);
	inv_M[0][2] = (M[0][1]*M[1][2] - M[1][1]*M[0][2]);

	inv_M[1][0] = (M[1][2]*M[2][0] - M[2][2]*M[1][0]);
	inv_M[1][1] = (M[0][0]*M[2][2] - M[2][0]*M[0][2]);
	inv_M[1][2] = (M[0][2]*M[1][0] - M[1][2]*M[0][0]);

	inv_M[2][0] = (M[1][0]*M[2][1] - M[2][0]*M[1][1]);
	inv_M[2][1] = (M[0][1]*M[2][0] - M[2][1]*M[0][0]);
	inv_M[2][2] = (M[0][0]*M[1][1] - M[1][0]*M[0][1]);

/*
	for (int i=0; i < 3; i++)
	{
		for (int j=0; j< 3; j++)
		{
			printf("%e ", inv_M[i][j]);
		}
		printf("\n");
	}
*/
	lc1 = (inv_M[0][0]*targetAng[0] + inv_M[0][1]*targetAng[1] + inv_M[0][2]*targetAng[2])/det;
	lc2 = (inv_M[1][0]*targetAng[0] + inv_M[1][1]*targetAng[1] + inv_M[1][2]*targetAng[2])/det;
	lc3 = (inv_M[2][0]*targetAng[0] + inv_M[2][1]*targetAng[1] + inv_M[2][2]*targetAng[2])/det;

	if (lc1 < 0) {lc1=0.;}
	if (lc2 < 0) {lc2=0.;}
	if (lc3 < 0) {lc3=0.;}



/*
	// Conserve Intensity
	interp_coeff[0]=lc1/(lc1+lc2+lc3);
	interp_coeff[1]=lc2/(lc1+lc2+lc3);
	interp_coeff[2]=lc3/(lc1+lc2+lc3);
*/




	// Conserve Flux
	interp_coeff[0]=lc1;
	interp_coeff[1]=lc2;
	interp_coeff[2]=lc3;



/*	if (*c1 + *c2 + *c3 > 1.01)
	{
		printf("ERROR: weight mismatch -- %e %e %e!\n", *c1, *c2, *c3);
	} 
*/
	return;
}









//	Initialization subrouting for constructing BSP trees
void initAngIndex(double angGridCoords[NUMANGLES][3], double angDualGridCoords[NUMDUALANGLES][3], int angGridIndexSort[NUMANGLES][3], int angDualGridIndexSort[NUMANGLES][3])
{
	int i,j,p;

// Generate a sorted angle grid locator array

	//initialize
	for (i=0; i < NUMANGLES; i++)
	{
		for (p=0; p < 3; p++)
		{
			angGridIndexSort[i][p] = i;
		}
	}


	//Sort angle index arrays to be increasing in the 3 xyz coordinates
	for (i=0; i < NUMANGLES; i++)
	{
		for (j=i+1; j < NUMANGLES; j++)
		{
			for (p=0; p < 3; p++)
			{
				// If coordinate is out of order
				if (angGridCoords[angGridIndexSort[i][p]][p] > angGridCoords[angGridIndexSort[j][p]][p])
				{
					//swap
					int temp = angGridIndexSort[j][p];
					angGridIndexSort[j][p] = angGridIndexSort[i][p];
					angGridIndexSort[i][p] = temp;
				}
			}
		}
	}	



	//Same exercise for Dual Angle Grid

	//initialize
	for (i=0; i < NUMDUALANGLES; i++)
	{
		for (p=0; p < 3; p++)
		{
			angDualGridIndexSort[i][p] = i;
		}
	}


	//Sort angle index arrays to be increasing in the 3 xyz coordinates
	for (i=0; i < NUMDUALANGLES; i++)
	{
		for (j=i+1; j < NUMDUALANGLES; j++)
		{
			for (p=0; p < 3; p++)
			{
				// If coordinate is out of order
				if (angDualGridCoords[angDualGridIndexSort[i][p]][p] > angDualGridCoords[angDualGridIndexSort[j][p]][p])
				{
					//swap
					int temp = angDualGridIndexSort[j][p];
					angDualGridIndexSort[j][p] = angDualGridIndexSort[i][p];
					angDualGridIndexSort[i][p] = temp;
				}
			}
		}
	}	

}


void transformI_stretch(int ix, int iy,int iz,double I_return[NUMANGLES], double M1_input[5])
{
  double F_final[3],Efinal,Efinalrad;

  F_final[0]=M1_input[1]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][0][0] 
    +M1_input[2]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][1][0]
    +M1_input[3]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][2][0];

  F_final[1]=M1_input[1]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][0][1] 
    +M1_input[2]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][1][1]
    +M1_input[3]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][2][1];

  F_final[2]=M1_input[1]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][0][2] 
    +M1_input[2]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][1][2]
    +M1_input[3]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][2][2];

  Efinal=M1_input[0];
  Efinalrad=M1_input[4];

  int i,j,p,l;
  double I_start[NUMANGLES];
  double res[3];

  double F_start[3], F_start_norm[3], F_final_norm[3], Fmag_start, Fmag_final, F_stretch;
  double fStart, F_Start_forward, F_Start_backward, fFinal, Estart;
  double stretchFactor;

  double rotM[3][3];
  double transformAng[NUMANGLES][3];


  for (i=0; i < NUMANGLES; i++)
    {
      I_start[i] = I_return[i];
      //printf("%e\n", I_start[i]);
      I_return[i] = 0.;
    }

  //printf("%e %e %e\n", F_final[0], F_final[1], F_final[2]);

  Estart=0.;
  //Calculate net flux direction
  for (l=0; l < 3; l++)
    {
      F_start[l] = 0.;
    }

  for (p=0; p < NUMANGLES; p++)
    {
      //if(p<10) printf("%e\n", I_start[p]);
      for (l=0; l < 3; l++)
	{
	  F_start[l] += I_start[p]*angGridCoords[p][l];
	}
      Estart += I_start[p];
    }

  //	printf("start: %e %e %e || %e\n", F_start[0], F_start[1], F_start[2], Estat);
  //	printf("F start: %e %e %e || %e\n", F_start[0], F_start[1], F_start[2], sqrt(F_start[0]*F_start[0] + F_start[1]*F_start[1] + F_start[2]*F_start[2]));




  Fmag_start = sqrt(F_start[0]*F_start[0] + F_start[1]*F_start[1] + F_start[2]*F_start[2]);
  Fmag_final = sqrt(F_final[0]*F_final[0] + F_final[1]*F_final[1] + F_final[2]*F_final[2]);

  if (Fmag_final/Efinal < 1.0e-10)
    {
      for (p=0; p<NUMANGLES; p++)
	{
	  I_return[p] = Efinal/NUMANGLES;
	}
      //	    printf("%e %e\n",I_start[0],I_return[0]);
      return;
    }

  for (l=0; l < 3; l++)
    {
      F_start_norm[l] = F_start[l]/Fmag_start;
      F_final_norm[l] = F_final[l]/Fmag_final;
    }

  /*
  if (Fmag_start/Estart < 1.0e-2)
    {
      Fmag_start = 1.0;   //Careful about renormalizing data when Fmag_start = 0

      F_start_norm[0]=1.0;
      F_start_norm[1]=0.0;
      F_start_norm[2]=0.0;
    }
  */

  if (Fmag_start/Estart < 1.0e-2)
    {
      //Calculate starting flux after reflecting half the rays about origin

      F_start_norm[0]=0.0;
      F_start_norm[1]=0.0;
      F_start_norm[2]=0.0;

      for (l=0; l < NUMANGLES; l++)
	{
	  if (angGridCoords[l][0] < 0) //reflect half the rays about origin
	    {
	      for (p=0; p < 3; p++)
		{
		  F_start_norm[p] += I_start[l]*(-angGridCoords[l][p]);
		}
	    }
	  else
	    {
	      for (p=0; p < 3; p++)
		{
		  F_start_norm[p] += I_start[l]*angGridCoords[l][p];
		}
	    }
	}

      double tempMag = sqrt(F_start_norm[0]*F_start_norm[0] +
			    F_start_norm[1]*F_start_norm[1] + F_start_norm[2]*F_start_norm[2]);

      for (p=0; p < 3; p++)
	{
	  F_start_norm[p] = F_start_norm[p]/tempMag;

	}

    }

  
  F_Start_forward = 0.;
  F_Start_backward = 0.;
  double cosang[NUMANGLES];

  for (p=0; p < NUMANGLES; p++)
    {


      cosang[p] = angGridCoords[p][0]*F_start_norm[0] + angGridCoords[p][1]*F_start_norm[1] + angGridCoords[p][2]*F_start_norm[2];

      if (cosang[p] > 0.)
	{
	  F_Start_forward += I_start[p]*cosang[p]; //positive fluxes in direction of net flux
	}
      else
	{
	 
	  F_Start_backward += I_start[p]*cosang[p]; //negative fluxes "" ""
	}
    }



  fStart = Fmag_start/Estart;
  fFinal = Fmag_final/Efinal;

  struct solverarg args;

  args.F_Start_forward=F_Start_forward;
  args.F_Start_backward=F_Start_backward;
  args.Intensities = &I_start[0];
  args.F_final_norm = &F_final_norm[0];
  args.fFinal = fFinal;
  args.cosang = &cosang[0];
  
  stretchFactor = calc_stretchFactor(&args);
  //exit(1);

  F_stretch = fabs(stretchFactor*F_Start_forward + F_Start_backward/stretchFactor);

  //printf("%e %e %e\n", stretchFactor, F_Start_forward, F_Start_backward);


  //stretchFactor = sqrt(fabs((fFinal*fFinal - fFinal*fFinal*fStart*fStart)/(fStart*fStart - fFinal*fFinal*fStart*fStart)));


  //printf("Fstart = %e %e %e\n", F_start[0], F_start[1], F_start[2]);

  //printf("SF = %e | Estart = %e, Fstart = %e, Ffinal = %e | F/E = %e\n", stretchFactor,Estart,  Fmag_start, Fmag_final, Fmag_start/Estart);

  if(!isfinite(stretchFactor) || stretchFactor<1.e-15) 
    {
	    
      for (i=0; i < NUMANGLES; i++)

	I_return[i] = I_start[i];
	

      return;
	   
    }

  //Calculate rotation matrix
  calc_rot_M(F_start_norm, F_final_norm, rotM);

  //	printf("start: %e %e %e\nfinish: %e %e %e\n", F_norm[0], F_norm[1], F_norm[2], F_final[0], F_final[1], F_final[2]);




  int probeAng;
  for (probeAng=0; probeAng < NUMANGLES; probeAng++)
    {

      //Apply Rotation Matrix to some initial intesity distribution

      for (i=0; i < 3; i++)
	{
	  transformAng[probeAng][i]=0;
	  for (j=0; j < 3; j++)
	    {
	      transformAng[probeAng][i] += rotM[i][j]*angGridCoords[probeAng][j];
	    }
	}




      //		printf("%e %e %e\n", angGridCoords[probeAng][0], angGridCoords[probeAng][1], angGridCoords[probeAng][2]);
      //		printf("%e %e %e\n", transformAng[probeAng][0], transformAng[probeAng][1], transformAng[probeAng][2]);


      //Keep track of component of E^2 that remains unchanged


      //Stretch initial intensity distrubution parallel to Ffinal

      //First, calculate parallel component of I
      double n_parallel[3], n_perp[3], n_final[3];
      double dotprod = 0., n_norm = 0.;


      for (l=0; l < 3; l++)
	{
	  dotprod += transformAng[probeAng][l]*F_final_norm[l];
	}


      for (l=0; l < 3; l++)
	{
	  n_parallel[l] = dotprod*F_final_norm[l];
	  n_perp[l] = transformAng[probeAng][l] - n_parallel[l];

	  if (dotprod > 0)
	    {
	      n_parallel[l] = n_parallel[l]*stretchFactor;
	    }
	  else
	    {
	      n_parallel[l] = n_parallel[l]/stretchFactor;
	    }
	  n_final[l] = n_perp[l] + n_parallel[l];

	}



      n_norm = sqrt(n_final[0]*n_final[0] + n_final[1]*n_final[1] + n_final[2]*n_final[2]);


      //		printf("%e\n", n_norm);

      for (l=0; l < 3; l++)
	{
	  n_final[l] = n_final[l]/n_norm;
	}


      //		printf("%e %e %e\n", n_final[0], n_final[1], n_final[2]);


      //		printf("n final = %e %e %e | %e\n", n_final[0], n_final[1], n_final[2], sqrt(n_final[0]*n_final[0] + n_final[1]*n_final[1] + n_final[2]*n_final[2]));





      //		struct bsptree *bspCurrentLoc = angDualGridRoot;
      double bestDistance = 1.0e10;
      int bestIndex = 0;


      bspGetNearestDualNeighbor(n_final, angDualGridCoords, angDualGridRoot, &bestDistance, &bestIndex);
      //bestIndex = get_angDualIndex(n_final, angDualGridCoords);

      int angNeighborIndex[3];
      double interpCoeffs[3];

      for (l=0; l < 3; l++)
	{
	  angNeighborIndex[l] = dualAdjacency[bestIndex][l];
	}


      //		printf("Best Angle: %d, %e %e %e\n", bestIndex, angDualGridCoords[bestIndex][0], angDualGridCoords[bestIndex][1], angDualGridCoords[bestIndex][2]);


      linComb(n_final, angGridCoords, angNeighborIndex, interpCoeffs);

      for (p=0; p < 3; p++)
	{
	  //			printf("%d %d | %e %e %e \n", bestIndex, angNeighborIndex[p], n_norm, I_start[probeAng], interpCoeffs[p]);

	  I_return[angNeighborIndex[p]] += n_norm*I_start[probeAng]*interpCoeffs[p];
	}

      //		printf("\n");

    }



  double rescaleFactor = Fmag_final/F_stretch;
  //Renormalize to get correct fluxes
  for (p=0; p < NUMANGLES; p++)
    {
      I_return[p] = I_return[p]*rescaleFactor;
    }
}


//decomposes M1 beam into intensities, uses energy densities and fluxes as the input
void ZERO_decomposeM1(int ix, int iy, int iz,double M1_Data[5], double I_return[NUMANGLES])
{
  int probeAng;
  if(M1_Data[0]<SMALL)
    {
      for (probeAng=0; probeAng < NUMANGLES; probeAng++)
	I_return[probeAng]=0.;
      return;
    }

  double F_final[3];

  F_final[0]=M1_Data[1]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][0][0] 
    +M1_Data[2]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][1][0]
    +M1_Data[3]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][2][0];

  F_final[1]=M1_Data[1]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][0][1] 
    +M1_Data[2]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][1][1]
    +M1_Data[3]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][2][1];

  F_final[2]=M1_Data[1]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][0][2] 
    +M1_Data[2]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][1][2]
    +M1_Data[3]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][2][2];


  double fmag = sqrt(F_final[0]*F_final[0] + 
		     F_final[1]*F_final[1] + 
		     F_final[2]*F_final[2]);

  double ff = fmag / M1_Data[0]; //F/Elab using input argument

  //printf("ff: %f\n",ff);
  double beta;

  if(ff<1.e-2) 
    beta=3.*ff/4.;
  else
    beta=(4.-sqrt(16.-12.*ff*ff))/2./ff;
	    
  double gamma2=1./(1.-beta*beta);

  double Erad=M1_Data[0]/(4./3.*gamma2-1./3.);

  if(fmag<SMALL)
    fmag=1.;

  double f_norm[3];
  f_norm[0] = F_final[0]/fmag;
  f_norm[1] = F_final[1]/fmag;
  f_norm[2] = F_final[2]/fmag;

 
  for (probeAng=0; probeAng < NUMANGLES; probeAng++)
    {

      double mu= f_norm[0]*angGridCoords[probeAng][0] + f_norm[1]*angGridCoords[probeAng][1] + f_norm[2]*angGridCoords[probeAng][2];
		      
      if(beta<1.e-5)
	mu=1.;

      double bm=1.-beta*mu;
     
      I_return[probeAng] = Erad/NUMANGLES/bm/bm/bm/bm/gamma2/gamma2;
    }

  return;
}




void ZERO_shortCharI(int ix, int iy, int iz,double delta_t, double I_Data[3][3][3][NUMANGLES], double source_Data[3][3][3][4],
		     double I_return[NUMANGLES],int verbose)
{
  double S[3][3][3];  //radiative source function

  double I_ray[NUMANGLES];  //intensity field evaluated at the center, for time independent problem
  double I_time[NUMANGLES];  //"" "", for time-dependent problem


  //Fill out source function for cube
  int i,j,k;
  for (i=0; i<3; i++)
    {
      for (j=0; j<3; j++)
	{
	  for (k=0; k<3; k++)
	    {
	      double alpha=source_Data[i][j][k][2], sigma=source_Data[i][j][k][3];
	      double eps;

	      if (alpha + sigma > 0)
		{
		  eps = alpha/(alpha+sigma);
		}
	      else
		{
		  eps = 1.0;
		}

	      S[i][j][k]=eps*pow(source_Data[i][j][k][0],4)*STEFAN_BOLTZMANN/PI + (1.0-eps)*source_Data[i][j][k][1]*LIGHT_C/4.0/PI;
	    }
	}
    }

  //Calculate solution to RT, given intensity information from previous iteration
  int probeAng;
  for (probeAng = 0; probeAng < NUMANGLES; probeAng++)
    {


      //First, get interpolated quantities along ray intersection boundary
      double interp_I[NUMANGLES];
      double interp_S = 0.;

      //initialize variables
      int p;
      for (p=0; p < NUMANGLES; p++)
	{
	  interp_I[p] = 0.;
	}



      for (p=0; p < 4; p++)
	{
	  int intersect_i = intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][0][p];
	  int intersect_j = intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][1][p];
	  int intersect_k = intersectGridIndices[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][2][p];


	  int q;
	  for (q=0 ; q<NUMANGLES; q++)
	    {
	      interp_I[q] += I_Data[intersect_i][intersect_j][intersect_k][q]*intersectGridWeights[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][p];

	      
	      if(interp_I[q]<0.) 
		{
		  printf("neg interp_I %e %e\n",
			 I_Data[intersect_i][intersect_j][intersect_k][q],
			 intersectGridWeights[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][p]);
		  //exit(1);
		}
	      
	            
	      //if(!isfinite(interp_I[q]))
	      //printf("%d %d > %d %d > %e %e\n",p,q,ix,iy,I_Data[intersect_i][intersect_j][intersect_k][q],intersectGridWeights[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][p]);

	    }
	  interp_S += S[intersect_i][intersect_j][intersect_k]*intersectGridWeights[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng][p];

	}
      double dtau = (source_Data[1][1][1][2]+source_Data[1][1][1][3]) * intersectDistances[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng];

	  double rotAng[3];
	  double rotPhi = intersectGridPhi[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng];
	  //printf("rotPhi %d : %f\n",probeAng,rotPhi);

	  rotAng[0] = cos(rotPhi)*angGridCoords[probeAng][0] - sin(rotPhi)*angGridCoords[probeAng][1];
	  rotAng[1] = sin(rotPhi)*angGridCoords[probeAng][0] + cos(rotPhi)*angGridCoords[probeAng][1];
	  rotAng[2] = angGridCoords[probeAng][2];

	  int bestAngIndex = get_angIndex(rotAng, angGridCoords);
	  

      I_ray[probeAng] = I_Solve(S[1][1][1], interp_S, interp_I[bestAngIndex], dtau);   //SOLVED BY SHORT CHARACTERISTICS!

      I_time[probeAng] = I_Data[1][1][1][probeAng] + (I_ray[probeAng] - I_Data[1][1][1][probeAng])/intersectDistances[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng] * LIGHT_C * delta_t; //apply time step

      
      if(delta_t > intersectDistances[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng])
	printf("dt larger than intersectDistances: %e %e %d %d. increase phi range?\n",intersectDistances[ix+NGCX][iy+NGCY][iz+NGCZ][probeAng],delta_t,ix,iy);
      

      I_return[probeAng] = I_time[probeAng];
      
      //		printf("%d - %e\n", probeAng, I_time[probeAng]);


    }
}

void
ZERO_calcVET(int ix, int iy, int iz,double I_time[NUMANGLES], double eddingtonFactor[3][3], double angGridCoords[NUMANGLES][3])
{

  //Calculate radiative moments using our RT solution to intensity field

  double dOmega = 4.0*PI/NUMANGLES;
  double targetDirection1[3], targetDirection2[3];
  double cos1, cos2;
  double P[3][3], E = 0;
  int probeAng;

  //intialize Pressure tensor P_ij
  int q,r;
  for (q=0; q < 3; q++)
    {
      for (r=0; r < 3; r++)
	{
	  P[q][r] = 0.;
	}
    }

  //Calculate P_ij by summing up contributions over all angles
  for (q=0; q < 3; q++)
    {
      targetDirection1[0] = carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][q][0];
      targetDirection1[1] = carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][q][1];
      targetDirection1[2] = carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][q][2];

      for (r=0; r < 3; r++)
	{
	  targetDirection2[0] = carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][r][0];
	  targetDirection2[1] = carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][r][1];
	  targetDirection2[2] = carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][r][2];
 
	  int probeAng;
	  for (probeAng = 0; probeAng < NUMANGLES; probeAng++)
	    {

	      cos1 = angGridCoords[probeAng][0]*targetDirection1[0] + angGridCoords[probeAng][1]*targetDirection1[1] + angGridCoords[probeAng][2]*targetDirection1[2];
	      cos2 = angGridCoords[probeAng][0]*targetDirection2[0] + angGridCoords[probeAng][1]*targetDirection2[1] + angGridCoords[probeAng][2]*targetDirection2[2];


	      P[q][r] += I_time[probeAng]*cos1*cos2;
	    }
	}
    }
	

  //Calculate E
  for (probeAng = 0; probeAng < NUMANGLES; probeAng++)
    {
      E += I_time[probeAng];
    }


  //	printf("E = %e  |  Pxx = %e\n", E, P[0][0]);


  //Set eddington tensor as P_ij/E
  for (q=0; q < 3; q++)
    {
      for (r=0; r < 3; r++)
	{
	  if (E > 0)
	    {
	      eddingtonFactor[q][r] = P[q][r]/E;
	    }
	  else
	    {
	      eddingtonFactor[q][r] = 0.;
	    }
	}
    }
return;
}

int zero_readangles()
{
  int angGridIndexSort[NUMANGLES][3];
  int angDualGridIndexSort[NUMDUALANGLES][3];


  int readStatus = readAngleFiles(angGridCoords, angDualGridCoords, dualAdjacency);

  if(readStatus==-1) 
    {
      my_err("zero_readangles() failed. missing files best-xyz.dat, best-dualtri-xyz.dat, or best-dualtri-adjIndex.dat.\n");
      exit(-1);
    }

  initAngIndex(angGridCoords, angDualGridCoords, angGridIndexSort, angDualGridIndexSort);
  //Calculate decision trees for BSP angle lookup
  splitAngGrid(NUMANGLES, angGridIndexSort, 0, angGridCoords, &angGridRoot);
  splitDualAngGrid(NUMDUALANGLES, angDualGridIndexSort, 0, angDualGridCoords, &angDualGridRoot);

  //Calculate interpolation weights

  //making backup acting as the previous time step
  int ii;
#pragma omp parallel for private(ii) schedule (static)
  for(ii=0;ii<Nloop_6;ii++) //domain + 1 layer 
    {
      int ix,iy,iz;
      ix=loop_6[ii][0];
      iy=loop_6[ii][1];
      iz=loop_6[ii][2];

      //tetrad
      if(RADCLOSURECOORDS==MINKCOORDS)
	{
 	  carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][0][0]=1.;
	  carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][0][1]=0.;
	  carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][0][2]=0.;
	 
	  carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][1][0]=0.;
	  carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][1][1]=1.;
	  carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][1][2]=0.;

	  carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][2][0]=0.;
	  carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][2][1]=0.;
	  carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][2][2]=1.;
	}

      if(RADCLOSURECOORDS==SPHCOORDS)
	{
	  double xxsph[4],th,ph;
	  get_xx_arb(ix,iy,iz,xxsph,SPHCOORDS);
	  th=xxsph[2];
	  ph=xxsph[3];

	  carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][0][0]=sin(th)*cos(ph);
	  carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][0][1]=sin(th)*sin(ph);
	  carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][0][2]=cos(th);
	 
	  carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][1][0]=cos(th)*cos(ph);
	  carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][1][1]=cos(th)*sin(ph);
	  carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][1][2]=-sin(th);

	  carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][2][0]=-sin(ph);
	  carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][2][1]=cos(ph);
	  carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][2][2]=0.;
	}

      //interpolation weights
      if(MYCOORDS==MINKCOORDS)
	{
	  if(TNZ==1)
	    setupInterpWeights_cart2D(ix,iy,iz,angGridCoords, intersectGridIndices, intersectGridWeights, intersectDistances);
	  else
	    setupInterpWeights_cart3D(ix,iy,iz,angGridCoords, intersectGridIndices, intersectGridWeights, intersectDistances);

	 
	}
      else if(MYCOORDS==SPHCOORDS)
	{
	   if(TNZ==1)
	     setupInterpWeights_sph2D(ix,iy,iz,angGridCoords, intersectGridIndices, intersectGridWeights, intersectDistances, intersectGridPhi);
	   else
	     setupInterpWeights_sph3D(ix,iy,iz,angGridCoords, intersectGridIndices, intersectGridWeights, intersectDistances);
	
	  
	}
      else
	{
	  my_err("Coordinate system not supported for VET\n"); exit(-1);
	}    
    }  


  return readStatus;
}




int ZEROtest_oldmain()

{



	double gridLoc[3][3][3][3];			//xyz locations of 3x3x3 cube
	double M1_Data[3][3][3][5];			//E,v0,v1,v2 output from M1 scheme for cube
	double source_Data[3][3][3][4];				//radiative source info for cube
	double I_Data[3][3][3][NUMANGLES];

	double eddingtonFactor[3][3];
	double I_start[NUMANGLES], I_return[NUMANGLES];

	





	clock_t begin, end;

	int i,j,k,l,p;






//------------------------------------------------------------------------------------------------

	for (p=0; p < NUMANGLES; p++)
	{
	  if (1 || p<25)
		{
			I_start[p]=1.;
		}
		else
		{
			I_start[p]=0.0;
		}

	}

	double M1data[5];

	M1data[0]=80.;
	M1data[1]=0.;
	M1data[2]=70.;
	M1data[3]=0.;

	//rotating, adjusting fluxes
	double fmag = sqrt(M1data[1]*M1data[1] + 
			   M1data[2]*M1data[2] + 
			   M1data[3]*M1data[3]);
	double ff = fmag / M1data[0];
	double beta;
	if(ff<1.e-2) 
	  beta=3.*ff/4.;
	else
	  beta=(4.-sqrt(16.-12.*ff*ff))/2./ff;	    
	
	double gamma=1./sqrt(1.-beta*beta);

	//Erad (Elab, beta)
	M1data[4]=M1data[0]/(4./3.*gamma*gamma-1./3.);
	

	printf("aim F/E: %e\n",ff);

	transformI(0,0,0,I_start, M1data);

	double Fnew[3] = {0.,0.,0.}, Efinal=0.;

	for (p=0; p < NUMANGLES; p++)
	  {
	    Efinal += I_start[p];

	    for (l=0; l < 3; l++)
	      {
		Fnew[l] += I_start[p]*angGridCoords[p][l];
	      }		

	    //printf("%e ",I_start[p]);
	  }
	printf("\n");

	double Fnorm = sqrt(Fnew[0]*Fnew[0] + Fnew[1]*Fnew[1] + Fnew[2]*Fnew[2]);

	printf("Efinal = %e, Ffinal =  %e - %e %e %e | F/E = %e\n", Efinal, Fnorm, Fnew[0], Fnew[1], Fnew[2], Fnorm/Efinal);

	FILE* fout1=fopen("beam.dat","w");
	FILE* fout2=fopen("beam.gp","w");
	fprintf(fout2,"set ylabel \"y\"\n",Fnew[0],Fnew[1],Fnew[2]);
	fprintf(fout2,"set xlabel \"x\"\n",Fnew[0],Fnew[1],Fnew[2]);
	fprintf(fout2,"set zlabel \"z\"\n",Fnew[0],Fnew[1],Fnew[2]);
	fprintf(fout2,"set view equal xyz\n",Fnew[0],Fnew[1],Fnew[2]);
	fprintf(fout2,"splot \"beam.dat\" u 1:2:3 w p pt 7 ps .1\n",Fnew[0],Fnew[1],Fnew[2]);
	int ifzero=0;
	double maxI=-1.;    
	for (p=0; p < NUMANGLES; p++)
	  if(maxI<I_start[p]) maxI=I_start[p];

	for (p=0; p < NUMANGLES; p++)
	  {
	    ifzero=0;
	    
	    if(I_start[p]<3.e-5*maxI)
	      {
		ifzero=1;
		I_start[p]=1.e-1*maxI;
	      }
	    

	    for (l=0; l < 3; l++)
	      {
		Fnew[l] = I_start[p]*angGridCoords[p][l];
	      }		

	    fprintf(fout1,"%f %f %f\n",Fnew[0],Fnew[1],Fnew[2]);

	    if(!ifzero)
	      fprintf(fout2,"set arrow from 0,0,0 to %f,%f,%f front lw 2 lc 2\n",Fnew[0],Fnew[1],Fnew[2]);
	    else	      
	      fprintf(fout2,"set arrow from 0,0,0 to %f,%f,%f nohead front lw 1 lc 3\n",Fnew[0],Fnew[1],Fnew[2]);

	    //fprintf(fout2,"set term pngcairo enhanced\n");
	    //fprintf(fout2,"set output \"beam.png\"\n");
	    //fprintf(fout2,"replot\n");
	
	  }
	fclose(fout1);
	fclose(fout2);
	


	exit(-1);

//------------------------------------------------------------------------------------------------

}



// Use method of lagrange multipliers to determine correct tranformation

int transformI_Lagrange(int ix,int iy,int iz,double I_return[NUMANGLES], double M1_input[5])
{
  double F_final[3],Efinal,Efinalrad;
 
  F_final[0]=M1_input[1]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][0][0] 
    +M1_input[2]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][1][0]
    +M1_input[3]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][2][0];

  F_final[1]=M1_input[1]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][0][1] 
    +M1_input[2]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][1][1]
    +M1_input[3]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][2][1];

  F_final[2]=M1_input[1]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][0][2] 
    +M1_input[2]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][1][2]
    +M1_input[3]*carttetrad[ix+NGCX][iy+NGCY][iz+NGCZ][2][2];
 
  Efinal=M1_input[0];
  Efinalrad=M1_input[4];

  int i,j,p,l;
  double I_start[NUMANGLES];

  double F_start[3], Fmag_start, Fmag_final;
  double Estart;

  double target_vec[4], Lagrange_Multiplier[4];
  double Lagrange_Matrix[4][4], inv_Lagrange_Matrix[4][4];



  for (i=0; i < NUMANGLES; i++)
    {
      I_start[i] = I_return[i];
      I_return[i] = 0.;
    }


  Estart=0.;
  //Calculate net flux direction
  for (l=0; l < 3; l++)
    {
      F_start[l] = 0.;
    }

  for (p=0; p < NUMANGLES; p++)
    {
      //if(p<10) printf("%e\n", I_start[p]);
      for (l=0; l < 3; l++)
	{
	  F_start[l] += I_start[p]*angGridCoords[p][l];
	}
      Estart += I_start[p];
    }




  Fmag_start = sqrt(F_start[0]*F_start[0] + F_start[1]*F_start[1] + F_start[2]*F_start[2]);
  Fmag_final = sqrt(F_final[0]*F_final[0] + F_final[1]*F_final[1] + F_final[2]*F_final[2]);


  target_vec[0]=Efinal-Estart;
  target_vec[1]=F_final[0]-F_start[0];
  target_vec[2]=F_final[1]-F_start[1];
  target_vec[3]=F_final[2]-F_start[2];


  //initialize various quantities to zero
  for (i=0; i < 4; i++)
    {
      Lagrange_Multiplier[i]=0.;
      for (j=0; j < 4; j++)
	{
	  Lagrange_Matrix[i][j]=0.;
	  inv_Lagrange_Matrix[i][j]=0.;
	}
    }

  //Calculate 4x4 matrix for Lagrange multiplier solution
  for (l=0; l<NUMANGLES; l++)
    {
      Lagrange_Matrix[0][0] += 1.0;
      Lagrange_Matrix[0][1] += angGridCoords[l][0];
      Lagrange_Matrix[0][2] += angGridCoords[l][1];
      Lagrange_Matrix[0][3] += angGridCoords[l][2];

      Lagrange_Matrix[1][0] += angGridCoords[l][0];
      Lagrange_Matrix[1][1] += angGridCoords[l][0]*angGridCoords[l][0];
      Lagrange_Matrix[1][2] += angGridCoords[l][0]*angGridCoords[l][1];
      Lagrange_Matrix[1][3] += angGridCoords[l][0]*angGridCoords[l][2];

      Lagrange_Matrix[2][0] += angGridCoords[l][1];
      Lagrange_Matrix[2][1] += angGridCoords[l][1]*angGridCoords[l][0];
      Lagrange_Matrix[2][2] += angGridCoords[l][1]*angGridCoords[l][1];
      Lagrange_Matrix[2][3] += angGridCoords[l][1]*angGridCoords[l][2];

      Lagrange_Matrix[3][0] += angGridCoords[l][2];
      Lagrange_Matrix[3][1] += angGridCoords[l][2]*angGridCoords[l][0];
      Lagrange_Matrix[3][2] += angGridCoords[l][2]*angGridCoords[l][1];
      Lagrange_Matrix[3][3] += angGridCoords[l][2]*angGridCoords[l][2];
    }


  //matrixInverse(Lagrange_Matrix, inv_Lagrange_Matrix);
  inverse_44matrix(Lagrange_Matrix, inv_Lagrange_Matrix);


  //Apply inverse matrix to get multipliers
  for (i=0; i < 4; i++)
    {
      for (j=0; j < 4; j++)
	{
	  Lagrange_Multiplier[i] += inv_Lagrange_Matrix[i][j]*target_vec[j];
	}
    }

  int isnegative=0;
  for(l=0; l < NUMANGLES; l++)
    {
      I_return[l] = I_start[l] + Lagrange_Multiplier[0] + angGridCoords[l][0]*Lagrange_Multiplier[1] + angGridCoords[l][1]*Lagrange_Multiplier[2] + angGridCoords[l][2]*Lagrange_Multiplier[3];

      if(I_return[l]<0.) 
	isnegative=1;	
    }

  return isnegative;
}

void transformI(int ix, int iy, int iz,double I0[NUMANGLES], double M1_input[5])
{
  transformI_stretch(ix,iy,iz,I0,M1_input);
  return;

  /*  
  double I1[NUMANGLES];
  int i,j;
  int ret,iter;
  
  for(i=0;i<NUMANGLES;i++)
    {
      I1[i]=I0[i];
    }
  
  for(iter=0;iter<1;iter++)
    {
      ret=transformI_Lagrange(ix,iy,iz,I0,M1_input);
  
      if(ret==1) //some intensities negative
	{	  
	  transformI_stretch(ix,iy,iz,I1,M1_input); //these are guaranteed to be positive
	  double minf=1.e100,f; //I = f I0 + (1-f) I1 
	  for(i=0;i<NUMANGLES;i++)
	    {
	      if(I0[i]<0.)
		{
		  f=I1[i]/(I1[i]-I0[i]);
		  if(f<minf) minf=f;
		}
	    }

	  if(minf<0.) minf=0.;
	  if(minf>1.) minf=1.;

	  // printf("combining with %f at %d %d\n",minf,ix,iy); 
	  for(i=0;i<NUMANGLES;i++)
	    {
	      I0[i]=minf*I0[i] + (1.-minf)*I1[i];
	      if(I0[i]<0.) I0[i]=0.;
	      I1[i]=I0[i];
	    }
	}
      else
	break;
    
    }
  

  return;
  */
}

