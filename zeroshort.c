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




void setupInterpWeights(double angGridCoords[NUMANGLES][3], int intersectGridIndices[NUMANGLES][3][4], double intersectGridWeights[NUMANGLES][4], double intersectDistances[NUMANGLES])
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
	      intersectGridIndices[probeAng][0][p] = realXIndex;
	    }
	  intersectGridIndices[probeAng][1][0] = lowerYIndex;
	  intersectGridIndices[probeAng][1][1] = lowerYIndex;
	  intersectGridIndices[probeAng][1][2] = upperYIndex;
	  intersectGridIndices[probeAng][1][3] = upperYIndex;

	  intersectGridIndices[probeAng][2][0] = lowerZIndex;
	  intersectGridIndices[probeAng][2][1] = upperZIndex;
	  intersectGridIndices[probeAng][2][2] = lowerZIndex;
	  intersectGridIndices[probeAng][2][3] = upperZIndex;

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
	      intersectGridIndices[probeAng][1][p] = realYIndex;
	    }
	  intersectGridIndices[probeAng][0][0] = lowerXIndex;
	  intersectGridIndices[probeAng][0][1] = lowerXIndex;
	  intersectGridIndices[probeAng][0][2] = upperXIndex;
	  intersectGridIndices[probeAng][0][3] = upperXIndex;

	  intersectGridIndices[probeAng][2][0] = lowerZIndex;
	  intersectGridIndices[probeAng][2][1] = upperZIndex;
	  intersectGridIndices[probeAng][2][2] = lowerZIndex;
	  intersectGridIndices[probeAng][2][3] = upperZIndex;

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
	      intersectGridIndices[probeAng][2][p] = realZIndex;
	    }
	  intersectGridIndices[probeAng][0][0] = lowerXIndex;
	  intersectGridIndices[probeAng][0][1] = lowerXIndex;
	  intersectGridIndices[probeAng][0][2] = upperXIndex;
	  intersectGridIndices[probeAng][0][3] = upperXIndex;

	  intersectGridIndices[probeAng][1][0] = lowerYIndex;
	  intersectGridIndices[probeAng][1][1] = upperYIndex;
	  intersectGridIndices[probeAng][1][2] = lowerYIndex;
	  intersectGridIndices[probeAng][1][3] = upperYIndex;

	  w0_high = posX - floor(posX);
	  w0_low = 1.0 - w0_high;
	  w1_high = posY - floor(posY);
	  w1_low = 1.0 - w1_high;
		
	}


      intersectGridWeights[probeAng][0] = w0_low*w1_low;
      intersectGridWeights[probeAng][1] = w0_low*w1_high;
      intersectGridWeights[probeAng][2] = w0_high*w1_low;
      intersectGridWeights[probeAng][3] = w0_high*w1_high;

      intersectDistances[probeAng] = maxL*GRID_SPACING;



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


void transformI(double I_return[NUMANGLES], double F_final[3], double fFinal, struct bsptree *angDualGridRoot, double angGridCoords[NUMANGLES][3], double angDualGridCoords[NUMDUALANGLES][3], int dualAdjacency[NUMDUALANGLES][3])
{
	int i,j,p,l;
	double I_start[NUMANGLES];
	double res[3];

	double F_start[3], F_start_norm[3], F_final_norm[3], Fmag_start, Fmag_final;
	double fStart, Estart;
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
	for (l=0; l < 3; l++)
	{
		F_start_norm[l] = F_start[l]/Fmag_start;
		F_final_norm[l] = F_final[l]/Fmag_final;
	}

	fStart = Fmag_start/Estart;

	stretchFactor = sqrt(fabs((fFinal*fFinal - fFinal*fFinal*fStart*fStart)/(fStart*fStart - fFinal*fFinal*fStart*fStart)));


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


			n_parallel[l] = n_parallel[l]*stretchFactor;
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
	//	bestIndex = get_angDualIndex(n_final, angDualGridCoords);

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



	double rescaleFactor = Fmag_final/Fmag_start/stretchFactor;
	//Renormalize to get correct fluxes
	for (p=0; p < NUMANGLES; p++)
	{
		I_return[p] = I_return[p]*rescaleFactor;
	}
}







void ZERO_shortChar(double delta_t, double M1_Data[3][3][3][5], double source_Data[3][3][3][4], double angGridCoords[NUMANGLES][3], int intersectGridIndices[NUMANGLES][3][4], double intersectGridWeights[NUMANGLES][4], double intersectDistances[NUMANGLES], double eddingtonFactor[3][3], double I_return[NUMANGLES], double F_return[3],int verbose)
{
  //Note:  M1 data has format:   E, v0, v1, v2

  double S[3][3][3];  //radiative source function
  double I_Data[3][3][3][NUMANGLES]; //Radiative intensity at all boundary points
  double f_norm[3][3][3][3]; //Normalized velocity 

  double I_ray[NUMANGLES];  //intensity field evaluated at the center, for time independent problem
  double I_time[NUMANGLES];  //"" "", for time-dependent problem


  int q,r;
  
  //Fill out source function and intensity field for cube

  double beta[3][3][3];
  double ffzero;
  int i,j,k;
  for (i=0; i<3; i++)
    {
      for (j=0; j<3; j++)
	{
	  for (k=0; k<3; k++)
	    {

	      // SOURCE FUNCTION
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

	      double temp=source_Data[i][j][k][0];
	      S[i][j][k]=eps*temp*temp*temp*temp*STEFAN_BOLTZMANN/PI + (1.0-eps)*source_Data[i][j][k][1]*LIGHT_C/4.0/PI;
		

	      double fmag = sqrt(M1_Data[i][j][k][1]*M1_Data[i][j][k][1] + 
				 M1_Data[i][j][k][2]*M1_Data[i][j][k][2] + 
				 M1_Data[i][j][k][3]*M1_Data[i][j][k][3]);

	      double ff = fmag / M1_Data[i][j][k][0]; //F/E using input argument
	      if(i==1 && j==1 && k==1) ffzero=ff;
	      if(ff<1.e-2) 
		beta[i][j][k]=3.*ff/4.;
	      else
		beta[i][j][k]=(4.-sqrt(16.-12.*ff*ff))/2./ff;	    

	      if(fmag<SMALL)
		fmag=1.;

	      f_norm[i][j][k][0] = M1_Data[i][j][k][1]/fmag;
	      f_norm[i][j][k][1] = M1_Data[i][j][k][2]/fmag;
	      f_norm[i][j][k][2] = M1_Data[i][j][k][3]/fmag;
	    }
	}
    }

  double Estart = 0., Fstart[3] = {0., 0., 0.};

  int probeAng;
 
  double gamma2 = 1.0/(1.0-beta[1][1][1]*beta[1][1][1]);
 


  for (probeAng=0; probeAng < NUMANGLES; probeAng++)
    {
      double mu= f_norm[1][1][1][0]*angGridCoords[probeAng][0] + f_norm[1][1][1][1]*angGridCoords[probeAng][1] + f_norm[1][1][1][2]*angGridCoords[probeAng][2];
		      
      if(beta[1][1][1]<SMALL)
	mu=1.;


      double mu_max = mu + 0.1;
      double mu_min = mu - 0.1;

      if (mu_max > 1.0)
	{
	  mu_max = 1.0;
	}
      if (mu_min < -1.0)
	{
	  mu_min = -1.0;
	}


      double bm=1-beta[1][1][1]*mu;
      double bmax=1-beta[1][1][1]*mu_max;
      double bmin=1-beta[1][1][1]*mu_min;

      //Factor 2 is to convert E_iso to I_iso, since I_iso = E_iso/4pi, where we also absorb *2pi factor for phi integral
      I_Data[1][1][1][probeAng] = M1_Data[1][1][1][4]/NUMANGLES/bm/bm/bm/bm/gamma2/gamma2;
      Estart += I_Data[1][1][1][probeAng];

      //if(probeAng<10) printf("%e\n", I_start[probeAng]);
      for (j=0; j<3; j++)
	{

	  Fstart[j]+=I_Data[1][1][1][probeAng]*angGridCoords[probeAng][j];
	}

      //I_Data[1][1][1][probeAng] = M1_Data[1][1][1][0]/gamma2/gamma2/beta/3.0*(1.0/bmax/bmax/bmax - 1.0/bmin/bmin/bmin)/(mu_max - mu_min);

    }

  double Ffinal[3];
  
  /*
  printf("INPT: Elab = (%e)  Erad = (%e)\n", M1_Data[1][1][1][0],M1_Data[1][1][1][4]);

  printf("INPT: F/E = (%e, %e, %e) beta1 = (%e)\n", M1_Data[1][1][1][1]/M1_Data[1][1][1][0], 
	 M1_Data[1][1][1][2]/M1_Data[1][1][1][0],
	 M1_Data[1][1][1][3]/M1_Data[1][1][1][0],beta[1][1][1]);
  printf("ZERO: F/E = (%e, %e, %e)\n", Fstart[0]/Estart, Fstart[1]/Estart, Fstart[2]/Estart);
  */

  transformI(&I_Data[1][1][1][0], &M1_Data[1][1][1][1], ffzero, angDualGridRoot, angGridCoords, angDualGridCoords, dualAdjacency);

  /*
  Estart=0.;Fstart[0]=Fstart[1]=Fstart[2]=0.;
  for (probeAng=0; probeAng < NUMANGLES; probeAng++)
    {
      Estart += I_Data[1][1][1][probeAng];
      for (j=0; j<3; j++)
	{
	  Fstart[j]+=I_Data[1][1][1][probeAng]*angGridCoords[probeAng][j];
	}
    }
  */
  //printf("ZERO: F/E = (%e, %e, %e)\n", Fstart[0]/Estart, Fstart[1]/Estart, Fstart[2]/Estart);
  //getch();
  

  




  //Calculate solution to RT, given M1 radiation shapes
  for (probeAng = 0; probeAng < NUMANGLES; probeAng++)
    {

      //	  printf("%d\n", probeAng);


      //First, get interpolated M1 quantities along ray intersection boundary
      double interp_I = 0.;
      double interp_S = 0.;

      int p;
      for (p=0; p < 4; p++)
	{
	  int intersect_i = intersectGridIndices[probeAng][0][p];
	  int intersect_j = intersectGridIndices[probeAng][1][p];
	  int intersect_k = intersectGridIndices[probeAng][2][p];




	  //INTENSITY FIELD2
	  //Calculate angle factor between ray angle and M1 velocity at boundary (beta mu)

	  double mu= f_norm[intersect_i][intersect_j][intersect_k][0]*angGridCoords[probeAng][0] + f_norm[intersect_i][intersect_j][intersect_k][1]*angGridCoords[probeAng][1] + f_norm[intersect_i][intersect_j][intersect_k][2]*angGridCoords[probeAng][2];

			
	  double betahere = beta[intersect_i][intersect_j][intersect_k];
	  double gamma2 = 1.0/(1.0-betahere*betahere);

	  if(betahere<SMALL)
	    mu=1.;

	  //Get radiative quantities

	 

	  double mu_max = mu + 0.1;
	  double mu_min = mu - 0.1;

	  if (mu_max > 1.0)
	    {
	      mu_max = 1.0;
	    }
	  if (mu_min < -1.0)
	    {
	      mu_min = -1.0;
	    }


	  double bm=1.-betahere*mu;
	  double bmax=1.-betahere*mu_max;
	  double bmin=1.-betahere*mu_min;


	  I_Data[intersect_i][intersect_j][intersect_k][probeAng] = M1_Data[intersect_i][intersect_j][intersect_k][4]/2.0/bm/bm/bm/bm/gamma2/gamma2;
	  //I_Data[intersect_i][intersect_j][intersect_k][probeAng] = M1_Data[intersect_i][intersect_j][intersect_k][0]/gamma2/gamma2/beta/3.0*(1.0/bmax/bmax/bmax - 1.0/bmin/bmin/bmin)/(mu_max - mu_min);
	   
	  interp_I += I_Data[intersect_i][intersect_j][intersect_k][probeAng]*intersectGridWeights[probeAng][p];
	  interp_S += S[intersect_i][intersect_j][intersect_k]*intersectGridWeights[probeAng][p];


	}






      

      double dtau = (source_Data[1][1][1][2]+source_Data[1][1][1][3]) * intersectDistances[probeAng];


      double iray=I_Solve(S[1][1][1], interp_S, interp_I, dtau);   //SOLVED BY SHORT CHARACTERISTICS!;
      I_ray[probeAng] = iray;
      

      double I_old = I_Data[1][1][1][probeAng];
      I_time[probeAng] = I_old + (I_ray[probeAng] - I_old)/intersectDistances[probeAng] * LIGHT_C * delta_t; //apply time step
      
      I_return[probeAng] = I_time[probeAng];
      
      /*
      if (verbose)
	{
	  printf("%d - %e %e\n", probeAng, interp_I, I_ray[probeAng]);
	}
      */
       
      
    }


  



  //Calculate radiative moments using our RT solution to intensity field

  double dOmega = 4.0*PI/NUMANGLES;
  double targetDirection1[3], targetDirection2[3];
  double cos1, cos2;
  double P[3][3], E = 0;

  //intialize Pressure tensor P_ij
  for (q=0; q < 3; q++)
    {
      F_return[q]=0;

      for (r=0; r < 3; r++)
	{
	  P[q][r] = 0.;
	}
    }

  //Calculate P_ij by summing up contributions over all angles
  for (q=0; q < 3; q++)
    {
      if (q == 0)
	{
	  targetDirection1[0] = 1.0;
	  targetDirection1[1] = 0.0;
	  targetDirection1[2] = 0.0;
	}
      if (q == 1)
	{
	  targetDirection1[0] = 0.0;
	  targetDirection1[1] = 1.0;
	  targetDirection1[2] = 0.0;
	}
      if (q == 2)
	{
	  targetDirection1[0] = 0.0;
	  targetDirection1[1] = 0.0;
	  targetDirection1[2] = 1.0;
	}


      for (r=0; r < 3; r++)
	{
	  if (r == 0)
	    {
	      targetDirection2[0] = 1.0;
	      targetDirection2[1] = 0.0;
	      targetDirection2[2] = 0.0;
	    }
	  if (r == 1)
	    {
	      targetDirection2[0] = 0.0;
	      targetDirection2[1] = 1.0;
	      targetDirection2[2] = 0.0;
	    }
	  if (r == 2)
	    {
	      targetDirection2[0] = 0.0;
	      targetDirection2[1] = 0.0;
	      targetDirection2[2] = 1.0;
	    }


	  int probeAng;
	  for (probeAng = 0; probeAng < NUMANGLES; probeAng++)
	    {

	      cos1 = angGridCoords[probeAng][0]*targetDirection1[0] + angGridCoords[probeAng][1]*targetDirection1[1] + angGridCoords[probeAng][2]*targetDirection1[2];
	      cos2 = angGridCoords[probeAng][0]*targetDirection2[0] + angGridCoords[probeAng][1]*targetDirection2[1] + angGridCoords[probeAng][2]*targetDirection2[2];


	      P[q][r] += I_time[probeAng]*cos1*cos2*dOmega;
	      if (r==0)
		{
		  F_return[q] += I_time[probeAng]*cos1*dOmega;
		}
	    }
	}
    }
	

  //Calculate E
  for (probeAng = 0; probeAng < NUMANGLES; probeAng++)
    {
      E += I_time[probeAng]*dOmega;
    }


  //	printf("E = %e  |  Pxx = %e\n", E, P[0][0]);


  //Set eddington tensor as P_ij/E
  // Normalize flux by E
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
	    

	  

	    /*
	  //test
	  if(q==r)
	    eddingtonFactor[q][r] = 1./3.;
	  else
	    eddingtonFactor[q][r] = 0.;
	  */
	}

      F_return[q] = F_return[q]/E;
    }

}


void ZERO_shortCharI(double delta_t, double I_Data[3][3][3][NUMANGLES], double source_Data[3][3][3][4], double angGridCoords[NUMANGLES][3], int intersectGridIndices[NUMANGLES][3][4], double intersectGridWeights[NUMANGLES][4], double intersectDistances[NUMANGLES], double eddingtonFactor[3][3], double I_return[NUMANGLES],int verbose)
{
  //Note:  M1 data has format:   E, v0, v1, v2

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
	  int intersect_i = intersectGridIndices[probeAng][0][p];
	  int intersect_j = intersectGridIndices[probeAng][1][p];
	  int intersect_k = intersectGridIndices[probeAng][2][p];

	  int q;
	  for (q=0 ; q<NUMANGLES; q++)
	    {
	      interp_I[q] += I_Data[intersect_i][intersect_j][intersect_k][q]*intersectGridWeights[probeAng][p];
	    }
	  interp_S += S[intersect_i][intersect_j][intersect_k]*intersectGridWeights[probeAng][p];

	}
      double dtau = (source_Data[1][1][1][2]+source_Data[1][1][1][3]) * intersectDistances[probeAng];

      I_ray[probeAng] = I_Solve(S[1][1][1], interp_S, interp_I[probeAng], dtau);   //SOLVED BY SHORT CHARACTERISTICS!
      I_time[probeAng] = I_Data[1][1][1][probeAng] + (I_ray[probeAng] - I_Data[1][1][1][probeAng])/intersectDistances[probeAng] * LIGHT_C * delta_t; //apply time step
      I_return[probeAng] = I_time[probeAng];


      //		printf("%d - %e\n", probeAng, I_time[probeAng]);


    }






  //Calculate radiative moments using our RT solution to intensity field

  double dOmega = 4.0*PI/NUMANGLES;
  double targetDirection1[3], targetDirection2[3];
  double cos1, cos2;
  double P[3][3], E = 0;

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
      if (q == 0)
	{
	  targetDirection1[0] = 1.0;
	  targetDirection1[1] = 0.0;
	  targetDirection1[2] = 0.0;
	}
      if (q == 1)
	{
	  targetDirection1[0] = 0.0;
	  targetDirection1[1] = 1.0;
	  targetDirection1[2] = 0.0;
	}
      if (q == 2)
	{
	  targetDirection1[0] = 0.0;
	  targetDirection1[1] = 0.0;
	  targetDirection1[2] = 1.0;
	}


      for (r=0; r < 3; r++)
	{
	  if (r == 0)
	    {
	      targetDirection2[0] = 1.0;
	      targetDirection2[1] = 0.0;
	      targetDirection2[2] = 0.0;
	    }
	  if (r == 1)
	    {
	      targetDirection2[0] = 0.0;
	      targetDirection2[1] = 1.0;
	      targetDirection2[2] = 0.0;
	    }
	  if (r == 2)
	    {
	      targetDirection2[0] = 0.0;
	      targetDirection2[1] = 0.0;
	      targetDirection2[2] = 1.0;
	    }


	  int probeAng;
	  for (probeAng = 0; probeAng < NUMANGLES; probeAng++)
	    {

	      cos1 = angGridCoords[probeAng][0]*targetDirection1[0] + angGridCoords[probeAng][1]*targetDirection1[1] + angGridCoords[probeAng][2]*targetDirection1[2];
	      cos2 = angGridCoords[probeAng][0]*targetDirection2[0] + angGridCoords[probeAng][1]*targetDirection2[1] + angGridCoords[probeAng][2]*targetDirection2[2];


	      P[q][r] += I_time[probeAng]*cos1*cos2*dOmega;
	    }
	}
    }
	

  //Calculate E
  for (probeAng = 0; probeAng < NUMANGLES; probeAng++)
    {
      E += I_time[probeAng]*dOmega;
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


}



int zero_readangles()
{
  int angGridIndexSort[NUMANGLES][3];
  int angDualGridIndexSort[NUMDUALANGLES][3];


  int readStatus = readAngleFiles(angGridCoords, angDualGridCoords, dualAdjacency);
  setupInterpWeights(angGridCoords, intersectGridIndices, intersectGridWeights, intersectDistances);
  initAngIndex(angGridCoords, angDualGridCoords, angGridIndexSort, angDualGridIndexSort);
  //Calculate decision trees for BSP angle lookup
  splitAngGrid(NUMANGLES, angGridIndexSort, 0, angGridCoords, &angGridRoot);
  splitDualAngGrid(NUMDUALANGLES, angDualGridIndexSort, 0, angDualGridCoords, &angDualGridRoot);



  if(readStatus==-1) 
    {
      my_err("zero_readangles() failed. exiting.\n");
      exit(-1);
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

	double Ffinal[3];






	clock_t begin, end;

	int i,j,k,l,p;






//------------------------------------------------------------------------------------------------

	for (p=0; p < NUMANGLES; p++)
	{
		if (p<20)
		{
			I_start[p]=1.0;
		}
		else
		{
			I_start[p]=0.0;
		}

	}

	Ffinal[0]=0.0;
	Ffinal[1]=4.0;
	Ffinal[2]=0.0;


	transformI(I_start, Ffinal, 0.7, angDualGridRoot, angGridCoords, angDualGridCoords, dualAdjacency);



	double Fnew[3] = {0.,0.,0.}, Efinal=0.;

	for (p=0; p < NUMANGLES; p++)
	  {
	    //Efinal += I_return[p];
	    Efinal += I_start[p];

	    for (l=0; l < 3; l++)
	      {
		//Fnew[l] += I_return[p]*angGridCoords[p][l];
		Fnew[l] += I_start[p]*angGridCoords[p][l];
	      }		
	  }

	double Fnorm = sqrt(Fnew[0]*Fnew[0] + Fnew[1]*Fnew[1] + Fnew[2]*Fnew[2]);

	printf("Efinal = %e, Ffinal =  %e - %e %e %e | F/E = %e\n", Efinal, Fnorm, Fnew[0], Fnew[1], Fnew[2], Fnorm/Efinal);


/*
	for (l=0; l < 3; l++)
	{
		Fnew[l] += Fnew[l]/Fnorm;
	}

	printf("%e\n", sqrt(F_end[0]*F_end[0] + F_end[1]*F_end[1] + F_end[2]*F_end[2])/Efinal);
*/




//	printf("%e %e\n", stretchFactor, Efinal);
//	printf("cos2avg = %e\n", cos2Avg/Itot);
	

//	printf("E static = %e, E stretch = %e\n", Estatic, Estretch);
//	printf("final: %e %e %e || %e \n", Fnew[0], Fnew[1], Fnew[2], Efinal);
//	printf("F end: %e %e %e || %e \n", F_end[0], F_end[1], F_end[2], sqrt(F_end[0]*F_end[0] + F_end[1]*F_end[1] + F_end[2]*F_end[2])/Efinal);






/*

	double interpAng[3] = {0., 0., 0.};


	for (l=0; l < 3; l++)
	{
		for (p=0; p < 3; p++)
		{
			interpAng[l] += angGridCoords[angNeighborIndex[p]][l]*interpCoeffs[p];
		}
	}

	printf("Interp Angle: %e %e %e\n", interpAng[0], interpAng[1], interpAng[2]);
*/


	exit(-1);

//------------------------------------------------------------------------------------------------


}
