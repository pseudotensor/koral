#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "ko.h"

#define GRID_SPACING 0.0588        //distance between adjacent grid cells (grid is assumed cubic with equal spacing in all directions)
#define LIGHT_C 1.		//Speed of light in code units
#define STEFAN_BOLTZMANN SIGMA_RAD //stefan boltzmann constant in code units
#define PI 3.14159265358979



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
		return ((exp(-dtau)*(S0 - S1*(1.0+dtau)) + S0*(dtau - 1.0) + S1)/dtau + I1*exp(-dtau));
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





void ZERO_shortChar(double delta_t, double M1_Data[3][3][3][4], double source_Data[3][3][3][4], double angGridCoords[NUMANGLES][3], int intersectGridIndices[NUMANGLES][3][4], double intersectGridWeights[NUMANGLES][4], double intersectDistances[NUMANGLES], double eddingtonFactor[3][3])
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




	//Calculate solution to RT, given M1 radiation shapes
	int probeAng;
	for (probeAng = 0; probeAng < NUMANGLES; probeAng++)
	{


		//First, get interpolated M1 quantities along ray intersection boundary
		double interp_M1_Data[4];
		double interp_S = 0.;

		//initialize variables
		int p;
		for (p=0; p < 4; p++)
		{
			interp_M1_Data[p] = 0.;
		}



		for (p=0; p < 4; p++)
		{
			int intersect_i = intersectGridIndices[probeAng][0][p];
			int intersect_j = intersectGridIndices[probeAng][1][p];
			int intersect_k = intersectGridIndices[probeAng][2][p];

			int q;
			for (q=0 ; q<4; q++)
			{
				interp_M1_Data[q] += M1_Data[intersect_i][intersect_j][intersect_k][q]*intersectGridWeights[probeAng][p];
			}
			interp_S += S[intersect_i][intersect_j][intersect_k]*intersectGridWeights[probeAng][p];

		}


		//Calculate angle factor between ray angle and M1 velocity at boundary (beta mu)
		double mu=0.;
		double v_ang[3], vmag = sqrt(interp_M1_Data[1]*interp_M1_Data[1] + interp_M1_Data[2]*interp_M1_Data[2] + interp_M1_Data[3]*interp_M1_Data[3]);

		if (vmag <= 0)
		{
			vmag= 1.0;
			mu = 1.0;
		}
		else
		{
		  int l;
			for (l=1; l < 4; l++)
			{
				v_ang[l-1]=interp_M1_Data[l]/vmag;
			}

			//n dot v/|v|
			mu = v_ang[0]*angGridCoords[probeAng][0] + v_ang[1]*angGridCoords[probeAng][1] + v_ang[2]*angGridCoords[probeAng][2];
		}

		double beta = sqrt(interp_M1_Data[1]*interp_M1_Data[1] + interp_M1_Data[2]*interp_M1_Data[2] + interp_M1_Data[3]*interp_M1_Data[3]);  //This assumes the velocities are already normalized by speed of light!!!
		double gamma = 1.0/sqrt(1.0-beta*beta);





		//Calculate same quantities as above, but for central cell
		double c_mu=0.;
		double c_v_ang[3], c_vmag = sqrt(M1_Data[1][1][1][1]*M1_Data[1][1][1][1] + M1_Data[1][1][1][2]*M1_Data[1][1][1][2] + M1_Data[1][1][1][3]*M1_Data[1][1][1][3]);
		if (c_vmag <= 0)
		{
			c_vmag= 1.0;
			c_mu = 1.0;

		}
		else
		{
		  int l;
			for (l=1; l < 4; l++)
			{
				c_v_ang[l-1]=M1_Data[1][1][1][l]/vmag;
			}



			//n dot v/|v|
			c_mu = c_v_ang[0]*angGridCoords[probeAng][0] + c_v_ang[1]*angGridCoords[probeAng][1] + c_v_ang[2]*angGridCoords[probeAng][2];
		}

		double c_beta = sqrt(M1_Data[1][1][1][1]*M1_Data[1][1][1][1] + M1_Data[1][1][1][2]*M1_Data[1][1][1][2] + M1_Data[1][1][1][3]*M1_Data[1][1][1][3]);  //This assumes the velocities are already normalized by speed of light!!!
		double c_gamma = 1.0/sqrt(1.0-beta*beta);




		//Get radiative quantities
		double I_boundary = interp_M1_Data[0]/pow(gamma*(1-beta*mu),4);
		double I_center = M1_Data[1][1][1][0]/pow(c_gamma*(1-c_beta*c_mu),4);
		double dtau = (source_Data[1][1][1][2]+source_Data[1][1][1][3]) * intersectDistances[probeAng];


		I_ray[probeAng] = I_Solve(S[1][1][1], interp_S, I_boundary, dtau);   //SOLVED BY SHORT CHARACTERISTICS!
		I_time[probeAng] = I_center + (I_ray[probeAng] - I_center)/intersectDistances[probeAng] * LIGHT_C * delta_t; //apply time step


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









int ZEROtest_oldmain()
{

	double angGridCoords[NUMANGLES][3];  		//Store xyz locations of angle grid
	double angDualGridCoords[NUMDUALANGLES][3]; 	//Store xyz locations of dual angle grid
	int dualAdjacency[NUMDUALANGLES][3]; 		//Store index information for adjacent angles
	int readStatus;

	double gridLoc[3][3][3][3];			//xyz locations of 3x3x3 cube
	double M1_Data[3][3][3][4];			//E,v0,v1,v2 output from M1 scheme for cube
	double source_Data[3][3][3][4];				//radiative source info for cube


	int intersectGridIndices[NUMANGLES][3][4];	//indices for gridLoc of 4 points bounding intersection
	double intersectGridWeights[NUMANGLES][4];	//weights corresponding to 4 points bounding intersection
	double intersectDistances[NUMANGLES];		//distance to intersection point from center

	double eddingtonFactor[3][3];



	readStatus = readAngleFiles(angGridCoords, angDualGridCoords, dualAdjacency);
	setupInterpWeights(angGridCoords, intersectGridIndices, intersectGridWeights, intersectDistances);

	int i,j,k;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < 3; k++)
			{
				M1_Data[i][j][k][0] = 1.0; //E
				M1_Data[i][j][k][1] = 0.0; //v0
				M1_Data[i][j][k][2] = 0.0; //v1
				M1_Data[i][j][k][3] = 0.0; //v2

				source_Data[i][j][k][0]=0.; // T_gas
				source_Data[i][j][k][1]=0.; // E_rad
				source_Data[i][j][k][2]=0.; // abs coeff: alpha (1/length)
				source_Data[i][j][k][3]=0.; // scat coeff: sigma "" ""

			}
		}
	}



	ZERO_shortChar(0.5, M1_Data, source_Data, angGridCoords, intersectGridIndices, intersectGridWeights, intersectDistances, eddingtonFactor);


	
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			printf("%e ", eddingtonFactor[i][j]);
		}
		printf("\n");
	}


	return 1;
}

int
zero_readangles()
{
  int readStatus = readAngleFiles(angGridCoords, angDualGridCoords, dualAdjacency);
  setupInterpWeights(angGridCoords, intersectGridIndices, intersectGridWeights, intersectDistances);

  if(readStatus==-1) 
    {
      my_err("zero_readangles() failed. exiting.\n");
      exit(-1);
    }

  return readStatus;
}
