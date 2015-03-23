#include "ko.h"

int
calc_avgs_throughout()
{
  /***************************/
  //scale height at each radius
#ifdef CALCHRONTHEGO
  int ix,iy,iz,gix,giy,giz; struct geometry geom, geomBL;
  ldouble sigma,scaleth,xxBL[4]; 
  for(gix=0;gix<TNX;gix++)
    sigma_otg_temp[gix]=scaleth_otg_temp[gix]=0.;
  for(ix=0;ix<NX;ix++)
    {
      sigma=scaleth=0.;
      for(iy=0;iy<NY;iy++)
	{
	  for(iz=0;iz<NZ;iz++)
	    {
	      fill_geometry(ix,iy,iz,&geom);
	      coco_N(geom.xxvec,xxBL,MYCOORDS,BLCOORDS);
	      sigma+=get_u(p,RHO,ix,iy,iz)*geom.gdet;
	      scaleth+=get_u(p,RHO,ix,iy,iz)*geom.gdet*(M_PI/2. - xxBL[2])*(M_PI/2. - xxBL[2]);
	    }
	}
      gix=ix+TOI;
      
      sigma_otg_temp[gix]=sigma;
      scaleth_otg_temp[gix]=scaleth;
    }

#ifdef MPI
  MPI_Allreduce(sigma_otg_temp, sigma_otg, TNX, MPI_LDOUBLE, MPI_SUM,
		MPI_COMM_WORLD);
  MPI_Allreduce(scaleth_otg_temp, scaleth_otg, TNX, MPI_LDOUBLE, MPI_SUM,
		MPI_COMM_WORLD);
#else
  for(ix=0;ix<NX;ix++)
    {
      gix=ix+TOI;
      sigma_otg[gix]=sigma_otg_temp[gix];
      scaleth_otg[gix]=scaleth_otg_temp[gix];
    }
#endif
  
for(ix=-NGCX;ix<NX+NGCX;ix++)
  {
    gix=ix+TOI;
    if(gix>0 && gix<TNX)
      scaleth_otg[gix]=sqrt(scaleth_otg[gix]/sigma_otg[gix]);
  }
#endif
/***************************/


  /***************************/
  //averaging velocities 
#ifdef CORRECT_POLARAXIS_3D
#ifdef POLARAXISAVGIN3D

 #ifdef PHIWEDGE
 if(PHIWEDGE<0.999*2.*M_PI)
   my_warning("POLARAXISAVGIN3D requires full 2 pi in azimuth!\n");
#endif
 
 int ix,iy,iz,gix,giy,giz,iv; 
 struct geometry geom, geomBL;
 ldouble v[4],ucon[4],r,th,ph;

 for(gix=0;gix<TNX;gix++)
   for(iv=0;iv<NV+2;iv++)
     axis1_primplus_temp[iv][gix]= axis1_primplus[iv][gix]=axis2_primplus_temp[iv][gix]= axis2_primplus[iv][gix]=0.;

 for(ix=0;ix<NX;ix++)
   {
     gix=ix+TOI;
     for(iz=0;iz<NZ;iz++)
       {
#ifdef MPI
	 if(TJ==0) //topmost tile
#endif
	   {
	     fill_geometry(ix,NCCORRECTPOLAR,iz,&geom);
	     fill_geometry_arb(ix,NCCORRECTPOLAR,iz,&geomBL,BLCOORDS);

	     axis1_primplus_temp[RHO][gix]+=get_u(p,RHO,ix,NCCORRECTPOLAR,iz);
	     axis1_primplus_temp[UU][gix]+=get_u(p,UU,ix,NCCORRECTPOLAR,iz);
	     //cartesian velocities in VX..VZ slots
	     decompose_vels(&get_u(p,0,ix,NCCORRECTPOLAR,iz),VX,v,&geom,&geomBL);
	     axis1_primplus_temp[VX][gix]+=v[1];
	     axis1_primplus_temp[VY][gix]+=v[2];
	     axis1_primplus_temp[VZ][gix]+=v[3];
	     //angular velocity in NV slot!!!
	     axis1_primplus_temp[NV][gix]+=get_u(p,VZ,ix,NCCORRECTPOLAR,iz);

#ifdef RADIATION
	     axis1_primplus_temp[EE][gix]+=get_u(p,EE,ix,NCCORRECTPOLAR,iz);
	     #ifdef NCOMPTONIZATION
	     axis1_primplus_temp[NF][gix]+=get_u(p,NF,ix,NCCORRECTPOLAR,iz);
	     #endif
	     //cartesian velocities in FX..FZ slots
	     decompose_vels(&get_u(p,0,ix,NCCORRECTPOLAR,iz),FX,v,&geom,&geomBL);
	     //if(ix==10) printf("1 %d > %e | %e %e %e\n",ix,get_u(p,EE,ix,NCCORRECTPOLAR,iz),v[1],v[2],v[3]);
	     axis1_primplus_temp[FX][gix]+=v[1];
	     axis1_primplus_temp[FY][gix]+=v[2];
	     axis1_primplus_temp[FZ][gix]+=v[3];
	     //rad angular velocity in NV+1 slot!!!
	     axis1_primplus_temp[NV+1][gix]+=get_u(p,FZ,ix,NCCORRECTPOLAR,iz);
#endif
	   }

#ifdef MPI
	 if(TJ==NTY-1) //bottommost tile
#endif
	   {
	     fill_geometry(ix,NY-NCCORRECTPOLAR-1,iz,&geom);
	     fill_geometry_arb(ix,NY-NCCORRECTPOLAR-1,iz,&geomBL,BLCOORDS);

	     axis2_primplus_temp[RHO][gix]+=get_u(p,RHO,ix,NY-NCCORRECTPOLAR-1,iz);

	     axis2_primplus_temp[UU][gix]+=get_u(p,UU,ix,NY-NCCORRECTPOLAR-1,iz);
	     //cartesian velocities in VX..VZ slots
	     decompose_vels(&get_u(p,0,ix,NY-NCCORRECTPOLAR-1,iz),VX,v,&geom,&geomBL);
	     axis2_primplus_temp[VX][gix]+=v[1];
	     axis2_primplus_temp[VY][gix]+=v[2];
	     axis2_primplus_temp[VZ][gix]+=v[3];
	     //angular velocity in NV slot!!!
	     axis2_primplus_temp[NV][gix]+=get_u(p,VZ,ix,NY-NCCORRECTPOLAR-1,iz);

#ifdef RADIATION
	     axis2_primplus_temp[EE][gix]+=get_u(p,EE,ix,NY-NCCORRECTPOLAR-1,iz);
	     #ifdef NCOMPTONIZATION
	     axis2_primplus_temp[NF][gix]+=get_u(p,NF,ix,NY-NCCORRECTPOLAR-1,iz);
	     #endif
	     //cartesian velocities in FX..FZ slots
	     decompose_vels(&get_u(p,0,ix,NY-NCCORRECTPOLAR-1,iz),FX,v,&geom,&geomBL);
	     //if(ix==10) printf("2 %d > %e | %e %e %e\n",ix,get_u(p,EE,ix,NY-NCCORRECTPOLAR-1,iz),v[1],v[2],v[3]);
	     axis2_primplus_temp[FX][gix]+=v[1];
	     axis2_primplus_temp[FY][gix]+=v[2];
	     axis2_primplus_temp[FZ][gix]+=v[3];
	     //angular velocity in NV+1 slot!!!
	     axis2_primplus_temp[NV+1][gix]+=get_u(p,FZ,ix,NY-NCCORRECTPOLAR-1,iz);
#endif
	   }
       }
     for(iv=0;iv<NV+2;iv++)
       {
	 axis1_primplus_temp[iv][gix]/=NZ;
	 axis2_primplus_temp[iv][gix]/=NZ;
       }
   }

#ifdef MPI
 MPI_Allreduce(&axis1_primplus_temp[0][0], &axis1_primplus[0][0], TNX*(NV+2), MPI_LDOUBLE, MPI_SUM,
		MPI_COMM_WORLD);
 MPI_Allreduce(&axis2_primplus_temp[0][0], &axis2_primplus[0][0], TNX*(NV+2), MPI_LDOUBLE, MPI_SUM,
		MPI_COMM_WORLD);
#else 
  for(ix=0;ix<NX;ix++)
    {
      gix=ix+TOI;
      for(iv=0;iv<NV+2;iv++)
      {
	axis1_primplus[iv][gix]=axis1_primplus_temp[iv][gix];
	axis2_primplus[iv][gix]=axis2_primplus_temp[iv][gix];
      }
    }
#endif
#endif
#endif

  
  /***************************/

return 0;
}


int 
mpi_exchangedata()
{
  //time mark
  struct timespec temp_clock;
  my_clock_gettime(&temp_clock);    
  mid1_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;

#ifdef MPI  
  MPI_Request reqs[MPIMSGBUFSIZE];
  int nreqs=0;
  mpi_senddata(reqs,&nreqs);
  mpi_recvdata(reqs,&nreqs);

  /****************************/
  /*
  // testing which messages are lost
  int count,out1[MPIMSGBUFSIZE],out2[MPIMSGBUFSIZE],tags[MPIMSGBUFSIZE];
  MPI_Status statuses[MPIMSGBUFSIZE];
  for(int i=0;i<MPIMSGBUFSIZE;i++)
    {
      out1[i]=out2[i]=0;
      tags[i]=-1;
    }
  for(;;)
    { 
      MPI_Barrier(MPI_COMM_WORLD);
      getch();
      MPI_Testsome(nreqs,reqs,&count,out1,statuses);
      if(PROCID==7) 
	{
	  printf("\n");
	  for(int i=0;i<count;i++)
	    {
	      out2[out1[i]]=1;
	      tags[out1[i]]=statuses[i].MPI_TAG;
	    }
	  printf("%2d > %3d (%3d) > ",PROCID,count ,nreqs);
	  for(int i=0;i<MPIMSGBUFSIZE;i++)
	      printf("%d ",out2[i]);
	  printf("\n");
	  printf("tags of completed msgs: \n");
	  for(int i=0;i<MPIMSGBUFSIZE;i++)
	    {	    
	      if(out2[i]==0 || 1) 
		{
		  printf("%d ",tags[i]);
		}
	    }
	  fflush(stdout);
	}     
    }
  // end of testing
  */
  /*****************************/

  MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
  mpi_savedata();  
#endif


  my_clock_gettime(&temp_clock);    
  mid2_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;

  return 0;
}

#ifdef MPI //thick bracket comments out recv and send so that non-mpi compiler can compile without MPI_Request
int
mpi_recvdata(MPI_Request *reqs, int *nreqs)
{
  int i,j,k,iv;
  MPI_Status status;
  double temp;
  int verbose=0;
  int tx,ty,tz;
  //upper x
  if(mpi_isitBC(XBCHI)==0)
    {
      tx=TI+1;
      ty=TJ;
      tz=TK;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      MPI_Irecv(msgbufs[6], NY*NZ*NV*NG, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
      if(verbose) printf("%d received MPI_MSG_XLO from %d\n",PROCID,mpi_tile2procid(TI+1,TJ,TK));
   }
  
  //lower x
  if(mpi_isitBC(XBCLO)==0)
    {
      tx=TI-1;
      ty=TJ;
      tz=TK;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      MPI_Irecv(msgbufs[7], NG*NY*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
      if(verbose) printf("%d received MPI_MSG_XHI from %d\n",PROCID,mpi_tile2procid(TI-1,TJ,TK));
    } 

  //upper y
  if(mpi_isitBC(YBCHI)==0)
    {
      tx=TI;
      ty=TJ+1;
      tz=TK;
      #ifdef PERIODIC_YBC
      if(ty>=NTY) ty-=NTY;
      #endif
      MPI_Irecv(msgbufs[8], NX*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
      if(verbose) printf("%d received MPI_MSG_YLO from %d\n",PROCID,mpi_tile2procid(TI,TJ+1,TK));
    }

  //lower y
  if(mpi_isitBC(YBCLO)==0)
    {
      tx=TI;
      ty=TJ-1;
      tz=TK;
      #ifdef PERIODIC_YBC
      if(ty<0) ty+=NTY;
      #endif
      MPI_Irecv(msgbufs[9], NX*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
      if(verbose) printf("%d received MPI_MSG_YHI from %d\n",PROCID,mpi_tile2procid(TI,TJ-1,TK));
    }
  //upper z
  if(mpi_isitBC(ZBCHI)==0)
    {							
      tx=TI;
      ty=TJ;
      tz=TK+1;
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif
      MPI_Irecv(msgbufs[10], NX*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_ZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }
  //lower z
  if(mpi_isitBC(ZBCLO)==0)
    {
      tx=TI;
      ty=TJ;
      tz=TK-1;
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif
      MPI_Irecv(msgbufs[11], NX*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_ZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
      if(verbose) printf("%d received MPI_MSG_ZHI from %d\n",PROCID,mpi_tile2procid(TI,TJ,TK-1));
    }

  //corners
#ifdef MPI4CORNERS
  //elongated along z
  //upper x upper y
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0)
    {
      tx=TI+1;
      ty=TJ+1;
      tz=TK;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty>=NTY) ty-=NTY;
      #endif
      MPI_Irecv(msgbufs[32], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }
  //upper x lower y
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0)
    {
      tx=TI+1;
      ty=TJ-1;
      tz=TK;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty<0) ty+=NTY;
      #endif
      MPI_Irecv(msgbufs[33], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }
  //lower x upper y
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0)
    {
      tx=TI-1;
      ty=TJ+1;
      tz=TK;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty>=NTY) ty-=NTY;
      #endif
      MPI_Irecv(msgbufs[34], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }
  //lower x lower y
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0)
    {
      tx=TI-1;
      ty=TJ-1;
      tz=TK;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty<0) ty+=NTY;
      #endif
      MPI_Irecv(msgbufs[35], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }
   //elongated along y
  //upper x upper z
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      tx=TI+1;
      ty=TJ;
      tz=TK+1;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif
      MPI_Irecv(msgbufs[36], NG*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }
  //upper x lower z
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      tx=TI+1;
      ty=TJ;
      tz=TK-1;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif
      MPI_Irecv(msgbufs[37], NG*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }
  //lower x upper z
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      tx=TI-1;
      ty=TJ;
      tz=TK+1;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif
      MPI_Irecv(msgbufs[38], NG*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }
  //lower x lower z
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      tx=TI-1;
      ty=TJ;
      tz=TK-1;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif
      MPI_Irecv(msgbufs[39], NG*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }
   //elongated along x
  //upper y upper z
  if(mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      tx=TI;
      ty=TJ+1;
      tz=TK+1;
      #ifdef PERIODIC_YBC
      if(ty>=NTY) ty-=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif
      MPI_Irecv(msgbufs[40], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YLOZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }
  //upper y lower z
  if(mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      tx=TI;
      ty=TJ+1;
      tz=TK-1;
      #ifdef PERIODIC_YBC
      if(ty>=NTY) ty-=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif
      MPI_Irecv(msgbufs[41], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YLOZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }
  //lower y upper z
  if(mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      tx=TI;
      ty=TJ-1;
      tz=TK+1;
      #ifdef PERIODIC_YBC
      if(ty<0) ty+=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif
      MPI_Irecv(msgbufs[42], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YHIZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }
  //lower x lower z
  if(mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      tx=TI;
      ty=TJ-1;
      tz=TK-1;
      #ifdef PERIODIC_YBC
      if(ty<0) ty+=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif
      MPI_Irecv(msgbufs[43], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YHIZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }

  /********** corners corners ************/
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      tx=TI+1;
      ty=TJ+1;
      tz=TK+1;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty>=NTY) ty-=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif
      MPI_Irecv(msgbufs[44], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYLOZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }

  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      tx=TI+1;
      ty=TJ+1;
      tz=TK-1;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty>=NTY) ty-=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif
      MPI_Irecv(msgbufs[45], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYLOZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }

 if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      tx=TI+1;
      ty=TJ-1;
      tz=TK+1;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty<0) ty+=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif
      MPI_Irecv(msgbufs[46], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYHIZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }

 if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      tx=TI+1;
      ty=TJ-1;
      tz=TK-1;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty<0) ty+=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif
      MPI_Irecv(msgbufs[47], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYHIZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }

  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      tx=TI-1;
      ty=TJ+1;
      tz=TK+1;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty>=NTY) ty-=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif
      MPI_Irecv(msgbufs[48], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYLOZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }

  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      tx=TI-1;
      ty=TJ+1;
      tz=TK-1;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty>=NTY) ty-=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif
      MPI_Irecv(msgbufs[49], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYLOZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }

 if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      tx=TI-1;
      ty=TJ-1;
      tz=TK+1;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty<0) ty+=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif
      MPI_Irecv(msgbufs[50], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYHIZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }

 if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      tx=TI-1;
      ty=TJ-1;
      tz=TK-1;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty<0) ty+=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif
      MPI_Irecv(msgbufs[51], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYHIZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }

#endif  

  return 0;
}


int
mpi_savedata()
{
  int i,j,k,iv;
  int verbose=0;
  //upper x
  if(mpi_isitBC(XBCHI)==0)
    {
      for(i=NX;i<NX+NG;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[6][(i-NX)*NY*NZ*NV + j*NZ*NV + k*NV + iv]);
	      //printf("r %d > %d %d > %e\n",PROCID,i,j,msgbufs[6][(i-NX)*NY*NZ*NV + j*NZ*NV + k*NV + 0]);
	    }
    }
  
  //lower x
  if(mpi_isitBC(XBCLO)==0)
    {
      for(i=-NG;i<0;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      set_u(p,iv,i,j,k,msgbufs[7][(i+NG)*NY*NZ*NV + j*NZ*NV + k*NV + iv]);
    } 

 //upper y
  if(mpi_isitBC(YBCHI)==0)
    {
      for(i=0;i<NX;i++)
	for(j=NY;j<NY+NG;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      set_u(p,iv,i,j,k,msgbufs[8][i*NG*NZ*NV + (j-NY)*NZ*NV + k*NV + iv]);
    }

  //lower y
  if(mpi_isitBC(YBCLO)==0)
    {
      for(i=0;i<NX;i++)
	for(j=-NG;j<0;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      set_u(p,iv,i,j,k,msgbufs[9][i*NG*NZ*NV + (j+NG)*NZ*NV + k*NV + iv]);
    }

  //upper z
  if(mpi_isitBC(ZBCHI)==0)
    {
      for(i=0;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=NZ;k<NZ+NG;k++)
	    for(iv=0;iv<NV;iv++)
	      set_u(p,iv,i,j,k,msgbufs[10][i*NY*NG*NV + j*NG*NV + (k-NZ)*NV + iv]);
    }
   //lower z
  if(mpi_isitBC(ZBCLO)==0)
    {
      for(i=0;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=-NG;k<0;k++)
	    for(iv=0;iv<NV;iv++)
	      set_u(p,iv,i,j,k,msgbufs[11][i*NY*NG*NV + j*NG*NV + (k+NG)*NV + iv]);
    }

  //corners
#ifdef MPI4CORNERS
  //elongated along z
  //upper x upper y
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0)
    {
      for(i=NX;i<NX+NG;i++)
	for(j=NY;j<NY+NG;j++)
	  for(k=0;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[32][(i-NX)*NG*NZ*NV + (j-NY)*NZ*NV + k*NV + iv]);
	    }
    }
  //upper x lower y
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0)
    {
      for(i=NX;i<NX+NG;i++)
	for(j=-NG;j<0;j++)
	  for(k=0;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[33][(i-NX)*NG*NZ*NV + (j+NG)*NZ*NV + k*NV + iv]);
	    }
    }
  //lower x upper y
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0)
    {
      for(i=-NG;i<0;i++)
	for(j=NY;j<NY+NG;j++)
	  for(k=0;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[34][(i+NG)*NG*NZ*NV + (j-NY)*NZ*NV + k*NV + iv]);
	    }
    }
  //lower x lower y
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0)
    {
      for(i=-NG;i<0;i++)
	for(j=-NG;j<0;j++)
	  for(k=0;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[35][(i+NG)*NG*NZ*NV + (j+NG)*NZ*NV + k*NV + iv]);
	    }
    }
  //elongated along y
  //upper x upper z
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      for(i=NX;i<NX+NG;i++)
	for(j=0;j<NY;j++)
	  for(k=NZ;k<NZ+NG;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[36][(i-NX)*NY*NG*NV + j*NG*NV + (k-NZ)*NV + iv]);
	    }
    }
  //upper x lower z
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      for(i=NX;i<NX+NG;i++)
	for(j=0;j<NY;j++)
	  for(k=-NG;k<0;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[37][(i-NX)*NY*NG*NV + j*NG*NV + (k+NG)*NV + iv]);
	    }
    }
  //lower x upper z
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      for(i=-NG;i<0;i++)
	for(j=0;j<NY;j++)
	  for(k=NZ;k<NZ+NG;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[38][(i+NG)*NY*NG*NV + j*NG*NV + (k-NZ)*NV + iv]);
	    }
    }
  //lower x lower z
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      for(i=-NG;i<0;i++)
	for(j=0;j<NY;j++)
	  for(k=-NG;k<0;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[39][(i+NG)*NY*NG*NV + j*NG*NV + (k+NG)*NV + iv]);
	    }
    }
  //elongated along x
  //upper y upper z
  if(mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      for(i=0;i<NX;i++)
	for(j=NY;j<NY+NG;j++)
	  for(k=NZ;k<NZ+NG;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[40][i*NG*NG*NV + (j-NY)*NG*NV + (k-NZ)*NV + iv]);
	    }
    }
  //upper y lower z
  if(mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      for(i=0;i<NX;i++)
	for(j=NY;j<NY+NG;j++)
	  for(k=-NG;k<0;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[41][i*NG*NG*NV + (j-NY)*NG*NV + (k+NG)*NV + iv]);
	    }
    }
  //lower y upper z
  if(mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      for(i=0;i<NX;i++)
	for(j=-NG;j<0;j++)
	  for(k=NZ;k<NZ+NG;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[42][i*NG*NG*NV + (j+NG)*NG*NV + (k-NZ)*NV + iv]);
	    }
    }
  //lower y lower z
  if(mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      for(i=0;i<NX;i++)
	for(j=-NG;j<0;j++)
	  for(k=-NG;k<0;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[43][i*NG*NG*NV + (j+NG)*NG*NV + (k+NG)*NV + iv]);
	    }
    }
  //corners corners

  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      for(i=NX;i<NX+NG;i++)
	for(j=NY;j<NY+NG;j++)
	  for(k=NZ;k<NZ+NG;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[44][(i-NX)*NG*NG*NV + (j-NY)*NG*NV + (k-NZ)*NV + iv]);
	    }
    }

  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      for(i=NX;i<NX+NG;i++)
	for(j=NY;j<NY+NG;j++)
	  for(k=-NG;k<0;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[45][(i-NX)*NG*NG*NV + (j-NY)*NG*NV + (k+NG)*NV + iv]);
	    }
    }


  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      for(i=NX;i<NX+NG;i++)
	for(j=-NG;j<0;j++)
	  for(k=NZ;k<NZ+NG;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[46][(i-NX)*NG*NG*NV + (j+NG)*NG*NV + (k-NZ)*NV + iv]);
	    }
    }

  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      for(i=NX;i<NX+NG;i++)
	for(j=-NG;j<0;j++)
	  for(k=-NG;k<0;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[47][(i-NX)*NG*NG*NV + (j+NG)*NG*NV + (k+NG)*NV + iv]);
	    }
    }

  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      for(i=-NG;i<0;i++)
	for(j=NY;j<NY+NG;j++)
	  for(k=NZ;k<NZ+NG;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[48][(i+NG)*NG*NG*NV + (j-NY)*NG*NV + (k-NZ)*NV + iv]);
	    }
    }

  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      for(i=-NG;i<0;i++)
	for(j=NY;j<NY+NG;j++)
	  for(k=-NG;k<0;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[49][(i+NG)*NG*NG*NV + (j-NY)*NG*NV + (k+NG)*NV + iv]);
	    }
    }


  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      for(i=-NG;i<0;i++)
	for(j=-NG;j<0;j++)
	  for(k=NZ;k<NZ+NG;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[50][(i+NG)*NG*NG*NV + (j+NG)*NG*NV + (k-NZ)*NV + iv]);
	    }
    }

  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      for(i=-NG;i<0;i++)
	for(j=-NG;j<0;j++)
	  for(k=-NG;k<0;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[51][(i+NG)*NG*NG*NV + (j+NG)*NG*NV + (k+NG)*NV + iv]);
	    }
    }


#endif


  return 0; 
}

int
mpi_senddata(MPI_Request *reqs, int *nreqs)
{  
  int i,j,k,iv;
  int tx,ty,tz;
  int verbose=0;
  ldouble temp;

  //lower x
  if(mpi_isitBC(XBCLO)==0)
    {
      tx=TI-1;
      ty=TJ;
      tz=TK;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif

      for(i=0;i<NG;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[0][(i)*NY*NZ*NV + j*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[0], NY*NZ*NV*NG, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
      if(verbose) printf("%d sent MPI_MSG_XLO to %d\n",PROCID,mpi_tile2procid(tx,ty,tz));
    }

  //upper x
  if(mpi_isitBC(XBCHI)==0)
    {
      tx=TI+1;
      ty=TJ;
      tz=TK;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      for(i=NX-NG;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      msgbufs[1][(i-NX+NG)*NY*NZ*NV + j*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
      MPI_Isend(msgbufs[1],NG*NY*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
      if(verbose) printf("%d sent MPI_MSG_XHI to %d\n",PROCID,mpi_tile2procid(tx,ty,tz));
    }

  //lower y
  if(mpi_isitBC(YBCLO)==0)
    {
      tx=TI;
      ty=TJ-1;
      tz=TK;
      #ifdef PERIODIC_YBC
      if(ty<0) ty+=NTY;
      #endif
      for(i=0;i<NX;i++)
	for(j=0;j<NG;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      msgbufs[2][i*NG*NZ*NV + (j)*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
      MPI_Isend(msgbufs[2], NX*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }

  //upper y
  if(mpi_isitBC(YBCHI)==0)
    {
      tx=TI;
      ty=TJ+1;
      tz=TK;
      #ifdef PERIODIC_YBC
      if(ty>=NTY) ty-=NTY;
      #endif
      for(i=0;i<NX;i++)
	for(j=NY-NG;j<NY;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      msgbufs[3][i*NG*NZ*NV + (j-NY+NG)*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
      MPI_Isend(msgbufs[3], NX*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }
  //lower z
  if(mpi_isitBC(ZBCLO)==0)
    {
      tx=TI;
      ty=TJ;
      tz=TK-1;
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif

      for(i=0;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NG;k++)
	    for(iv=0;iv<NV;iv++)
	      msgbufs[4][i*NY*NG*NV + j*NG*NV + (k)*NV + iv]=get_u(p,iv,i,j,k);
      MPI_Isend(msgbufs[4], NX*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_ZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }
  //upper z
  if(mpi_isitBC(ZBCHI)==0)
    {
      tx=TI;
      ty=TJ;
      tz=TK+1;
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif
      for(i=0;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=NZ-NG;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      msgbufs[5][i*NY*NG*NV + j*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
      MPI_Isend(msgbufs[5], NX*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_ZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }

#ifdef MPI4CORNERS

  /***************************/
  //elongated corners - along z
  //lower x lower y
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0)
    {
      tx=TI-1;
      ty=TJ-1;
      tz=TK;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty<0) ty+=NTY;
      #endif

      for(i=0;i<NG;i++)
	for(j=0;j<NG;j++)
	  for(k=0;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[12][i*NG*NZ*NV + j*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[12], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }  

  //lower x higher y
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0)
    {
      tx=TI-1;
      ty=TJ+1;
      tz=TK;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty>=NTY) ty-=NTY;
      #endif

      for(i=0;i<NG;i++)
	for(j=NY-NG;j<NY;j++)
	  for(k=0;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[13][i*NG*NZ*NV + (j-NY+NG)*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[13], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }
  
  //higher x lower y
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0)
    {
      tx=TI+1;
      ty=TJ-1;
      tz=TK;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty<0) ty+=NTY;
      #endif

      for(i=NX-NG;i<NX;i++)
	for(j=0;j<NG;j++)
	  for(k=0;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[14][(i-NX+NG)*NG*NZ*NV + j*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[14], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }  

  //higher x higher y
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0)
    {
      tx=TI+1;
      ty=TJ+1;
      tz=TK;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty>=NTY) ty-=NTY;
      #endif

      for(i=NX-NG;i<NX;i++)
	for(j=NY-NG;j<NY;j++)
	  for(k=0;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[15][(i-NX+NG)*NG*NZ*NV + (j-NY+NG)*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[15], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }  

  /***************************/
  //elongated corners - along y
  //lower x lower z
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      tx=TI-1;
      ty=TJ;
      tz=TK-1;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif

      for(i=0;i<NG;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NG;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[16][i*NY*NG*NV + j*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[16], NG*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }  

  //lower x higher z
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      tx=TI-1;
      ty=TJ;
      tz=TK+1;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif

      for(i=0;i<NG;i++)
	for(j=0;j<NY;j++)
	  for(k=NZ-NG;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[17][i*NY*NG*NV + j*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[17], NG*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }
  
  //higher x lower z
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      tx=TI+1;
      ty=TJ;
      tz=TK-1;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif

      for(i=NX-NG;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NG;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[18][(i-NX+NG)*NY*NG*NV + j*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[18], NG*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }  

  //higher x higher z
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      tx=TI+1;
      ty=TJ;
      tz=TK+1;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif

      for(i=NX-NG;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=NZ-NG;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[19][(i-NX+NG)*NY*NG*NV + j*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[19], NG*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }  

   /***************************/
  //elongated corners - along x
  //lower y lower z
  if(mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      tx=TI;
      ty=TJ-1;
      tz=TK-1;
      #ifdef PERIODIC_YBC
      if(ty<0) ty+=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif

      for(i=0;i<NX;i++)
	for(j=0;j<NG;j++)
	  for(k=0;k<NG;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[20][i*NG*NG*NV + j*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[20], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YLOZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }  

  //lower y higher z
  if(mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      tx=TI;
      ty=TJ-1;
      tz=TK+1;
      #ifdef PERIODIC_YBC
      if(ty<0) ty+=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif

      for(i=0;i<NX;i++)
	for(j=0;j<NG;j++)
	  for(k=NZ-NG;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[21][i*NG*NG*NV + j*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[21], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YLOZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }
  
  //higher y lower z
  if(mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      tx=TI;
      ty=TJ+1;
      tz=TK-1;
      #ifdef PERIODIC_YBC
      if(ty>=NTY) ty-=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif

      for(i=0;i<NX;i++)
	for(j=NY-NG;j<NY;j++)
	  for(k=0;k<NG;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[22][i*NG*NG*NV + (j-NY+NG)*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[22], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YHIZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }  

  //higher y higher z
  if(mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      tx=TI;
      ty=TJ+1;
      tz=TK+1;
      #ifdef PERIODIC_YBC
      if(ty>=NTY) ty-=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif

      for(i=0;i<NX;i++)
	for(j=NY-NG;j<NY;j++)
	  for(k=NZ-NG;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[23][i*NG*NG*NV + (j-NY+NG)*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[23], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YHIZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }  

  /***************************/
  //corners corners 
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      tx=TI-1;
      ty=TJ-1;
      tz=TK-1;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty<0) ty+=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif

      for(i=0;i<NG;i++)
	for(j=0;j<NG;j++)
	  for(k=0;k<NG;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[24][i*NG*NG*NV + j*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[24], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYLOZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }  
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      tx=TI-1;
      ty=TJ-1;
      tz=TK+1;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty<0) ty+=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif

      for(i=0;i<NG;i++)
	for(j=0;j<NG;j++)
	  for(k=NZ-NG;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[25][i*NG*NG*NV + j*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[25], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYLOZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }  
 if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      tx=TI-1;
      ty=TJ+1;
      tz=TK-1;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty>=NTY) ty-=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif

      for(i=0;i<NG;i++)
	for(j=NY-NG;j<NY;j++)
	  for(k=0;k<NG;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[26][i*NG*NG*NV + (j-NY+NG)*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[26], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYHIZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }  
 if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      tx=TI-1;
      ty=TJ+1;
      tz=TK+1;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty>=NTY) ty-=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif

      for(i=0;i<NG;i++)
	for(j=NY-NG;j<NY;j++)
	  for(k=NZ-NG;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[27][i*NG*NG*NV + (j-NY+NG)*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[27], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYHIZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }  
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      tx=TI+1;
      ty=TJ-1;
      tz=TK-1;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty<0) ty+=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif

      for(i=NX-NG;i<NX;i++)
	for(j=0;j<NG;j++)
	  for(k=0;k<NG;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[28][(i-NX+NG)*NG*NG*NV + j*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[28], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYLOZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }  
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      tx=TI+1;
      ty=TJ-1;
      tz=TK+1;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty<0) ty+=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif

      for(i=NX-NG;i<NX;i++)
	for(j=0;j<NG;j++)
	  for(k=NZ-NG;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[29][(i-NX+NG)*NG*NG*NV + j*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[29], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYLOZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }  
 if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      tx=TI+1;
      ty=TJ+1;
      tz=TK-1;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty>=NTY) ty-=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif

      for(i=NX-NG;i<NX;i++)
	for(j=NY-NG;j<NY;j++)
	  for(k=0;k<NG;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[30][(i-NX+NG)*NG*NG*NV + (j-NY+NG)*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[30], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYHIZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }  
 if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      tx=TI+1;
      ty=TJ+1;
      tz=TK+1;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      #ifdef PERIODIC_YBC
      if(ty>=NTY) ty-=NTY;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif

      for(i=NX-NG;i<NX;i++)
	for(j=NY-NG;j<NY;j++)
	  for(k=NZ-NG;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[31][(i-NX+NG)*NG*NG*NV + (j-NY+NG)*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[31], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYHIZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }  
#endif



  return 0;
}
#endif

//verify if there is real BC at all, if set_bc() needed
int
mpi_hasBC()
{
#ifndef MPI
  return 1;
#else
  if(TI==0 || TI==NTX-1 || TJ==0 || TJ==NTY-1 || TK==0 || TK==NTZ-1)
    return 1; //this cell has some real BC
  else
    return 0;
#endif
}

//verify if given cell from set_bc() falls into a real BC
int
mpi_isitBC(int BCtype)
{
#ifndef MPI
  return 1;
#else //check here if real BC
  int perx,pery,perz; //is periodic in x,y,z?
  perx=pery=perz=0; 
  #ifdef PERIODIC_XBC
  perx=1;
  #endif
  #ifdef PERIODIC_YBC
  pery=1;
  #endif
  #ifdef PERIODIC_ZBC
  perz=1;
  #endif

  if(BCtype==XBCLO && TI==0 && perx==0)
    return 1;
  if(BCtype==XBCHI && TI==NTX-1 && perx==0)
    return 1;
  if(BCtype==YBCLO && TJ==0 && pery==0)
    return 1;
  if(BCtype==YBCHI && TJ==NTY-1 && pery==0)
    return 1;
  if(BCtype==ZBCLO && TK==0 && perz==0)
    return 1;
  if(BCtype==ZBCHI && TK==NTZ-1 && perz==0)
    return 1;

  return 0; 
#endif
}

void
mpi_synchtiming(ldouble *time)
{
#ifdef MPI

  MPI_Barrier(MPI_COMM_WORLD);
  
  global_time=*time;
  
  MPI_Allreduce(&tstepdenmax, &global_tstepdenmax, 1, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);

  #ifdef SELFTIMESTEP
  MPI_Allreduce(&tstepdenmin, &global_tstepdenmin, 1, MPI_DOUBLE, MPI_MIN,
                MPI_COMM_WORLD);
  #endif

  //maximal time taken by information exchange
  ldouble localmp_time=mid2_time-mid1_time;
  //MPI_Allreduce(&localmp_time, &maxmp_time, 1, MPI_DOUBLE, MPI_MAX,        MPI_COMM_WORLD);   
  //only PROCID==0 writes to stdout
  MPI_Reduce(&localmp_time, &maxmp_time, 1, MPI_DOUBLE, MPI_MAX,     0,   MPI_COMM_WORLD);   

  //total operation time
  ldouble local_u2ptime=end_u2ptime-start_u2ptime;
  struct {
    double time;
    int ti;
  } in,outmin,outmax;
  in.time=local_u2ptime;
  in.ti=TI;

  MPI_Reduce(&in, &outmax, 1, MPI_DOUBLE_INT, MPI_MAXLOC,0,
		 MPI_COMM_WORLD);  
  MPI_Reduce(&in, &outmin, 1, MPI_DOUBLE_INT, MPI_MINLOC,0,
		 MPI_COMM_WORLD);  

  max_u2ptime=outmax.time;
  max_u2ptime_loc=outmax.ti;
  min_u2ptime=outmin.time;
  min_u2ptime_loc=outmin.ti;  

  MPI_Bcast(&global_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  MPI_Barrier(MPI_COMM_WORLD);

  *time=global_time;
  tstepdenmax=global_tstepdenmax;
  tstepdenmin=global_tstepdenmin;
 #endif
}

void
mpi_myinit(int argc, char *argv[])
{
#ifdef MPI
  int i,j,k;

  //check for conflicts in declarations
  #ifndef OUTPUTPERCORE
  #ifdef RESOUTPUT_ASCII
  my_err("RESOUTPUT_ASCII requires MPI_OUTPUTPERCORE\n");exit(-1);
  #endif
  #ifdef AVGOUTPUT_ASCII
  my_err("RESOUTPUT_ASCII requires MPI_OUTPUTPERCORE\n");exit(-1);
  #endif
  #endif
  #ifdef CALCHRONTHEGO
  if(TNZ>1) //3D?
    {
      if(PROCID==0) printf("CALCHRONTHEGO not implemented for 3D.\n");
      exit(-1);
    }
  #endif

  //initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &PROCID);
  MPI_Comm_size(MPI_COMM_WORLD, &NPROCS);

  if(NPROCS!=NTX*NTY*NTZ)
    {
      if(PROCID==0) printf("Wrong number of processes. Problem set up for: %d x %d x %d = %d processes.\n",NTX,NTY,NTZ,NTX*NTY*NTZ);
      exit(-1);
    }
  
  mpi_procid2tile(PROCID,&TI,&TJ,&TK);
  mpi_tileorigin(TI,TJ,TK,&TOI,&TOJ,&TOK);

  if(PROCID==0) printf("pid: %d/%d; tot.res: %dx%dx%d; tile.res:  %dx%dx%d\n"
	 "tile: %d,%d,%d; tile orig.: %d,%d,%d\n",PROCID,NPROCS,TNX,TNY,TNZ,NX,NY,NZ,TI,TJ,TK,TOI,TOJ,TOK);

#else
  TI=TJ=TK=0;
  TOI=TOJ=TOK=0;
  PROCID=0;
  NPROCS=1;
#endif
}

void
mpi_myfinalize()
{
#ifdef MPI
  MPI_Finalize();
#endif
}

void
mpi_tileorigin(int ti, int tj, int tk, int* toi, int* toj, int* tok)
{
  *toi = ti * (TNX/NTX);
  *toj = tj * (TNY/NTY);
  *tok = tk * (TNZ/NTZ);
}

void
mpi_global2localidx(int gix,int giy, int giz, int *lix, int *liy, int *liz)
{
  #ifdef MPI
  int tilei,tilej,tilek;
  int toi,toj,tok;
  mpi_procid2tile(PROCID,&tilei,&tilej,&tilek);
  mpi_tileorigin(tilei,tilej,tilek,&toi,&toj,&tok);
  *lix = gix - toi;
  *liy = giy - toj;
  *liz = giz - tok;
  #else
  *lix = gix;
  *liy = giy;
  *liz = giz;
  #endif
}


void
mpi_local2globalidx(int lix,int liy, int liz, int *gix, int *giy, int *giz)
{
  #ifdef MPI
  int tilei,tilej,tilek;
  int toi,toj,tok;
  mpi_procid2tile(PROCID,&tilei,&tilej,&tilek);
  mpi_tileorigin(tilei,tilej,tilek,&toi,&toj,&tok);
  *gix = lix + toi;
  *giy = liy + toj;
  *giz = liz + tok;
  #else
  *gix = lix;
  *giy = liy;
  *giz = liz;
  #endif
}

void
mpi_procid2tile(int procid, int* tilei, int* tilej, int* tilek)
{
  *tilek = floor(procid / (NTX * NTY));
  *tilej = floor((procid - (*tilek) * NTX * NTY) / NTX);
  *tilei = procid - NTX * NTY * (*tilek) - NTX * (*tilej);
}

int
mpi_tile2procid(int tilei, int tilej, int tilek)
{
  return tilek * NTX * NTY + tilej * NTX + tilei;
}

//parasitic openMP 

int
omp_myinit()
{
#ifdef OMP
  printf("MPI does not work with OMP.\n"); exit(-1);
#ifdef MPI
  printf("MPI does not work with OMP.\n"); exit(-1);
#endif

  //omp_set_dynamic(0);
  //omp_set_num_threads(NTX*NTY*NTZ);

  NPROCS=omp_get_num_threads();
  PROCID=0;
  TI=TJ=TK=0;
  TOI=TOJ=TOK=0;
  PROCID=0;  

#else
  NPROCS=1;
  TI=TJ=TK=0;
  TOI=TOJ=TOK=0;
  PROCID=0;  
#endif

  return 0;
}

