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
      MPI_Irecv(msgbufs[16], NG*NG*NZ*NV, MPI_DOUBLE,
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
      MPI_Irecv(msgbufs[17], NG*NG*NZ*NV, MPI_DOUBLE,
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
      MPI_Irecv(msgbufs[18], NG*NG*NZ*NV, MPI_DOUBLE,
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
      MPI_Irecv(msgbufs[19], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
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
  //upper x upper y
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0)
    {
      for(i=NX;i<NX+NG;i++)
	for(j=NY;j<NY+NG;j++)
	  for(k=0;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[16][(i-NX)*NG*NZ*NV + (j-NY)*NZ*NV + k*NV + iv]);
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
		set_u(p,iv,i,j,k,msgbufs[17][(i-NX)*NG*NZ*NV + (j+NG)*NZ*NV + k*NV + iv]);
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
		set_u(p,iv,i,j,k,msgbufs[18][(i+NG)*NG*NZ*NV + (j-NY)*NZ*NV + k*NV + iv]);
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
		set_u(p,iv,i,j,k,msgbufs[19][(i+NG)*NG*NZ*NV + (j+NG)*NZ*NV + k*NV + iv]);
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
  if(TNX>1 && TNY>1 && TNZ>1)
    {
      printf("message passing of corners for flux_ct not implemented for 3D yet. sorry.\n");
      exit(-1);
    }

  if(TNZ>1)
    {
      printf("message passing of corners for flux_ct not implemented for TNZ>1 yet. sorry.\n");
      exit(-1);
    }

  //corners
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
		mpi_tile2procid(TI+1,TJ+1,TK), MPI_MSG_XHIYHI, MPI_COMM_WORLD,&reqs[*nreqs]);
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

  //maximal time taken by information exchange
  ldouble localmp_time=mid2_time-mid1_time;
  MPI_Allreduce(&localmp_time, &maxmp_time, 1, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);   

  MPI_Bcast(&global_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  MPI_Barrier(MPI_COMM_WORLD);

  *time=global_time;
  tstepdenmax=global_tstepdenmax;
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
  #ifdef SUBZONES
  printf("SUBZONES do not work with OMP.\n"); exit(-1);
  #endif
  #ifdef MPI
  printf("MPI does not work with OMP.\n"); exit(-1);
  #endif

  omp_set_dynamic(0);
  omp_set_num_threads(NTX*NTY*NTZ);

  #pragma omp parallel
  {
    NPROCS=omp_get_num_threads();
    PROCID=omp_get_thread_num();
    if(NPROCS!=NTX*NTY*NTZ)
      {
	if(PROCID==0) 
	  {
	    printf("Wrong number of threads. Problem set up for: %d x %d x %d = %d threads (openMP).\n",NTX,NTY,NTZ,NTX*NTY*NTZ);
	    exit(-1);
	  }
      } 
    mpi_procid2tile(PROCID,&TI,&TJ,&TK);
    mpi_tileorigin(TI,TJ,TK,&TOI,&TOJ,&TOK);
    if(PROCID==0) printf("tid: %d/%d; tot.res: %dx%dx%d; tile.res:  %dx%dx%d\n"
			 "tile: %d,%d,%d; tile orig.: %d,%d,%d\n",PROCID,NPROCS,TNX,TNY,TNZ,NX,NY,NZ,TI,TJ,TK,TOI,TOJ,TOK);
  }
#else
  NPROCS=1;
  TI=TJ=TK=0;
  TOI=TOJ=TOK=0;
  PROCID=0;  
#endif

  return 0;
}

