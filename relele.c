//KORAL - relele.c
//some relativistic routines

#include "ko.h"

//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts velocities u1->u2
int
conv_vels(ldouble *u1,ldouble *u2,int which1,int which2,ldouble gg[][5],ldouble GG[][5])
{
  int i,j;
  ldouble ut[4];

  /*************** VEL3 -> VEL3 ***************/
  if(which1==VEL3 && which2==VEL3)
    {
      for(i=0;i<4;i++) ut[i]=u1[i];      
    }
  /*************** VEL4 -> VEL4 ***************/
  else if(which1==VEL4 && which2==VEL4)
    {
      //assumes u^t unknown
      ldouble a,b,c;
      a=gg[0][0];
      b=0.;
      c=1.;
      for(i=1;i<4;i++)
	{
	  b+=2.*u1[i]*gg[0][i];
	  for(j=1;j<4;j++)
	    {
	      c+=u1[i]*u1[j]*gg[i][j];
	    }
	}
      ldouble delta=b*b-4.*a*c;
      if(delta<0.) my_err("delta.lt.0 in VEL4->VEL4\n");
      ut[0]=(-b-sqrt(delta))/2./a;
      //TODO: more strict criterion
      if(ut[0]<0.) ut[0]=(-b+sqrt(delta))/2./a;
     
      for(i=1;i<4;i++) ut[i]=u1[i];      
    }
  /*************** VELR -> VELR ***************/
  else if(which1==VELR && which2==VELR)
    {
      for(i=0;i<4;i++) ut[i]=u1[i];      
    }
  /*************** VEL4 -> VEL3 ***************/
  else if(which1==VEL4 && which2==VEL3)
    {
      //assumes u^t unknown
      ldouble a,b,c;
      a=gg[0][0];
      b=0.;
      c=1.;
      for(i=1;i<4;i++)
	{
	  b+=2.*u1[i]*gg[0][i];
	  for(j=1;j<4;j++)
	    {
	      c+=u1[i]*u1[j]*gg[i][j];
	    }
	}
      ldouble delta=b*b-4.*a*c;
      if(delta<0.) my_err("delta.lt.0 in VEL4->VEL4\n");
      ut[0]=(-b-sqrt(delta))/2./a;
      if(ut[0]<1.) ut[0]=(-b+sqrt(delta))/2./a;

      for(i=1;i<4;i++) ut[i]=u1[i]/ut[0];      
    }
  /*************** VEL3 -> VEL4 ***************/
  else if(which1==VEL3 && which2==VEL4)
    {
      ldouble a,b,c;
      a=b=0.;
      for(i=1;i<4;i++)
	{
	  a+=2.*u1[i]*gg[0][i];
	  for(j=1;j<4;j++)
	    {
	      b+=u1[i]*u1[j]*gg[i][j];
	    }
	}
						
      ut[0]=sqrt(-1./(gg[0][0]+a+b));
      if(ut[0]<1. || isnan(ut[0]))
	{
	  print_4vector(u1);
	  printf("ut.nan in conv_vels(%d,%d)\n",which1,which2); getchar();
	}
      ut[1]=u1[1]*ut[0];
      ut[2]=u1[2]*ut[0];
      ut[3]=u1[3]*ut[0];      
    }
  /*************** VEL3 -> VELR ***************/
  else if(which1==VEL3 && which2==VELR)
    {
      ldouble a,b,c;
      a=b=0.;
      for(i=1;i<4;i++)
	{
	  a+=2.*u1[i]*gg[0][i];
	  for(j=1;j<4;j++)
	    {
	      b+=u1[i]*u1[j]*gg[i][j];
	    }
	}
						
      ut[0]=sqrt(-1./(gg[0][0]+a+b));
      if(ut[0]<1. || isnan(ut[0]))
	{
	  printf("ut.nan in conv_vels(%d,%d)\n",which1,which2); getchar();
	}
      //to 4-velocity
      ut[1]=u1[1]*ut[0];
      ut[2]=u1[2]*ut[0];
      ut[3]=u1[3]*ut[0];   
      //to relative velocity
      for(i=1;i<4;i++)
	ut[i]=ut[i]-ut[0]*GG[0][i]/GG[0][0];
    }
  /*************** VEL4 -> VELR ***************/
  else if (which1==VEL4 && which2==VELR)
    {
      //assumes u^t unknown
      ldouble a,b,c;
      a=gg[0][0];
      b=0.;
      c=1.;
      for(i=1;i<4;i++)
	{
	  b+=2.*u1[i]*gg[0][i];
	  for(j=1;j<4;j++)
	    {
	      c+=u1[i]*u1[j]*gg[i][j];
	    }
	}
      ldouble delta=b*b-4.*a*c;
      if(delta<0.) my_err("delta.lt.0 in VEL4->VEL4\n");
      ut[0]=(-b-sqrt(delta))/2./a;
      if(ut[0]<1.) ut[0]=(-b+sqrt(delta))/2./a;

      for(i=1;i<4;i++)
	ut[i]=u1[i]-ut[0]*GG[0][i]/GG[0][0];
    }

  /*************** VELR -> VEL4 ***************/
  else if (which1==VELR && which2==VEL4)
    {
      ldouble qsq=0.;
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  qsq+=u1[i]*u1[j]*gg[i][j];
      ldouble gamma2=1.+qsq;
      ldouble alpha2=-1./GG[0][0];
      ut[0]=sqrt(gamma2/alpha2);
      for(i=1;i<4;i++)
	ut[i]=u1[i]+ut[0]*GG[0][i]/GG[0][0];
    }
  /*************** VELR -> VEL3 ***************/
  else if (which1==VELR && which2==VEL3)
    {
      ldouble qsq=0.;
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  qsq+=u1[i]*u1[j]*gg[i][j];
      ldouble gamma2=1.+qsq;
      ldouble alpha2=-1./GG[0][0];
      ut[0]=sqrt(gamma2/alpha2);
      for(i=1;i<4;i++)
	ut[i]=u1[i]+ut[0]*GG[0][i]/GG[0][0];
      ut[1]/=ut[0];
      ut[2]/=ut[0];
      ut[3]/=ut[0];
    }

  /*************** not supported  ***************/
  else
    {
      my_err("velocity conversion not supported.\n");
      return -1;
    }

  //write to index 0 only when needed
  if(which2==VEL4)
    u2[0]=ut[0];
  u2[1]=ut[1];
  u2[2]=ut[2];
  u2[3]=ut[3];
  
  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts hydro velocities as above but takes full vector of primitives
//as an input and outputs to primitives the spatial components;
int
conv_velsinprims(ldouble *pp,int which1, int which2,ldouble gg[][5],ldouble GG[][5])
{
  ldouble v1[4],v2[4];
  v1[0]=0.; //not considered
  v1[1]=pp[2];
  v1[2]=pp[3];
  v1[3]=pp[4];
  conv_vels(v1,v2,which1,which2,gg,GG);
  pp[2]=v2[1];
  pp[3]=v2[2];
  pp[4]=v2[3];
  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//returns contravariant four-velocity of a normal observer
//n_mu=(-alp,0,0,0)
int
calc_normalobs_4vel(ldouble GG[][5], ldouble *ncon)
{
  ldouble alp=1./sqrt(-GG[0][0]);
  ldouble ncov[4]={-alp,0.,0.,0.};
  indices_12(ncov,ncon,GG);
  return 0.;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//returns hydro primitives for an atmosphere
//velocities already in VELPRIM
int
set_hdatmosphere(ldouble *pp,ldouble *xx,ldouble gg[][5],ldouble GG[][5],int atmtype)
{
  if(atmtype==0)
    {
      //normal observer
      ldouble ucon[4];
      calc_normalobs_4vel(GG,ucon);
      conv_vels(ucon,ucon,VEL4,VELPRIM,gg,GG);
      pp[2]=ucon[1];
      pp[3]=ucon[2];
      pp[4]=ucon[3];
      //something, to be changed by user in problem specific file
      pp[0]=1.;
      pp[1]=0.001;
      return 0;
    }
  else
    my_err("atmtype value not handled in set_atmosphere()\n");
  return 0.;
}





//**********************************************************************
//**********************************************************************
//**********************************************************************
//picks a tensor from cell face arr at ix,iy,iz
int
pick_Tb(ldouble* arr,int ix,int iy,int iz,int idim,ldouble T[][4])
{
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      T[i][j]=get_Tb(arr,i,j,ix,iy,iz,idim);
  return 0;
}

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
//picks a tensor from arr at ix,iy,iz
int
pick_T(ldouble* arr,int ix,int iy,int iz,ldouble T[][4])
{
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      T[i][j]=get_T(arr,i,j,ix,iy,iz);
  return 0;
}

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
//picks a metric at ix,iy,iz
int
pick_g(int ix,int iy,int iz,ldouble gg[][5])
{
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<5;j++)
      gg[i][j]=get_g(g,i,j,ix,iy,iz);

  return 0;
}

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
//picks an inversed metric at ix,iy,iz
int
pick_G(int ix,int iy,int iz,ldouble GG[][5])
{
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<5;j++)
      GG[i][j]=get_g(G,i,j,ix,iy,iz);

  return 0;
}

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
//picks a metric at cell faces
int
pick_gb(int ix,int iy,int iz,int idim,ldouble gg[][5])
{
  ldouble g00,g03,g11,g22,g33,gdet,dlgdet0,dlgdet1,dlgdet2;

  //ix,iy,iz correspond to indices in cell-faces arrays
  if(idim==0)
    {
      gg[0][0]=get_gb(gbx,0,0,ix,iy,iz,0);
      gg[0][1]=get_gb(gbx,0,1,ix,iy,iz,0);
      gg[0][2]=get_gb(gbx,0,2,ix,iy,iz,0);
      gg[0][3]=get_gb(gbx,0,3,ix,iy,iz,0);
      gg[1][0]=get_gb(gbx,1,0,ix,iy,iz,0);
      gg[1][1]=get_gb(gbx,1,1,ix,iy,iz,0);
      gg[1][2]=get_gb(gbx,1,2,ix,iy,iz,0);
      gg[1][3]=get_gb(gbx,1,3,ix,iy,iz,0);
      gg[2][0]=get_gb(gbx,2,0,ix,iy,iz,0);  
      gg[2][1]=get_gb(gbx,2,1,ix,iy,iz,0);  
      gg[2][2]=get_gb(gbx,2,2,ix,iy,iz,0);  
      gg[2][3]=get_gb(gbx,2,3,ix,iy,iz,0);  
      gg[3][0]=get_gb(gbx,3,0,ix,iy,iz,0);
      gg[3][1]=get_gb(gbx,3,1,ix,iy,iz,0);
      gg[3][2]=get_gb(gbx,3,2,ix,iy,iz,0);
      gg[3][3]=get_gb(gbx,3,3,ix,iy,iz,0);
      gg[0][4]=get_gb(gbx,0,4,ix,iy,iz,0);
      gg[1][4]=get_gb(gbx,1,4,ix,iy,iz,0);
      gg[2][4]=get_gb(gbx,2,4,ix,iy,iz,0);
      gg[3][4]=get_gb(gbx,3,4,ix,iy,iz,0);      
    }
  if(idim==1)
    {
      gg[0][0]=get_gb(gby,0,0,ix,iy,iz,1);
      gg[0][1]=get_gb(gby,0,1,ix,iy,iz,1);
      gg[0][2]=get_gb(gby,0,2,ix,iy,iz,1);
      gg[0][3]=get_gb(gby,0,3,ix,iy,iz,1);
      gg[1][0]=get_gb(gby,1,0,ix,iy,iz,1);
      gg[1][1]=get_gb(gby,1,1,ix,iy,iz,1);
      gg[1][2]=get_gb(gby,1,2,ix,iy,iz,1);
      gg[1][3]=get_gb(gby,1,3,ix,iy,iz,1);
      gg[2][0]=get_gb(gby,2,0,ix,iy,iz,1);  
      gg[2][1]=get_gb(gby,2,1,ix,iy,iz,1);  
      gg[2][2]=get_gb(gby,2,2,ix,iy,iz,1);  
      gg[2][3]=get_gb(gby,2,3,ix,iy,iz,1);  
      gg[3][0]=get_gb(gby,3,0,ix,iy,iz,1);
      gg[3][1]=get_gb(gby,3,1,ix,iy,iz,1);
      gg[3][2]=get_gb(gby,3,2,ix,iy,iz,1);
      gg[3][3]=get_gb(gby,3,3,ix,iy,iz,1);
      gg[0][4]=get_gb(gby,0,4,ix,iy,iz,1);
      gg[1][4]=get_gb(gby,1,4,ix,iy,iz,1);
      gg[2][4]=get_gb(gby,2,4,ix,iy,iz,1);
      gg[3][4]=get_gb(gby,3,4,ix,iy,iz,1);      
    }
  if(idim==2)
    {
      gg[0][0]=get_gb(gbz,0,0,ix,iy,iz,2);
      gg[0][1]=get_gb(gbz,0,1,ix,iy,iz,2);
      gg[0][2]=get_gb(gbz,0,2,ix,iy,iz,2);
      gg[0][3]=get_gb(gbz,0,3,ix,iy,iz,2);
      gg[1][0]=get_gb(gbz,1,0,ix,iy,iz,2);
      gg[1][1]=get_gb(gbz,1,1,ix,iy,iz,2);
      gg[1][2]=get_gb(gbz,1,2,ix,iy,iz,2);
      gg[1][3]=get_gb(gbz,1,3,ix,iy,iz,2);
      gg[2][0]=get_gb(gbz,2,0,ix,iy,iz,2);  
      gg[2][1]=get_gb(gbz,2,1,ix,iy,iz,2);  
      gg[2][2]=get_gb(gbz,2,2,ix,iy,iz,2);  
      gg[2][3]=get_gb(gbz,2,3,ix,iy,iz,2);  
      gg[3][0]=get_gb(gbz,3,0,ix,iy,iz,2);
      gg[3][1]=get_gb(gbz,3,1,ix,iy,iz,2);
      gg[3][2]=get_gb(gbz,3,2,ix,iy,iz,2);
      gg[3][3]=get_gb(gbz,3,3,ix,iy,iz,2);
      gg[0][4]=get_gb(gbz,0,4,ix,iy,iz,2);
      gg[1][4]=get_gb(gbz,1,4,ix,iy,iz,2);
      gg[2][4]=get_gb(gbz,2,4,ix,iy,iz,2);
      gg[3][4]=get_gb(gbz,3,4,ix,iy,iz,2);      
  }  

  return 0;

}

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
//picks an inversed metric at cell faces
int
pick_Gb(int ix,int iy,int iz,int idim,ldouble gg[][5])
{
  ldouble g00,g03,g11,g22,g33,gdet,dlgdet0,dlgdet1,dlgdet2;

  //ix,iy,iz correspond to indices in cell-faces arrays
  if(idim==0)
    {
      gg[0][0]=get_gb(Gbx,0,0,ix,iy,iz,0);
      gg[0][1]=get_gb(Gbx,0,1,ix,iy,iz,0);
      gg[0][2]=get_gb(Gbx,0,2,ix,iy,iz,0);
      gg[0][3]=get_gb(Gbx,0,3,ix,iy,iz,0);
      gg[1][0]=get_gb(Gbx,1,0,ix,iy,iz,0);
      gg[1][1]=get_gb(Gbx,1,1,ix,iy,iz,0);
      gg[1][2]=get_gb(Gbx,1,2,ix,iy,iz,0);
      gg[1][3]=get_gb(Gbx,1,3,ix,iy,iz,0);
      gg[2][0]=get_gb(Gbx,2,0,ix,iy,iz,0);  
      gg[2][1]=get_gb(Gbx,2,1,ix,iy,iz,0);  
      gg[2][2]=get_gb(Gbx,2,2,ix,iy,iz,0);  
      gg[2][3]=get_gb(Gbx,2,3,ix,iy,iz,0);  
      gg[3][0]=get_gb(Gbx,3,0,ix,iy,iz,0);
      gg[3][1]=get_gb(Gbx,3,1,ix,iy,iz,0);
      gg[3][2]=get_gb(Gbx,3,2,ix,iy,iz,0);
      gg[3][3]=get_gb(Gbx,3,3,ix,iy,iz,0);
      gg[0][4]=get_gb(Gbx,0,4,ix,iy,iz,0);
      gg[1][4]=get_gb(Gbx,1,4,ix,iy,iz,0);
      gg[2][4]=get_gb(Gbx,2,4,ix,iy,iz,0);
      gg[3][4]=get_gb(Gbx,3,4,ix,iy,iz,0);      
    }
  if(idim==1)
    {
      gg[0][0]=get_gb(Gby,0,0,ix,iy,iz,1);
      gg[0][1]=get_gb(Gby,0,1,ix,iy,iz,1);
      gg[0][2]=get_gb(Gby,0,2,ix,iy,iz,1);
      gg[0][3]=get_gb(Gby,0,3,ix,iy,iz,1);
      gg[1][0]=get_gb(Gby,1,0,ix,iy,iz,1);
      gg[1][1]=get_gb(Gby,1,1,ix,iy,iz,1);
      gg[1][2]=get_gb(Gby,1,2,ix,iy,iz,1);
      gg[1][3]=get_gb(Gby,1,3,ix,iy,iz,1);
      gg[2][0]=get_gb(Gby,2,0,ix,iy,iz,1);  
      gg[2][1]=get_gb(Gby,2,1,ix,iy,iz,1);  
      gg[2][2]=get_gb(Gby,2,2,ix,iy,iz,1);  
      gg[2][3]=get_gb(Gby,2,3,ix,iy,iz,1);  
      gg[3][0]=get_gb(Gby,3,0,ix,iy,iz,1);
      gg[3][1]=get_gb(Gby,3,1,ix,iy,iz,1);
      gg[3][2]=get_gb(Gby,3,2,ix,iy,iz,1);
      gg[3][3]=get_gb(Gby,3,3,ix,iy,iz,1);
      gg[0][4]=get_gb(Gby,0,4,ix,iy,iz,1);
      gg[1][4]=get_gb(Gby,1,4,ix,iy,iz,1);
      gg[2][4]=get_gb(Gby,2,4,ix,iy,iz,1);
      gg[3][4]=get_gb(Gby,3,4,ix,iy,iz,1);      
    }
  if(idim==2)
    {
      gg[0][0]=get_gb(Gbz,0,0,ix,iy,iz,2);
      gg[0][1]=get_gb(Gbz,0,1,ix,iy,iz,2);
      gg[0][2]=get_gb(Gbz,0,2,ix,iy,iz,2);
      gg[0][3]=get_gb(Gbz,0,3,ix,iy,iz,2);
      gg[1][0]=get_gb(Gbz,1,0,ix,iy,iz,2);
      gg[1][1]=get_gb(Gbz,1,1,ix,iy,iz,2);
      gg[1][2]=get_gb(Gbz,1,2,ix,iy,iz,2);
      gg[1][3]=get_gb(Gbz,1,3,ix,iy,iz,2);
      gg[2][0]=get_gb(Gbz,2,0,ix,iy,iz,2);  
      gg[2][1]=get_gb(Gbz,2,1,ix,iy,iz,2);  
      gg[2][2]=get_gb(Gbz,2,2,ix,iy,iz,2);  
      gg[2][3]=get_gb(Gbz,2,3,ix,iy,iz,2);  
      gg[3][0]=get_gb(Gbz,3,0,ix,iy,iz,2);
      gg[3][1]=get_gb(Gbz,3,1,ix,iy,iz,2);
      gg[3][2]=get_gb(Gbz,3,2,ix,iy,iz,2);
      gg[3][3]=get_gb(Gbz,3,3,ix,iy,iz,2);
      gg[0][4]=get_gb(Gbz,0,4,ix,iy,iz,2);
      gg[1][4]=get_gb(Gbz,1,4,ix,iy,iz,2);
      gg[2][4]=get_gb(Gbz,2,4,ix,iy,iz,2);
      gg[3][4]=get_gb(Gbz,3,4,ix,iy,iz,2);      
  }  

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//prints primitives
int
print_p(ldouble *p)
{
  printf("rho:   %10e\nuu:    %10e\nvr:    %10e\nvth:   %10e\nvph:   %10e\nS:     %10e\n",p[0],p[1],p[2],p[3],p[4],p[5]);
#ifdef RADIATION
  printf("E:     %10e\nFx:    %10e\nFy:    %10e\nFz:    %10e\n\n",p[6],p[7],p[8],p[9]);
#endif
  return 0;
}

//prints conserved
int
print_u(ldouble *u)
{
  printf("rhout: %10e\nTtt:   %10e\nTtr:   %10e\nTtth:  %10e\nTtph:  %10e\nSut:   %10e\n",u[0],u[1]-u[0],u[2],u[3],u[4],u[5]);
#ifdef RADIATION
  printf("Rtt:   %10e\nRt1:   %10e\nRt2:   %10e\nRt3:   %10e\n\n",u[6],u[7],u[8],u[9]);
#endif
  return 0;
}


