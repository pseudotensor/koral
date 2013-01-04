//KORAL - relele.c
//some relativistic routines

#include "ko.h"

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates conserved in given cell using global array p[]
int
calc_conserved(int ix,int iy,int iz)
{
  int iv;
  ldouble uu[NV],pp[NV];
  ldouble gg[4][5],tlo[4][4],tup[4][4];
  
  pick_g(ix,iy,iz,gg);
  pick_T(tmuup,ix,iy,iz,tup);
  pick_T(tmulo,ix,iy,iz,tlo);

  for(iv=0;iv<NV;iv++)
    {
      pp[iv]=get_u(p,iv,ix,iy,iz);
    }

  p2u(pp,uu,gg,tup,tlo);

  for(iv=0;iv<NV;iv++)
    {
      set_u(u,iv,ix,iy,iz,uu[iv]);
    }

  return 0;
}
 

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates primitives in given cell basing on global array u[]
int
calc_primitives(int ix,int iy,int iz)
{
  int iv,u2pret,u2pretav;
  ldouble uu[NV],uuav[NV],pp[NV],ppav[NV];
  ldouble gg[4][5], GG[4][5],tlo[4][4],tup[4][4];

  pick_g(ix,iy,iz,gg);
  pick_G(ix,iy,iz,GG);
  pick_T(tmuup,ix,iy,iz,tup);
  pick_T(tmulo,ix,iy,iz,tlo);

  for(iv=0;iv<NV;iv++)
    {
      uu[iv]=get_u(u,iv,ix,iy,iz);
      pp[iv]=get_u(p,iv,ix,iy,iz);
    }

  //converting to primitives
  u2pret=u2p(uu,pp,gg,GG,tup,tlo);

  //sets the flag to mark if hot conversion did not succeed - the entropy will not be updated
  set_cflag(0,ix,iy,iz,u2pret); 

  ldouble vr=pp[2];
  ldouble vth=pp[3];
  ldouble vph=pp[4];

  ldouble ut2=-1./(gg[0][0] + 2.*vph*gg[0][3] + vr*vr*gg[1][1] + vth*vth*gg[2][2] + vph*vph*gg[3][3] );
  
  //checking on floors & ceilings
  int updateu=0;
 
  if(ut2<0.)
    {
      my_err("ut2.lt.0 in calc_primitives!\n"); ut2=0.;
    }
  
  for(iv=0;iv<NV;iv++)    
    set_u(p,iv,ix,iy,iz,pp[iv]);	      

  return 0;
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
//precalculates metric etc. and saves it to arrays
int
calc_metric()
{
  ldouble xx[4];
  int ix,iy,iz,i,j,k;
  ldouble gloc[4][5];
  ldouble Kr[4][4][4];
  ldouble eup[4][4],elo[4][4];
  ldouble tup[4][4],tlo[4][4];

  printf("Precalculating metrics... ");
  
  for(ix=-NG;ix<NX+NG;ix++)
    {
      for(iy=-NG;iy<NY+NG;iy++)
	{
	  for(iz=-NG;iz<NZ+NG;iz++)
	    {
	      //cell centers
	      xx[0]=0.;
	      xx[1]=get_x(ix,0);
	      xx[2]=get_x(iy,1);
	      xx[3]=get_x(iz,2);
	      calc_g(xx,gloc);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		    set_g(g,i,j,ix,iy,iz,gloc[i][j]);
	      for(j=0;j<3;j++)
		set_g(g,j,4,ix,iy,iz,calc_dlgdet(xx,j));
	      set_g(g,3,4,ix,iy,iz,calc_gdet(xx));

	      calc_LNRFes(gloc,eup,elo);
	      calc_tetrades(gloc,tup,tlo);
	      
	      /*
	      printf("Bardeen eup elo\n");
	      print_tensor(eup);
	      print_tensor(elo);

	      printf("Tetrad tup tlo\n");
	      print_tensor(tup);
	      print_tensor(tlo);

	      getchar();
	      */

	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  {
		    set_T(emuup,i,j,ix,iy,iz,eup[i][j]);
		    set_T(emulo,i,j,ix,iy,iz,elo[i][j]);
		    set_T(tmuup,i,j,ix,iy,iz,tup[i][j]);
		    set_T(tmulo,i,j,ix,iy,iz,tlo[i][j]);
		  }	      

	      calc_G(xx,gloc);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		    set_g(G,i,j,ix,iy,iz,gloc[i][j]);
	      
	      calc_Krzysie(xx,Kr);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  for(k=0;k<4;k++)
		    set_gKr(i,j,k,ix,iy,iz,Kr[i][j][k]);


	      	      
	      //x-faces
	      if(ix==-NG)
		{
		  xx[0]=0.;
		  xx[1]=get_xb(ix,0);
		  xx[2]=get_x(iy,1);
		  xx[3]=get_x(iz,2);
		  calc_g(xx,gloc);
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      set_gb(gbx,i,j,ix,iy,iz,gloc[i][j],0);

		  calc_LNRFes(gloc,eup,elo);
		  calc_tetrades(gloc,tup,tlo);

		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      {
			set_Tb(emuupbx,i,j,ix,iy,iz,eup[i][j],0);
			set_Tb(emulobx,i,j,ix,iy,iz,elo[i][j],0);
			set_Tb(tmuupbx,i,j,ix,iy,iz,tup[i][j],0);
			set_Tb(tmulobx,i,j,ix,iy,iz,tlo[i][j],0);
		      }	      


		  calc_G(xx,gloc);
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      set_gb(Gbx,i,j,ix,iy,iz,gloc[i][j],0);
		  for(j=0;j<3;j++)
		    set_gb(gbx,j,4,ix,iy,iz,calc_dlgdet(xx,j),0);
		  set_gb(gbx,3,4,ix,iy,iz,calc_gdet(xx),0);

		  calc_Krzysie(xx,Kr);
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      for(k=0;k<4;k++)
			set_gKrb(i,j,k,ix,iy,iz,Kr[i][j][k],0);


		}
	      xx[0]=0.;
	      xx[1]=get_xb(ix+1,0);
	      xx[2]=get_x(iy,1);
	      xx[3]=get_x(iz,2);
	      calc_g(xx,gloc);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  set_gb(gbx,i,j,ix+1,iy,iz,gloc[i][j],0);

	      calc_LNRFes(gloc,eup,elo);
	      calc_tetrades(gloc,tup,tlo);

	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  {
		    set_Tb(emuupbx,i,j,ix+1,iy,iz,eup[i][j],0);
		    set_Tb(emulobx,i,j,ix+1,iy,iz,elo[i][j],0);
		    set_Tb(tmuupbx,i,j,ix+1,iy,iz,tup[i][j],0);
		    set_Tb(tmulobx,i,j,ix+1,iy,iz,tlo[i][j],0);
		  }	      

	      calc_G(xx,gloc);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  set_gb(Gbx,i,j,ix+1,iy,iz,gloc[i][j],0);
	      for(j=0;j<3;j++)
		set_gb(gbx,j,4,ix+1,iy,iz,calc_dlgdet(xx,j),0);
	      set_gb(gbx,3,4,ix+1,iy,iz,calc_gdet(xx),0);

	      calc_Krzysie(xx,Kr);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  for(k=0;k<4;k++)
		    set_gKrb(i,j,k,ix+1,iy,iz,Kr[i][j][k],0);

		  
	      //y-faces
	      if(iy==-NG)
		{
		  xx[0]=0.;
		  xx[1]=get_x(ix,0);
		  xx[2]=get_xb(iy,1);
		  xx[3]=get_x(iz,2);
		  calc_g(xx,gloc);
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      set_gb(gby,i,j,ix,iy,iz,gloc[i][j],1);

		  calc_LNRFes(gloc,eup,elo);
		  calc_tetrades(gloc,tup,tlo);

		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      {
			set_Tb(emuupby,i,j,ix,iy,iz,eup[i][j],1);
			set_Tb(emuloby,i,j,ix,iy,iz,elo[i][j],1);
			set_Tb(tmuupby,i,j,ix,iy,iz,tup[i][j],1);
			set_Tb(tmuloby,i,j,ix,iy,iz,tlo[i][j],1);
		      }	      

		  calc_G(xx,gloc);
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      set_gb(Gby,i,j,ix,iy,iz,gloc[i][j],1);
		  for(j=0;j<3;j++)
		    set_gb(gby,j,4,ix,iy,iz,calc_dlgdet(xx,j),1);
		  set_gb(gby,3,4,ix,iy,iz,calc_gdet(xx),1);

		  calc_Krzysie(xx,Kr);
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      for(k=0;k<4;k++)
			set_gKrb(i,j,k,ix,iy,iz,Kr[i][j][k],1);

		}
	      xx[0]=0.;
	      xx[1]=get_x(ix,0);
	      xx[2]=get_xb(iy+1,1);
	      xx[3]=get_x(iz,2);
	      calc_g(xx,gloc);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  set_gb(gby,i,j,ix,iy+1,iz,gloc[i][j],1);

	      calc_LNRFes(gloc,eup,elo);
	      calc_tetrades(gloc,tup,tlo);

	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  {
		    set_Tb(tmuupby,i,j,ix,iy+1,iz,tup[i][j],1);
		    set_Tb(tmuloby,i,j,ix,iy+1,iz,tlo[i][j],1);
		  }	      

	      calc_G(xx,gloc);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  set_gb(Gby,i,j,ix,iy+1,iz,gloc[i][j],1);
	      for(j=0;j<3;j++)
		set_gb(gby,j,4,ix,iy+1,iz,calc_dlgdet(xx,j),1);
	      set_gb(gby,3,4,ix,iy+1,iz,calc_gdet(xx),1);
		  
	      calc_Krzysie(xx,Kr);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  for(k=0;k<4;k++)
		    set_gKrb(i,j,k,ix,iy+1,iz,Kr[i][j][k],1);

	      //z-faces
	      if(iz==-NG)
		{
		  xx[0]=0.;
		  xx[1]=get_x(ix,0);
		  xx[2]=get_x(iy,1);
		  xx[3]=get_xb(iz,2);
		  calc_g(xx,gloc);
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      set_gb(gbz,i,j,ix,iy,iz,gloc[i][j],2);

		  calc_LNRFes(gloc,eup,elo);
		  calc_tetrades(gloc,tup,tlo);

		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      {
			set_Tb(emuupbz,i,j,ix,iy,iz,eup[i][j],2);
			set_Tb(emulobz,i,j,ix,iy,iz,elo[i][j],2);
			set_Tb(tmuupbz,i,j,ix,iy,iz,tup[i][j],2);
			set_Tb(tmulobz,i,j,ix,iy,iz,tlo[i][j],2);
		      }	      

		  calc_G(xx,gloc);
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      set_gb(Gbz,i,j,ix,iy,iz,gloc[i][j],2);
		  for(j=0;j<3;j++)
		    set_gb(gbz,j,4,ix,iy,iz,calc_dlgdet(xx,j),2);
		  set_gb(gbz,3,4,ix,iy,iz,calc_gdet(xx),2);

		  calc_Krzysie(xx,Kr);
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      for(k=0;k<4;k++)
			set_gKrb(i,j,k,ix,iy,iz,Kr[i][j][k],2);

		}
	      xx[0]=0.;
	      xx[1]=get_x(ix,0);
	      xx[2]=get_x(iy,1);
	      xx[3]=get_xb(iz+1,2);
	      calc_g(xx,gloc);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  set_gb(gbz,i,j,ix,iy,iz+1,gloc[i][j],2);	  

	      calc_LNRFes(gloc,eup,elo);
	      calc_tetrades(gloc,tup,tlo);

	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  {
		    set_Tb(emuupbz,i,j,ix,iy,iz+1,eup[i][j],2);
		    set_Tb(emulobz,i,j,ix,iy,iz+1,elo[i][j],2);
		    set_Tb(tmuupbz,i,j,ix,iy,iz+1,tup[i][j],2);
		    set_Tb(tmulobz,i,j,ix,iy,iz+1,tlo[i][j],2);
		  }	      

	      calc_G(xx,gloc);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  set_gb(Gbz,i,j,ix,iy,iz+1,gloc[i][j],2);	  
	      for(j=0;j<3;j++)
		set_gb(gbz,j,4,ix,iy,iz+1,calc_dlgdet(xx,j),2);
	      set_gb(gbz,3,4,ix,iy,iz+1,calc_gdet(xx),2);

	      calc_Krzysie(xx,Kr);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  for(k=0;k<4;k++)
		    set_gKrb(i,j,k,ix,iy,iz+1,Kr[i][j][k],2);

	      
	    }
	}
    }

  printf("done!\n");

  return 0;
}

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
//calculates energy-momentum tensor components basing on vector of primitivies p and given metric g
int
calc_Tmunu( ldouble *p, ldouble gg[][5], ldouble T[][4], ldouble *ut_ret)
{
  ldouble rho=p[0];
  ldouble uu=p[1];
  ldouble vr=p[2];
  ldouble vth=p[3];
  ldouble vph=p[4];

  ldouble gtt=gg[0][0];
  ldouble gtph=gg[0][3];
  ldouble grr=gg[1][1];
  ldouble gthth=gg[2][2];
  ldouble gphph=gg[3][3];


  ldouble ut2=-1./(gtt + 2.*vph*gtph + vr*vr*grr + vth*vth*gthth + vph*vph*gphph );
  if(ut2<0.)
    {
      my_err("ut2.lt.0 in calc_Tmunu\n"); ut2=0.;
     }

  ldouble ut=sqrtl(ut2);

  //passes it up to avoid doubling sqrtl
  *ut_ret=ut;
  
  ldouble ur=vr*ut;
  ldouble uth=vth*ut;
  ldouble uph=vph*ut;

  ldouble w=rho+GAMMA*uu;

  ldouble Ttt=w*ut*(ut*gtt+uph*gtph)+(GAMMA-1.)*uu;
  ldouble Ttr=w*ut*ur*grr;
  ldouble Ttth=w*ut*uth*gthth;
  ldouble Ttph=w*ut*(uph*gphph+ut*gtph);

  ldouble Trt=w*ur*(ut*gtt+uph*gtph);
  ldouble Trr=w*ur*ur*grr+(GAMMA-1.)*uu;
  ldouble Trth=w*ur*uth*gthth;
  ldouble Trph=w*ur*(uph*gphph+ut*gtph);

  ldouble Ttht=w*uth*(ut*gtt+uph*gtph);
  ldouble Tthr=w*uth*ur*grr;
  ldouble Tthth=w*uth*uth*gthth+(GAMMA-1.)*uu;
  ldouble Tthph=w*uth*(uph*gphph+ut*gtph);

  ldouble Tpht=w*uph*(ut*gtt+uph*gtph);
  ldouble Tphr=w*uph*ur*grr;
  ldouble Tphth=w*uph*uth*gthth;
  ldouble Tphph=w*uph*(uph*gphph+ut*gtph)+(GAMMA-1.)*uu;

  T[0][0]=Ttt;
  T[0][1]=Ttr;
  T[0][2]=Ttth;
  T[0][3]=Ttph;

  T[1][0]=Trt;
  T[1][1]=Trr;
  T[1][2]=Trth;
  T[1][3]=Trph;

  T[2][0]=Ttht;
  T[2][1]=Tthr;
  T[2][2]=Tthth;
  T[2][3]=Tthph;

  T[3][0]=Tpht;
  T[3][1]=Tphr;
  T[3][2]=Tphth;
  T[3][3]=Tphph;
  
  return 0;
}

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
//updates entropy (p[5]) basing on new primitives or stays with the old one if entropy u2p solver was involved
int
update_entropy(int ix,int iy,int iz,int u2pflag)
{
  ldouble gg[4][5];
  pick_g(ix,iy,iz,gg);

  ldouble gtt=gg[0][0];
  ldouble gtph=gg[0][3];
  ldouble grr=gg[1][1];
  ldouble gthth=gg[2][2];
  ldouble gphph=gg[3][3];

  ldouble rho=get_u(p,0,ix,iy,iz);
  ldouble uu=get_u(p,1,ix,iy,iz);
  ldouble vr=get_u(p,2,ix,iy,iz);
  ldouble vth=get_u(p,3,ix,iy,iz);
  ldouble vph=get_u(p,4,ix,iy,iz);

  ldouble ut2=-1./(gtt + 2.*vph*gtph + vr*vr*grr + vph*vph*gphph + vth*vth*gthth);
  ldouble ut=sqrtl(ut2);

  ldouble S,Sut;

  //u2p_hot worked
  if(u2pflag==0 && uu>0. && rho>0.)
    {
      S=calc_Sfromu(rho,uu);      
      set_u(p,5,ix,iy,iz,S);
      set_u(u,5,ix,iy,iz,S*ut); 
    }
  //u2p_entropy worked
  else if(u2pflag==-1  && uu>0. && rho>0.)
    {
      Sut=get_u(u,5,ix,iy,iz);
      S=Sut/ut;
      set_u(p,5,ix,iy,iz,S);
    }
  //somnething else - leave entropy as it was
  else
    {
      //nothing
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
  printf("rho:   %10Le\nuu:    %10Le\nvr:    %10Le\nvth:   %10Le\nvph:   %10Le\nS:     %10Le\n",p[0],p[1],p[2],p[3],p[4],p[5]);
#ifdef RADIATION
  printf("E:     %10Le\nFx:    %10Le\nFy:    %10Le\nFz:    %10Le\n\n",p[6],p[7],p[8],p[9]);
#endif
  return 0;
}

//prints conserved
int
print_u(ldouble *u)
{
  printf("rhout: %10Le\nTtt:   %10Le\nTtr:   %10Le\nTtth:  %10Le\nTtph:  %10Le\nSut:   %10Le\n",u[0],u[1]-u[0],u[2],u[3],u[4],u[5]);
#ifdef RADIATION
  printf("Rtt:   %10Le\nRt1:   %10Le\nRt2:   %10Le\nRt3:   %10Le\n\n",u[6],u[7],u[8],u[9]);
#endif
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
