
/*
int
set_initial_profile()
{
  int ix,iy,iz;
  for(iz=0;iz<NZ;iz++)
    {
      for(iy=0;iy<NY;iy++)
	{
	  for(ix=0;ix<NX;ix++)
	    {
*/

ldouble rho,mx,my,mz,m,E,uint,E0,Fx,Fy,Fz,pLTE,vx;  
ldouble xx,yy,zz;
ldouble uu[NV];

xx=get_x(ix,0);
yy=get_x(iy,1);
zz=get_x(iz,2);
ldouble gg[4][5],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
calc_LNRFes(gg,eup,elo);

ldouble pp[NV],T;

/************************/

ldouble t=0.;
rho=RHO;
uint=UINT;
vx=0.;
E=ERAD*4.*SIGMA_RAD*TEMP*TEMP*TEMP*TEMP;
Fx=0.;


      Fz=Fy=0.;
      pp[0]=rho;
      pp[1]=uint;
pp[2]=vx;
      pp[3]=0.;
      pp[4]=0.;
      pp[5]=calc_Sfromu(rho,uint);
#ifdef RADIATION
      pp[6]=E;
      pp[7]=Fx;
      pp[8]=Fy;
      pp[9]=Fz; 
#endif

//print_Nvector(pp,NV); getchar();
p2u(pp,uu,gg,eup,elo);	 


/***********************************************/

	      int iv;

	      for(iv=0;iv<NV;iv++)
		{
		  set_u(u,iv,ix,iy,iz,uu[iv]);
		  set_u(p,iv,ix,iy,iz,pp[iv]);
		}

	      //entropy
	      update_entropy(ix,iy,iz,0);

	      //if(isnan(get_u(p,5,ix,iy,iz))) {printf("pr: %d %d %d S: %Le\n",ix,iy,iz,0.);getchar();}

	      //mark initialy succesfull u2p_hot step
	      set_cflag(0,ix,iy,iz,0);
