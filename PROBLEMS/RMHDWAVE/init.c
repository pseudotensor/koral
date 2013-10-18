

ldouble rho,mx,my,mz,m,E,uint,E0,Fx,Fy,Fz,pLTE,vx,vy,Bx,By;
ldouble xx,yy,zz;
ldouble uu[NV];

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

xx=get_x(ix,0);
yy=get_x(iy,1);
zz=get_x(iz,2);
ldouble gg[4][5],GG[4][5],tup[4][4],tlo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);

ldouble pp[NV],T;

/************************/

ldouble t=0.;


rho=RHOZERO*(1+DRRE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx)-DRIM/DRRE*sin(OMRE*t-KK*xx)));
uint=UZERO*(1.+DURE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx)-DUIM/DURE*sin(OMRE*t-KK*xx))) ;
vx=0.+DV1RE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx)-DV1IM/DV1RE*sin(OMRE*t-KK*xx)) ; 
vy=0.+DV2RE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx)-DV2IM/DV2RE*sin(OMRE*t-KK*xx)) ; 
Bx=B1ZERO;
By=B2ZERO+DB2RE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx)-DB2IM/DB2RE*sin(OMRE*t-KK*xx)) ; 
#ifdef RADIATION
E=EEZERO*(1+DEERE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx)-DEEIM/DEERE*sin(OMRE*t-KK*xx)));
Fx=0.+ERAD*DF1RE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx)-DF1IM/DF1RE*sin(OMRE*t-KK*xx));
Fy=0.+ERAD*DF2RE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx)-DF2IM/DF2RE*sin(OMRE*t-KK*xx));
#endif

pp[0]=rho;
pp[1]=uint;
pp[2]=vx;
pp[3]=vy;
pp[4]=0.;
pp[5]=calc_Sfromu(rho,uint);
pp[B1]=Bx;
pp[B2]=By;
#ifdef RADIATION
pp[6]=E;
pp[7]=Fx;
pp[8]=Fy;
pp[9]=0.; 
prad_ff2lab(pp,pp,&geom);
#endif

//print_Nvector(pp,NV); getchar();
p2u(pp,uu,&geom);	 


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
