//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp) {

ldouble gdet_src,gdet_bc;
int iix,iiy,iiz,iv;  	  

gdet_bc=get_g(g,3,4,ix,iy,iz);  
ldouble gg[4][5],GG[4][5],tlo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
pick_T(tmulo,ix,iy,iz,tlo);
ldouble xx=get_x(ix,0);
ldouble yy=get_x(iy,1);

//printf("aa\n");

//source of light
//x<0 && yy>.5)
if (iy>=NY || (ix<0 && yy>.3 ))
  {
    ldouble Fx,Fy,Fz,rho,E,uint,vx;
    iix=ix;
    iiy=iy;
    iiz=0;
    rho=1.;

    Fx=0.;
    Fy=0.;
    Fz=0.;	      	      

    
    uint=1.;

    
    vx=0.;

    //directly
    E=pow(100,4.);
    Fx=URFX;
    Fy=URFY;
      
    pp[0]=rho;
    pp[1]=uint;
    pp[2]=vx;
    pp[3]=0.;
    pp[4]=0.;
    pp[5]=calc_Sfromu(rho,uint);
    pp[6]=E;
    pp[7]=Fx;
    pp[8]=Fy;
    pp[9]=Fz;

    //prad_ff2lab(pp,pp,gg,GG,tlo);

    p2u(pp,uu,gg,GG);
    return 0.;
  }

if (ix<0.)
  {
    ldouble Fx,Fy,Fz,rho,E,uint,vx;
    iix=ix;
    iiy=iy;
    iiz=0;
    rho=1.;

    Fx=0.;
    Fy=0.;
    Fz=0.;	      	      

    
    uint=1.;

    
    vx=0.;

    //directly
    E=1.;
    Fx=0.;
    Fy=0.;
      
    pp[0]=rho;
    pp[1]=uint;
    pp[2]=vx;
    pp[3]=0.;
    pp[4]=0.;
    pp[5]=calc_Sfromu(rho,uint);
    pp[6]=E;
    pp[7]=Fx;
    pp[8]=Fy;
    pp[9]=Fz;

    //prad_ff2lab(pp,pp,gg,GG,tlo);

    p2u(pp,uu,gg,GG);
    return 0.;
  }

if(ix>=NX && 0) //x wall
  {
    iix=NX-1-(ix-NX);
    iiy=iy;
    iiz=iz;

    for(iv=0;iv<NV;iv++)
      {      
	pp[iv]=get_u(p,iv,iix,iiy,iiz);      
      }

    pp[7]*=-1.;

    p2u(pp,uu,gg,GG);

    return 0;
  }

 if(iy<0 && 1) //y wall
  {
    iix=ix;
    iiy=-iy-1;
    iiz=iz;

    for(iv=0;iv<NV;iv++)
      {      
	pp[iv]=get_u(p,iv,iix,iiy,iiz);      
      }

    pp[8]*=-1.;

    p2u(pp,uu,gg,GG);

    return 0;
  }
 

iix=ix;
iiz=iz;
iiy=iy;

//outflow
while(iix<0)    iix=0 ;
while(iix>=NX)    iix=NX-1;
//outflow
while(iiy<0)    iiy=0 ;
while(iiy>=NY)    iiy=NY-1;
//periodic
while(iiz<0)    iiz+=NZ;
while(iiz>=NZ)    iiz-=NZ; 
 
for(iv=0;iv<NV;iv++)
  {      
    pp[iv]=get_u(p,iv,iix,iiy,iiz);      
  }
p2u(pp,uu,gg,GG);
return 0;
  
