
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

ldouble rho,mx,my,mz,m,E,uint,E0,Fx,Fy,Fz,pLTE;  
ldouble xx,yy,zz;
ldouble uu[NV],xxvec[4];

get_xx(ix,iy,iz,xxvec);
coco_N(xxvec,xxvec,MYCOORDS,BLCOORDS);
xx=xxvec[1];
yy=xxvec[2];
zz=xxvec[3];


ldouble gg[4][5],GG[4][5],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
calc_ZAMOes(gg,eup,elo,MYCOORDS);

ldouble pp[NV],T;

//working in BL
ldouble ggBL[4][5],GGBL[4][5];
calc_g_arb(xxvec,ggBL,KERRCOORDS);
calc_G_arb(xxvec,GGBL,KERRCOORDS);
ldouble eupBL[4][4],eloBL[4][4];
ldouble tupBL[4][4],tloBL[4][4];
calc_tetrades(ggBL,tupBL,tloBL,KERRCOORDS);
calc_ZAMOes(ggBL,eupBL,eloBL,KERRCOORDS);

	 
ldouble podpierd=-(GGBL[0][0]-2.*ELL*GGBL[0][3]+ELL*ELL*GGBL[3][3]);
ldouble ut=-1./sqrt(podpierd);



ut/=UTPOT; //rescales rin
//printf("%Lf %Lf ut: %Lf\n",xx,yy,podpierd);getchar();
ldouble Vphi,Vr;
ldouble D,W,eps,uT,uphi,uPhi;
if(ut<-1 || podpierd<0. || xx<3. || NODONUT || INFLOWING)
  {
    ldouble r=xx;
    //		  printf("%Lf\n",xx); getchar();
    rho=RHO_AMB*powl(xx/2.,-1.5);


    uint=U_AMB*powl(xx/2.,-5./2.);
    Vphi=0.;
    rho=PAR_D/(r*r*sqrtl(2./r));  
    uPhi=GGBL[3][3]*uphi+GGBL[0][3]*ut;
    //		  Vphi=0.*uPhi/uT;

    //zaczynam jednak od profilu analitycznego:   
		  
    ldouble mD=PAR_D/(r*r*sqrtl(2./r*(1.-2./r)));
    ldouble mE=PAR_E/(powl(r*r*sqrtl(2./r),GAMMA)*powl(1.-2./r,(GAMMA+1.)/4.));
    ldouble V=sqrtl(2./r)*(1.-2./r);
    W=1./sqrtl(1.-V*V*ggBL[1][1]);
    ldouble mrho=mD/W;
    ldouble muint=mE/W;	      
    //corrected rho:
    rho=PAR_D/(r*r*sqrtl(2./r));    
    Vr=V;
    uint=muint;

  }
 else
   {
     /*
       ldouble phi=zz;
       uphi=-ELL*ut;
       uT=GG[0][0]*ut+GG[0][3]*uphi;
       eps=1./GAMMA*(-1./ut-1.);
       W=uT/sqrt(-gg[0][0]);
       D=powl(eps*(GAMMA-1.)/KKK,1./(GAMMA-1.))*W;
       E=eps*D*W;
       rho=D/W;
       uint=E/W;
       uPhi=GG[3][3]*uphi+GG[0][3]*ut;
       Vphi=uPhi/uT;
       #ifdef DONUT3D
       ldouble ran;
       ran=random()/(ldouble)RAND_MAX;
       //		  printf("%Le\n",ran); getchar();
       rho*=1.+(PERTP*(ran-.5));
       uint*=1.+PERTP*(ran-.5);		  
       #endif
     */
     ldouble h=-1./ut;
     ldouble eps=(h-1.)/GAMMA;
     rho=powl(eps*(GAMMA-1.)/KKK,1./(GAMMA-1.));
     uint=rho*eps;
     //		  uint=KKK*powl(rho,GAMMA)/(GAMMA-1.);
     uphi=-ELL*ut;
     uT=GGBL[0][0]*ut+GGBL[0][3]*uphi;
     uPhi=GGBL[3][3]*uphi+GGBL[0][3]*ut;
     Vphi=uPhi/uT;
     Vr=0.;
   }       
//	      if(rho<RHO_AMB) rho=RHO_AMB;
//	      if(uint<U_AMB) uint=U_AMB;

//4-velocity in BL
ldouble ucon[4]={0.,-Vr,0.,Vphi};


conv_vels(ucon,ucon,VEL3,VEL4,ggBL,GGBL);

//converting to KS
trans2_coco(xxvec,ucon,ucon,BLCOORDS,MYCOORDS);

//to be written in primitives
conv_vels(ucon,ucon,VEL4,VELPRIM,gg,GG);


pp[2]=ucon[1]; 
pp[3]=ucon[2];
pp[4]=ucon[3];

//density etc.
pp[0]=rho; pp[1]=uint; 
pp[5]=calc_Sfromu(pp[0],pp[1]);

//testing if interpolated primitives make sense
check_floors_hd(pp,VELPRIM,gg,GG);
//end of floor section

p2u(pp,uu,gg,GG);

/*
  print_Nvector(pp,NV);
  print_Nvector(uu,NV);
  getchar();
*/

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
