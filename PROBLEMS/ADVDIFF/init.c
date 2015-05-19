

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

ldouble pp[NV],uu[NV];

/************************/

ldouble Tgas=TAMB*(1.+PULSEMAG*exp(-(geom.xx-0.)*(geom.xx-0.)/5./5.));


//if(ix==NX/2 ) {printf("%d %Le %Le %Le %Le %Le %e\n",ix,xx,yy,zz,Tgas,(ldouble)4.*SIGMA_RAD*Tgas*Tgas*Tgas*Tgas,SIGMA_RAD); getchar();}

pp[0]=RHOAMB;//*(1.+PULSEMAG*exp(-(geom.xx-0.5)*(geom.xx-0.5)/0.01));
pp[1]=calc_PEQ_ufromTrho(TAMB,RHOAMB);
pp[2]=VELX;
pp[3]=0.;
pp[4]=0.;
pp[5]=calc_Sfromu(pp[0],pp[1]);
#ifdef RADIATION
pp[6]=calc_LTE_EfromT(Tgas);//UUAMB*(1.+PULSEMAG*exp(-(geom.xx-0.)*(geom.xx-0.)/5./5.));
pp[7]=VELX;
pp[8]=0.;
pp[9]=0.;
#endif	 

p2u(pp,uu,&geom);	 


/***********************************************/

int iv;

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }

