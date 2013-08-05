//v1 - 24
//..
//v7 - 30

  ldouble ucond[4],ucovd[4];
  ldouble bcond[4],bcovd[4],bsqd;

  //**********************************************************************
  //***** four velocity **************************************************
  //**********************************************************************

  for(iv=1;iv<4;iv++)
    ucond[iv]=pp[1+iv];
  ucond[0]=0.;
  conv_vels(ucond,ucond,VELPRIM,VEL4,gg,GG);
  indices_21(ucond,ucovd,gg);

#ifdef MAGNFIELD
  bcon_calc(pp,ucond,ucovd,bcond);
  indices_21(bcond,bcovd,gg); 
  bsqd = dot(bcond,bcovd);
  v1=bsqd/2.;
#endif
