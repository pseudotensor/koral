//user output - called from fileop.c

//v1 - 24
//..
//v7 - 30

v1=calc_PEQ_Tfromurho(pp[UU],pp[RHO]);
v2=calc_LTE_TfromE(pp[EE0]);

ldouble urad00[4],Rtt;
calc_ff_Rtt(pp,&Rtt,urad00,&geom);
v3=-Rtt;
