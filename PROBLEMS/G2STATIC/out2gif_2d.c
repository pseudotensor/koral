//int
//convert_out2gif_2d(char *fname,char *fname2,int niter,ldouble t)
//{
//  FILE *fgnu=fopen("plot.gp","w");
//  char bufor[50];

ldouble minx,miny,maxx,maxy;

#ifdef EQPLANEOUTPUT
minx= -1.1*(exp(get_xb(NX,0))+MKS1R0);
maxx= 1.1*(exp(get_xb(NX,0))+MKS1R0);
miny= -1.1*(exp(get_xb(NX,0))+MKS1R0);
maxy= 1.1*(exp(get_xb(NX,0))+MKS1R0);
#endif
#ifdef VERTPLANEOUTPUT
minx= -0.1*(exp(get_xb(NX,0))+MKS1R0);
maxx= 1.1*(exp(get_xb(NX,0))+MKS1R0);
miny= -1.1*(exp(get_xb(NX,0))+MKS1R0);
maxy= 1.1*(exp(get_xb(NX,0))+MKS1R0);
#endif

  fprintf(fgnu,
	  "set table \"table.gp\"\n"
	  "set contour base\n"
	  "unset surface\n"
	  "set log z\n"
	  "set cntrparam levels discrete %f,.1\n"
#ifdef EQPLANEOUTPUT
	  "splot \"%s\" u (($1)*sin($3)):(($1)*cos($3)):24 w l\n"
#endif
#ifdef VERTPLANEOUTPUT
	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):24 w l\n"
#endif
	  "unset dgrid3d\n"
	  "unset log z\n"
	  "unset table\n"
	  "unset contour\n"
	  "unset surface\n"

#ifdef EQPLANEOUTPUT
	  "set term gif large size 700,600\n"
#endif
#ifdef VERTPLANEOUTPUT
	  "set term gif large size 600,900\n"
#endif
	  "set output \"%s\"\n"
	  "set size 1,1\n"
	  "set origin 0,0\n"
	  "set multiplot\n"

	  "set view map\n"
	  "set pm3d\n"
	  "unset surface\n"
	  "unset key\n"
	  "set style line 1 lt 1 lw 2 lc 3\n"
	  "set style line 11 lt 1 lw 2 lc 6\n"
	  "set style line 2 lt 1 lw 2 lc 2\n"
	  "set style line 3 lt 1 lw 2 lc 2\n"
	  "set style line 21 lt 3 lw 1 lc -1\n"
	 
	  "set style arrow 1 ls 2\n"
	  
	  "set palette model RGB rgbformulae 7,5,15\n"
	
	  "set autoscale\n"
	  "set xrange [%f:%f]\n"
	  "set yrange [%f:%f]\n"

	  "set lmargin at screen 0.1\n"
	  "set rmargin at screen 0.81\n"
	  "set bmargin at screen .10\n"
	  "set tmargin at screen .95\n"
	  "set ylabel \"z\"\n"
	  "set xlabel \"x\" offset 0,1\n"
	  "set cblabel \"\"\n"
	  "set title \"rho\" offset 0,-1\n"
	  "set format cb \"%%.1e\"\n"

	  "set log cb\n"
	  "set autoscale cb\n"
	  "set cbrange [1.e-12:1e-10]\n"
	  
#ifdef EQPLANEOUTPUT
	  "splot \"%s\" u (($1)*sin($3)):(($1)*cos($3)):($14) ti \"\" w l ls 1\n"
#endif
#ifdef VERTPLANEOUTPUT
	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):($14) ti \"\" w l ls 1\n"
#endif

	  "set isosam 10,10\n"
	  "set ylabel \"\"\n"
	  "set xlabel \"\"\n"
	  "unset tics\n"
	  "unset title\n"
	  "unset border\n"
	  "unset log cb\n"

	  /*
	  "plot \"%s\" u (($1)*sin($3)):(($1)*cos($3)):"
	  "(($16*sin($3)+$18*cos($3)))*%f:"
	  "((-$18*sin($3)+$16*cos($3)))*%f every %d:%d w vectors arrowstyle 1 ti \"\"\n"
	  */

	  /*	    
	  "plot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):"
	  "(($16*sin($2)+$17*cos($2)))*%f:"
	  "((-$17*sin($2)+$16*cos($2)))*%f every %d:%d w vectors arrowstyle 1 ti \"\"\n"
	  */

#ifdef VERTPLANEOUTPUT
	  "plot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):"
	  "(($16*sin($2)+$17*cos($2))/(($16*sin($2)+$17*cos($2))**2+(-$17*sin($2)+$16*cos($2))**2)**.5)*%f:"
	  "((-$17*sin($2)+$16*cos($2))/(($16*sin($2)+$17*cos($2))**2+(-$17*sin($2)+$16*cos($2))**2)**.5)*%f every %d:%d w vectors arrowstyle 1 ti \"\"\n"
#endif
#ifdef EQPLANEOUTPUT
	  "plot \"%s\" u (($1)*sin($3)):(($1)*cos($3)):"
	  "(($16*sin($3)+$18*cos($3))/(($16*sin($3)+$18*cos($3))**2+(-$18*sin($3)+$16*cos($3))**2)**.5)*%f:"
	  "((-$18*sin($3)+$16*cos($3))/(($16*sin($3)+$18*cos($3))**2+(-$18*sin($3)+$16*cos($3))**2)**.5)*%f every %d:%d w vectors arrowstyle 1 ti \"\"\n"
#endif

	  "unset tics\n"
	  "unset border\n"
	  "unset pm3d\n"
	  "unset surface\n"
	  "plot \"table.gp\" w l ls 1\n"
	  "set pm3d\n"
	  "set tics\n"
	  "set border\n"

	 
 	  ,
	  MINTRACE,
	  fname,
	  fname2,
	  minx,
	  maxx,
	  miny,
	  maxy,
	  fname,
	  fname,

	  //maxx/21.*100.,maxx/21.*100.,
	  maxx/21/2,maxx/21/2,


#ifdef VERTPLANEOUTPUT
	  NX/21+1,NY/21+1
#endif
#ifdef EQPLANEOUTPUT
	  NX/21+1,NZ/21+1
#endif
	 
	  );  
//#endif	    

/*
  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  int i=system("gnuplot plot.gp ");
  return 0;
*/
