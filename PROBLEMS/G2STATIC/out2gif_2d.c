//int
//convert_out2gif_2d(char *fname,char *fname2,int niter,ldouble t)
//{
//  FILE *fgnu=fopen("plot.gp","w");
//  char bufor[50];

ldouble minx,miny,maxx,maxy;

#if (OUTCOORDS==KERRCOORDS)
minx= -.02*(exp(get_xb(-NG,0))+MKS1R0);
maxx= 1.01*(exp(get_xb(NX,0))+MKS1R0);
miny= -1.01*(exp(get_xb(NX,0))+MKS1R0);
maxy= 1.01*(exp(get_xb(NX,0))+MKS1R0);
#else

minx= -.02*get_xb(NX,0);
maxx= 1.02*get_xb(NX,0);
miny= -1.1*get_xb(NX,0);
maxy= 1.1*get_xb(NX,0);
#endif

  fprintf(fgnu,
	  "set term gif large size 600,900\n"
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
	  //"set log cb\n"
	  "set ylabel \"z\"\n"
	  "set xlabel \"x\" offset 0,1\n"
	  "set cblabel \"\"\n"
	  "set title \"rho\" offset 0,-1\n"
	  "set format cb \"%%.1e\"\n"

	  "set autoscale cb\n"
	  //"set cbrange [0:3e-20]\n"
	  
	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):($14) ti \"\" w l ls 1\n"

	  "set isosam 10,10\n"
	  "set ylabel \"\"\n"
	  "set xlabel \"\"\n"
	  "unset tics\n"
	  "unset title\n"
	  "unset border\n"
	  "unset log cb\n"

	  //	  "plot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):"
	  //	  "(($16*sin($2)+$17*cos($2))*%f):"
	  //	  "((-$17*sin($2)+$16*cos($2))*%f) every %d:%d w vectors arrowstyle 1 ti \"\"\n"
	  "plot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):"
	  "(($16*sin($2)+$17*cos($2))/(($16*sin($2)+$17*cos($2))**2+(-$17*sin($2)+$16*cos($2))**2)**.5*%f):"
	  "((-$17*sin($2)+$16*cos($2))/(($16*sin($2)+$17*cos($2))**2+(-$17*sin($2)+$16*cos($2))**2)**.5*%f) every %d:%d w vectors arrowstyle 1 ti \"\"\n"

	 
	 
 	  ,fname2,
	  minx,
	  maxx,
	  miny,
	  maxy,
	  fname,
	  fname,.95,.95,NX/21+1,NY/21+1
	 
	  );  
//#endif	    

/*
  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  int i=system("gnuplot plot.gp ");
  return 0;
*/
