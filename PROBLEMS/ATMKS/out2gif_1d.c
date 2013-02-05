//int!
//convert_out2gif_1d(char *fname,char *fname2,int niter,ldouble t)
//{
//  FILE *fgnu=fopen("plot.gp","w");
//  char bufor[50];

/****************************************/
/****************************************/
/****************************************/
   fprintf(fgnu,
	  "set style line 1 lw 2 lc 1 lt 1\n"
	  "set style line 2 lw 2 lc 2 lt 1\n"
	  "set style line 3 lw 2 lc 3 lt 1\n"
	  "set style line 4 lw 2 lc 4 lt 1\n"
	  "set style line 10 lw 2 lc 3 lt 3\n"
	  "set style line 11 lw 2 lc 4 lt 3\n"
	  "set style line 12 lw 2 lc 5 lt 3\n"
	  "set style line 13 lw 2 lc 0 lt 3\n"
	  "set style line 14 lw 2 lc 7 lt 3\n"
	  "set term gif large size 800,600\n"
	  "set output \"%s\"\n"
	  "set size 1.,1.\n"
	  "set origin 0,0\n"
	  "set multiplot\n"
	  "set autoscale\n" 
	   //	   "set log \n"
	  "set label \"t=%.2e (%.2e s)\" at screen .48, .98\n"

	  "set autoscale\n"
	  "set xrange [%f:%f]\n"

	  "set lmargin at screen 0.1\n"
	  "set rmargin at screen 0.45\n"
	  "set bmargin at screen .5\n"
	  "set tmargin at screen .9\n"
	  //	  "unset log y\n"
	  "set format x \"\"\n"
	  "set format y \"%%.1e\"\n" 
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"
	  //	  "set log y\n"
	  "plot \"%s\" u 1:14 w lp ls 3 pt 7 ps .5 ti \"rho\"\n"
	  //	  "unset log y\n"


	  "set lmargin at screen 0.1\n"
	  "set rmargin at screen 0.45\n"
	  "set bmargin at screen .07\n"
	  "set tmargin at screen .42\n"
	  //	  "unset log y\n"
	  "set format x \"%%.1f\"\n"
	  "set format y \"%%.1e\"\n" 
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"
	  "plot \"%s\" u 1:($15) w lp ls 1 pt 7 ps .5 ti \"u\"\n"

	  "set lmargin at screen 0.60\n"
	  "set rmargin at screen 0.95\n"
	  "set bmargin at screen .07\n"
	  "set tmargin at screen .42\n"
	  //	  "unset log y\n"
	  "set format x \"%%.1f\"\n"
	  "set format y \"%%.1e\"\n" 
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"
	  "plot \"%s\" u 1:($16) w lp ls 2 pt 7 ps .5  ti \"u^r/u^t\"\n"


	  "set lmargin at screen 0.60\n"
	  "set rmargin at screen 0.95\n"
	  "set bmargin at screen .5\n"
	  "set tmargin at screen .9\n"
	  "set format x \"\"\n"
	  "set format y \"%%.1e\"\n" 
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"
	  "plot \"%s\" u 1:($6) w lp ls 4 pt 7 ti \"utilde^r\""
	  ,fname2,t,t/CCC,get_xb(-NG,0),get_xb(NX+NG,0),fname,fname,fname,fname,fname);



/****************************************/
/****************************************/
/****************************************/

//  fprintf(fgnu,"\n");
//  fclose(fgnu);   
//  
//  int i=system("gnuplot plot.gp ");
//  return 0;
//}
