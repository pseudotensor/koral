//int
//convert_out2gif_1d(char *fname,char *fname2,int niter,ldouble t)
//{
//  FILE *fgnu=fopen("plot.gp","w");
//  char bufor[50];

/****************************************/
/****************************************/
/****************************************/
   fprintf(fgnu,
	  "set style line 1 lw 2 lc 1 lt 1\n"
	  "set style line 2 lw 3 lc 2 lt 1\n"
	  "set style line 3 lw 3 lc 3 lt 1\n"
	  "set style line 4 lw 2 lc 4 lt 1\n"
	  "set style line 10 lw 2 lc 3 lt 3\n"
	  "set style line 11 lw 2 lc 4 lt 3\n"
	  "set style line 12 lw 2 lc 5 lt 3\n"
	  "set style line 13 lw 2 lc 0 lt 3\n"
	  "set style line 14 lw 2 lc 7 lt 3\n"
	  "set term gif large size 500,600\n"
	  "set output \"%s\"\n"
	  "set size 1,1\n"
	  "set origin 0,0\n"
	  "set multiplot\n"
	  "set autoscale\n" 
	  "set label \"t=%.2f (%.2f s)\" at screen .48, .98\n"
	   "set xrange [%f:%f]\n"

	  "set lmargin at screen 0.15\n"
	  "set rmargin at screen 0.95\n"
	  "set bmargin at screen .55\n"
	  "set tmargin at screen .95\n"
	  "set format x \"%%.1f\"\n"
	  "set format y \"%%.1e\"\n" 
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"
	  "set yrange [0:4.]\n"
#ifdef SKIP_MULTIRADFLUID
	  "plot \"%s\" u 1:($20+$24) w l ls 4 ti \"E_total\", \"%s\" u 1:20 w l ls 2 ti \"E_1\", \"%s\" u 1:24 w l ls 3  ti \"E_2\" \n"
#else
	  "plot \"%s\" u 1:20 w l ls 2 ti \"E\"\n"
#endif
	  "set lmargin at screen 0.15\n"
	  "set rmargin at screen 0.95\n"
	  "set bmargin at screen .05\n"
	  "set tmargin at screen .45\n"
	  //	  "unset log y\n"
	  "set format x \"%%.1f\"\n"
	  "set format y \"%%.1e\"\n" 
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"
	  "set yrange [-1.5:1.5]\n"
#ifdef SKIP_MULTIRADFLUID
	   "plot  \"%s\" u 1:($21+$25) w l ls 4  ti \"F_total\",\"%s\" u 1:21 w l ls 2  ti \"F_1\", \"%s\" u 1:25 w l ls 3  ti \"F_2\"\n"
#else
	  "plot \"%s\" u 1:($21+1.e-80) w l ls 2  ti \"F\"\n"
#endif

	  


#ifdef SKIP_MULTIRADFLUID
	   ,fname2,t,t/CCC,get_xb(-NG,0),get_xb(NX+NG,0),fname,fname,fname,fname,fname,fname);
#else
	  ,fname2,t,t/CCC,get_xb(-NG,0),get_xb(NX+NG,0),fname,fname);
#endif



/****************************************/
/****************************************/
/****************************************/

//  fprintf(fgnu,"\n");
//  fclose(fgnu);   
//  
//  int i=system("gnuplot plot.gp ");
//  return 0;
//}
