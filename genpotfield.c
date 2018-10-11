/* genpotfield.c

Program to generate simple artificial force fields using a Gaussian
mixture, around a number of hierarchically generated centers 
where the argument would be the radius r from this centers, 
for each pixel x,y.

Generating parameters: 

  genpotfield -generate [-s seedval|clock] [-h <r>] [-varfrac r] \
                         -p lev1_nmeans lev1_sd lev1_ampl \
                         -p lev2_nmeans lev2_sd lev2_ampl \
			   .
			   .
		         >example.par

Sampling a distribution:

  genpotfield -sample -parfile example.par -nsamples 1000 > example.samples

Creating a regular grid:

  genpotfield -sample -parfile example.par -nsamples grid > example.samples
  
  Options in -generate:
  
       -h (hetereogeneity) is a random variation in npeaks (default = 0)
            as a proportion of a given nmean. The number of peaks
            will then uniformly vary from nmean -/+ heterogeneity*nmean,
            with the provision that the minimum number of members is 2.
            (all branches have the same length). Default value 0.0: 
            the given number of peaks for the level without random variation.
	    
       -varfrac Proportion for sampling the sd value for a peak. It is
                based on varfrac times the sd of the mother peak in
                the hierarchical scheme. Default value 0.1

   Options in -sample
       -nsamples i or 'grid'   Number of randomly drawn samples with
                               three columns in the output file:
			        x y f(x,y)
                                 
                               If instead of an integer the keyword 
                               is given, full grid is separated, with
                               empty lines between rows such that 
			       gnuplot immediately can use splot to
			       show the field in 3D.
       -ngrid i   Number of grid points, one sided, per axis. 
                  The range for x and y will be [-1:+1] in 2*ngrid points, 
		  for both axes x and y, yielding (2*ngrid)^2 points.
		  Example: ngrid=30 will generate 3600 points and the
		  output file will have 3660 records, i.e., including
		  the empty separator lines intended for gnuplot.

Author: Lambert Schomaker
Revised: Oct 2018
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>


#define MAXLEV 5000 /* Actually: MAXLEV on generate, MAXPEAK on sample */

/* #define DEBUG */

/* void srand48();
   double drand48(); */

#define VARFRAC 0.1
#define NGRID 30

double rangau(double sd)
{
	int i;
	double r;
	
	r = -6.0;
	for(i = 0; i < 12; ++i) {
		r += drand48();
	}
	return(r*sd);
}

double ranbipolar()
{
	double r;
	
	r = 2.0 * (drand48() - 0.5);
	return(r);
}

double perturb(double x, double varfrac)
{
        double r;
        
        r = x + varfrac * ranbipolar() * x;
        
        if(r <= 0.0) r = 0.001; /* sd and ampl cannot be zero or negative */
        
        return(r);
}

void expect(char item[], char key[])
{
      if(strcmp(item,key) != 0) {
         fprintf(stderr,"Expecting '%s' not [%s] in input file\n", key,item);
         exit(1);
      }
}

void read_field(char fname[],double x[], double y[], double sd[], double ampl[], int *npeaks)
{
   FILE *fp;
   int itag;
   char item[30], str[30];
   
   fp = fopen(fname,"r");
   if(fp == NULL) {
      fprintf(stderr,"Error opening input file [%s]\n", fname);
      exit(1);
   }
   *npeaks = 0;
   while(fscanf(fp,"%s %d",item,&itag) != EOF) {
      expect(item,"itag");
      fscanf(fp,"%s %lf",item,&x[*npeaks]);
      expect(item,"xpk");
      fscanf(fp,"%s %lf",item,&y[*npeaks]);
      expect(item,"ypk");
      fscanf(fp,"%s %lf",item,&sd[*npeaks]);
      expect(item,"sd");
      fscanf(fp,"%s %lf",item,&ampl[*npeaks]);
      expect(item,"ampl");
      fscanf(fp,"%s %s",item,str);
      expect(item,"pk");
      fscanf(fp,"%s",item);
      if(strcmp(item,"-->") !=0 && strcmp(item,"*") != 0) {
         fprintf(stderr,"Expecting terminator --> or *, not [%s]\n", item);
      }
#ifdef DEBUG
      fprintf(stderr,"itag %d xpk %f ypk %f sd %f ampl %f \n"
                              , itag, x[*npeaks], y[*npeaks], sd[*npeaks], ampl[*npeaks]);
#endif
      *npeaks += 1;
   }
   fclose(fp);
}


void gen_field(int itag,int nmeans[MAXLEV],double sd[MAXLEV],double peakampl[MAXLEV]
               ,int nlev,int ilev,int imean
               ,double xpk, double ypk, double het, double varfrac)
{
	int i;
	double xxpk, yypk;
	double rmn, v;
	
        v = varfrac;
	
	if(ilev >= nlev) {
#ifdef DEBUG
                fprintf(stderr,"ilev=%d >= nlev=%d\n", ilev, nlev);
#endif
                fprintf(stdout,"itag %d xpk %f ypk %f sd %f ampl %f pk %d/%d *\n"
                              , itag, xpk, ypk, perturb(sd[ilev-1],v), perturb(peakampl[ilev-1],v), itag % 10,nmeans[ilev-1]);

	} else {
	        if(ilev > 0) {
                  fprintf(stdout,"itag %d xpk %f ypk %f sd %f ampl %f pk %d/%d -->\n"
                              , itag, xpk, ypk, perturb(sd[ilev-1],v), perturb(peakampl[ilev-1],v), itag % 10, nmeans[ilev-1]);
                } else {
                  i = 0;
#ifdef DEBUG                  
                  fprintf(stdout,"itag=%d xpk=%f ypk=%f sd=%f ampl=%f %d/%d [root]\n"
                              , itag, xpk, ypk, 1.0, 1.0, 1, 1);
#endif
                }
#ifdef DEBUG
                fprintf(stderr,"ilev=%d < nlev=%d, itag=%d\n", ilev, nlev, itag);
#endif
	        xxpk = xpk;
	        yypk = ypk;
		
		rmn =  nmeans[ilev];	
		if(het != 0.0) {
			rmn = rmn + het * rmn * (drand48() - 0.5) * 2.0;
			if (rmn < 2.0) rmn = 2.0;
		}

		for(i = 0; i < (int) rmn; ++i) {
			xxpk = xpk + rangau(sd[ilev]);
			yypk = ypk + rangau(sd[ilev]);
			gen_field(10*itag+i+1,nmeans,sd,peakampl,nlev,ilev+1,imean
			                                ,xxpk,yypk,het,varfrac);
		}
	}
}

void usage(char reason[])
{
printf("?? %s\n", reason);
printf("genpotfield() (Lambert Schomaker/RuG, 2018)\n\n");
printf("Usage: genpotfield [-id] [-s seedval|clock] [-repl previous.log] [-h <r>] L1_nmeans L1_sd   L2_nmeans L2_sd ...   >out.dat\n");
printf("\n");
printf("       Give flag -help for more information\n");

}

int getseed_clock()
{
    struct timeval ct;

    gettimeofday(&ct, NULL); 
    return ct.tv_usec;   
}

void parse_opt(int argc, char *argv[], char mode[], double *het, double *varfrac
                                     , int *nsamples, int *ngrid, char parfile[], int *iarg0)
{

   int seed;
   char *eop; 
   	
	*iarg0 = 1;
	*het = 0.0;
	*varfrac = VARFRAC;
	*ngrid = NGRID;
	*nsamples = 10;
	strcpy(mode,"-generate");
	
	while (*iarg0 < argc && argv[*iarg0][0]=='-') {
	        /* Flag -p ends the preliminary flag sequence */
	        if(strcmp(argv[*iarg0], "-p") == 0) {
	           break;
	        }
	        
	        if(strcmp(argv[*iarg0], "-generate") == 0) {
	        	strcpy(mode,argv[*iarg0]);
	        	*iarg0 += 1;
	        	
	        } else if(strcmp(argv[*iarg0], "-sample") == 0 ) {
	        	strcpy(mode,argv[*iarg0]);
	        	*iarg0 += 1;
	        	
	        } else if(strstr(argv[*iarg0], "-help") != NULL) {
	        	usage("-help"); /* was help(); */
	        	exit(1);
	        	
	        } else if(strcmp(argv[*iarg0], "-h") == 0) {
	                if(strcmp(mode,"-h") == 0) {
	                   fprintf(stderr,"-h (heterogeneity of npeaks/level) has no meaning in mode -sample, only in mode -generate\n");
	                   exit(1);
	                }
			*het = strtod(argv[*iarg0+1], &eop);
			if(*het < 0.0 || *het > 1.0) {
				fprintf(stderr,"Expecting -h r with r being a proportion, i.e., r = [0.0,1.0]  Given was: %s\n", argv[*iarg0+1]);
				exit(1);
			}
			if(eop == NULL) {
				fprintf(stderr,"Expecting a number at [%s]\n", argv[*iarg0+1]);
				exit(1);
			}
			*iarg0 += 2;
	        	
	        } else if(strcmp(argv[*iarg0], "-varfrac") == 0) {
	                if(strcmp(mode,"-sample") == 0) {
	                   fprintf(stderr,"-varfrac has no meaning in mode -sample, only in mode -generate\n");
	                   exit(1);
	                }
			*varfrac = strtod(argv[*iarg0+1], &eop);
			if(*varfrac < 0.0) {
				fprintf(stderr,"Expecting -varfrac r with r being a proportion, i.e., r = [0.0,1.0]  Given was: %s\n", argv[*iarg0+1]);
				exit(1);
			}
			if(eop == NULL) {
				fprintf(stderr,"Expecting a number at [%s]\n", argv[*iarg0+1]);
				exit(1);
			}
			*iarg0 += 2;
	        	
	        } else if(strcmp(argv[*iarg0], "-parfile") == 0) {
			strcpy(parfile,argv[*iarg0+1]);
			*iarg0 += 2;
			fprintf(stderr,"Parameter input file [%s]\n", parfile);
	        	
	        } else if(strcmp(argv[*iarg0], "-nsamples") == 0) {
	                if(strcmp(argv[*iarg0+1],"grid") == 0) {
	                   *nsamples = -999;
	                } else {
   			   *nsamples = strtol(argv[*iarg0+1], &eop, 10);
			   if(*nsamples < 1) {
				fprintf(stderr,"Expecting -nsamples i Given was: %s\n", argv[*iarg0+1]);
				exit(1);
			   }
			   if(eop == NULL) {
				fprintf(stderr,"Expecting a number at [%s]\n", argv[*iarg0+1]);
				exit(1);
			   }
			}
			*iarg0 += 2;
			fprintf(stderr,"Nsamples %d\n", *nsamples);
	        	
	        } else if(strcmp(argv[*iarg0], "-ngrid") == 0) {
	                if(strcmp(mode,"-generate") == 0) {
	                   fprintf(stderr,"-ngrid has no meaning in mode -generate, only in mode -sample\n");
	                   exit(1);
	                }
 			*ngrid = strtol(argv[*iarg0+1], &eop, 10);
			if(*ngrid < 1) {
			    fprintf(stderr,"Expecting -ngrid i Given was: %s\n", argv[*iarg0+1]);
			    exit(1);
			}
			if(eop == NULL) {
			    fprintf(stderr,"Expecting a number at [%s]\n", argv[*iarg0+1]);
			    exit(1);
			}
			*iarg0 += 2;
			fprintf(stderr,"Ngrid %d\n", *ngrid);
	        	
		} else {
			if(strcmp(argv[*iarg0], "-s") == 0) {
			     if(strcmp(argv[*iarg0+1],"clock") == 0) {
			        seed = getseed_clock();
			     }  else {
				seed = atoi(argv[*iarg0+1]);
				if (seed<0) {
					fprintf(stderr,"Expecting a number at [%s] or 'clock'\n", argv[*iarg0+1]);
					exit(1);
				}
			     }
			     srand48(seed);
			     fprintf(stderr,"Seed %d\n", seed);
			     *iarg0 += 2;
			} else {
				fprintf(stderr,"Expecting a known flag at [%s]\n"
				              ,argv[*iarg0]);
				usage("unexpected flag");
				exit(1);
			}
		}
	}

	if(strcmp(mode,"-generate") == 0) {
	   if(argc - *iarg0 + 1 < 2) {
		usage("no relevant flags found");
		exit(1);
	   }
	}

#ifdef AAP	
	if(((argc - *iarg0 + 1) % 2) != 1) {
                usage("Error: Give tuples -p nmeans sd ampl -p nmeans sd ampl ...");
		exit(1);
	}
#endif

	if( ((argc - *iarg0)/2) > MAXLEV) {
		fprintf(stderr,"Maximum number of levels=%d\n", MAXLEV);
		exit(1);
	}
}

void all_gen_field(int argc, char *argv[], int iarg0, double het, double varfrac)
{
   int ilev, iarg, imean, itag, nlev;
   double xpk, ypk;
   int nmeans[MAXLEV];
   double sd[MAXLEV];	
   double peakampl[MAXLEV];	
   char *eop;
   
/* For all levels: enter number of centroids and their sd */

	ilev = 0;
	for(iarg = iarg0; iarg < argc; iarg += 4) {
	        if(strcmp(argv[iarg],"-p") != 0) {
		  fprintf(stderr,"Expecting '-p' at [%s]\n", argv[iarg]);
		  exit(1);
	        }
	        
		nmeans[ilev] = strtol(argv[iarg+1], &eop, 10);
		if(eop == NULL) {
		  fprintf(stderr,"Expecting nmeans at [%s]\n", argv[iarg]);
		  exit(1);
		}
		
		sd[ilev] = strtod(argv[iarg+2], &eop);
		if(eop == NULL) {
		  fprintf(stderr,"Expecting sd at [%s]\n", argv[iarg+1]);
		  exit(1);
		}
		
		peakampl[ilev] = strtod(argv[iarg+3], &eop);
		if(eop == NULL) {
		  fprintf(stderr,"Expecting peak amplitude at [%s]\n", argv[iarg+2]);
		  exit(1);
		}
#ifdef DEBUG
		fprintf(stderr,"ilev=%d nm=%d sd=%f amp=%f\n"
		              , ilev, nmeans[ilev], sd[ilev], peakampl[ilev]);
#endif
		++ilev;
	}
	nlev = ilev;
#ifdef DEBUG
        fprintf(stderr,"nlev=%d (inferred from argv)\n", nlev);
#endif
	
	
	ilev =  0;
	imean = 0;
	xpk = 0.0;
	ypk = 0.0;
	itag = 0;

	gen_field(itag,nmeans,sd,peakampl,nlev,ilev,imean,xpk,ypk,het,varfrac);

}

double truncated_rangau(double sd, double rmin, double rmax)
{
   double r;
   
   r = rmin - 1.0;
   while(r < rmin || r > rmax) {
      r = rangau(sd);
   }
   return(r);
}


#define PI ((double) 3.141592653589793)

double gauss(double x, double mu, double sd)
{
   double f,q,g;
   
   f = 1.0 / sqrt(2. * PI * sd * sd);
   
   q = - (x - mu) * (x - mu) / (2. * sd * sd);
   
   g = f * exp(q);
   
   return(g);
}

double sample(double xx, double yy, double x[], double y[], double sd[], double ampl[], int npeaks)
{
   int i;
   double p,r,dx,dy;
   
   p = 0.0;
   for(i = 0; i < npeaks; i++) {
      dx = xx - x[i];
      dy = yy - y[i];
      r = sqrt(dx*dx + dy*dy);
      p = p + ampl[i] * gauss(r,0.0,sd[i]);
   }
   
   return(p);
}


void all_sample_field(int argc, char *argv[], int nsamples, int ngrid, char fname[])
{
   int npeaks, isample, i, j;
   double x[MAXLEV], y[MAXLEV], sd[MAXLEV], ampl[MAXLEV];
   double xx, yy;
   
   read_field(fname,x,y,sd,ampl,&npeaks);
   
   if(nsamples != -999) {
      
     /* Random samples x,y from the field */
   
     for(isample = 0; isample < nsamples; ++isample) {
       xx = truncated_rangau(1.0,-1.0,1.0);
       yy = truncated_rangau(1.0,-1.0,1.0);
      
       fprintf(stdout,"%f %f %f\n", xx,yy,sample(xx,yy,x,y,sd,ampl,npeaks));
     }
     
   } else {
      
     /* Grid: full field, regularly sampled */
   
      for(i = -ngrid; i < ngrid; ++i) {
         for(j = -ngrid; j < ngrid; ++j) {
            xx = (double) i / (double) (2*ngrid);
            yy = (double) j / (double) (2*ngrid);
            fprintf(stdout,"%f %f %f\n", xx,yy,sample(xx,yy,x,y,sd,ampl,npeaks));
         }
         fprintf(stdout,"\n");
      }
   }
   
}

int main(int argc, char *argv[])
{
	int iarg0, nsamples, ngrid;
	double het, varfrac;
	char mode[100], parfile[1000];

	if(argc < 2) {
		usage("not enough arguments");
		exit(1);
	}
	
        parse_opt(argc, argv, mode, &het, &varfrac, &nsamples, &ngrid, parfile, &iarg0);
        
        if(strcmp(mode,"-generate") == 0) {
          fprintf(stderr,"Mode: -generate\n");
          all_gen_field(argc,argv,iarg0,het,varfrac);
          
        } else if(strcmp(mode,"-sample") == 0) {
          fprintf(stderr,"Mode: -sample Input file [%s]\n", parfile);

          all_sample_field(argc,argv,nsamples,ngrid,parfile); 
          
        } else {
          fprintf(stderr,"genpotfield() Unknown mode [%s] Should be -generate or -sample\n", mode);
          exit(1); 
        }
   return(0);
}
/* End genpotfield.c */


