# GenPotField
This is a small package, mainly intended for machine-learning researchers who want to take a look at potential-field estimation as is done in Molecular Dynamics (MD) or molecular simulations. The tool only generates a potential fields. None of the realistic geometric problems of molecular symmetry etc. are handled, just a 2D mountain landcape is generated using a randomly drawn Gaussian mixture. However, some structuring is present. A main peak is drawn randomly and is given the largest amplitude. Then, in its vicinity, secondary peaks are drawn, and so on, recursively. Standard deviations, number of levels and number of peaks per hierarchical level can be specified.

The Bash script 'demo' will generate 100 potential fields (parameters, samples, grid and plot), with each set
or experiment containing 1000 samples. The data sets have a random name of the form potfield-[RANDNUM].par etc. Data are organized in four directories: ./Par ./Data ./Grid and ./Plot
	
Working environment is assumed to be Linux with bash, sed and gnuplot.	

See: potfield-example.png

![Example image](potfield-example.png)

![Example image](potfield-example2.png)

genpotfield.c

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
  
       -s <integer> or 'clock' Seed value for drand48()
       
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
			        x y f(x,y)  like so:
		
                                0.286470 0.483686 1.247119 
                               -0.660440 0.307425 0.042538 
                               -0.068826 0.930838 0.046059 
                                0.513680 0.822562 0.027226 

                        If instead of an integer the keyword 'grid'
                        is given, a full grid is generated, with
                        empty lines between rows such that 
                        gnuplot immediately can use splot to
                        show the field in 3D. See -ngrid.
			       
       -ngrid <i> Number of grid points, one sided, per axis. 
                  The range for x and y will be [-1:+1] in 2*ngrid points, 
                  for both axes x and y, yielding (2*ngrid)^2 points.
                  Example: ngrid=30 will generate 3600 points and the
                  output file will have 3660 records, i.e., including
                  the empty separator lines intended for gnuplot.

Author: Lambert Schomaker
Revised: Oct 2018

