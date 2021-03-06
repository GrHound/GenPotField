genpotfield.c

This is a small package, mainly intended for machine-learning researchers who
want to take a look at potential-field estimation as is done in Molecular Dynamics
(MD) or molecular simulations. The tool only generates a potential fields.
None of the realistic geometric problems of molecular symmetry etc. are
handled, just a 2D mountain landcape is generated using a randomly drawn
Gaussian mixture. However, some structuring is present. A main peak is drawn
randomly and is given the largest amplitude. Then, in its vicinity,
secondary peaks are drawn, and so on, recursively. Standard deviations,
number of levels and number of peaks per hierarchical level can be
specified.

This tool was generated at the Lorentz Center, Leiden, 2018, during the
workshop Molecular Simulations Meets Machine Learning and Artificial Intelligence 
from 8 Oct 2018 through 12 Oct 2018.

Lambert Schomaker - Thu Oct 11 10:35:10 CEST 2018

Package description

README.txt    This text

potfield-example.png  Example image, created by gnuplot

Makefile      Do: make
makit         or start the script ./makit to generate the binary

genpotfield.c ANSI C source
genpotfield.o Object file
genpotfield   Executable binary

demo          Bash script that generates 100 data sets. Per set:

              First it randomly draws the parameters for the hierarchical Gaussian
              mixture, then it samples 1000 points from the distribution, writing
              x y and f(x,y), it produces a regular grid for gnuplot and creates
              an SVG plot.

do-clean      Bash script to delete all the data files

do-plot       Bash script invoking gnuplot

gnuplot.gnu   Gnuplot template for interactive plotting as opposed to SVG generation.


Directories

Each data set or 'experiment' obtains a random integer ID and is put in the directories,
with four file types: parameters, sampled data, a grid for visualisation, and the .svg
image for the field. Example: potfield-10092.par, potfield-10092.samples etc.

./Par/*.par      Parameter files containing the x,y,sd and amplitudes of the peaks in the mixture
./Data/*.samples Randomly sampled data
./Grid/*.grid    Regular grid for plotting
./Plot/*.svg     Plots



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
