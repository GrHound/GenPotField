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

