# README #

This is an R package that aims to develop a framework for using transmission models for TB burden estimation. Further details on use can be browsed in the markdown vignette in the vignettes folder.

If you have an R installation with a working development environment and devtools installed, you should be able to install the package by typing:

devtools::install_github('petedodd/esteban',dependencies=TRUE,build_vignettes=FALSE)

The fortran random number generation routines are linked to R's random number stream modifying and extending some fortran routines by Alan Miller.

## Acknowledgements ##

* Particular thanks to TB MAC (http://www.tb-mac.org) and also to WHO for funding the development of this work.
* Many thanks to Philippe Glaziou and Carel Pretorius for discussions informing the development of this work.


For any queries, feedback or bugs, please contact Pete Dodd: p.j.dodd+esteban@sheffield.ac.uk
