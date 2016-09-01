.onAttach <- function(libname='esteban',pkgname='esteban'){
    packageStartupMessage("Welcome to esteban!\nAn experimental package for using transmission models to estimate TB burden.\n(c) P.J. Dodd 2016: released under CC BY-SA 4.0 licence (see licence or https://creativecommons.org/licenses/by-sa/4.0/).\nSee GitHub README and vignette for more details on use (https://github.com/petedodd/esteban).")
    data(rngs)                          #load up
}

## This is almost certainlty not best practice and would be good to do away with
.onLoad <- function(libname='esteban',pkgname='esteban'){

}
