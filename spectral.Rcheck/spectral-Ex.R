pkgname <- "spectral"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('spectral')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("decibels")
### * decibels

flush(stderr()); flush(stdout())

### Name: decibels
### Title: Convert to decibels
### Aliases: decibels dB

### ** Examples

decibels(1)  # return 0
decibels(10) # returns 10
decibels(0)  # return -Inf



cleanEx()
nameEx("idecibels")
### * idecibels

flush(stderr()); flush(stdout())

### Name: idecibels
### Title: Convert from decibels
### Aliases: idecibels idB

### ** Examples

idecibels(0)  # return 1
idecibels(10) # returns 10
idecibels(-Inf)  # return 0



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
