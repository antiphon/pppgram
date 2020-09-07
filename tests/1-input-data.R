# input test

devtools::load_all()

library(spatstat)
xl <- list(bei, ants)


check_pp(bei) # ok

check_pp(ants) # fail

check_pp(gordon) # fail

