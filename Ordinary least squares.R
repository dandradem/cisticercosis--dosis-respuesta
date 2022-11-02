rm(list=ls())
#OBSERVED DATA

##Oral

oral <- data.frame(x1 = c(0, 10, 28, 100, 103, 113, 244, 1000, 
                          10000, 10168, 18189, 20000, 39800, 100000),
                   x2 = c(0, 10, 28, 100, 103, 113, 244, 1000, 
                          10000, 10168, 18189, 20000, 37582, 100000),
                   any = c(0, 0.75, 0.833, 0.667, 1, 0.75, 0.8, 0.909,
                           1, 1, 0.875, 1, 0.933, 1), 
                   viable = c(0, 0.25, 0.667, 0.222, 1, 0.25, 0.767,
                              0.273, 1, 0.875, 0.625, 1, 0.9, 1), 
                   brain = c(0, 0, 0.167, 0, 0.167, 0, 0.067, 0, 0.222,
                             0.125, 0.375, 0.8, 0.5, 0.6))

oral$any_y <- mean(oral$any)
oral$viable_y <- mean(oral$viable)
oral$brain_y <- mean(oral$brain)

##Proglottids

prog <- data.frame(x1 = c(0, 10168, 18189, 39800),
                   x2 = c(0, 10168, 18189, 37582),
                   any = c(0, 1, .875, .933), 
                   viable = c(0, .875, .625, .9), 
                   brain = c(0, .125, .375, 0.5))

prog$any_y <- mean(prog$any)
prog$viable_y <- mean(prog$viable)
prog$brain_y <- mean(prog$brain)

##Eggs

egg <- data.frame(x = c(0, 10, 100, 1000, 10000, 20000, 100000), 
                  any = c(0, .75, .667, .909, 1, 1, 1), 
                  viable = c(0, .25, .222, .273, 1, 1, 1), 
                  brain = c(0, 0, 0, 0, .222, .8, .6))

egg$any_y <- mean(egg$any)
egg$viable_y <- mean(egg$viable)
egg$brain_y <- mean(egg$brain)

##Beetles

bee <- data.frame(x = c(0, 28, 103, 113, 244), 
                  any = c(0, .833, 1, .75, .8), 
                  viable = c(0, .667, 1, .25, .767),
                  brain = c(0, .167, .167, 0, .067))

bee$any_y <- mean(bee$any)
bee$viable_y <- mean(bee$viable)
bee$brain_y <- mean(bee$brain)

##Carotid

car <- data.frame(x = c(0, 2500, 5000, 10000, 45000, 50000),
                  any = c(0, 1, 1, .818, 1, 1), 
                  viable = c(0, 1, 1, .818, 1, 1), 
                  brain = c(0, .6, 1, .636, 1, .8))

car$any_y <- mean(car$any)
car$viable_y <- mean(car$viable)
car$brain_y <- mean(car$brain)

#TWO-PARAMETER LOG-LOGISTIC

##Parameters
predict.ll2 <- function(dose, bslope, bid50) {
  result <- 1/(1 + exp(bslope * (log(dose) - log(bid50))))
  return(result)
}

##BRAIN CYSTS

###Eggs
bslope_egg_brain <- -0.79
bid50_egg_brain <- 5.74e04

egg$ll2_brain <- predict.ll2(egg$x, bslope_egg_brain, bid50_egg_brain)





#SIMPLE LOGISTIC REGRESSION

##Parameter
predict.slr <- function(dose, binter, bslope) {
  result <- 1/(1 + exp(-(binter + bslope * dose)))
  return(result)
}

##ANY CYST

###Proglottids
binter_any_prog <- 3.339e-01
bslope_any_prog <- 6.926e-05

prog$slr_any <- predict.slr(prog$x1, binter_any_prog, bslope_any_prog)

###Eggs
binter_any_egg <- -0.374904
bslope_any_egg <- 0.002897
  
egg$slr_any <- predict.slr(egg$x, binter_any_egg, bslope_any_egg)

###Beetles
binter_any_bee <- 1.892815
bslope_any_bee <- -0.002040

bee$slr_any <- predict.slr(bee$x, binter_any_bee, bslope_any_bee)

###Carotid
binter_any_car <- 2.196e+00
bslope_any_car <- 2.886e-05

car$slr_any <- predict.slr(car$x, binter_any_car, bslope_any_car)


##VIABLE CYSTS

###Proglottids
binter_via_prog <- -2.737e-01
bslope_via_prog <- 6.458e-05
  
prog$slr_via <- predict.slr(prog$x1, binter_via_prog, bslope_via_prog)

###Eggs
binter_via_egg <- -1.8287013
bslope_via_egg <- 0.0009082
  
egg$slr_via <- predict.slr(egg$x, binter_via_egg, bslope_via_egg)

###Beetles
binter_via_bee <- 0.208357
bslope_via_bee <- 0.003639
  
bee$slr_via <- predict.slr(bee$x, binter_via_bee, bslope_via_bee)

###Carotid
binter_via_car <- 2.196e+00
bslope_via_car <- 2.886e-05
  
car$slr_via <- predict.slr(car$x, binter_via_car, bslope_via_car)


##BRAIN CYSTS

###Proglottids
binter_brain_prog <- -2.256e+00
bslope_brain_prog <- 6.309e-05
  
prog$slr_brain <- predict.slr(prog$x1, binter_brain_prog, bslope_brain_prog)

###Eggs
binter_brain_egg <- -2.173e+00
bslope_brain_egg <- 3.206e-05
  
egg$slr_brain <- predict.slr(egg$x, binter_brain_egg, bslope_brain_egg)

###Beetles
binter_brain_bee <- -1.786151
bslope_brain_bee <- -0.003885
  
bee$slr_brain <- predict.slr(bee$x, binter_brain_bee, bslope_brain_bee)

###Carotid
binter_brain_car <- 9.296e-01
bslope_brain_car <- 1.119e-05
  
car$slr_brain <- predict.slr(car$x, binter_brain_car, bslope_brain_car)





#EXPONENTIAL REGRESSION

##Parameter
predict.er <- function(dose, r) {
  result <- 1 - exp(-r * dose)
  return(result)
}

##ANY CYST

###Eggs
r_any_egg <- 0.015

egg$er_any <- predict.er(egg$x, r_any_egg)


##VIABLE CYSTS

###Proglottids
r_via_prog <- 1.01e-04

prog$er_via <- predict.er(prog$x1, r_via_prog)

###Eggs
r_via_egg <- 3.62e-04

egg$er_via <- predict.er(egg$x, r_via_egg)

###Beetles
r_via_bee <- 0.007

bee$er_via <- predict.er(bee$x, r_via_bee)


##BRAIN CYSTS

###Proglottids
r_brain_prog <- 1.94e-05

prog$er_brain <- predict.er(prog$x2, r_brain_prog)

###Eggs
r_brain_egg <- 4.126e-05

egg$er_brain <- predict.er(egg$x, r_brain_egg)

###Carotid
r_brain_car <- 3.28e-04

car$er_brain <- predict.er(car$x, r_brain_car)




#APPROXIMATE BETA-POISSON

##ANY CYST

###Proglottids
prog$abp_any <- c(1, 1, 1, 1)

###Eggs
egg$abp_any <- c(0.003, 0.58, 0.843, 0.9381, 0.9752, 0.9805, 0.9907)

###Beetles
bee$abp_any <- c(0.8037, 0.969, 0.9777, 0.979, 0.985)

###Carotid
car$abp_any <- c(1, 1, 1, 1, 1, 1)


##VIABLE CYSTS

###Proglottids
prog$abp_via <- c(0, 0.77, 0.82, 0.88)

###Eggs
egg$abp_via <- c(0, 0.014, 0.1205, 0.5727, 0.9212, 0.9583, 0.9911)

###Beetles
bee$abp_via <- c(0.08, 0.706, 0.786, 0.79, 0.83)

###Carotid
car$abp_via <- c(1, 1, 1, 1, 1, 1)


##BRAIN CYSTS

###Proglottids
prog$abp_brain <- c(0, 0.282, 0.345, 0.42)

###Eggs
egg$abp_brain <- c(0, 0.0024, 0.0219, 0.1265, 0.3259, 0.3903, 0.5373)

###Beetles
bee$abp_brain <- c(0.2032, 0.424, 0.472, 0.474, 0.51)

###Carotid
car$abp_brain <- c(0.3425, 0.871, 0.8830, 0.8954, 0.919, 0.921)





#EXACT BETA-POISSON

##ANY CYST

###Oral
oral$ebp_any <- c(0, 0.6665, 0.7414, 0.8098, 0.8118, 0.8158, 0.8469, 0.8901, 
                  0.9338, 0.9340, 0.9425, 0.9437, 0.9518, 0.9611)

###Proglottids
prog$ebp_any <- c(0, 0.9347, 0.9348, 0.9413)

###Eggs
egg$ebp_any <- c(0, 0.532, 0.8157, 0.927, 0.9705, 0.9768, 0.9878)

###Beetles
bee$ebp_any <- c(0, 0.7877, 0.82, 0.82, 0.84)

###Carotid
car$ebp_any <- c(0, 0.9285, 0.9285, 0.9325, 0.9558, 0.9571)


##VIABLE CYSTS

###Oral
oral$ebp_via <- c(0, 0.2918, 0.4177, 0.5522, 0.5567, 0.5652, 0.629, 0.7281,
                  0.8333, 0.835, 0.8535, 0.8564, 0.8765, 0.899)

###Proglottids
prog$ebp_via <- c(0, 0.77, 0.824, 0.8881)

###Eggs
egg$ebp_via <- c(0, 0.013, 0.1205, 0.5718, 0.9216, 0.9583, 0.9913)

###Beetles
bee$ebp_via <- c(0, 0.596, 0.6798, 0.686, 0.739)

###Carotid
car$ebp_via <- c(0, 0.9285, 0.9285, 0.9325, 0.9558, 0.9571)


##BRAIN CYSTS

###Oral
oral$ebp_brain <- c(0, 0.0038, 0.0102, 0.032, 0.0333, 0.0361, 0.0647, 0.1551,
                    0.3366, 0.34, 0.3853, 0.3926, 0.4385, 0.5065)

###Proglottids
prog$ebp_brain <- c(0, 0.2837, 0.347, 0.4156)

###Eggs
egg$ebp_brain <- c(0, 0.0024, 0.0218, 0.1265, 0.3254, 0.3901, 0.5372)

###Beetles
bee$ebp_brain <- c(0, 0.06, 0.075, 0.0762, 0.08)

###Carotid
car$ebp_brain <- c(0, 0.713, 0.7482, 0.7592, 0.7933, 0.7961)





#SUM OF SQUARED ERRORS OF PREDICTION (SSE)

##Two-parameter log-logistic
sum((egg$brain - egg$ll2_brain)^2)


##Simple logistic regression
sum((prog$any - prog$slr_any)^2)
sum((egg$any - egg$slr_any)^2)
sum((bee$any - bee$slr_any)^2)
sum((car$any - car$slr_any)^2)

sum((prog$viable - prog$slr_via)^2)
sum((egg$viable - egg$slr_via)^2)
sum((bee$viable - bee$slr_via)^2)
sum((car$viable - car$slr_via)^2)


sum((prog$brain - prog$slr_brain)^2)
sum((egg$brain - egg$slr_brain)^2)
sum((bee$brain - bee$slr_brain)^2)
sum((car$brain - car$slr_brain)^2)


##Exponential regression
sum((egg$any - egg$er_any)^2)

sum((prog$viable - prog$er_via)^2)
sum((egg$viable - egg$er_via)^2)
sum((bee$viable - bee$er_via)^2)

sum((prog$brain - prog$er_brain)^2)
sum((egg$brain - egg$er_brain)^2)
sum((car$brain - car$er_brain)^2)


##Approximate beta-poisson
sum((prog$any - prog$abp_any)^2)
sum((egg$any - egg$abp_any)^2)
sum((bee$any - bee$abp_any)^2)
sum((car$any - car$abp_any)^2)

sum((prog$viable - prog$abp_via)^2)
sum((egg$viable - egg$abp_via)^2)
sum((bee$viable - bee$abp_via)^2)
sum((car$viable - car$abp_via)^2)


sum((prog$brain - prog$abp_brain)^2)
sum((egg$brain - egg$abp_brain)^2)
sum((bee$brain - bee$abp_brain)^2)
sum((car$brain - car$abp_brain)^2)


##Exact beta-poisson
sum((oral$any - oral$ebp_any)^2)
sum((prog$any - prog$ebp_any)^2)
sum((egg$any - egg$ebp_any)^2)
sum((bee$any - bee$ebp_any)^2)
sum((car$any - car$ebp_any)^2)

sum((oral$viable - oral$ebp_via)^2)
sum((prog$viable - prog$ebp_via)^2)
sum((egg$viable - egg$ebp_via)^2)
sum((bee$viable - bee$ebp_via)^2)
sum((car$viable - car$ebp_via)^2)


sum((oral$brain - oral$ebp_brain)^2)
sum((prog$brain - prog$ebp_brain)^2)
sum((egg$brain - egg$ebp_brain)^2)
sum((bee$brain - bee$ebp_brain)^2)
sum((car$brain - car$ebp_brain)^2)





#COEFFICIENT OF DETERMINATION (R2)

##Two-parameter log-logistic
(sum((egg$ll2_brain - egg$brain_y)^2) / (sum((egg$ll2_brain - egg$brain_y)^2) + sum((egg$brain - egg$ll2_brain)^2)))


##Simple logistic regression
(sum((prog$slr_any - prog$any_y)^2) / (sum((prog$slr_any - prog$any_y)^2) + sum((prog$any - prog$slr_any)^2)))
(sum((egg$slr_any - egg$any_y)^2) / (sum((egg$slr_any - egg$any_y)^2) + sum((egg$any - egg$slr_any)^2)))
(sum((bee$slr_any - bee$any_y)^2) / (sum((bee$slr_any - bee$any_y)^2) + sum((bee$any - bee$slr_any)^2)))
(sum((car$slr_any - car$any_y)^2) / (sum((car$slr_any - car$any_y)^2) + sum((car$any - car$slr_any)^2)))

(sum((prog$slr_via - prog$viable_y)^2) / (sum((prog$slr_via - prog$viable_y)^2) + sum((prog$viable - prog$slr_via)^2)))
(sum((egg$slr_via - egg$viable_y)^2) / (sum((egg$slr_via - egg$viable_y)^2) + sum((egg$viable - egg$slr_via)^2)))
(sum((bee$slr_via - bee$viable_y)^2) / (sum((bee$slr_via - bee$viable_y)^2) + sum((bee$viable - bee$slr_via)^2)))
(sum((car$slr_via - car$viable_y)^2) / (sum((car$slr_via - car$viable_y)^2) + sum((car$viable - car$slr_via)^2)))


(sum((prog$slr_brain - prog$brain_y)^2) / (sum((prog$slr_brain - prog$brain_y)^2) + sum((prog$brain - prog$slr_brain)^2)))
(sum((egg$slr_brain - egg$brain_y)^2) / (sum((egg$slr_brain - egg$brain_y)^2) + sum((egg$brain - egg$slr_brain)^2)))
(sum((bee$slr_brain - bee$brain_y)^2) / (sum((bee$slr_brain - bee$brain_y)^2) + sum((bee$brain - bee$slr_brain)^2)))
(sum((car$slr_brain - car$brain_y)^2) / (sum((car$slr_brain - car$brain_y)^2) + sum((car$brain - car$slr_brain)^2)))


##Exponential regression
(sum((egg$er_any - egg$any_y)^2) / (sum((egg$er_any - egg$any_y)^2) + sum((egg$any - egg$er_any)^2)))

(sum((prog$er_via - prog$viable_y)^2) / (sum((prog$er_via - prog$viable_y)^2) + sum((prog$viable - prog$er_via)^2)))
(sum((egg$er_via - egg$viable_y)^2) / (sum((egg$er_via - egg$viable_y)^2) + sum((egg$viable - egg$er_via)^2)))
(sum((bee$er_via - bee$viable_y)^2) / (sum((bee$er_via - bee$viable_y)^2) + sum((bee$viable - bee$er_via)^2)))

(sum((prog$er_brain - prog$brain_y)^2) / (sum((prog$er_brain - prog$brain_y)^2) + sum((prog$brain - prog$er_brain)^2)))
(sum((egg$er_brain - egg$brain_y)^2) / (sum((egg$er_brain - egg$brain_y)^2) + sum((egg$brain - egg$er_brain)^2)))
(sum((car$er_brain - car$brain_y)^2) / (sum((car$er_brain - car$brain_y)^2) + sum((car$brain - car$er_brain)^2)))


##Approximate beta-poisson
(sum((prog$abp_any - prog$any_y)^2) / (sum((prog$abp_any - prog$any_y)^2) + sum((prog$any - prog$abp_any)^2)))
(sum((egg$abp_any - egg$any_y)^2) / (sum((egg$abp_any - egg$any_y)^2) + sum((egg$any - egg$abp_any)^2)))
(sum((bee$abp_any - bee$any_y)^2) / (sum((bee$abp_any - bee$any_y)^2) + sum((bee$any - bee$abp_any)^2)))
(sum((car$abp_any - car$any_y)^2) / (sum((car$abp_any - car$any_y)^2) + sum((car$any - car$abp_any)^2)))

(sum((prog$abp_via - prog$viable_y)^2) / (sum((prog$abp_via - prog$viable_y)^2) + sum((prog$viable - prog$abp_via)^2)))
(sum((egg$abp_via - egg$viable_y)^2) / (sum((egg$abp_via - egg$viable_y)^2) + sum((egg$viable - egg$abp_via)^2)))
(sum((bee$abp_via - bee$viable_y)^2) / (sum((bee$abp_via - bee$viable_y)^2) + sum((bee$viable - bee$abp_via)^2)))
(sum((car$abp_via - car$viable_y)^2) / (sum((car$abp_via - car$viable_y)^2) + sum((car$viable - car$abp_via)^2)))


(sum((prog$abp_brain - prog$brain_y)^2) / (sum((prog$abp_brain - prog$brain_y)^2) + sum((prog$brain - prog$abp_brain)^2)))
(sum((egg$abp_brain - egg$brain_y)^2) / (sum((egg$abp_brain - egg$brain_y)^2) + sum((egg$brain - egg$abp_brain)^2)))
(sum((bee$abp_brain - bee$brain_y)^2) / (sum((bee$abp_brain - bee$brain_y)^2) + sum((bee$brain - bee$abp_brain)^2)))
(sum((car$abp_brain - car$brain_y)^2) / (sum((car$abp_brain - car$brain_y)^2) + sum((car$brain - car$abp_brain)^2)))


##Exact beta-poisson
(sum((oral$ebp_any - oral$any_y)^2) / (sum((oral$ebp_any - oral$any_y)^2) + sum((oral$any - oral$ebp_any)^2)))
(sum((prog$ebp_any - prog$any_y)^2) / (sum((prog$ebp_any - prog$any_y)^2) + sum((prog$any - prog$ebp_any)^2)))
(sum((egg$ebp_any - egg$any_y)^2) / (sum((egg$ebp_any - egg$any_y)^2) + sum((egg$any - egg$ebp_any)^2)))
(sum((bee$ebp_any - bee$any_y)^2) / (sum((bee$ebp_any - bee$any_y)^2) + sum((bee$any - bee$ebp_any)^2)))
(sum((car$ebp_any - car$any_y)^2) / (sum((car$ebp_any - car$any_y)^2) + sum((car$any - car$ebp_any)^2)))

(sum((oral$ebp_via - oral$viable_y)^2) / (sum((oral$ebp_via - oral$viable_y)^2) + sum((oral$viable - oral$ebp_via)^2)))
(sum((prog$ebp_via - prog$viable_y)^2) / (sum((prog$ebp_via - prog$viable_y)^2) + sum((prog$viable - prog$ebp_via)^2)))
(sum((egg$ebp_via - egg$viable_y)^2) / (sum((egg$ebp_via - egg$viable_y)^2) + sum((egg$viable - egg$ebp_via)^2)))
(sum((bee$ebp_via - bee$viable_y)^2) / (sum((bee$ebp_via - bee$viable_y)^2) + sum((bee$viable - bee$ebp_via)^2)))
(sum((car$ebp_via - car$viable_y)^2) / (sum((car$ebp_via - car$viable_y)^2) + sum((car$viable - car$ebp_via)^2)))

(sum((oral$ebp_brain - oral$brain_y)^2) / (sum((oral$ebp_brain - oral$brain_y)^2) + sum((oral$brain - oral$ebp_brain)^2)))
(sum((prog$ebp_brain - prog$brain_y)^2) / (sum((prog$ebp_brain - prog$brain_y)^2) + sum((prog$brain - prog$ebp_brain)^2)))
(sum((egg$ebp_brain - egg$brain_y)^2) / (sum((egg$ebp_brain - egg$brain_y)^2) + sum((egg$brain - egg$ebp_brain)^2)))
(sum((bee$ebp_brain - bee$brain_y)^2) / (sum((bee$ebp_brain - bee$brain_y)^2) + sum((bee$brain - bee$ebp_brain)^2)))
(sum((car$ebp_brain - car$brain_y)^2) / (sum((car$ebp_brain - car$brain_y)^2) + sum((car$brain - car$ebp_brain)^2)))
