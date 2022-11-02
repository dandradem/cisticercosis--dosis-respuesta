###############################################################################
#######################NEGATIVE BINOMIAL REGRESSION############################
###############################################################################

#Clean the R environment
rm(list=ls())

#Establish the Working Directory
this.dir<-dirname(rstudioapi::getSourceEditorContext()$path)   
setwd(this.dir)

#Require necessary packages
if("MASS" %in% rownames(installed.packages()) == FALSE) +
  {install.packages("MASS"); require(MASS)} else{require(MASS)}
if("ggplot2" %in% rownames(installed.packages()) == FALSE) +
  {install.packages("ggplot2"); require(ggplot2)} else{require(ggplot2)}
if("scales" %in% rownames(installed.packages()) == FALSE) +
  {install.packages("scales"); require(scales)} else{require(scales)}

#Source the dataset script 
source('dataset.R')

db$exp_route_b <- ifelse(db$exp_route%in%c("proglottids","eggs","beetles"), "oral", "carotid")

#############
#TOTAL CYSTS#
#############

#PARAMETERS ESTIMATION

##Oral
total_oral <- glm.nb(total_cyst ~ dose_a, data = db[db$exp_route_b=="oral",])
summary(total_oral)

##Proglottids
total_prog <- glm.nb(total_cyst ~ dose_a, data = db[db$exp_route=="proglottids",])
summary(total_prog)

##Eggs
total_egg <- glm.nb(total_cyst ~ dose_a, data = db[db$exp_route=="eggs",])
summary(total_egg)

##Beetles
total_bee <- glm.nb(total_cyst ~ dose_a, data = db[db$exp_route=="beetles",])
summary(total_bee)

##Carotid
total_car <- glm.nb(total_cyst ~ dose_a, data = db[db$exp_route=="carotid",])
summary(total_car)


#PLOT DESIGN

##PROGLOTTIDS

###Extract the data for graphing the curve
curve.total.prog <- data.frame(
  dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "proglottids"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.total.prog <- cbind(curve.total.prog, predict(total_prog, curve.total.prog, type = "link",
                                                    se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.total.prog <- within(curve.total.prog, {
  total_cyst <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})


##EGGS

###Extract the data for graphing the curve
curve.total.egg <- data.frame(
  dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "eggs"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.total.egg <- cbind(curve.total.egg, predict(total_egg, curve.total.egg, type = "link",
                                                  se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.total.egg <- within(curve.total.egg, {
  total_cyst <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})


##BEETLES

###Extract the data for graphing the curve
curve.total.bee <- data.frame(
  dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "beetles"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.total.bee <- cbind(curve.total.bee, predict(total_bee, curve.total.bee, type = "link",
                                                  se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.total.bee <- within(curve.total.bee, {
  total_cyst <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})


##CAROTID

###Extract the data for graphing the curve
curve.total.car <- data.frame(
  dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "carotid"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.total.car <- cbind(curve.total.car, predict(total_car, curve.total.car, type = "link",
                                                  se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.total.car <- within(curve.total.car, {
  total_cyst <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})


##MERGING THE DATA OF THE ROUTES OF EXPOSURE
curve.total <- rbind(curve.total.prog, curve.total.egg, curve.total.bee, curve.total.car)
curve.total$exp_route <- factor(curve.total$exp_route, 
                                levels = c("proglottids","eggs","beetles","carotid"))


#PLOT
exp_route.labs <- c("PROGLOTIS", "HUEVOS", "ESCARABAJOS","CARÓTIDA")
names(exp_route.labs) <- c("proglottids", "eggs", "beetles","carotid")

total_labels <- data.frame(exp_route = c("proglottids","eggs","beetles","carotid"),
                           label = c("A","B","C","D"))
total_labels$exp_route <- factor(total_labels$exp_route,
                                 levels = c("proglottids","eggs","beetles","carotid"))

total.plot <- ggplot(curve.total, aes(x = dose_a, y = total_cyst)) +
  geom_line(aes(y = LL), size = 0.5, linetype = 5) +
  geom_line(aes(y = UL), size = 0.5, linetype = 5) +  
  geom_line(size = 1)+
  coord_cartesian(xlim = c(1e-02,1e06), ylim = c(0,7000)) +
  xlab("Dosis (número de huevos)") +
  ylab("Número de quistes desarrollados") +
  scale_x_log10(breaks = c(1e-02, 1e00, 1e02, 1e04, 1e06), 
                label = trans_format('log10',math_format(10^.x)))+
  geom_point(data=db, size = 2) +
  facet_wrap(~exp_route,
             labeller = labeller(exp_route = exp_route.labs),
             nrow = 4) +
  theme(strip.text.x = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
  geom_text(x = -2, y = 7000, aes(label = label), data = total_labels)

total.plot



##############
#VIABLE CYSTS#
##############

#PARAMETERS ESTIMATION

##Oral
ves_oral <- glm.nb(ves_cyst ~ dose_a, data = db[db$exp_route_b=="oral",])
summary(ves_oral)

##Proglottids
ves_prog <- glm.nb(ves_cyst ~ dose_a, data = db[db$exp_route=="proglottids",])
summary(ves_prog)

##Eggs
ves_egg <- glm.nb(ves_cyst ~ dose_a, data = db[db$exp_route=="eggs",])
summary(ves_egg)

##Beetles
ves_bee <- glm.nb(ves_cyst ~ dose_a, data = db[db$exp_route=="beetles",])
summary(ves_bee)

##Carotid
ves_car <- glm.nb(ves_cyst ~ dose_a, data = db[db$exp_route=="carotid",])
summary(ves_car)


#PLOT DESIGN

##PROGLOTTIDS

###Extract the data for graphing the curve
curve.ves.prog <- data.frame(
  dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "proglottids"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.ves.prog <- cbind(curve.ves.prog, predict(ves_prog, curve.ves.prog, type = "link",
                                                    se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.ves.prog <- within(curve.ves.prog, {
  ves_cyst <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})


##EGGS

###Extract the data for graphing the curve
curve.ves.egg <- data.frame(
  dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "eggs"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.ves.egg <- cbind(curve.ves.egg, predict(ves_egg, curve.ves.egg, type = "link",
                                                  se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.ves.egg <- within(curve.ves.egg, {
  ves_cyst <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})


##BEETLES

###Extract the data for graphing the curve
curve.ves.bee <- data.frame(
  dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "beetles"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.ves.bee <- cbind(curve.ves.bee, predict(ves_bee, curve.ves.bee, type = "link",
                                                  se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.ves.bee <- within(curve.ves.bee, {
  ves_cyst <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})


##CAROTID

###Extract the data for graphing the curve
curve.ves.car <- data.frame(
  dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "carotid"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.ves.car <- cbind(curve.ves.car, predict(ves_car, curve.ves.car, type = "link",
                                                  se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.ves.car <- within(curve.ves.car, {
  ves_cyst <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})


##MERGING THE DATA OF THE ROUTES OF EXPOSURE
curve.ves <- rbind(curve.ves.prog, curve.ves.egg, curve.ves.bee, curve.ves.car)
curve.ves$exp_route <- factor(curve.ves$exp_route, 
                                levels = c("proglottids","eggs","beetles","carotid"))


#PLOT
exp_route.labs <- c("PROGLOTIS", "HUEVOS", "ESCARABAJOS","CARÓTIDA")
names(exp_route.labs) <- c("proglottids", "eggs", "beetles","carotid")

ves_labels <- data.frame(exp_route = c("proglottids","eggs","beetles","carotid"),
                           label = c("A","B","C","D"))
ves_labels$exp_route <- factor(ves_labels$exp_route,
                                 levels = c("proglottids","eggs","beetles","carotid"))

ves.plot <- ggplot(curve.ves, aes(x = dose_a, y = ves_cyst)) +
  geom_line(aes(y = LL), size = 0.5, linetype = 5) +
  geom_line(aes(y = UL), size = 0.5, linetype = 5) +  
  geom_line(size = 1)+
  coord_cartesian(xlim = c(1e-02,1e06), ylim = c(0,7000)) +
  xlab("Dosis (número de huevos)") +
  ylab("Número de quistes desarrollados") +
  scale_x_log10(breaks = c(1e-02, 1e00, 1e02, 1e04, 1e06), 
                label = trans_format('log10',math_format(10^.x)))+
  geom_point(data=db, size = 2) +
  facet_wrap(~exp_route,
             labeller = labeller(exp_route = exp_route.labs),
             nrow = 4) +
  theme(strip.text.x = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
  geom_text(x = -2, y = 7000, aes(label = label), data = ves_labels)

ves.plot



#############
#BRAIN CYSTS#
#############

#PARAMETERS ESTIMATION

##Oral
brain_oral <- glm.nb(brain_cyst ~ dose_b, data = db[db$exp_route_b=="oral",])
summary(brain_oral)

##Proglottids
brain_prog <- glm.nb(brain_cyst ~ dose_b, data = db[db$exp_route=="proglottids",])
summary(brain_prog)

##Eggs
brain_egg <- glm.nb(brain_cyst ~ dose_b, data = db[db$exp_route=="eggs",])
summary(brain_egg)

##Beetles
brain_bee <- glm.nb(brain_cyst ~ dose_b, data = db[db$exp_route=="beetles",])
summary(brain_bee)

##Carotid
brain_car <- glm.nb(brain_cyst ~ dose_b, data = db[db$exp_route=="carotid",])
summary(brain_car)


#PLOT DESIGN

##PROGLOTTIDS

###Extract the data for graphing the curve
curve.brain.prog <- data.frame(
  dose_b = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "proglottids"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.brain.prog <- cbind(curve.brain.prog, predict(brain_prog, curve.brain.prog, type = "link",
                                                se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.brain.prog <- within(curve.brain.prog, {
  brain_cyst <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})


##EGGS

###Extract the data for graphing the curve
curve.brain.egg <- data.frame(
  dose_b = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "eggs"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.brain.egg <- cbind(curve.brain.egg, predict(brain_egg, curve.brain.egg, type = "link",
                                              se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.brain.egg <- within(curve.brain.egg, {
  brain_cyst <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})


##BEETLES

###Extract the data for graphing the curve
curve.brain.bee <- data.frame(
  dose_b = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "beetles"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.brain.bee <- cbind(curve.brain.bee, predict(brain_bee, curve.brain.bee, type = "link",
                                              se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.brain.bee <- within(curve.brain.bee, {
  brain_cyst <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})


##CAROTID

###Extract the data for graphing the curve
curve.brain.car <- data.frame(
  dose_b = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  exp_route = "carotid"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.brain.car <- cbind(curve.brain.car, predict(brain_car, curve.brain.car, type = "link",
                                              se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.brain.car <- within(curve.brain.car, {
  brain_cyst <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})


##MERGING THE DATA OF THE ROUTES OF EXPOSURE
curve.brain <- rbind(curve.brain.prog, curve.brain.egg, curve.brain.bee, curve.brain.car)
curve.brain$exp_route <- factor(curve.brain$exp_route, 
                              levels = c("proglottids","eggs","beetles","carotid"))


#PLOT
exp_route.labs <- c("PROGLOTIS", "HUEVOS", "ESCARABAJOS","CARÓTIDA")
names(exp_route.labs) <- c("proglottids", "eggs", "beetles","carotid")

brain_labels <- data.frame(exp_route = c("proglottids","eggs","beetles","carotid"),
                         label = c("A","B","C","D"))
brain_labels$exp_route <- factor(brain_labels$exp_route,
                               levels = c("proglottids","eggs","beetles","carotid"))

brain.plot <- ggplot(curve.brain, aes(x = dose_b, y = brain_cyst)) +
  geom_line(aes(y = LL), size = 0.5, linetype = 5) +
  geom_line(aes(y = UL), size = 0.5, linetype = 5) +  
  geom_line(size = 1)+
  coord_cartesian(xlim = c(1e-02,1e06), ylim = c(0,7000)) +
  xlab("Dosis (número de huevos)") +
  ylab("Número de quistes desarrollados") +
  scale_x_log10(breaks = c(1e-02, 1e00, 1e02, 1e04, 1e06), 
                label = trans_format('log10',math_format(10^.x)))+
  geom_point(data=db, size = 2) +
  facet_wrap(~exp_route,
             labeller = labeller(exp_route = exp_route.labs),
             nrow = 4) +
  theme(strip.text.x = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
  geom_text(x = -2, y = 7000, aes(label = label), data = brain_labels)

brain.plot


#ORAL PLOT

##Total cysts

###Extract the data for graphing the curve
curve.total.oral <- data.frame(
  dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  cyst_type = "total"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.total.oral <- cbind(curve.total.oral, predict(total_oral, curve.total.oral, type = "link",
                                                    se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.total.oral <- within(curve.total.oral, {
  cysts <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

curve.total.oral <- curve.total.oral %>%
  rename(
    dose = dose_a
  )


##Viable cysts

###Extract the data for graphing the curve
curve.ves.oral <- data.frame(
  dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  cyst_type = "vesicular"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.ves.oral <- cbind(curve.ves.oral, predict(ves_oral, curve.ves.oral, type = "link",
                                                    se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.ves.oral <- within(curve.ves.oral, {
  cysts <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

curve.ves.oral <- curve.ves.oral %>%
  rename(
    dose = dose_a
  )


##Brain cysts

###Extract the data for graphing the curve
curve.brain.oral <- data.frame(
  dose_b = exp(seq(log(1.00e06), log(1.00e-02), length=10000)),
  cyst_type = "brain"
)

###Predict the regression value by the dose and route of exposure
###with the created model
curve.brain.oral <- cbind(curve.brain.oral, predict(brain_oral, curve.brain.oral, type = "link",
                                                se.fit = TRUE))

###Calculate the probability of infection, the lower and upper limits 
##with the predicted values
curve.brain.oral <- within(curve.brain.oral, {
  cysts <- exp(fit)
  LL <- exp(fit - 1.96 * se.fit)
  UL <- exp(fit + 1.96 * se.fit)
})

curve.brain.oral <- curve.brain.oral %>%
  rename(
    dose = dose_b
  )

##MERGING THE DATA OF THE TYPES OF CYST
curve.oral <- rbind(curve.total.oral, curve.ves.oral, curve.brain.oral)
curve.oral$cyst_type <- factor(curve.oral$cyst_type,
                               levels = c("total","vesicular","brain"))

##EXTRACTING CYSTS DATA FOR DOT PLOT
geom_oral_total <- db[db$exp_route_b=="oral",c(1,2,6,10)]
geom_oral_total$cyst_type <- "total"
geom_oral_total <- geom_oral_total %>%
  rename(
    dose = dose_a,
    cysts = total_cyst
  )

geom_oral_ves <- db[db$exp_route_b=="oral",c(1,2,4,10)]
geom_oral_ves$cyst_type <- "vesicular"
geom_oral_ves <- geom_oral_ves %>%
  rename(
    dose = dose_a,
    cysts = ves_cyst
  )

geom_oral_brain <- db[db$exp_route_b=="oral",c(1,7,9,10)]
geom_oral_brain$cyst_type <- "brain"
geom_oral_brain <- geom_oral_brain %>%
  rename(
    dose = dose_b,
    cysts = brain_cyst
  )

geom_oral_cyst <- rbind(geom_oral_total, geom_oral_ves, geom_oral_brain)
geom_oral_cyst <- geom_oral_cyst[geom_oral_cyst$dose>0,]
geom_oral_cyst$cyst_type <- factor(geom_oral_cyst$cyst_type,
                                   levels = c("total","vesicular","brain"))

##PLOT
cyst_type.labs <- c("QUISTES EN GENERAL", "QUISTES VIABLES", "QUISTES CEREBRALES")
names(cyst_type.labs) <- c("total", "vesicular", "brain")

oral_labels <- data.frame(cyst_type = c("total","vesicular","brain"),
                         label = c("A","B","C"))
oral_labels$cyst_type <- factor(oral_labels$cyst_type,
                               levels = c("total","vesicular","brain"))

oral_plot <- ggplot(curve.oral, aes(x = dose, y = cysts)) +
  geom_line(aes(y = LL), size = 0.5, linetype = 5) +
  geom_line(aes(y = UL), size = 0.5, linetype = 5) +  
  geom_line(size = 1)+
  coord_cartesian(xlim = c(1e-02,1e06), ylim = c(0,7000)) +
  xlab("Dosis (número de huevos)") +
  ylab("Número de quistes desarrollados") +
  scale_x_log10(breaks = c(1e-02, 1e00, 1e02, 1e04, 1e06), 
                label = trans_format('log10',math_format(10^.x))) +
  geom_point(aes(shape = exp_route), data=geom_oral_cyst, size = 2) +
  facet_wrap(~cyst_type,
             labeller = labeller(cyst_type = cyst_type.labs),
             nrow = 4) +
  theme(strip.text.x = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20)) +
  geom_text(x = -2, y = 7000, aes(label = label), data = oral_labels)

oral_plot
