###############################################################################
#########################TWO-PARAMETER LOG-LOGISTIC############################
###############################################################################

#Clean the R environment
rm(list=ls())

#Establish the Working Directory
this.dir<-dirname(rstudioapi::getSourceEditorContext()$path)   
setwd(this.dir)

#Require necessary packages
if("rlang" %in% rownames(installed.packages()) == FALSE) +
  {install.packages("rlang"); require(rlang)} else{require(rlang)}
if("drc" %in% rownames(installed.packages()) == FALSE) +
  {install.packages("drc"); require(drc)} else{require(drc)}
if("ggplot2" %in% rownames(installed.packages()) == FALSE) +
  {install.packages("ggplot2"); require(ggplot2)} else{require(ggplot2)}
if("scales" %in% rownames(installed.packages()) == FALSE) +
  {install.packages("scales"); require(scales)} else{require(scales)}

#Source the dataset script 
source('dataset_15nc.R')



#############
#TOTAL CYSTS#
#############

#PARAMETERS ESTIMATION

##Proglottids
total.prog <- drm(nInf/nExp~dose_a, 
                  data = db_exp_total[db_exp_total$exp_route=="proglottids",],
                  fct = LL.2())
summary(total.prog)

##Eggs
total.egg <- drm(nInf/nExp~dose_a, 
                 data = db_exp_total[db_exp_total$exp_route=="eggs",],
                 fct = LL.2())
summary(total.egg)

##Beetles
total.bee <- drm(nInf/nExp~dose_a, 
                 data = db_exp_total[db_exp_total$exp_route=="beetles",],
                 fct = LL.2())
summary(total.bee)

##Carotid
total.car <- drm(nInf/nExp~dose_a, 
                 data = db_exp_total[db_exp_total$exp_route=="carotid",],
                 fct = LL.2())
summary(total.car)

#Extract the calculated mean infective dose (ID50)
ID50.total.prog <- total.prog$coefficients[2]
ID50.total.egg <- total.egg$coefficients[2]

ID50.total <- data.frame(exp_route = c("proglottids","eggs","beetles","carotid"))
ID50.total$exp_route <- factor(ID50.total$exp_route, 
                               levels = c("proglottids", "eggs", "beetles", "carotid"))

##Create a data.frame to storage the calculated ID50
attach(ID50.total)
ID50.total$ID50 <- ifelse(exp_route=="proglottids",ID50.total.prog,
                          ifelse(exp_route=="eggs", ID50.total.egg, NA))
detach(ID50.total)


#PLOT DESIGN

#Extract the data for graphing the curve

##Proglottids
curve.total.prog <- data.frame(dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.total.prog$exp_route <- "proglottids"

##Predict the probability of infection, upper limit and by the dose and route of exposure
##with the created model
pred.total.prog <- predict(total.prog, newdata=curve.total.prog, interval="confidence") 

curve.total.prog$p <- pred.total.prog[,1]
curve.total.prog$LL <- pred.total.prog[,2]
curve.total.prog$UL <- pred.total.prog[,3]

##Eggs
curve.total.egg <- data.frame(dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.total.egg$exp_route <- "eggs"

##Predict the probability of infection, upper limit and by the dose and route of exposure
##with the created model
pred.total.egg <- predict(total.egg, newdata=curve.total.egg, interval="confidence") 

curve.total.egg$p <- pred.total.egg[,1]
curve.total.egg$LL <- pred.total.egg[,2]
curve.total.egg$UL <- pred.total.egg[,3]

##Beetles
curve.total.bee <- data.frame(dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.total.bee$exp_route <- "beetles"

curve.total.bee$p <- NA
curve.total.bee$LL <- NA
curve.total.bee$UL <- NA

##Carotid
curve.total.car <- data.frame(dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.total.car$exp_route <- "carotid"

curve.total.car$p <- NA
curve.total.car$LL <- NA
curve.total.car$UL <- NA

#Merge the predicted values from each route of exposure
curve.total <- rbind(curve.total.prog,curve.total.car,curve.total.egg,
                     curve.total.bee)
curve.total$exp_route <- factor(curve.total$exp_route,
                                levels = c("proglottids", "eggs", "beetles", "carotid"))

###PLOT 1
exp_route.labs <- c("PROGLOTTIDS", "EGGS", "BEETLES","CAROTID")
names(exp_route.labs) <- c("proglottids", "eggs", "beetles","carotid")

letter_labels <- data.frame(exp_route = c("proglottids","eggs","beetles","carotid"),
                            label = c("A","B","C","D"))
letter_labels$exp_route <- factor(letter_labels$exp_route, 
                                  levels = c("proglottids","eggs","beetles","carotid"))

total.plot <- ggplot(curve.total, aes(x = dose_a, y = p)) +
  geom_ribbon(aes(ymin=LL, ymax=UL), alpha=0.2) +
  geom_vline(data = ID50.total, aes(xintercept = ID50), linetype = 5) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(1e-02,1e06), ylim = c(0,1)) +
  xlab("Dose (number of eggs)") +
  ylab("P(inf)") +
  scale_y_continuous(labels = percent) +
  scale_x_log10(breaks = c(1e-02, 1e00, 1e02, 1e04, 1e06), 
                label = trans_format('log10',math_format(10^.x))) +
  geom_point(data=geom_total, aes(x = dose_a, y = nInf/nExp), size = 2) +
  geom_errorbarh(data=geom_total, 
                 aes(x = dose_a, y = nInf/nExp, xmax = dose_a + dose_sd_a, 
                     xmin = dose_a - dose_sd_a, height = .05)) +
  facet_wrap(~exp_route,
             labeller = labeller(exp_route = exp_route.labs),
             nrow = 4) +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  geom_text(x = -2, y = 1, aes(label = label), data = letter_labels)

total.plot



#################
#VESICULAR CYSTS#
#################

#PARAMETERS ESTIMATION

##Proglottids
ves.prog <- drm(nInf/nExp~dose_a, 
                  data = db_exp_ves[db_exp_ves$exp_route=="proglottids",],
                  fct = LL.2())
summary(ves.prog)

##Eggs
ves.egg <- drm(nInf/nExp~dose_a, 
                 data = db_exp_ves[db_exp_ves$exp_route=="eggs",],
                 fct = LL.2())
summary(ves.egg)

##Beetles
ves.bee <- drm(nInf/nExp~dose_a, 
                 data = db_exp_ves[db_exp_ves$exp_route=="beetles",],
                 fct = LL.2())
summary(ves.bee)

##Carotid
ves.car <- drm(nInf/nExp~dose_a, 
                 data = db_exp_ves[db_exp_ves$exp_route=="carotid",],
                 fct = LL.2())
summary(ves.car)

#Extract the calculated mean infective dose (ID50)
ID50.ves.prog <- ves.prog$coefficients[2]
ID50.ves.egg <- ves.egg$coefficients[2]

ID50.ves <- data.frame(exp_route = c("proglottids","eggs","beetles","carotid"))
ID50.ves$exp_route <- factor(ID50.ves$exp_route, 
                               levels = c("proglottids", "eggs", "beetles", "carotid"))

##Create a data.frame to storage the calculated ID50
attach(ID50.ves)
ID50.ves$ID50 <- ifelse(exp_route=="proglottids", ID50.ves.prog, 
                        ifelse(exp_route=="eggs", ID50.ves.egg, NA))
detach(ID50.ves)


#PLOT DESIGN

#Extract the data for graphing the curve

##Proglottids
curve.ves.prog <- data.frame(dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.ves.prog$exp_route <- "proglottids"

##Predict the probability of infection, upper limit and by the dose and route of exposure
##with the created model
pred.ves.prog <- predict(ves.prog, newdata=curve.ves.prog, interval="confidence") 

curve.ves.prog$p <- pred.ves.prog[,1]
curve.ves.prog$LL <- pred.ves.prog[,2]
curve.ves.prog$UL <- pred.ves.prog[,3]

##Eggs
curve.ves.egg <- data.frame(dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.ves.egg$exp_route <- "eggs"

##Predict the probability of infection, upper limit and by the dose and route of exposure
##with the created model
pred.ves.egg <- predict(ves.egg, newdata=curve.ves.egg, interval="confidence") 

curve.ves.egg$p <- pred.ves.egg[,1]
curve.ves.egg$LL <- pred.ves.egg[,2]
curve.ves.egg$UL <- pred.ves.egg[,3]

##Beetles
curve.ves.bee <- data.frame(dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.ves.bee$exp_route <- "beetles"

curve.ves.bee$p <- NA
curve.ves.bee$LL <- NA
curve.ves.bee$UL <- NA

##Carotid
curve.ves.car <- data.frame(dose_a = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.ves.car$exp_route <- "carotid"

curve.ves.car$p <- NA
curve.ves.car$LL <- NA
curve.ves.car$UL <- NA

#Merge the predicted values from each route of exposure
curve.ves <- rbind(curve.ves.prog,curve.ves.car,curve.ves.egg,
                     curve.ves.bee)
curve.ves$exp_route <- factor(curve.ves$exp_route,
                                levels = c("proglottids", "eggs", "beetles", "carotid"))

###PLOT 1
exp_route.labs <- c("PROGLOTTIDS", "EGGS", "BEETLES","CAROTID")
names(exp_route.labs) <- c("proglottids", "eggs", "beetles","carotid")

ves.plot <- ggplot(curve.ves, aes(x = dose_a, y = p)) +
  geom_ribbon(aes(ymin=LL, ymax=UL), alpha=0.2) +
  geom_vline(data = ID50.ves, aes(xintercept = ID50), linetype = 5) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(1e-02,1e06), ylim = c(0,1)) +
  xlab("Dose (number of eggs)") +
  ylab("P(inf)") +
  scale_y_continuous(labels = percent) +
  scale_x_log10(breaks = c(1e-02, 1e00, 1e02, 1e04, 1e06), 
                label = trans_format('log10',math_format(10^.x))) +
  geom_point(data=geom_ves, aes(x = dose_a, y = nInf/nExp), size = 2) +
  geom_errorbarh(data=geom_ves, 
                 aes(x = dose_a, y = nInf/nExp, xmax = dose_a + dose_sd_a, 
                     xmin = dose_a - dose_sd_a, height = .05)) +
  facet_wrap(~exp_route,
             labeller = labeller(exp_route = exp_route.labs),
             nrow = 4) +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  geom_text(x = -2, y = 1, aes(label = label), data = letter_labels)

ves.plot



#############
#BRAIN CYSTS#
#############

#PARAMETERS ESTIMATION

##Proglottids
brain.prog <- drm(nInf/nExp~dose_b, 
                  data = db_exp_brain[db_exp_brain$exp_route=="proglottids",],
                  fct = LL.2())
summary(brain.prog)

##Eggs
brain.egg <- drm(nInf/nExp~dose_b, 
                 data = db_exp_brain[db_exp_brain$exp_route=="eggs",],
                 fct = LL.2())
summary(brain.egg)

##Beetles
brain.bee <- drm(nInf/nExp~dose_b, 
                 data = db_exp_brain[db_exp_brain$exp_route=="beetles",],
                 fct = LL.2())
summary(brain.bee)

##Carotid
brain.car <- drm(nInf/nExp~dose_b, 
                 data = db_exp_brain[db_exp_brain$exp_route=="carotid",],
                 fct = LL.2())
summary(brain.car)

#Extract the calculated mean infective dose (ID50)
ID50.brain.prog <- brain.prog$coefficients[2]
ID50.brain.egg <- brain.egg$coefficients[2]
ID50.brain.bee <- brain.bee$coefficients[2]
ID50.brain.car <- brain.car$coefficients[2]

ID50.brain <- data.frame(exp_route = c("proglottids","eggs","beetles","carotid"))
ID50.brain$exp_route <- factor(ID50.brain$exp_route, 
                               levels = c("proglottids", "eggs", "beetles", "carotid"))

##Create a data.frame to storage the calculated ID50
attach(ID50.brain)
ID50.brain$ID50 <- ifelse(exp_route=="proglottids",ID50.brain.prog,
                          ifelse(exp_route=="eggs", ID50.brain.egg, 
                                ifelse(exp_route=="beetles", ID50.brain.bee, ID50.brain.car)))
detach(ID50.brain)


#PLOT DESIGN

#Extract the data for graphing the curve

##Proglottids
curve.brain.prog <- data.frame(dose_b = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.brain.prog$exp_route <- "proglottids"

##Predict the probability of infection, upper limit and by the dose and route of exposure
##with the created model
pred.brain.prog <- predict(brain.prog, newdata=curve.brain.prog, interval="confidence") 

curve.brain.prog$p <- pred.brain.prog[,1]
curve.brain.prog$LL <- pred.brain.prog[,2]
curve.brain.prog$UL <- pred.brain.prog[,3]

##Eggs
curve.brain.egg <- data.frame(dose_b = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.brain.egg$exp_route <- "eggs"

##Predict the probability of infection, upper limit and by the dose and route of exposure
##with the created model
pred.brain.egg <- predict(brain.egg, newdata=curve.brain.egg, interval="confidence") 

curve.brain.egg$p <- pred.brain.egg[,1]
curve.brain.egg$LL <- pred.brain.egg[,2]
curve.brain.egg$UL <- pred.brain.egg[,3]

##Beetles
curve.brain.bee <- data.frame(dose_b = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.brain.bee$exp_route <- "beetles"

##Predict the probability of infection, upper limit and by the dose and route of exposure
##with the created model
pred.brain.bee <- predict(brain.bee, newdata=curve.brain.bee, interval="confidence") 

curve.brain.bee$p <- pred.brain.bee[,1]
curve.brain.bee$LL <- pred.brain.bee[,2]
curve.brain.bee$UL <- pred.brain.bee[,3]

##Carotid
curve.brain.car <- data.frame(dose_b = exp(seq(log(1.00e06), log(1.00e-02), length=10000)))
curve.brain.car$exp_route <- "carotid"

##Predict the probability of infection, upper limit and by the dose and route of exposure
##with the created model
pred.brain.car <- predict(brain.car, newdata=curve.brain.car, interval="confidence") 

curve.brain.car$p <- pred.brain.car[,1]
curve.brain.car$LL <- pred.brain.car[,2]
curve.brain.car$UL <- pred.brain.car[,3]

#Merge the predicted values from each route of exposure
curve.brain <- rbind(curve.brain.prog,curve.brain.car,curve.brain.egg,
                     curve.brain.bee)
curve.brain$exp_route <- factor(curve.brain$exp_route,
                                levels = c("proglottids", "eggs", "beetles", "carotid"))

###PLOT 1
exp_route.labs <- c("PROGLOTTIDS", "EGGS", "BEETLES","CAROTID")
names(exp_route.labs) <- c("proglottids", "eggs", "beetles","carotid")

brain.plot <- ggplot(curve.brain, aes(x = dose_b, y = p)) +
  geom_ribbon(aes(ymin=LL, ymax=UL), alpha=0.2) +
  geom_vline(data = ID50.brain, aes(xintercept = ID50), linetype = 5) +
  geom_line(size = 1) +
  coord_cartesian(xlim = c(1e-02,1e06), ylim = c(0,1)) +
  xlab("Dose (number of eggs)") +
  ylab("P(inf)") +
  scale_y_continuous(labels = percent) +
  scale_x_log10(breaks = c(1e-02, 1e00, 1e02, 1e04, 1e06), 
                label = trans_format('log10',math_format(10^.x))) +
  geom_point(data=geom_brain, aes(x = dose_b, y = nInf/nExp), size = 2) +
  geom_errorbarh(data=geom_brain, 
                 aes(x = dose_b, y = nInf/nExp, xmax = dose_b + dose_sd_b, 
                     xmin = dose_b - dose_sd_b, height = .05)) +
  facet_wrap(~exp_route,
             labeller = labeller(exp_route = exp_route.labs),
             nrow = 4) +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  geom_text(x = -2, y = 1, aes(label = label), data = letter_labels)

brain.plot

