## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(DImodels)

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("DImodels")
#  library("DImodels")

## ---- eval = FALSE------------------------------------------------------------
#  ?DImodels

## ---- eval = FALSE------------------------------------------------------------
#  ?sim3

## ---- echo = FALSE, results='asis'--------------------------------------------
library(DImodels)
data("design_a")
knitr::kable(head(design_a))

## ---- echo = TRUE, results='asis'---------------------------------------------
data("sim3")
knitr::kable(head(sim3, 10))

## ---- echo = TRUE-------------------------------------------------------------
hist(sim3$response, xlab = "Response", main = "")
# Similar graphs can also be generated for the other species proportions.
plot(sim3$p1, sim3$response, xlab = "Proportion of species 1", ylab = "Response")
summary(sim3$response)

## ---- echo = TRUE-------------------------------------------------------------
auto1 <- autoDI(y = "response", prop = 4:12, treat = "treatment", 
                FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"), data = sim3, 
                selection = "Ftest")

## ---- eval = FALSE------------------------------------------------------------
#  ?autoDI

## ---- echo = TRUE-------------------------------------------------------------
summary(auto1)

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  theta_CI(auto1, conf = .95)

## ---- echo = TRUE-------------------------------------------------------------
m1 <- DI(y = "response", prop = 4:12, 
         FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"), treat = "treatment", 
         DImodel = "FG", data = sim3)
summary(m1)

## ---- echo = TRUE-------------------------------------------------------------
m1_theta <- update_DI(object = m1, estimate_theta = TRUE)
coef(m1_theta)

## ---- echo = TRUE-------------------------------------------------------------
m1_group <- update_DI(object = m1_theta, 
                      ID = c("ID1", "ID1", "ID1", "ID1", "ID1",
                             "ID1", "ID1", "ID1", "ID1"))
coef(m1_group)

## ---- echo = TRUE-------------------------------------------------------------
m1_group2 <- update_DI(object = m1_theta, 
                       ID = c("ID1", "ID1", "ID1", 
                              "ID2", "ID2", "ID2", 
                              "ID3", "ID3", "ID3"))
coef(m1_group2)

## ---- eval = FALSE------------------------------------------------------------
#  ?DI
#  ?autoDI

## ---- echo = TRUE-------------------------------------------------------------
m2 <- DI(y = "response", prop = 4:12, 
         FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"), treat = "treatment", 
         DImodel = "FG", extra_formula = ~ (p1 + p2 + p3 + p4):treatment,
         data = sim3)
summary(m2)

## ---- echo = TRUE-------------------------------------------------------------
FG_matrix <- DI_data(prop = 4:12, FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"), 
                     data = sim3, what = "FG")
sim3a <- data.frame(sim3, FG_matrix)

## ---- echo = TRUE-------------------------------------------------------------
m3 <- DI(y = "response", prop = 4:12, 
         FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"),
         treat = "treatment", DImodel = "FG", 
         extra_formula = ~ (bfg_FG1_FG2 + bfg_FG1_FG3 + bfg_FG2_FG3 +
                              wfg_FG1 + wfg_FG2 + wfg_FG3) : treatment, data = sim3a)
summary(m3)

## ---- echo = TRUE-------------------------------------------------------------
sim3a$treatmentA <- as.numeric(sim3a$treatment == "A")

## ---- echo = TRUE-------------------------------------------------------------
m3 <- DI(y = "response",
         custom_formula = response ~ 0 + p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 +
           treatmentA + bfg_FG1_FG2 + bfg_FG1_FG3 + bfg_FG2_FG3, data = sim3a)
summary(m3)

## ---- echo = TRUE-------------------------------------------------------------
# Fit model
m3 <- DI(y = "response", prop = 4:12, 
         treat = "treatment", DImodel = "AV", 
         extra_formula = ~ (AV) : treatment, data = sim3a)

predict_data <- sim3[c(1, 79, 352), 3:12]
# Only species proportions and treatment is needed
print(predict_data)
# Make prediction
predict(m3, newdata = predict_data)

## -----------------------------------------------------------------------------
# The interval and level parameters can be used to calculate the 
# uncertainty around the predictions

# Get confidence interval around prediction
predict(m3, newdata = predict_data, interval = "confidence")

# Get prediction interval around prediction
predict(m3, newdata = predict_data, interval = "prediction")

# The function returns a 95% interval by default, 
# this can be changed using the level argument
predict(m3, newdata = predict_data, 
        interval = "prediction", level = 0.9)

## ---- echo = TRUE-------------------------------------------------------------
contr <- list("p1vsp2" = c(1, -1, 0, 0,  0,  0, 0, 0,  0, 0, 0, 0),
              "p3vsp5" = c(0,  0, 1, 0, -1,  0, 0, 0,  0, 0, 0, 0),
              "p4vsp6" = c(0,  0, 0, 1,  0, -1, 0, 0,  0, 0, 0, 0),
              "p7vsp9" = c(0,  0, 0, 0,  0,  0, 1, 0, -1, 0, 0, 0))
the_C <- contrasts_DI(m3, contrast = contr)
summary(the_C)

## ---- echo = TRUE-------------------------------------------------------------
contr <- list("treatAvsB" = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0))
the_C <- contrasts_DI(m3, contrast = contr)
summary(the_C)

## ---- echo = TRUE-------------------------------------------------------------
mixA <- c(0.25, 0,      0.25, 0,      0.25, 0,      0.25, 0, 0, 0, 0, 0)
mixB <- c(0,    0.3333, 0,    0.3333, 0,    0.3333, 0,    0, 0, 0, 0, 0)

# We have the proportions of the individual species in the mixtures, however
# we still need to calculate the interaction effect for these communities
contr_data <- data.frame(rbind(mixA, mixB))
colnames(contr_data) <- names(coef(m3))

# Adding the interaction effect of the two mixtures
contr_data$AV <- DI_data_E_AV(prop = 1:9, data = contr_data)$AV
print(contr_data)

# We can now subtract the respective values in each column of the two 
# mixtures and get our contrast
my_contrast <- as.matrix(contr_data[1, ] - contr_data[2, ])
rownames(my_contrast) <- "mixAvsB"

the_C <- contrasts_DI(m3, contrast = my_contrast)
summary(the_C)

