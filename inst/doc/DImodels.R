## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(DImodels)

## ----eval = FALSE-------------------------------------------------------------
# install.packages("DImodels")
# library("DImodels")

## ----eval = FALSE-------------------------------------------------------------
# ?DImodels

## ----eval = FALSE-------------------------------------------------------------
# ?sim3

## ----echo = FALSE, results='asis'---------------------------------------------
library(DImodels)
data("design_a")
knitr::kable(head(design_a))

## ----echo = TRUE, results='asis'----------------------------------------------
data("sim3")
knitr::kable(head(sim3, 10))

## ----echo = TRUE--------------------------------------------------------------
hist(sim3$response, xlab = "Response", main = "")
# Similar graphs can also be generated for the other species proportions.
plot(sim3$p1, sim3$response, xlab = "Proportion of species 1", ylab = "Response")
summary(sim3$response)

## ----echo = TRUE--------------------------------------------------------------
auto1 <- autoDI(y = "response", prop = 4:12, treat = "treatment", 
                FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"), data = sim3, 
                selection = "Ftest")

## ----eval = FALSE-------------------------------------------------------------
# ?autoDI

## ----echo = TRUE--------------------------------------------------------------
summary(auto1)

## ----eval = FALSE, echo = TRUE------------------------------------------------
# theta_CI(auto1, conf = .95)

## ----echo = TRUE--------------------------------------------------------------
m1 <- DI(y = "response", prop = 4:12, 
         FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"), treat = "treatment", 
         DImodel = "FG", data = sim3)
summary(m1)

## ----echo = TRUE--------------------------------------------------------------
m1_theta <- update_DI(object = m1, estimate_theta = TRUE)
coef(m1_theta)

## ----echo = TRUE--------------------------------------------------------------
m1_group <- update_DI(object = m1_theta, 
                      ID = c("ID1", "ID1", "ID1", "ID1", "ID1",
                             "ID1", "ID1", "ID1", "ID1"))
coef(m1_group)

## ----echo = TRUE--------------------------------------------------------------
m1_group2 <- update_DI(object = m1_theta, 
                       ID = c("ID1", "ID1", "ID1", 
                              "ID2", "ID2", "ID2", 
                              "ID3", "ID3", "ID3"))
coef(m1_group2)

## ----eval = FALSE-------------------------------------------------------------
# ?DI
# ?autoDI

## ----echo = TRUE--------------------------------------------------------------
m2 <- DI(y = "response", prop = 4:12, 
         FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"), treat = "treatment", 
         DImodel = "FG", extra_formula = ~ (p1 + p2 + p3 + p4):treatment,
         data = sim3)
summary(m2)

## ----echo = TRUE--------------------------------------------------------------
FG_matrix <- DI_data(prop = 4:12, FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"), 
                     data = sim3, what = "FG")
sim3a <- data.frame(sim3, FG_matrix)

## ----echo = TRUE--------------------------------------------------------------
m3 <- DI(y = "response", prop = 4:12, 
         FG = c("FG1","FG1","FG1","FG1","FG1","FG2","FG2","FG3","FG3"),
         treat = "treatment", DImodel = "FG", 
         extra_formula = ~ (bfg_FG1_FG2 + bfg_FG1_FG3 + bfg_FG2_FG3 +
                              wfg_FG1 + wfg_FG2 + wfg_FG3) : treatment, data = sim3a)
summary(m3)

## ----echo = TRUE--------------------------------------------------------------
sim3a$treatmentA <- as.numeric(sim3a$treatment == "A")

## ----echo = TRUE--------------------------------------------------------------
m3 <- DI(y = "response",
         custom_formula = response ~ 0 + p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 +
           treatmentA + bfg_FG1_FG2 + bfg_FG1_FG3 + bfg_FG2_FG3, data = sim3a)
summary(m3)

## ----echo = TRUE--------------------------------------------------------------
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

## ----echo = TRUE--------------------------------------------------------------
contr <- data.frame(p1 = c(1,  0),
                    p2 = c(-1, 0),
                    p7 = c(0,  1),
                    p9 = c(0, -1))
rownames(contr) <- c("p1_vs_p2", "p7_vs_p9")
  
the_C <- contrasts_DI(m3, contrast_vars = contr)
summary(the_C)

## ----echo = TRUE--------------------------------------------------------------
contr <- data.frame("treatmentA" = 1)
rownames(contr) <- "p1_TreatmentAvsB"
the_C <- contrasts_DI(m3, contrast_vars = contr)
summary(the_C)

## ----echo = TRUE--------------------------------------------------------------
# Suppose these are species communities we wish to compare
mixA <- c(0.25, 0,      0.25, 0,      0.25, 0,      0.25, 0, 0)
mixB <- c(0,    0.3333, 0,    0.3333, 0,    0.3333, 0,    0, 0)

# The contrast can be created by subtracting the species proportions
contr <- matrix(mixA - mixB, nrow = 1)
colnames(contr) <- paste0("p", 1:9)
print(contr)

# The values for the interaction terms will be calculated 
# automatically 
the_C <- contrasts_DI(m3, contrast_vars = contr)
summary(the_C)

## ----echo = TRUE--------------------------------------------------------------
# p1, p2, p5 and p7 equi-proportional mixture at treatment A
mixA <- contrast_matrix(object = m3, 
                        contrast_vars = data.frame("p1" = 0.25, 
                                                   "p2" = 0.25,                                                          "p5" = 0.25,
                                                   "p7" = 0.25,
                                                   "treatmentA" = 1))
# p2, p4, and p6 equi-proportional mixture at treatment A
mixB <- contrast_matrix(object = m3, 
                        contrast_vars = data.frame("p2" = 1/3, 
                                                   "p4" = 1/3,                                                           "p6" = 1/3,
                                                   "treatmentA" = 1))
# Subtracting these two values would give us the contrast for 
# comparing these mixtures
my_contrast <- mixA - mixB
rownames(my_contrast) <- "4_sp_mix vs 3_sp_mix"

# This contrast can be passed to the `contrast` parameter in `contrasts_DI`
the_C <- contrasts_DI(m3, contrast = my_contrast)
summary(the_C)

## -----------------------------------------------------------------------------
comms_to_compare <- sim3[c(36, 79, 352), 3:12]
compare_communities(m3, data = comms_to_compare)

## -----------------------------------------------------------------------------
# ref accepts either a row-number of rowname of the community to be used as reference
# adjust should be one of the following
# c("single-step", "Shaffer", "Westfall", "free", "holm",
#   "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none") 
compare_communities(m3, data = comms_to_compare,
                    ref = 3, adjust = "bonferroni")

