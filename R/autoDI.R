# y = response variable index or name
# block = block variable index or name
# density = density variable index or name
# prop = vector of indices or names of the species proportions
# treat = environmental covariate index or name
# data = dataset in data.frame format
# selection = "Ftest" to perform F (or X2) tests to select best model,
#             "AIC" to use select model with smallest AIC

autoDI <- function(y, block, density, prop, treat, FG = NULL, data,
                   selection = c("Ftest","AIC","AICc","BIC","BICc"),
                   step0 = FALSE, step4 = TRUE) {
  if(missing(y)) stop("You must supply a response variable name or column index through the argument 'y'.\n")
  # family / binomial denominator lock
  # include family and total as function arguments to lift the lock
  family <- "gaussian"
  total <- NULL
  # family lock
  if(!(family %in% c("gaussian","normal")))
    stop("As of version ", packageVersion("DImodels"),
         " DI models are implemented for family = 'gaussian' (= 'normal') only")
  # set theta flag
  estimate_theta <- TRUE
  # check if user chose "Ftest", "AIC", "AICc" or "BIC"
  selection <- match.arg(selection)
  # default family set as normal
  if(missing(family) || family == "normal") family <- "gaussian"
  # checks for binomial
  if(family %in% c("binomial","quasibinomial")) {
    if(missing(total)) {
      if(any(!(data[,y] %in% c(0,1)))) {
        stop("total must be informed for non-binary discrete proportion data")
      } else total <- 1
    } else total <- data[,total]
  } else total <- NULL
  # checks for quasi and information criteria
  if(family %in% c("quasipoisson","quasibinomial") & selection %in% c("AIC","AICc","BIC","BICc"))
    stop("cannot compute information criteria for quasi models, use selection = 'Ftest'")
  # flags if block/density and/or treat are missing
  if(missing(block)) block <- NA
  if(missing(density)) density <- NA
  treat_flag <- FALSE
  if(missing(treat)) {
    treat <- NA
    treat_flag <- TRUE
  }
  
  # preparing new data object
  data_obj <- DI_data_prepare(y = y, block = block, density = density, prop = prop, treat = treat, FG = FG, data = data)
  newdata <- data_obj$newdata
  nSpecies <- data_obj$nSpecies
  ## fitting the DI models
  all_models <- allmodels_autoDI(y = data_obj$y, block = data_obj$block, density = data_obj$density,
                                 prop = data_obj$prop, FG = data_obj$FG, treat = data_obj$treat, data = newdata,
                                 family = family, total = total, estimate_theta = estimate_theta,
                                 nSpecies = nSpecies, treat_flag = treat_flag,
                                 even_flag = data_obj$even_flag, P_int_flag = data_obj$P_int_flag,
                                 FGnames = FG)
  
  ## Step 0
  if(step0) {
    autoDI_step0(y = y, block = block, density = density, treat = treat, family = family, data = data)
  }
  
  ## first do selection
  if(treat_flag) {
    ## without treatment covariate -- single step: full selection
    message("\n", strrep("-", getOption("width")))
    message("Step 1: Investigating the diversity effect\n")
    if(!is.na(block) & is.na(density)) message("All models include block\n")
    if(is.na(block) & !is.na(density)) message("All models include density\n")
    if(!is.na(block) & !is.na(density)) message("All models include block and density\n")
    # loglikelihoods, information criteria and no. of parameters
    llik <- sapply(all_models$model_list, function(x) as.numeric(logLik(x)))
    mAIC <- sapply(all_models$model_list, AIC2)
    mAICc <- sapply(all_models$model_list, AICc)
    mBIC <- sapply(all_models$model_list, BIC2)
    mBICc <- sapply(all_models$model_list, BICc)
    ndf <- sapply(all_models$model_list, function(x) length(coef(x)))
    # model selection
    selected_model <- switch(selection,
                             Ftest = test_autoDI(model_list = all_models$model_list, family = family, treat = treat),
                             AIC = AICsel_autoDI(model_list = all_models$model_list, mAIC = mAIC, treat = treat),
                             AICc = AICcsel_autoDI(model_list = all_models$model_list, mAICc = mAICc, treat = treat),
                             BIC = BICsel_autoDI(model_list = all_models$model_list, mBIC = mBIC, treat = treat),
                             BICc = BICcsel_autoDI(model_list = all_models$model_list, mBICc = mBICc, treat = treat)
                             )
    if(is.null(FG)) {
      message("\nFunctional groups (argument 'FG') were not specified, and therefore not investigated.") 
    }
    # get model formula
    #fmla_char <- as.character(formula(all_models$model_list[[selected_model]]))[c(2,1,3)]
    # 'pseudo-algebraic' formula specification
    fmla_full <- paste(y, " = ",
                       paste0(names(coef(all_models$model_list[[selected_model]])), collapse = " + "))
    fmla_full <- gsub(":", "*", fmla_full)
    # print info
    message("\nSelected model: ",
        namesub_autoDI(selected_model), "\n",
        "Formula: ",
        sep = "")
    message(fmla_full, sep = " ")
    message("\n", strrep("-", getOption("width")))
    message("Step 2: No investigation of treatment effect included, since no treatment was specified
        (argument 'treat' omitted)")
  } else {
    ## with treatment covariate, two steps
    ## step 1: full model selection
    message("\n", strrep("-", getOption("width")))
    message("Step 1: Investigating the diversity effect\n")
    if(!is.na(block) & is.na(density)) message("All models include block and treatment\n")
    if(is.na(block) & !is.na(density)) message("All models include density and treatment\n")
    if(!is.na(block) & !is.na(density)) message("All models include block, density and treatment\n")
    # loglikelihoods, information criteria and no. of parameters
    llik <- sapply(all_models$model_list_treat, function(x) as.numeric(logLik(x)))
    mAIC <- sapply(all_models$model_list_treat, AIC2)
    mAICc <- sapply(all_models$model_list_treat, AICc)
    mBIC <- sapply(all_models$model_list_treat, BIC2)
    mBICc <- sapply(all_models$model_list_treat, BICc)
    ndf <- sapply(all_models$model_list_treat, function(x) length(coef(x)))
    # model selection
    selected_model <- switch(selection,
                             Ftest = test_autoDI(model_list = all_models$model_list_treat, family = family, treat = treat),
                             AIC = AICsel_autoDI(model_list = all_models$model_list_treat, mAIC = mAIC, treat = treat),
                             AICc = AICcsel_autoDI(model_list = all_models$model_list_treat, mAICc = mAICc, treat = treat),
                             BIC = BICsel_autoDI(model_list = all_models$model_list_treat, mBIC = mBIC, treat = treat),
                             BICc = BICcsel_autoDI(model_list = all_models$model_list_treat, mBICc = mBICc, treat = treat)
                             )
    if(is.null(FG)) {
      message("\nFunctional groups (argument 'FG') were not specified, and therefore not investigated.") 
    }
    # get model formula
    fmla_char <- as.character(formula(all_models$model_list_treat[[selected_model]]))[c(2,1,3)]
    fmla_full <- paste(y, " = ",
                       paste0(names(coef(all_models$model_list_treat[[selected_model]])), collapse = " + "))
    fmla_full <- gsub(":", "*", fmla_full)
    # print info
    message("\nSelected model: ",
        namesub_autoDI(selected_model), "\n",
        "Formula: ",
        sep = "")
    message(fmla_full, sep = " ")
    
    ## step 2 -- test the effect of treatment for the selected model
    message("\n", strrep("-", getOption("width")))
    message("Step 2: Investigating the treatment effect\n")
    string_id <- substr(selected_model, 1, 4)
    candidate_model_id <- grep(string_id, names(all_models$model_list))
    new_model_list <- list(all_models$model_list[[candidate_model_id]],
                           all_models$model_list_treat[[selected_model]])
    names(new_model_list) <- c(names(all_models$model_list)[candidate_model_id],
                               names(all_models$model_list_treat)[candidate_model_id])
    mAIC_new <- sapply(new_model_list, AIC2)
    mAICc_new <- sapply(new_model_list, AICc)
    mBIC_new <- sapply(new_model_list, BIC2)
    mBICc_new <- sapply(new_model_list, BICc)
    selected_model_new <- switch(selection,
                                 Ftest = test_autoDI(model_list = new_model_list, family = family, treat = treat),
                                 AIC = AICsel_autoDI(model_list = new_model_list, mAIC = mAIC_new, treat = treat),
                                 AICc = AICcsel_autoDI(model_list = new_model_list, mAICc = mAICc_new, treat = treat),
                                 BIC = BICsel_autoDI(model_list = new_model_list, mBIC = mBIC_new, treat = treat),
                                 BICc = BICcsel_autoDI(model_list = new_model_list, mBICc = mBICc_new, treat = treat)
                                 )
    fmla_char <- as.character(formula(new_model_list[[selected_model_new]]))[c(2,1,3)]
    fmla_full <- paste(y, " = ",
                       paste0(names(coef(new_model_list[[selected_model_new]])), collapse = " + "))
    fmla_full <- gsub(":", "*", fmla_full)
    message("\nSelected model: ",
        namesub_autoDI(selected_model_new), "\n",
        "Formula: ",
        sep = "")
    message(fmla_full, sep = " ")
  }

  ## finally, test for theta
  ## check first if selected model is either STR or ID
  if(treat_flag) {
    STR_flag <- length(grep("STR", selected_model) == 0)
    ID_flag <- length(grep("ID", selected_model) == 0)
  } else {
    STR_flag <- length(grep("STR", selected_model_new) == 0)
    ID_flag <- length(grep("ID", selected_model_new) == 0)
  }
  STR_ID_flag <- STR_flag | ID_flag
 
  if(!treat_flag) {
    select_from_treat <- grep("treat", selected_model_new)
  } else {
    select_from_treat <- numeric(0)
    selected_model_new <- selected_model
  }
  if(length(select_from_treat) == 0) {
    theta_models_list <- all_models$model_list_theta
    selected_model_new_object <- all_models$model_list[[selected_model_new]]
  } else {
    theta_models_list <- all_models$model_list_theta_treat
    selected_model_new_object <- all_models$model_list_treat[[selected_model_new]]
  }
  ## if selected model is STR or ID, do not test for theta
  if(STR_ID_flag) {
    message("\n", strrep("-", getOption("width")))
    message("No investigation for theta available, since theta cannot be estimated for the selected model,\nbecause there is no diversity effect")
    final_model_list <- list(selected_model_new_object)
    names(final_model_list) <- c(selected_model_new)
    selected_model_final <- selected_model_new
  } else {
    ## otherwise, proceed with the test for theta
    message("\n", strrep("-", getOption("width")))
    message("Step 3: Investigating whether theta is equal to 1 or not")
    theta_candidate_index <- grep(selected_model_new, names(theta_models_list))
    theta_candidate_name <- names(theta_models_list)[theta_candidate_index]
    theta_candidate_model <- theta_models_list[[theta_candidate_name]]
    th <- theta_candidate_model$coef["theta"]
    message("\nTheta estimate: ", round(th, 4), "\n", sep = "")
    final_model_list <- list(selected_model_new_object,
                             theta_candidate_model)
    names(final_model_list) <- c(selected_model_new,
                                 paste(selected_model_new, "_theta", sep = ""))
  
    mAIC_final <- sapply(final_model_list, AIC2)
    mAICc_final <- sapply(final_model_list, AICc)
    mBIC_final <- sapply(final_model_list, BIC2)
    mBICc_final <- sapply(final_model_list, BICc)
    selected_model_final <- switch(selection,
                                 Ftest = test_autoDI(model_list = final_model_list, family = family, treat = treat),
                                 AIC = AICsel_autoDI(model_list = final_model_list, mAIC = mAIC_final, treat = treat),
                                 AICc = AICcsel_autoDI(model_list = final_model_list, mAICc = mAICc_final, treat = treat),
                                 BIC = BICsel_autoDI(model_list = final_model_list, mBIC = mBIC_final, treat = treat),
                                 BICc = BICcsel_autoDI(model_list = final_model_list, mBICc = mBICc_final, treat = treat)
                                 )
    #fmla_char <- as.character(formula(final_model_list[[selected_model_final]]))[c(2,1,3)]
    message("\nSelected model: ",
        namesub_autoDI(selected_model_final), "\n",
        #"Formula: ",
        sep = "")
    #message(fmla_char, "\n", sep = " ")
  }
  
  if(step4) {
    ## step 4: lack-of-fit test (community model)
    message("\n", strrep("-", getOption("width")))
    message("Step 4: Comparing the final selected model with the reference (community) model")
    community <- get_community(prop = prop, data = data)
    model_to_compare <- final_model_list[[selected_model_final]]
    model_to_compare_data <- model_to_compare$data
    model_to_compare_data$community <- community
    ref_model <- update(model_to_compare, . ~ . + community,
                        data = model_to_compare_data)
    reftest_autoDI(model_to_compare, ref_model, family)
  }
  
  # returned object
  ret <- list("model_list" = all_models$model_list,
              "model_list_theta" = all_models$model_list_theta,
              "model_list_treat" = all_models$model_list_treat,
              "model_list_theta_treat" = all_models$model_list_theta_treat,
              "logLik" = llik,
              "AIC" = mAIC,
              "AICc" = mAICc,
              "BIC" = mBIC,
              "BICc" = mBICc,
              "df" = ndf,
              "selected_model_code" = selected_model_final,
              "selected_model" = namesub_autoDI(selected_model_final),
              "selected_model_obj" = final_model_list[[selected_model_final]],
              "data" = newdata,
              "family" = family)
  class(ret) <- "autoDI"
  message("\n", strrep("-", getOption("width")))
  message("autoDI is limited in terms of model selection. Exercise caution when choosing your final model.")
  message(strrep("-", getOption("width")))
  return(ret)
}

AIC2 <- function(obj) {
  n <- length(na.omit(obj$y))
  p <- length(na.omit(obj$coef))
  #n <- nrow(obj$data)
  #p <- length(obj$coef)
  mu_hat <- fitted(obj)
  sigma_hat <- sqrt(sum(obj$residuals^2)/(obj$df.residual) * (n - p)/n)
  ll <- sum(dnorm(obj$y, mu_hat, sigma_hat, log = TRUE))
  np <- p + 1
  aic <- - 2*ll + 2*np
  return(aic)
}

AICc.default <- function(obj) {
  aic <- AIC2(obj)
  np <- length(obj$coef) + 1
  n <- nrow(obj$data)
  aicc <- aic + (2*np^2 + 2*np)/(n - np - 1)
  return(aicc)
}

BIC2 <- function(obj) {
  n <- length(na.omit(obj$y))
  p <- length(na.omit(obj$coef))
  #n <- nrow(obj$data)
  #p <- length(obj$coef)
  mu_hat <- fitted(obj)
  sigma_hat <- sqrt(sum(obj$residuals^2)/(obj$df.residual) * (n - p)/n)
  ll <- sum(dnorm(obj$y, mu_hat, sigma_hat, log = TRUE))
  np <- p + 1
  bic <- - 2*ll + log(n)*np
  return(bic)
}

BICc.default <- function(obj) {
  bic <- BIC2(obj)
  np <- length(obj$coef) + 1
  n <- nrow(obj$data)
  bicc <- bic + (log(n)*(np+1)*np)/(n - np - 1)
  return(bicc)
}

#################################
## S3 methods for autoDI class ##
#################################

summary.autoDI <- function(object, ...) {
  summary_results <- summary(object$selected_model_obj, ...)
  return(summary_results)
}

print.autoDI <- function(x, ...) {
  print(x$selected_model_obj, ...)
}

anova.autoDI <- function(object, ...) {
  anova_results <- anova(object$selected_model_obj, ...)
  return(anova_results)
}

coef.autoDI <- function(object, ...) {
  coef_results <- coef(object$selected_model_obj, ...)
  return(coef_results)
}

model.matrix.autoDI <- function(object, ...) {
  mm_results <- model.matrix(object$selected_model_obj, ...)
  return(mm_results)
}

formula.autoDI <- function(x, ...) {
  f_results <- formula(x$selected_model_obj, ...)
  return(f_results)
}

hnp.autoDI <- function(object, ...) {
  hnp_results <- hnp(object$selected_model_obj, ...)
  return(invisible(hnp_results))
}

plot.autoDI <- function(x, ...) {
  plot(x$selected_model_obj, ...)
}

logLik.autoDI <- function(object, ...) {
  llik <- object$logLik
  npar <- object$df
  models <- namesub_autoDI(names(object$model_list))
  fam <- family(object$model_list[[1]])$family
  if (fam %in% c("gaussian", "Gamma", "inverse.gaussian")) 
    npar <- npar + 1
  for(i in 1:length(models)) {
    cat(models[i], ": ", llik[i], " (df = ", npar[i], ")\n", sep = "") 
  }
  return(invisible(llik))
}

AIC.autoDI <- function(object, ...) {
  mAIC <- sapply(object$model_list, AIC2)
  npar <- object$df
  models <- namesub_autoDI(names(object$model_list))
  for(i in 1:length(models)) {
    cat(models[i], ": ", mAIC[i], " (df = ", npar[i], ")\n", sep = "") 
  }
  return(invisible(mAIC))
}

AICc.autoDI <- function(obj) {
  mAICc <- sapply(obj$model_list, AICc)
  npar <- obj$df
  models <- namesub_autoDI(names(obj$model_list))
  for(i in 1:length(models)) {
    cat(models[i], ": ", mAICc[i], " (df = ", npar[i], ")\n", sep = "") 
  }
  return(invisible(mAICc))
}

BIC.autoDI <- function(object, ...) {
  mBIC <- sapply(object$model_list, BIC2)
  npar <- object$df
  models <- namesub_autoDI(names(object$model_list))
  for(i in 1:length(models)) {
    cat(models[i], ": ", mBIC[i], " (df = ", npar[i], ")\n", sep = "") 
  }
  return(invisible(mBIC))
}

BICc.autoDI <- function(obj) {
  mBICc <- sapply(obj$model_list, BICc)
  npar <- obj$df
  models <- namesub_autoDI(names(obj$model_list))
  for(i in 1:length(models)) {
    cat(models[i], ": ", mBICc[i], " (df = ", npar[i], ")\n", sep = "") 
  }
  return(invisible(mBICc))
}

AICc <- function(obj) {
  UseMethod("AICc")
}

BICc <- function(obj) {
  UseMethod("BICc")
}