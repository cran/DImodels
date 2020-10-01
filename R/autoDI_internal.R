allmodels_autoDI <- function(y, block, density, prop, treat, FG, data, family, total, estimate_theta,
                             nSpecies, treat_flag, even_flag, P_int_flag, FGnames) {
  ## detecting pathway
  ## if FG declared => use DI_FG, otherwise use DI_ADD
  if(is.null(FG)) use_FG <- FALSE else use_FG <- TRUE
  ## if treatment covariate is included, do extra test on treat removal from selected model
  ## treat_flag = TRUE when treat is missing
  
  ### standard pathway
  # STR model
  mod_STR <- DI_STR(y = y, block = block, density = density, data = data, family = family, total = total)$model
  # species identity model
  mod_ID <- DI_ID(y = y, block = block, density = density, prop = prop, data = data, family = family, total = total)$model
  # evenness model
  if(!even_flag) {
    mod_AV_both <- DI_AV(y = y, block = block, density = density, prop = prop, data = data, family = family,
                       estimate_theta = estimate_theta, nSpecies = nSpecies, total = total)
    mod_AV <- mod_AV_both$model
  }
  if(use_FG) {
    # functional groups model
    mod_FG_both <- DI_FG(y = y, block = block, density = density, prop = prop, FG = FG, data = data, family = family,
                         estimate_theta = estimate_theta, nSpecies = nSpecies, total = total, FGnames = FGnames)
    mod_FG <- mod_FG_both$model
  } else {
    # additive species contributions to interactions model
    if(!P_int_flag) {
      mod_ADD_both <- DI_ADD(y = y, block = block, density = density, prop = prop, data = data, family = family,
                             estimate_theta = estimate_theta, nSpecies = nSpecies, total = total)
      mod_ADD <- mod_ADD_both$model
    }
  }
  # separate pairwise interactions model
  fmla_FULL <- as.formula(paste("~", block, "+", density, " + (", paste(prop, collapse = "+"),")^2"))
  X_FULL <- model.matrix(fmla_FULL, data = data)
  X_pairwise <- X_FULL[,grep(":", colnames(X_FULL))]
  FULL_flag <- nrow(X_FULL) > ncol(X_FULL) # check if we have enough data points to estimate FULL model; if TRUE, we're good to use this model
  if(FULL_flag) {
    FULL_flag <- DI_matrix_check(X_pairwise) # double-check if pairwise interactions matrix is of full rank; if TRUE, we're good to use this model
  }
  if(!FULL_flag) {
    message("Not all pairwise interactions can be estimated.\nTherefore, the FULL model is not included in the selection process.\n") 
  }
  if(FULL_flag) {
    mod_FULL_both <- DI_FULL(y = y, block = block, density = density, prop = prop, data = data, family = family,
                             estimate_theta = estimate_theta, nSpecies = nSpecies, total = total)
    mod_FULL <- mod_FULL_both$model
    if(even_flag) {
      model_list_theta <- list("FULL_model" = mod_FULL_both$theta)
    } else if(use_FG) {
      model_list_theta <- list("AV_model" = mod_AV_both$theta,
                               "FG_model" = mod_FG_both$theta,
                               "FULL_model" = mod_FULL_both$theta)
    } else if(P_int_flag) {
      model_list_theta <- list("AV_model" = mod_AV_both$theta,
                               "FULL_model" = mod_FULL_both$theta)
    } else {
      model_list_theta <- list("AV_model" = mod_AV_both$theta,
                               "ADD_model" = mod_ADD_both$theta,
                               "FULL_model" = mod_FULL_both$theta) 
    }
  } else {
    if(even_flag) {
      model_list_theta <- list()
    } else if(use_FG) {
      model_list_theta <- list("AV_model" = mod_AV_both$theta,
                               "FG_model" = mod_FG_both$theta)
    } else if(P_int_flag) {
      model_list_theta <- list("AV_model" = mod_AV_both$theta)
    } else {
      model_list_theta <- list("AV_model" = mod_AV_both$theta,
                               "ADD_model" = mod_ADD_both$theta) 
    }
  }
    
  if(!treat_flag) {
    ## include treatment covariate
    # STR model
    mod_STR_treat <- DI_STR_treat(y = y, block = block, density = density, treat = treat, data = data, family = family, total = total)$model
    # species identity model
    mod_ID_treat <- DI_ID_treat(y = y, block = block, density = density, prop = prop, treat = treat, data = data, family = family, total = total)$model
    # evenness model
    if(!even_flag) {
      mod_AV_both_treat <- DI_AV_treat(y = y, block = block, density = density, prop = prop, treat = treat, data = data, family = family,
                           estimate_theta = estimate_theta, nSpecies = nSpecies, total = total)
      mod_AV_treat <- mod_AV_both_treat$model
    }
    if(use_FG) {
      # functional groups model
      mod_FG_both_treat <- DI_FG_treat(y = y, block = block, density = density, prop = prop, FG = FG, treat = treat, 
                                 data = data, family = family, estimate_theta = estimate_theta, nSpecies = nSpecies, total = total, FGnames = FGnames)
      mod_FG_treat <- mod_FG_both_treat$model
    } else {
      # additive species contributions to interactions model
      if(!P_int_flag) {
        mod_ADD_both_treat <- DI_ADD_treat(y = y, block = block, density = density, prop = prop, treat = treat, data = data, family = family,
                               estimate_theta = estimate_theta, nSpecies = nSpecies, total = total)
        mod_ADD_treat <- mod_ADD_both_treat$model
      }
    }
    # separate pairwise interactions model
    fmla_FULL_treat <- as.formula(paste("~", block, "+", density, "+", treat, " + (", paste(prop, collapse = "+"),")^2"))
    X_FULL_treat <- model.matrix(fmla_FULL_treat, data = data)
    X_pairwise_treat <- X_FULL_treat[,grep(":", colnames(X_FULL_treat))]
    FULL_flag <- nrow(X_FULL_treat) > ncol(X_FULL_treat) # check if we have enough data points to estimate FULL model; if TRUE, we're good to use this model
    if(FULL_flag) {
      FULL_flag <- DI_matrix_check(X_pairwise_treat) # double-check if pairwise interactions matrix is of full rank; if TRUE, we're good to use this model
    }
    if(!FULL_flag) {
      message("Not all pairwise interactions can be estimated.\nTherefore, the FULL model is not included in the selection process.\n") 
    }
    if(FULL_flag) {
      mod_FULL_both_treat <- DI_FULL_treat(y = y, block = block, density = density, prop = prop, treat = treat, data = data, family = family,
                                           estimate_theta = estimate_theta, nSpecies = nSpecies, total = total)
      mod_FULL_treat <- mod_FULL_both_treat$model
      if(even_flag) {
        model_list_theta_treat <- list("FULL_model_treat" = mod_FULL_both_treat$theta)
      } else if(use_FG) {
        model_list_theta_treat <- list("AV_model_treat" = mod_AV_both_treat$theta,
                                 "FG_model_treat" = mod_FG_both_treat$theta,
                                 "FULL_model_treat" = mod_FULL_both_treat$theta)
      } else if(P_int_flag) {
        model_list_theta_treat <- list("AV_model_treat" = mod_AV_both_treat$theta,
                                 "FULL_model_treat" = mod_FULL_both_treat$theta)
      } else {
        model_list_theta_treat <- list("AV_model_treat" = mod_AV_both_treat$theta,
                                 "ADD_model_treat" = mod_ADD_both_treat$theta,
                                 "FULL_model_treat" = mod_FULL_both_treat$theta) 
      }
    } else {
      if(even_flag) {
        model_list_theta_treat <- list()
      } else if(use_FG) {
        model_list_theta_treat <- list("AV_model_treat" = mod_AV_both_treat$theta,
                                       "FG_model_treat" = mod_FG_both_treat$theta)
      } else if(P_int_flag) {
        model_list_theta_treat <- list("AV_model_treat" = mod_AV_both_treat$theta)
      } else {
        model_list_theta_treat <- list("AV_model_treat" = mod_AV_both_treat$theta,
                                       "ADD_model_treat" = mod_ADD_both_treat$theta) 
      }
    }
    
  }
  
  ## obtaining the four model lists
  ## list 1: standard models, no treatment covariate, theta = 1
  model_list <- list("STR_model" = mod_STR,
                     "ID_model" = mod_ID)
  if(!even_flag) {
    model_list$AV_model <- mod_AV
  }
  if(use_FG) {
    model_list$FG_model <- mod_FG 
  } else if(!P_int_flag) {
    model_list$ADD_model <- mod_ADD 
  }
  if(FULL_flag) {
    model_list$FULL_model <- mod_FULL
  }
  
  ## list 2: models with treatment covariate, theta = 1
  if(!treat_flag) {
    model_list_treat <- list("STR_model_treat" = mod_STR_treat,
                       "ID_model_treat" = mod_ID_treat)
    if(!even_flag) {
      model_list_treat$AV_model_treat <- mod_AV_treat
    }
    if(use_FG) {
      model_list_treat$FG_model_treat <- mod_FG_treat
    } else if(!P_int_flag) {
      model_list_treat$ADD_model_treat <- mod_ADD_treat
    }
    if(FULL_flag) {
      model_list_treat$FULL_model_treat <- mod_FULL_treat
    }
  }
  
  if(treat_flag) {
    return(list("model_list" = model_list, "model_list_theta" = model_list_theta,
                "model_list_treat" = list(), "model_list_theta_treat" = list()))
  } else {
    return(list("model_list" = model_list, "model_list_theta" = model_list_theta,
                "model_list_treat" = model_list_treat, "model_list_theta_treat" = model_list_theta_treat))
  }
}

namesub_autoDI <- Vectorize(function(name) {
  thename <- switch(name,
                    "STR_model" = "Structural 'STR' DImodel",
                    "ID_model" = "Species identity 'ID' DImodel",
                    "FULL_model" = "Separate pairwise interactions 'FULL' DImodel",
                    "AV_model" = "Average interactions 'AV' DImodel",
                    "E_model" = "Evenness 'E' DImodel",
                    "ADD_model" = 
                      "Additive species contributions to interactions 'ADD' DImodel",
                    "FG_model" = "Functional group effects 'FG' DImodel",
                    "STR_model_treat" = "Structural 'STR' DImodel with treatment",
                    "ID_model_treat" = "Species identity 'ID' DImodel with treatment",
                    "FULL_model_treat" = "Separate pairwise interactions 'FULL' DImodel with treatment",
                    "AV_model_treat" = "Average interactions 'AV' DImodel with treatment",
                    "E_model_treat" = "Evenness 'E' DImodel with treatment",
                    "ADD_model_treat" = 
                      "Additive species contributions to interactions 'ADD' DImodel with treatment",
                    "FG_model_treat" = "Functional group effects 'FG' DImodel with treatment",
                    "FULL_model_theta" = "Separate pairwise interactions 'FULL' DImodel, estimating theta",
                    "AV_model_theta" = "Average interactions 'AV' DImodel, estimating theta",
                    "E_model_theta" = "Evenness 'E' DImodel, estimating theta",
                    "ADD_model_theta" = 
                      "Additive species contributions to interactions 'ADD' DImodel, estimating theta",
                    "FG_model_theta" = "Functional group effects 'FG' DImodel, estimating theta",
                    "FULL_model_treat_theta" = "Separate pairwise interactions 'FULL' DImodel with treatment, estimating theta",
                    "AV_model_treat_theta" = "Average interactions 'AV' DImodel with treatment, estimating theta",
                    "E_model_treat_theta" = "Evenness 'E' DImodel with treatment, estimating theta",
                    "ADD_model_treat_theta" = 
                      "Additive species contributions to interactions 'ADD' DImodel with treatment, estimating theta",
                    "FG_model_treat_theta" = "Functional group effects 'FG' DImodel with treatment, estimating theta",
                    stop("not yet implemented"))
  return(thename)
}, "name")

reftest_autoDI <- function(model_to_compare, ref_model, family) {
  if(family %in% c("poisson","binomial")) {
    ref_test <- "Chisq"
  } else {
    ref_test <- "F"
  }
  anovas <- anova(model_to_compare, ref_model, test = ref_test)
  ## formatting the output
  anovas_format <- as.data.frame(anovas)
  anovas_format <- round(anovas_format, 4)
  anovas_format[is.na(anovas_format)] <- ""
  names(anovas_format)[c(2,4)] <- c("Resid. SSq","SSq")
  anovas_format$"Resid. MSq" <- round(anovas_format[,2]/anovas_format[,1], 4)
  anovas_format <- anovas_format[,c(1,2,7,3:6)]
  model_tokens <- c("Selected","Reference")
  anovas_format$model <- model_tokens
  anovas_format <- anovas_format[,c(8,1:7)]
  row.names(anovas_format) <- paste("DI Model", row.names(anovas_format))
  anovas_format[,8][anovas_format[,8] == 0] <- "<0.0001"
  message("\n")
  old <- options()
  options(scipen = 999)
  on.exit(options(old))
  print(anovas_format)
  #options(scipen = 0)
  #message("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
}

test_autoDI <- function(model_list, family, treat) {
  if(family %in% c("poisson","binomial")) {
    message("Selection using X2 tests", "\n")
    Test <- "Chisq"
  } else {
    message("Selection using F tests", "\n")
    Test <- "F"
  }
  
  model_names <- names(model_list)
  treat_flag <- grep("treatment", namesub_autoDI(model_names))
  theta_flag <- grep("theta", namesub_autoDI(model_names))
  
  treat_output <- rep("none", length(model_names))
  treat_output[treat_flag] <- paste("'", treat, "'", sep = "")
  theta_output <- rep(FALSE, length(model_names))
  theta_output[theta_flag] <- TRUE
  
  model_tokens <- gsub("_", "", gsub("[:a-z:]", "", names(model_list)))
  names(model_list) <- NULL
  #anovas <- do.call(anova, model_list) ## if using lm objects
  anovas <- eval(parse(text = paste("anova(",
                                    paste("model_list[[", 1:length(model_list), "]]",
                                          sep = "", collapse = ","),
                                    ",test ='", Test, "')", sep = "")
                                   ))
  if(family %in% c("poisson","binomial")) {
    p_values <- anovas$"Pr(>Chi)"
  } else {
    p_values <- anovas$"Pr(>F)"
  }
  p_less <- which(p_values < .05)
  p_value_selected <- ifelse(length(p_less) == 0, 1, max(p_less))
  selected <- model_names[p_value_selected]
  ## formatting the output
  anovas_format <- as.data.frame(anovas)
  anovas_format <- round(anovas_format, 4)
  anovas_format[is.na(anovas_format)] <- ""
  names(anovas_format)[c(2,4)] <- c("Resid. SSq","SSq")
  anovas_format$"Resid. MSq" <- round(anovas_format[,2]/anovas_format[,1], 4)
  anovas_format <- anovas_format[,c(1,2,7,3:6)]
  anovas_format$model <- model_tokens
  anovas_format$treat <- treat_output
  anovas_format$theta <- theta_output
  anovas_format <- anovas_format[,c(8:10,1:7)]
  names(anovas_format)[1] <- "DI_model"
  names(anovas_format)[3] <- "estimate_theta"
  row.names(anovas_format) <- paste("DI Model", row.names(anovas_format))
  anovas_format[,10][anovas_format[,10] == 0] <- "<0.0001"
  desc_table <- data.frame("Description" = namesub_autoDI(model_names))
  row.names(desc_table) <- paste("DI Model", 1:nrow(anovas_format))
  print(desc_table, right = FALSE)
  message("\n")
  old <- options()
  options(scipen = 999)
  on.exit(options(old))
  print(anovas_format)
  #options(scipen = 0)
  #message("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  return(selected)
}

AICsel_autoDI <- function(model_list, mAIC, treat) {
  message("Selection by AIC\nWarning: DI Model with the lowest AIC will be selected, even if the difference is very small.\nPlease inspect other models to see differences in AIC.", "\n\n", sep = "")
  
  model_names <- names(model_list)
  treat_flag <- grep("treatment", namesub_autoDI(model_names))
  theta_flag <- grep("theta", namesub_autoDI(model_names))
  
  treat_output <- rep("none", length(model_names))
  treat_output[treat_flag] <- paste("'", treat, "'", sep = "")
  theta_output <- rep(FALSE, length(model_names))
  theta_output[theta_flag] <- TRUE
  
  model_tokens <- gsub("_", "", gsub("[:a-z:]", "", model_names))
  model_descriptions <- namesub_autoDI(model_names)
  the_table <- data.frame("AIC" = mAIC,
                          "DI_Model" = model_tokens,
                          "treat" = treat_output,
                          "theta" = theta_output,
                          "Description" = model_descriptions,
                          row.names = NULL)
  print(the_table, right = FALSE)
  selected <- names(model_list)[which.min(mAIC)]
  return(selected)
}

AICcsel_autoDI <- function(model_list, mAICc, treat) {
  message("Selection by AICc\nWarning: DI Model with the lowest AICc will be selected, even if the difference is very small.\nPlease inspect other models to see differences in AICc.", "\n\n", sep = "")
  
  model_names <- names(model_list)
  treat_flag <- grep("treatment", namesub_autoDI(model_names))
  theta_flag <- grep("theta", namesub_autoDI(model_names))
  
  treat_output <- rep("none", length(model_names))
  treat_output[treat_flag] <- paste("'", treat, "'", sep = "")
  theta_output <- rep(FALSE, length(model_names))
  theta_output[theta_flag] <- TRUE
  
  model_tokens <- gsub("_", "", gsub("[:a-z:]", "", model_names))
  model_descriptions <- namesub_autoDI(model_names)
  the_table <- data.frame("AICc" = mAICc,
                          "DI_Model" = model_tokens,
                          "treat" = treat_output,
                          "theta" = theta_output,
                          "Description" = model_descriptions,
                          row.names = NULL)
  print(the_table, right = FALSE)
  selected <- names(model_list)[which.min(mAICc)]
  return(selected)
}

BICsel_autoDI <- function(model_list, mBIC, treat) {
  message("Selection by BIC\nWarning: DI Model with the lowest BIC will be selected, even if the difference is very small.\nPlease inspect other models to see differences in BIC.", "\n\n", sep = "")
  
  model_names <- names(model_list)
  treat_flag <- grep("treatment", namesub_autoDI(model_names))
  theta_flag <- grep("theta", namesub_autoDI(model_names))
  
  treat_output <- rep("none", length(model_names))
  treat_output[treat_flag] <- paste("'", treat, "'", sep = "")
  theta_output <- rep(FALSE, length(model_names))
  theta_output[theta_flag] <- TRUE
  
  model_tokens <- gsub("_", "", gsub("[:a-z:]", "", model_names))
  model_descriptions <- namesub_autoDI(model_names)
  the_table <- data.frame("BIC" = mBIC,
                          "DI_Model" = model_tokens,
                          "treat" = treat_output,
                          "theta" = theta_output,
                          "Description" = model_descriptions,
                          row.names = NULL)
  print(the_table, right = FALSE)
  selected <- names(model_list)[which.min(mBIC)]
  return(selected)
}

BICcsel_autoDI <- function(model_list, mBICc, treat) {
  message("Selection by BICc\nWarning: DI Model with the lowest BICc will be selected, even if the difference is very small.\nPlease inspect other models to see differences in BICc.", "\n\n", sep = "")
  
  model_names <- names(model_list)
  treat_flag <- grep("treatment", namesub_autoDI(model_names))
  theta_flag <- grep("theta", namesub_autoDI(model_names))
  
  treat_output <- rep("none", length(model_names))
  treat_output[treat_flag] <- paste("'", treat, "'", sep = "")
  theta_output <- rep(FALSE, length(model_names))
  theta_output[theta_flag] <- TRUE
  
  model_tokens <- gsub("_", "", gsub("[:a-z:]", "", model_names))
  model_descriptions <- namesub_autoDI(model_names)
  the_table <- data.frame("BICc" = mBICc,
                          "DI_Model" = model_tokens,
                          "treat" = treat_output,
                          "theta" = theta_output,
                          "Description" = model_descriptions,
                          row.names = NULL)
  print(the_table, right = FALSE)
  selected <- names(model_list)[which.min(mBICc)]
  return(selected)
}

DI_matrix_check <- function(model_matrix) {
  n_parms <- ncol(model_matrix)
  matrix_rank <- qr(model_matrix)$rank
  if(matrix_rank == n_parms) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

autoDI_step0 <- function(y, block, density, treat, family, data) {
  if(is.na(block) & is.na(density) & is.na(treat)) {
    return(invisible())
  }
  
  message("\n", strrep("-", getOption("width")))
  message("\nSequential analysis: Investigating only non-diversity experimental design structures\n")
  fmla1 <- paste(y, "~", 1)
  fit1 <- glm(fmla1, family = family, data = data)
  
  if(!is.na(block) & !is.na(density) & !is.na(treat)) {
    fmla2 <- paste(y, "~", block)
    fmla3 <- paste(y, "~", block, "+", density)
    fmla4 <- paste(y, "~", block, "+", density, "+", treat)
    fit2 <- glm(fmla2, family = family, data = data)
    fit3 <- glm(fmla3, family = family, data = data)
    fit4 <- glm(fmla4, family = family, data = data)
    model_list <- list(fit1, fit2, fit3, fit4)
    model_tokens <- c("Intercept only","block","block + density","block + density + treat")
  }
  
  if(!is.na(block) & !is.na(density) & is.na(treat)) {
    fmla2 <- paste(y, "~", block)
    fmla3 <- paste(y, "~", block, "+", density)
    fit2 <- glm(fmla2, family = family, data = data)
    fit3 <- glm(fmla3, family = family, data = data)
    model_list <- list(fit1, fit2, fit3)
    model_tokens <- c("Intercept only","block","block + density")
  }
  
  if(!is.na(block) & is.na(density) & !is.na(treat)) {
    fmla2 <- paste(y, "~", block)
    fmla3 <- paste(y, "~", block, "+", treat)
    fit2 <- glm(fmla2, family = family, data = data)
    fit3 <- glm(fmla3, family = family, data = data)
    model_list <- list(fit1, fit2, fit3)
    model_tokens <- c("Intercept only","block","block + treat")
  }
  
  if(is.na(block) & !is.na(density) & !is.na(treat)) {
    fmla2 <- paste(y, "~", density)
    fmla3 <- paste(y, "~", density, "+", treat)
    fit2 <- glm(fmla2, family = family, data = data)
    fit3 <- glm(fmla3, family = family, data = data)
    model_list <- list(fit1, fit2, fit3)
    model_tokens <- c("Intercept only","density","density + treat")
  }
  
  if(!is.na(block) & is.na(density) & is.na(treat)) {
    fmla2 <- paste(y, "~", block)
    fit2 <- glm(fmla2, family = family, data = data)
    model_list <- list(fit1, fit2)
    model_tokens <- c("Intercept only","block")
  }
  
  if(is.na(block) & !is.na(density) & is.na(treat)) {
    fmla2 <- paste(y, "~", density)
    fit2 <- glm(fmla2, family = family, data = data)
    model_list <- list(fit1, fit2)
    model_tokens <- c("Intercept only","density")
  }
  
  if(is.na(block) & is.na(density) & !is.na(treat)) {
    fmla2 <- paste(y, "~", treat)
    fit2 <- glm(fmla2, family = family, data = data)
    model_list <- list(fit1, fit2)
    model_tokens <- c("Intercept only","treat")
  }
  
  if(family %in% c("poisson","binomial")) {
    Test <- "Chisq"
  } else {
    Test <- "F"
  }
  
  anovas <- eval(parse(text = paste("anova(",
                                    paste("model_list[[", 1:length(model_list), "]]",
                                          sep = "", collapse = ","),
                                    ",test ='", Test, "')", sep = "")
  ))
  
  ## formatting the output
  anovas_format <- as.data.frame(anovas)
  anovas_format <- round(anovas_format, 4)
  anovas_format[is.na(anovas_format)] <- ""
  names(anovas_format)[c(2,4)] <- c("Resid. SSq","SSq")
  anovas_format$"Resid. MSq" <- round(anovas_format[,2]/anovas_format[,1], 4)
  anovas_format <- anovas_format[,c(1,2,7,3:6)]
  anovas_format$model <- model_tokens
  anovas_format <- anovas_format[,c(8,1:7)]
  anovas_format[,8][anovas_format[,8] == 0] <- "<0.0001"
  message("\n")
  old <- options()
  options(scipen = 999)
  on.exit(options(old))
  print(anovas_format)
  #options(scipen = 0)
}