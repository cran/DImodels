DI_data_prepare <- function(y, block, density, prop, treat, FG = NULL, data, theta = 1) {
  if(any(is.na(data))) {
    stop("The dataset contains missing values. Please remove them prior to analysis.") 
  }
  if(missing(y)) {
    y <- rep(0, nrow(data))
    warning("y was not supplied, so a vector of zeros was created for the response\n")
  }
  if(missing(block)) block <- NA
  if(missing(density)) density <- NA
  if(missing(treat)) treat <- NA
  # variables directly obtained from dataset
  if(is.na(block)) {
    block <- "block_zero"
    data$block_zero <- 0
  } else {
    if(!is.character(block)) block <- names(data)[block] # variable block (name)
    ## warning for continuous block variables
    if(!is.factor(data[,block])) {
      #data[,block] <- as.factor(data[,block])
      warning("'", block, "' has been declared as a block variable, but is not a factor.")#,
              #" It has been converted to a factor.")
    }
  }
  if(is.na(density)) {
    density <- "density_zero"
    data$density_zero <- 0
  } else {
    if(!is.character(density)) density <- names(data)[density] # variable density (name)
    ## warning for continuous density variables
    if(!is.factor(data[,density])) {
      #data[,density] <- as.factor(data[,density])
      warning("'", density, "' has been declared as a density variable, but is not a factor.")#,
      #" It has been converted to a factor.")
    }
  }
  if(is.na(treat)) {
    treat <- "treat_zero"
    data$treat_zero <- 0
  }
  
  if(!is.character(y)) y <- names(data)[y] # response variable (name)
  
  ## getting the indices for the columns corresponding to the species proportions (Pind)
  Pind_and_prop <- get_P_indices(prop = prop, data = data)
  Pind <- Pind_and_prop$Pind
  
  ## checking if the Pi's sum to 1
  prop_check <- DI_prop_check(Pind = Pind, data = data)
  if(prop_check == "error") {
    stop("One or more rows have species proportions that do not sum to 1. This must be corrected prior to analysis.\n")
  }
  if(prop_check == "minor") {
    Pi_sums <- apply(data[,Pind], 1, sum)
    data[,Pind] <- data[,Pind]/Pi_sums
    warning("One or more rows have species proportions that sum to approximately 1, but not exactly 1. This is typically a rounding issue, and has been corrected internally prior to analysis.\n")
  }
  
  ## species proportions column names (prop)
  prop <- Pind_and_prop$prop

  if(!is.character(treat)) treat <- names(data)[treat] # environmental covariate (name)
  ## calculating the E and AV variables
  E_AV <- DI_data_E_AV(prop = prop, data = data, theta = theta)
  newdata <- data
  newdata$AV <- E_AV$AV
  newdata$E <- E_AV$E
  even_flag <- E_AV$even_flag
  ## calculating the P*(1-P) (P_add) variables
  ADD <- DI_data_ADD(prop = prop, data = data, theta = theta)
  newdata <- data.frame(newdata, ADD$ADD)
  P_int_flag <- ADD$P_int_flag
  ## calculating FG variables
  if(!is.null(FG)) FG <- DI_data_FG(prop = prop, FG = FG, data = data, theta = theta)$FG
  ## return object
  return(list("newdata" = newdata, "y" = y, "block" = block, density = density,
              "prop" = prop, "treat" = treat, "FG" = FG,
              "P_int_flag" = P_int_flag, "even_flag" = even_flag, "nSpecies" = length(prop)))
}

get_P_indices <- function(prop, data) {
  if(!is.character(prop)) {
    Pind <- prop # indices for the columns with species proportions P_i
    prop <- names(data[, prop]) # species proportions P_i (names)
  } else {
    vec_grep <- Vectorize(grep, "pattern")
    Pind <- as.numeric(vec_grep(paste("\\<", prop, "\\>", sep = ""),
                                names(data)))
  }
  return(list(Pind = Pind, prop = prop))
}

DI_data_E_AV_internal <- function(prop, data) {
  if(!is.character(prop)) {
    prop <- get_P_indices(prop = prop, data = data)$prop
  }
  nSpecies <- length(prop) # number of species
  if(nSpecies <= 1) stop("must have at least 2 species to fit DI models")
  nComb <- choose(nSpecies, 2) # number of pairwise combinations
  pairwise_fmla <- as.formula(paste("~", "0+", "(",
                                    paste(prop, collapse = "+"), ")^2"))
  normE <- 2 * nSpecies/(nSpecies - 1)
  Ematrix <- model.matrix(pairwise_fmla, data = data)
  # evenness and "AV" variable
  if(nSpecies > 2) {
    AV <- rowSums(Ematrix[,(nSpecies+1):ncol(Ematrix)])
    E <- normE * AV
    even_flag <- FALSE
    return(list("AV" = AV, "E" = E, "even_flag" = even_flag))
  } else {
    even_flag <- TRUE
    return(list("even_flag" = even_flag))
  }
}

DI_data_E_AV <- function(prop, data, theta = 1) {
  if(theta != 1) {
    data[,prop] <- data[,prop]^theta 
  }
  result <- DI_data_E_AV_internal(prop = prop, data = data)
  return(result)
}

DI_data_ADD_internal <- function(prop, data) {
  Pind <- get_P_indices(prop = prop, data = data)$Pind
  nSpecies <- length(prop)
  if(nSpecies > 3) {
    ADD_vars <- as.data.frame(apply(data[, Pind], 2, function(x) x*(1-x)))
    prop_names <- names(data)[Pind]
    names(ADD_vars) <- paste(prop_names, "_add", sep = "")
    P_int_flag <- FALSE
  } else {
    ADD_vars <- NULL
    P_int_flag <- TRUE
  }
  return(list("ADD" = ADD_vars, "P_int_flag" = P_int_flag))
}

DI_data_ADD_theta <- function(prop, data, theta) {
  Pind <- get_P_indices(prop = prop, data = data)$Pind
  nSpecies <- length(prop)
  P_matrix <- data[,Pind]
  Pi_theta <- P_matrix^theta
  ADD_vars_theta <- matrix(NA, ncol = nSpecies, nrow = nrow(data))
  for(i in 1:nSpecies) {
    sum_Pj_theta <- apply(Pi_theta[,-i], 1, sum) 
    ADD_vars_theta[,i] <- Pi_theta[,i] * sum_Pj_theta
  }
  ADD_vars_theta <- as.data.frame(ADD_vars_theta)
  prop_names <- names(data)[Pind]
  names(ADD_vars_theta) <- paste(prop_names, "_add", sep = "")
  return(list("ADD_theta" = ADD_vars_theta))
}

DI_data_ADD <- function(prop, data, theta = 1) {
  if(theta == 1) {
    result <- DI_data_ADD_internal(prop = prop, data = data)
  } else {
    result <- DI_data_ADD_theta(prop = prop, data = data, theta = theta) 
  }
  return(result)
}

DI_data_FG_internal <- function(prop, FG, data) {
  # number of functional groups
  nfg <- length(unique(FG))
  # n for checking at the end
  n_check <- choose(nfg, 2) + nfg
  
  Pind <- get_P_indices(prop = prop, data = data)$Pind
  if(any(!is.character(FG))) stop("FG argument takes character strings with functional",
                                  " group names referring to each species, in order")
  fg_index <- FG
  fg_index_names <- levels(as.factor(fg_index))
  testdata <- data[,Pind]
  prop <- colnames(testdata)
  prop <- paste(prop, fg_index, sep = "_infg_")
  testdata2 <- testdata
  colnames(testdata2) <- prop
  nSpecies <- length(prop)
  pairwise_fmla <- as.formula(paste("~", "0+", "(",
                                    paste(prop, collapse = "+"), ")^2"))
  Pmatrix <- model.matrix(pairwise_fmla, data = testdata2)
  FGmatrix_raw <- Pmatrix[,(nSpecies+1):ncol(Pmatrix)]
  FGrawnames <- colnames(FGmatrix_raw)
  ## matching between FG effects
  #nfg <- length(fg_index_names)
  FGmatch <- list()
  for(i in 1:nfg) {
    FGmatch[[i]] <- grep(fg_index_names[i], FGrawnames)
  }
  FGnames <- rep(NA, length(FGrawnames))
  for(i in 1:nfg) {
    for(j in 1:nfg) {
      if(i < j) {
        FGnames[FGmatch[[i]][FGmatch[[i]] %in% FGmatch[[j]]]] <- 
          paste("bfg_", fg_index_names[i], ":", fg_index_names[j], sep = "")
      }
    }
  }
  ## matching within FG effects
  FGmatch2 <- list()
  for(i in 1:nfg) {
    FGmatch2[[i]] <- grep(paste(fg_index_names[i], ":", sep = ""), FGrawnames)
  }
  for(i in 1:nfg) {
    helper_index <- FGmatch[[i]][FGmatch[[i]] %in% FGmatch2[[i]]]
    FGnames[helper_index[is.na(FGnames[helper_index])]] <- 
      paste("wfg_", fg_index_names[i], sep = "")
  }
  ## summing variables belonging to same FG effect
  FGeffects <- levels(as.factor(FGnames))
  FG <- matrix(NA, nrow = nrow(FGmatrix_raw), ncol = length(FGeffects))
  for(i in 1:length(FGeffects)) {
    FG[,i] <- rowSums(as.matrix(FGmatrix_raw[,FGnames == FGeffects[i]]))
  }
  ## names of FG variables
  FGeffects <- gsub(":","_", FGeffects)
  colnames(FG) <- FGeffects
  
  ## final checks
  if(n_check != ncol(FG)) {
    stop("Expected ", n_check, " terms, but have ", ncol(FG),
         ". Please give your functional groups a different name.",
         " One or more of the options are being used internally.",
         " Perhaps use upper-case names?")
  }
  
  return(list("FG" = FG))
}

DI_data_FG <- function(prop, FG, data, theta = 1) {
  if(theta != 1) {
    data[,prop] <- data[,prop]^theta 
  }
  result <- DI_data_FG_internal(prop = prop, FG = FG, data = data)
  return(result)
}

DI_data_fullpairwise <- function(prop, data, theta = 1) {
  prop <- get_P_indices(prop = prop, data = data)$prop
  fmla <- as.formula(paste("~ 0 + ", "(", paste(prop, collapse = "+"), ")^2"))
  all_pairwise <- model.matrix(fmla, data = data)
  n_species <- length(prop)
  obj <- all_pairwise[,-c(1:n_species)]^theta
  return(obj)
}

DI_prop_check <- function(Pind, data) {
  props <- data[,Pind]
  pi_sums <- apply(props, 1, sum)
  if(any(pi_sums > 1.0001) | any(pi_sums < .9999)) {
    return("error")
  } else if(any(pi_sums < 1 & pi_sums > .9999) | any(pi_sums > 1 & pi_sums < 1.0001)) {
    return("minor") 
  } else return("ok")
}