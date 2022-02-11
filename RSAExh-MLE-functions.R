
# Function to fit a pair of production/comprehension data sets given corresponding likelihood functions and an init.
# Returns maximum likelihood estimates for all parameters, the max log-likelihood, and number of dropped data points if any (should be 0)

TobitModelFit <- function(
  Production_LogLik_FUN, # assumed to take the same number of arguments as provided in init, in the same order.
  Comprehension_LogLik_FUN, # ditto
  init, # parameters assumed to be lambda, deltaAB, deltaAnB, any number of optional parameters in [0,1], sigma_a, sigma_ab, in this order. eps is always initialized as -4
  prod_data, # assumed to have columns Prior, World, and Message
  comp_data,
  ignore_prod=F # whether to ignore the production data (only used for separate parameters model fits)
){
  # Initialize
  init <- unlist(init) # In case it was given as a list
  Np <- length(init) # number of parameters
  init[2:3] <- if_else(init[2:3]>0,init[2:3]*init[1],0) # used scaled costs lambda*Delta for fitting
  # Transform init to the correct scale (logit for optional, log for all others)
  par=suppressWarnings(c(qlogis(init),eps=-4))
  par[c(1:3,(Np-1):Np)] <- log(init[c(1:3,(Np-1):Np)])
  boundary_values <- which(!is.finite(par))
  negloglik <- function(par){
    # Transform the parameters back to feed likelihood functions
    trans_par <- plogis(par)
    trans_par[c(1:3,(Np-1):(Np+1))] <- exp(par[c(1:3,(Np-1):(Np+1))])
    trans_par[2:3] <- trans_par[2:3]/trans_par[1] # retrieve the actual cost, as required by the likelihood FUN
    # Leave the boundary values untouched:
    trans_par[boundary_values] <- init[boundary_values]
    # Compute the log-likelihood for prod and comp data:
    tmp_prod_data <- prod_data
    tmp_prod_data$loglik=do.call(Production_LogLik_FUN, unname(c(list(tmp_prod_data$Message,tmp_prod_data$World,tmp_prod_data$Prior),as.list(trans_par)[3:Np-2])))
    if(!ignore_prod){
      tmp_prod_data$loglik <- log((exp(tmp_prod_data$loglik)+trans_par[Np+1])/(1+4*trans_par[Np+1]))
    }
    # Remove data points that have 0 likelihood (just in case, shouldn't happen)
    tmp_prod_data <- tmp_prod_data[is.finite(tmp_prod_data$loglik),]
    tmp_comp_data <- comp_data
    tmp_comp_data$loglik <-do.call(Comprehension_LogLik_FUN, unname(c(list(tmp_comp_data$Posterior,tmp_comp_data$Utterance,tmp_comp_data$Prior),as.list(trans_par)[1:Np])))
    tmp_comp_data <- tmp_comp_data[is.finite(tmp_comp_data$loglik),] # Remove data points that have 0 likelihood (just in case, shouldn't happen)
    negloglik=-sum(tmp_prod_data$loglik)-sum(tmp_comp_data$loglik)
    return(negloglik)
  }
  result <- optim(par,negloglik,method="L-BFGS-B",lower=-7,upper=7)
  trans_par <- plogis(result$par)
  trans_par[c(1:3,(Np-1):(Np+1))] <- exp(result$par[c(1:3,(Np-1):(Np+1))])
  trans_par[2:3] <- trans_par[2:3]/trans_par[1]
  trans_par[boundary_values] <- init[boundary_values]
  tmp_prod_data <- prod_data %>%
    mutate(
      loglik=do.call(Production_LogLik_FUN, unname(c(list(Message,World,Prior),as.list(trans_par)[3:Np-2]))))
  if(!ignore_prod){
    tmp_prod_data <- mutate(tmp_prod_data,loglik=log((exp(loglik)+trans_par[Np+1])/(1+4*trans_par[Np+1])))
  }
  tmp_comp_data <- comp_data %>%
    mutate(
      loglik=do.call(Comprehension_LogLik_FUN, unname(c(list(Posterior,Utterance,Prior),as.list(trans_par)[1:Np])))
    )
  dropped=sum(!is.finite(tmp_prod_data$loglik))+sum(!is.finite(tmp_comp_data$loglik))
  output <- c(trans_par,-result$value,dropped)
  names(output) <- c(names(init),"eps","logLik","dropped")
  return(output)
}

# Function to apply TobitModelFit twice in order to fit parameters independently for production and comprehension (as promised in the pre-registration)
TobitSeparateModelFit <- function(
  Production_LogLik_FUN,
  Comprehension_LogLik_FUN,
  init_comp,
  init_prod=init_comp,
  prod_data, # assumed to have columns Prior, World, and Message
  comp_data # assumed to have columns Prior, Utterance, and Posterior
){
  Np <- length(init_comp)
  # Define fake likelihood functions which will essentially ignore the production and comprehension data respectively
  emptyProdFUN <- function(Message,World,Prior,...){0}
  emptyCompFUN <- function(Posterior,Utterance,Prior,...){0}
  # Run the main function but replacing either the production or the comprehension likelihood function by the constant.
  CompResults <- TobitModelFit(emptyProdFUN,Comprehension_LogLik_FUN,init_comp,prod_data,comp_data,ignore_prod = T)
  ProdResults <- TobitModelFit(Production_LogLik_FUN,emptyCompFUN,init_prod,prod_data,comp_data)
  # Take the name of the parameters from the init that was provided for comprehension
  Names <- names(CompResults)
  # Add _comp and _prod suffixes to the parameters
  names(CompResults) <- paste(Names,"_comp",sep="")
  names(ProdResults) <- paste(Names,"_prod",sep="")
  # Order everything nicely and return the log-likelihood on each data set, plus the sum.
  Out <- c(CompResults[1:(Np)],ProdResults[c(1:(Np-2),Np+1)],CompResults[2:3+Np],ProdResults[2:3+Np])
  Out[["logLik"]] <- Out[["logLik_comp"]]+Out[["logLik_prod"]]
  Out[["dropped"]] <- Out[["dropped_comp"]]+Out[["dropped_prod"]]
  return(Out)
}

# Does a complete fit given a list of initialization (in parallel using futures).
# Systematically tests for extreme values in case a parameter gets close to boundary or +Inf (except for lambda).
# Computes AIC and BIC for joint and separate fits using the best results if different initializations gave divergent results.
complete_model_fit <- function(
  Production_LogLik_FUN,
  Comprehension_LogLik_FUN,
  init_list,
  prod_data, # assumed to have columns Prior, World, and Message
  comp_data,
  N_param=length(init_list[[1]]),
  model_name="Model"
){
  param_names=names(init_list[[1]])
  # Fit for each initialization
  Results_joint <- bind_rows(future_lapply(
    init_list,
    function(init){TobitModelFit(Production_LogLik_FUN,Comprehension_LogLik_FUN,init,prod_data,comp_data)}
  ))
  # Keep only the best:
  parms_joint <- Results_joint %>% filter(dropped==min(dropped)) %>% filter(logLik==max(logLik,na.rm=T)) %>% head(1)
  # Find candidate limit parameters:
  par <- as_vector(parms_joint[1:N_param])
  par[2:3] <- par[2:3]*par[1]
  trans_par <- suppressWarnings(qlogis(par))
  trans_par[c(1:3,(N_param-1):N_param)] <- log(par[c(1:3,(N_param-1):N_param)])
  low_param <- which(round(trans_par) == -7)
  high_param <- which(round(trans_par[2:N_param]) == 7)+1 # don't test infinite lambda
  init_extreme <- NULL
  k=1
  # Define all possible combination of extreme parameters:
  for(i in 1:2^length(low_param)){
    for(j in 1:2^length(high_param)){
      trans_init=trans_par
      trans_init[low_param[as.logical(intToBits(i-1))[1:length(low_param)]]] <- -Inf
      trans_init[high_param[as.logical(intToBits(j-1))[1:length(high_param)]]] <- +Inf
      init <- plogis(trans_init)
      init[c(1:3,(N_param-1):N_param)] <- exp(trans_init[c(1:3,(N_param-1):N_param)])
      if(init[1]>0){init[2:3] <- init[2:3]/init[1]}
      init_extreme[[k]] <- init
      k <- k+1
    }
  }
  Results_joint_extreme <- bind_rows(future_lapply(
    init_extreme,
    function(inits){TobitModelFit(Production_LogLik_FUN,Comprehension_LogLik_FUN,inits,prod_data=prod_data,comp_data=comp_data)}
  ))
  # See whether the extreme values actually improves log-likelihood, and if so overwrite the old parameters.
  # Remove ties, keeping the version with more boundary values
  Results_joint_extreme <- Results_joint_extreme %>% filter(dropped==min(dropped)) %>% filter(logLik==max(logLik,na.rm=T)) %>% tail(1)
  if(Results_joint_extreme$logLik>parms_joint$logLik&Results_joint_extreme$dropped<=parms_joint$dropped){
    parms_joint <- Results_joint_extreme
  }
  # Add empty prod parameters for compatibility with the separate parameters version:
  parms_joint[,paste0(param_names[1:(N_param-2)],"_prod")] <- NA
  parms_joint <- parms_joint %>%
    select(-logLik,-dropped,logLik,dropped)
  # Test separate parameters for comprehension and production:
  Results_separate <- bind_rows(future_lapply(
    init_list,
    function(inits){TobitSeparateModelFit(Production_LogLik_FUN,Comprehension_LogLik_FUN,init_comp=inits,prod_data=prod_data,comp_data=comp_data)}
  ))
  # Keep only the best:
  parms_separate <- Results_separate %>% filter(dropped==min(dropped))%>% filter(logLik==max(logLik,na.rm=T)) %>% head(1)
  # Find candidate limit parameters for comprehension:
  comp_par <- as_vector(parms_separate[1:N_param])
  comp_par[2:3] <- comp_par[2:3]*comp_par[1]
  trans_par <- suppressWarnings(qlogis(comp_par))
  trans_par[c(1:3,(N_param-1):N_param)] <- log(comp_par[c(1:3,(N_param-1):N_param)])
  low_param <- which(round(trans_par) == -7)
  high_param <- which(round(trans_par[2:N_param]) == 7)+1 # don't test infinite lambda
  init_comp_extreme <- NULL
  k=1
  # Define all possible combination of extreme parameters:
  for(i in 1:2^length(low_param)){
    for(j in 1:2^length(high_param)){
      trans_init=trans_par
      trans_init[low_param[as.logical(intToBits(i-1))[1:length(low_param)]]] <- -Inf
      trans_init[high_param[as.logical(intToBits(j-1))[1:length(high_param)]]] <- +Inf
      init <- plogis(trans_init)
      init[c(1:3,(N_param-1):N_param)] <- exp(trans_init[c(1:3,(N_param-1):N_param)])
      if(init[1]>0){init[2:3] <- init[2:3]/init[1]}
      names(init) <- names(parms_joint)[1:N_param]
      init_comp_extreme[[k]] <- init
      k <- k+1
    }
  }
  # Find candidate limit parameters for production:
  prod_par <- as_vector(parms_separate[1:(N_param-2)+N_param])
  prod_par[2:3] <- prod_par[2:3]*prod_par[1]
  trans_par <- suppressWarnings(qlogis(prod_par))
  trans_par[1:3] <- log(prod_par[1:3])
  low_param <- which(round(trans_par) == -7)
  high_param <- which(round(trans_par[2:(N_param-2)]) == 7)+1 # don't test infinite lambda
  init_prod_extreme <- NULL
  k=1
  # Define all possible combinations of extreme parameters:
  for(i in 1:2^length(low_param)){
    for(j in 1:2^length(high_param)){
      trans_init=trans_par
      trans_init[low_param[as.logical(intToBits(i-1))[1:length(low_param)]]] <- -Inf
      trans_init[high_param[as.logical(intToBits(j-1))[1:length(high_param)]]] <- +Inf
      init <- plogis(trans_init)
      init[1:3] <- exp(trans_init[1:3])
      if(init[1]>0){init[2:3] <- init[2:3]/init[1]}
      init <- c(init,as_vector(parms_separate[(N_param-1):N_param]))
      init_prod_extreme[[k]] <- init
      k <- k+1
    }
  }
  # Combine all possible extreme values for production and comprehension:
  init_extreme <- NULL
  k=1
  for(i in 1:length(init_comp_extreme)){
    for(j in 1:length(init_prod_extreme)){
      init_extreme[[k]] <- list(init_comp_extreme[[i]],init_prod_extreme[[j]])
      k=k+1
    }
  }
  # Fit each:
  Results_separate_extreme <- bind_rows(future_lapply(
    init_extreme,
    function(inits){TobitSeparateModelFit(Production_LogLik_FUN,Comprehension_LogLik_FUN,init_comp=inits[[1]],init_prod = inits[[2]],prod_data=prod_data,comp_data=comp_data)}
  ))
  # See whether the extreme values actually improves log-likelihood, and if so overwrite the old parameters.
  Results_separate_extreme <- Results_separate_extreme %>% filter(dropped==min(dropped)) %>% filter(logLik==max(logLik,na.rm=T)) %>% tail(1)
  if(Results_separate_extreme$logLik>parms_separate$logLik&Results_separate_extreme$dropped<=parms_separate$dropped){
    parms_separate <- Results_separate_extreme
  }
  # Adjust names:
  names(parms_separate)[1:N_param] <- names(parms_joint)[1:N_param]
  parms_separate <- parms_separate %>%
    rename(eps=eps_prod) %>%
    select(-ends_with("prod"),-starts_with("logLik"),-starts_with("dropped"),ends_with("prod")&!matches("dropped_prod|logLik_prod"),logLik,dropped) %>%
    rename_with(~str_remove(., '_comp'))
  # Combine the results of the two submodels, compute AIC and order
  parms <- bind_rows(parms_joint,parms_separate,.id = "submodel") %>%
    mutate(model=model_name,
           submodel=c("joint","separate"),
           AIC = 2*dropped+2*rowSums(!is.na(.[2:(2*N_param)]))-2*logLik, # All and only meaningful parameters are non-na.
           BIC = log(nrow(prod_data)+nrow(comp_data)-dropped)*(dropped+rowSums(!is.na(.[2:(2*N_param)])))-2*logLik,
    ) %>%
    arrange(BIC) %>%
    select(model,everything())
  return(parms)
}
