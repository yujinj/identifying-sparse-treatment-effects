
library(readr)
library(glmnet)
library(Hotelling)

cgm_avg <- read_csv("/Users/yujin/Downloads/csv_data/cgm_avg.csv")

get_tir_windows <- function(glucose, window=24, gmax=180, gmin=100, fortau=FALSE){
  m = length(glucose)/window
  tir = rep(0, 2*m)
  for(i in 1:m){
    tir[(2*i-1)] = mean(glucose[ceiling((1:length(glucose))/window) == i] >= gmax, na.rm = TRUE) # hyperglycemia
    tir[(2*i)] = mean(glucose[ceiling((1:length(glucose))/window) == i] <= gmin, na.rm = TRUE) # hypoglycemia
  }
  if(fortau){return(tir[(m+1):(2*m)])}
  return(tir)
}

# Matching with two-hour time intervals

tir = t(apply(cgm_avg, 1, function(glucose, window = 24){
                                      m = (ncol(cgm_avg)/window)
                                      tir = rep(0, m)
                                      for(i in 1:m){
                                        tir[i] = mean((glucose[ceiling((1:ncol(cgm_avg))/window) == i] <= 180) &
                                                        (glucose[ceiling((1:ncol(cgm_avg))/window) == i] >= 100) , na.rm = TRUE) 
                                      }
                                      return(tir)
                                    }))

library(nbpMatching)
df <- data.frame(tir[,1:(ncol(tir)/2)])
df.dist <- gendistance(df)
df.mdm <- distancematrix(df.dist)
df.match <- nonbimatch(df.mdm)

##########################################################################
#######################           LASSO           ########################
##########################################################################


single_split <- function(tir_synth, treatment){
  
  n = nrow(tir_synth)
  m = as.integer(ncol(tir_synth)/2)
  
  # Fit LASSO
  X_centered = scale(tir_synth[1:n, 1:m] - tir_synth[1:n, (m+1):(2*m)], center = TRUE, scale = FALSE)
  y_centered = treatment-mean(treatment)
  glmfit = glmnet(x = X_centered, 
                  y = y_centered, 
                  family = "gaussian",
                  alpha = 1, 
                  nlambda = 200,
                  intercept = FALSE, 
                  standardize = FALSE)
  coef_list = predict(glmfit, type = "nonzero")
  for (i in 1:length(coef_list)){
    selected_outcomes = coef_list[[i]]
    if (length(selected_outcomes) > 0) {
      break
    }
  }
  lm.fit = lm(y_centered~X_centered[1:n, selected_outcomes]+0)
  rss = mean(summary(lm.fit)$residuals^2)
  return(list('selected_outcomes'=selected_outcomes, 'score'=rss))
}

multi_split <- function(tir_synth_list, treatment, B=100, etol=2e-4) {
  
  n = nrow(tir_synth_list[[1]])
  p.vals_list = list() # for each level
  for (i in 1:length(tir_synth_list)){
    p.vals_list[[i]] = rep(1, B)
  }
  
  for (i_b in 1:B){
    set.seed(i_b)
    w = (sample.int(n, size = n, replace = FALSE) %% 10 >= 4)
    
    # First Stage: Performing Lasso
    min_level = -1; min_score = 10000
    for (i in 1:length(tir_synth_list)){
      tir_synth = tir_synth_list[[i]]
      subset_res = single_split(tir_synth[w==0, ], treatment[w==0])
      if (subset_res$score < (min_score+etol)) {
        min_level = i
        min_score = subset_res$score
        selected_outcomes = subset_res$selected_outcomes
      }
    }
    tir_synth = tir_synth_list[[min_level]]
    
    # Second Stage: Calculating Treatment Effect
    m = as.integer(ncol(tir_synth)/2)
    y1 = tir_synth[w==1, ][treatment[w==1] == 1, (selected_outcomes+m)] - tir_synth[w==1, ][treatment[w==1] == 1, (selected_outcomes)]
    y0 = tir_synth[w==1, ][treatment[w==1] == 0, (selected_outcomes+m)] - tir_synth[w==1, ][treatment[w==1] == 0, (selected_outcomes)]
    if (length(selected_outcomes)==1){
      p.vals_list[[min_level]][i_b] = t.test(y1, y0)$p.value
    } else {
      p.vals_list[[min_level]][i_b] = Hotelling::hotelling.test(y1, y0)$pval
    }
  }
  adjusted_p.val <- min(sapply(p.vals_list, function(pv) {min(quantile(pv, 0.1) * 10, 1)}))
  
  return(adjusted_p.val)
}


##########################################################################
#################           Multiple Testing           ###################
##########################################################################


multiple_testing <- function(tir_synth, treatment){
  m = as.integer(ncol(tir_synth)/2)
  pvals = rep(0, m)
  for(i in 1:m){
    d1 = tir_synth[treatment == 1,(m+i)]-tir_synth[treatment == 1,i]
    d0 = tir_synth[treatment == 0,(m+i)]-tir_synth[treatment == 0,i]
    if (var(d1) > 0 && var(d0) > 0){
      pvals[i] = min(t.test(d1, d0)$p.value*m, 1) 
    } else {
      pvals[i] = 1
    }
  }
  return(pvals)
}


##########################################################################
#################           Simultaneous Inf           ###################
##########################################################################


psi <- function(tir_synth, treatment, B=1000, seed=0) {
  
  m <- as.integer(ncol(tir_synth)/2)
  n <- nrow(tir_synth)
  Y_tilde = tir_synth[,((m+1):(2*m))]-tir_synth[,(1:m)]
  psi = treatment * Y_tilde / mean(treatment) - (1-treatment) * Y_tilde / mean(1-treatment)
  
  tau_hat <- colMeans(psi)                    
  sd_psi  <- apply(psi, 2, sd)
  sd_psi[!is.finite(sd_psi) | sd_psi == 0] <- Inf
  se_hat  <- sd_psi / sqrt(n)
  T_obs   <- tau_hat / se_hat 
  
  psi_centered <- scale(psi, center = tau_hat, scale = FALSE) 
  
  # multiplier bootstrap 
  set.seed(seed)
  E <- matrix(rnorm(n * B), nrow = n, ncol = B)
  boot_raw <- (t(psi_centered) %*% E) / sqrt(n)
  boot_t <- t( scale(t(boot_raw), center = FALSE, scale = sd_psi))
  T_star <- apply(abs(boot_t), 2, max)
  c_alpha <- as.numeric(quantile(T_star, 0.95))
  
  return(which(abs(T_obs) > c_alpha))
}

##########################################################################
####################           Simulations          ######################
##########################################################################

n <- nrow(cgm_avg)
alpha_seq = seq(from = 5, to = 20, by = 0.5)
power_results = matrix(0, nrow = length(alpha_seq), ncol = 3) 

tau_size=6 # treatment effect window size: 30 minute

for (seed in 1:500){
  
  print(seed)
  set.seed(seed)
  # random treatment assignments
  treatment = as.numeric(assign.grp(df.match, seed = seed)[,6] == "A") #rbinom(n, size = 1, prob = 0.5)
  # random locations of treatment effects
  tau_start = sample.int(288-tau_size-1, size = 1)
  skip1 = FALSE; skip2 = FALSE; skip3 = FALSE
  
  for(a in 1:length(alpha_seq)){ # magnitudes of treatment effects
    
    if (!skip1|!skip2 |!skip3){
    
      # Create Synthetic Treatment Effect
      
      tau = rep(0, 576)
      noise = rep(0, 576)
      
      tau_idx = ((288 + tau_start + 1):(288 + tau_start + tau_size))
      tau[(288+tau_start+1):(288+tau_start+tau_size)] = -alpha_seq[a] # constant treatment effect
      noise = c(rep(0, 288), alpha_seq[a] * runif(288, 0.5, 1) * sample(c(-1, 1), size = 288, replace = TRUE)) # random noise
      
      cgm_synth = cgm_avg
      cgm_synth = t(apply(cgm_synth, 1, function(x){return(x+noise)}))
      cgm_synth[treatment == 1, ] = t(apply(cgm_synth[treatment == 1, ], 1, function(x){return(x+tau)}))
      
      tir_synth_s = t(apply(cgm_synth, 1, function(x){get_tir_windows(x, window = 2)}))
      tir_synth_m = t(apply(cgm_synth, 1, function(x){get_tir_windows(x, window = 6)}))
      tir_synth_l = t(apply(cgm_synth, 1, function(x){get_tir_windows(x, window = 12)}))
      
      sig_dim_s = which(get_tir_windows(-tau, window=2, gmax = 1, gmin=-1, fortau=TRUE)>0)
      sig_dim_m = which(get_tir_windows(-tau, window=6, gmax = 1, gmin=-1, fortau=TRUE)>0)
      sig_dim_l = which(get_tir_windows(-tau, window=12, gmax = 1, gmin=-1, fortau=TRUE)>0)
      
      # 1) Multi-sample splitting with our subset selection method
      
      if(!skip1){
        tir_synth_list = list(tir_synth_l, tir_synth_m, tir_synth_s)
        pval1 = multi_split(tir_synth_list, treatment)
        if (pval1 < 0.05) {
          power_results[a:length(alpha_seq), 1] = power_results[a:length(alpha_seq), 1] + 1
          skip1=TRUE
        }
      }
      
      # 2) Multiple testing with window size 30 min
      
      if (!skip2){
        pvals = multiple_testing(tir_synth_m, treatment)
        ms_selected_set = which(pvals < 0.05)
        if (length(intersect(ms_selected_set, sig_dim_m)) > 0) {
          power_results[a:length(alpha_seq), 2] = power_results[a:length(alpha_seq), 2] + 1
          skip2=TRUE
        }
      } 
      
      # 3) Simultaneous inference with window size 30 min
      
      if (!skip3){
        si_selected_set = psi(tir_synth_m, treatment)
        if (length(intersect(si_selected_set, sig_dim_m)) > 0) {
          power_results[a:length(alpha_seq), 3] = power_results[a:length(alpha_seq), 3] + 1
          skip3=TRUE
        }
      }
      
    }
  } 
}