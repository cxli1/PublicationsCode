ODYSSEY_ppfxn <- function(boot=1000, #Number of bootstraps
                          seed=252, #To set seed for the bootstraps replication
                          st, #time of run for saving data
                          data, #Data set to use
                          knots1 = c(0.33, 0.67), #knots for continuous vars
                          knots2 = c(0.5), #knots for VL
                          n_mod, #Covariates from data to put in numerator model for stablized weights
                          g_mod, #Covariates from data to put in PS model
                          truncate = FALSE, #if TRUE, truncate IPCW to truncate_max
                          truncate_max = 25,
                          estimator #Specify either "gcomp", "IPW", "TMLE", or "all"
) {
  set.seed(seed)
  if (boot==0) {B=1}
  else {B = boot}
  b_rand_seed <- sample(1:100000000, B, replace=FALSE)
  #Create data frames to store estimates in
  estimates <- data.frame(y_a1c0 = rep(NA, B),
                          se1 = rep(NA, B),
                          y_a0c0 = rep(NA, B),
                          se0 = rep(NA, B),
                          ate = rep(NA, B),
                          y_a1 = rep(NA, B),
                          y_a0 = rep(NA, B),
                          itt = rep(NA, B))
  estimates_gcomp <- estimates_ipw <- estimates_aipw <- estimates_tmle <- estimates
  wt_dist <- list(ipcw = data.frame(matrix(NA, nrow=B, ncol=6)),
                  s_ipcw = data.frame(matrix(NA, nrow=B, ncol=6))
  )
  #1. Bootstrap data B times
  for (b in 1:B) {
    print(paste0("Starting Bootstrap...", b))
    set.seed(b_rand_seed[b]) #3379785
    print(paste0("Seed set to...", b_rand_seed[b]))

    #Bootstrap: clustered bootstrap by trialno, sample the indices 1 to n with replacement, stratified by arm
    data_a1 <- data[data$tot_rx=="DTG",] #subset data to A=1
    data_a0 <- data[data$tot_rx=="SOC",] #subset data to A=0
    n_a1 <- data_a1 %>% select(TRIALNO) %>% distinct() %>% nrow() #n unique IDs for A=DTG
    n_a0 <- data_a0 %>% select(TRIALNO) %>% distinct() %>% nrow() #n unique IDs for A=SOC

    #create dataset of selected individuals from boosttrap resample b for A=1
    nested_data_a1 <- data_a1 %>% group_by(TRIALNO) %>% nest()
    bootIndices_a1 <- sample(1:n_a1, replace=T)
    bootA1 <- nested_data_a1[bootIndices_a1,]
    bootA1 <- bootA1 %>% ungroup(TRIALNO) %>% mutate(ID = seq(1:n_a1))
    bootA1 <- bootA1 %>% group_by(TRIALNO)%>% unnest(cols = c(data))

    #create dataset of selected individuals from boosttrap resample b for A=0
    nested_data_a0 <- data_a0 %>% group_by(TRIALNO) %>% nest()
    bootIndices_a0 <- sample(1:n_a0, replace=T)
    bootA0 <- nested_data_a0[bootIndices_a0,]
    bootA0 <- bootA0 %>% ungroup(TRIALNO) %>% mutate(ID = seq(1:n_a0)+n_a1)
    bootA0 <- bootA0 %>% group_by(TRIALNO)%>% unnest(cols = c(data))

    longBoot <- rbind(bootA1, bootA0) #bind together bootstrapped dataset

    if (boot==0) {
      bootA1 <- data[data$tot_rx=="DTG",] %>% mutate(ID = TRIALNO)
      bootA0 <- data[data$tot_rx=="SOC",] %>% mutate(ID = TRIALNO)
      longBoot <- data %>% mutate(ID = TRIALNO)
    }

    bootNaive <- longBoot %>% filter(visweek==12)
    km_naive <- survfit(Surv(t_wks96, y) ~ tot_rx, data = bootNaive)
    km_naive_dat <- fortify(km_naive)
    km_naive_dat$risk <- 1 - km_naive_dat$surv
    r_itt_DTG <- max(km_naive_dat$risk[km_naive_dat$strata == "DTG"])
    r_itt_SOC <- max(km_naive_dat$risk[km_naive_dat$strata == "SOC"])
    itt <- r_itt_DTG - r_itt_SOC

    if (estimator != "IPW"){
      #We need to convert longBoot to wideBoot
      dat <- longBoot %>% ungroup() %>% group_by(ID)

      #We need to handcode the interaction term with Cohort A/B (for baseline cov only here)
      dat <- dat %>% mutate(cohort_sex = interaction(tot_ab, sex),
                            cohort_mode = interaction(tot_ab, mode_3cat),
                            cohort_country = interaction(tot_ab, country_5cat))

      # Make confounders wide
      conf_cd4 <- dat %>%
        select(ID, t_end, lag1_cd4) %>%
        pivot_wider(names_from=t_end, names_prefix="Z1_", values_from=lag1_cd4)
      #Create splines for tv cd4
      conf_cd4_spline <- conf_cd4
      for(i in seq(from = 12, to = 96, by = 12)) {
        conf_cd4_spline <- add_spline(data = conf_cd4_spline,
                                      var = paste("Z1", i, sep = "_"),
                                      names = c(paste("Z1", i, 1, sep = "_"),
                                                paste("Z1", i, 2, sep = "_")),
                                      df = 2)}
      #Add cohort interaction terms for tv cd4
      cohort <- dat %>% filter(t_end==12) %>% select(ID, tot_ab)
      conf_cd4_spline_int <- conf_cd4_spline %>% left_join(cohort, by = "ID") %>%
        mutate(across(starts_with('Z'), ~case_when(tot_ab == "A" ~ 0, TRUE ~ .))) %>% select(-tot_ab) %>%
        rename_with(~paste0(.x, "_cohort"), .cols = starts_with("Z"))

      conf_vl <- dat %>%
        select(ID, t_end, log10vl) %>%
        pivot_wider(names_from=t_end, names_prefix="Z2_", values_from=log10vl)
      #Create splines for tv vl
      conf_vl_spline <- conf_vl
      for(i in seq(from = 12, to = 96, by = 12)) {
        conf_vl_spline <- add_spline(data = conf_vl_spline,
                                     var = paste("Z2", i, sep = "_"),
                                     names = c(paste("Z2", i, 1, sep = "_")
                                               #paste("Z2", i, 2, sep = "_")
                                               ),
                                     df = 1)}
      #Add cohort interaction terms for tv vl
      conf_vl_spline_int <- conf_vl_spline %>% left_join(cohort, by = "ID") %>%
        mutate(across(starts_with('Z'), ~case_when(tot_ab == "A" ~ 0, TRUE ~ .))) %>% select(-tot_ab) %>%
        rename_with(~paste0(.x, "_cohort"), .cols = starts_with("Z"))

      # Make exposure wide
      exposure <- dat %>%
        select(ID, t_end, tot_rx) %>%
        mutate(tot_rx01 = ifelse(tot_rx=="DTG", 1, 0)) %>%
        pivot_wider(names_from=t_end, names_prefix="A", values_from=tot_rx01) %>% select(-tot_rx)

      # Make censoring wide (for lmtp c=0 is censored...)
      censor <- dat %>%
        mutate(C = abs(c-1)) %>%
        select(ID, t_end, C) %>%
        pivot_wider(names_from=t_end, names_prefix="C_", values_from=C)

      # Make outcome wide (once Y=1, needs to be carried forward)
      outcome <- dat %>%
        select(ID, t_end, y_obs) %>%
        pivot_wider(names_from=t_end, names_prefix="Y", values_from=y_obs) %>%
        mutate(Y24 = ifelse(Y12==1, 1, Y24),
               Y36 = ifelse(Y24==1, 1, Y36),
               Y48 = ifelse(Y36==1, 1, Y48),
               Y60 = ifelse(Y48==1, 1, Y60),
               Y72 = ifelse(Y60==1, 1, Y72),
               Y84 = ifelse(Y72==1, 1, Y84),
               Y96 = ifelse(Y84==1, 1, Y96))
      # Interleave together in correct order (baseline, Z1, Z2, X, D, Y)
      baseline <- dat %>% filter(t_end==12) %>% select(ID, tot_ab, sex, cohort_sex, mode_3cat, cohort_mode, country_5cat, cohort_country,
                                                       bage_v2, bweight, bheight, bmcd4, log10_bvl_final)
      wide <-baseline %>%
        mutate(across(!starts_with("b"), as.factor))

      #Make splines for baseline age
      baseline <- add_spline(data = baseline,
                             var = "bage_v2",
                             names = c("bage_1", "bage_2"),
                             df = 2)
      #Make interaction terms with cohort for baseline age splines
      baseline <- baseline %>% mutate(cohort_bage_1 = ifelse(tot_ab=="A", 0, bage_1),
                                      cohort_bage_2 = ifelse(tot_ab=="A", 0, bage_2))

      #Make splines for baseline height
      baseline <- add_spline(data = baseline,
                             var = "bheight",
                             names = c("bheight_1", "bheight_2"),
                             df = 2)
      #Make interaction terms with cohort for baseline height splines
      baseline <- baseline %>% mutate(cohort_bheight_1 = ifelse(tot_ab=="A", 0, bheight_1),
                                      cohort_bheight_2 = ifelse(tot_ab=="A", 0, bheight_2))

      #Make splines for baseline weight
      baseline <- add_spline(data = baseline,
                             var = "bweight",
                             names = c("bweight_1", "bweight_2"),
                             df = 2)
      #Make interaction terms with cohort for baseline weight splines
      baseline <- baseline %>% mutate(cohort_bweight_1 = ifelse(tot_ab=="A", 0, bweight_1),
                                      cohort_bweight_2 = ifelse(tot_ab=="A", 0, bweight_2))

      #Make splines for baseline cd4 splines
      baseline <- add_spline(data = baseline,
                             var = "bmcd4",
                             names = c("bmcd4_1", "bmcd4_2"),
                             df = 2)
      #Make interaction terms with cohort for baseline cd4 splines
      baseline <- baseline %>% mutate(cohort_bmcd4_1 = ifelse(tot_ab=="A", 0, bmcd4_1),
                                      cohort_bmcd4_2 = ifelse(tot_ab=="A", 0, bmcd4_2))

      #Make splines for baseline log10 vl splines
      baseline <- add_spline(data = baseline,
                             var = "log10_bvl_final",
                             names = c("bmvl_1", "bmvl_2"),
                             df = 2)
      #Make interaction terms with cohort for baseline log10vl splines
      baseline <- baseline %>% mutate(cohort_bmvl_1 = ifelse(tot_ab=="A", 0, bmvl_1),
                                      cohort_bmvl_2 = ifelse(tot_ab=="A", 0, bmvl_2))

      wideBoot <-baseline %>%
        mutate(across(!starts_with("b"), as.factor))

      for (i in 1:8){
        wideBoot <- merge(wideBoot, exposure[ , c(1, i+1)], by="ID")
        wideBoot <- merge(wideBoot, conf_cd4_spline[ , c(1, i*2, i*2+1)], by="ID")
        wideBoot <- merge(wideBoot, conf_cd4_spline_int[ , c(1, i*2, i*2+1)], by="ID")
        wideBoot <- merge(wideBoot, conf_vl_spline[ , c(1, i+1)], by="ID")
        wideBoot <- merge(wideBoot, conf_vl_spline_int[ , c(1, i+1)], by="ID")
        wideBoot <- merge(wideBoot, censor[ , c(1, i+1)], by="ID")
        wideBoot <- merge(wideBoot, outcome[ , c(1, i+1)], by="ID")
      }
      message(c("rows wide:   ", nrow(wideBoot)))
      message(c("event:  ", sum(wideBoot$Y96, na.rm = T)))
      message(c("protocol deviation:   ", sum(colSums(abs(censor[,2:9]-1), na.rm = T))))
    }

    message(c("rows long:   ", nrow(longBoot)))
    message(c("event:  ", sum(longBoot$delta_t)))
    message(c("protocol deviation:   ", sum(1-longBoot$nltfu_t)))


    if (estimator=="IPW"|estimator=="all") {
      #first save naive
      estimates_ipw[b,6] <- r_itt_DTG
      estimates_ipw[b,7] <- r_itt_SOC
      estimates_ipw[b,8] <- itt

      # Fitting pooled logistic model, stratified by DTG and SOC
      #For DTG:
      #Specifying knot points for splines
      knots_age <- quantile(bootA1$bage_v2[bootA1$y_obs==0], probs = knots1, na.rm=T)
      knots_age_boundary <- quantile(bootA1$bage_v2[bootA1$y_obs==0], probs = c(0.05, 0.95), na.rm=T)
      knots_weight <- quantile(bootA1$bweight[bootA1$y_obs==0], probs = knots1, na.rm=T)
      knots_weight_boundary <- quantile(bootA1$bweight[bootA1$y_obs==0], probs = c(0.05, 0.95), na.rm=T)
      knots_height <- quantile(bootA1$bheight[bootA1$y_obs==0], probs = knots1, na.rm=T)
      knots_height_boundary <- quantile(bootA1$bheight[bootA1$y_obs==0], probs = c(0.05, 0.95), na.rm=T)

      knots_vl_1 <- quantile(bootA1$lag1_log10vl[bootA1$y_obs==0], probs = knots2, na.rm=T)
      knots_vl_1_boundary <- quantile(bootA1$lag1_log10vl[bootA1$y_obs==0], probs = c(0.05, 0.95), na.rm=T)

      knots_vl_2 <- quantile(bootA1$lag2_log10vl[bootA1$y_obs==0], probs = knots2, na.rm=T)
      knots_vl_2_boundary <- quantile(bootA1$lag2_log10vl[bootA1$y_obs==0], probs = c(0.05, 0.95), na.rm=T)

      knots_cd4 <- quantile(bootA1$cd4[bootA1$y_obs==0], probs = knots1, na.rm=T)
      knots_cd4_boundary <- quantile(bootA1$cd4[bootA1$y_obs==0], probs = c(0.05, 0.95), na.rm=T)
      knots_cd4_1 <- quantile(bootA1$lag1_cd4[bootA1$y_obs==0], probs = knots1, na.rm=T)
      knots_cd4_1_boundary <- quantile(bootA1$lag1_cd4[bootA1$y_obs==0], probs = c(0.05, 0.95), na.rm=T)
      knots_cd4_2 <- quantile(bootA1$lag2_cd4[bootA1$y_obs==0], probs = knots1, na.rm=T)
      knots_cd4_2_boundary <- quantile(bootA1$lag2_cd4[bootA1$y_obs==0], probs = c(0.05, 0.95), na.rm=T)

      assign("knots_age", knots_age, envir = .GlobalEnv)
      assign("knots_weight", knots_weight, envir = .GlobalEnv)
      assign("knots_height", knots_height, envir = .GlobalEnv)
      assign("knots_vl_1", knots_vl_1, envir = .GlobalEnv)
      assign("knots_vl_2", knots_vl_2, envir = .GlobalEnv)
      assign("knots_cd4", knots_cd4, envir = .GlobalEnv)
      assign("knots_cd4_1", knots_cd4_1, envir = .GlobalEnv)
      assign("knots_cd4_2", knots_cd4_2, envir = .GlobalEnv)

      assign("knots_age_boundary", knots_age_boundary, envir = .GlobalEnv)
      assign("knots_weight_boundary", knots_weight_boundary, envir = .GlobalEnv)
      assign("knots_height_boundary", knots_height_boundary, envir = .GlobalEnv)
      assign("knots_vl_1_boundary", knots_vl_1_boundary, envir = .GlobalEnv)
      assign("knots_vl_2_boundary", knots_vl_2_boundary, envir = .GlobalEnv)
      assign("knots_cd4_boundary", knots_cd4_boundary, envir = .GlobalEnv)
      assign("knots_cd4_1_boundary", knots_cd4_1_boundary, envir = .GlobalEnv)
      assign("knots_cd4_2_boundary", knots_cd4_2_boundary, envir = .GlobalEnv)


      logit_den_a1 <- glm(g_mod,
                          data=bootA1,
                          #data=bootA1 %>% filter(y_obs==0),
                          family='binomial'(link = "logit"))
      logit_num_a1 <- glm(n_mod,
                          data=bootA1,
                          #data=bootA1 %>% filter(y_obs==0),
                          family='binomial'(link = "logit"))
      bootA1$pr_c_den <- predict(logit_den_a1, bootA1, type='response') # predicted probabilities
      bootA1$pr_c_num <- predict(logit_num_a1, bootA1, type='response') # predicted probabilities

      #For SOC:
      #Specifying knot points for splines
      knots_age <- quantile(bootA0$bage_v2[bootA0$y_obs==0], probs = knots1, na.rm=T)
      knots_age_boundary <- quantile(bootA0$bage_v2[bootA0$y_obs==0], probs = c(0.05, 0.95), na.rm=T)
      knots_weight <- quantile(bootA0$bweight[bootA0$y_obs==0], probs = knots1, na.rm=T)
      knots_weight_boundary <- quantile(bootA0$bweight[bootA0$y_obs==0], probs = c(0.05, 0.95), na.rm=T)
      knots_height <- quantile(bootA0$bheight[bootA0$y_obs==0], probs = knots1, na.rm=T)
      knots_height_boundary <- quantile(bootA0$bheight[bootA0$y_obs==0], probs = c(0.05, 0.95), na.rm=T)

      knots_vl_1 <- quantile(bootA0$lag1_log10vl[bootA0$y_obs==0], probs = knots2, na.rm=T)
      knots_vl_1_boundary <- quantile(bootA0$lag1_log10vl[bootA0$y_obs==0], probs = c(0.05, 0.95), na.rm=T)

      knots_vl_2 <- quantile(bootA0$lag2_log10vl[bootA0$y_obs==0], probs = knots2, na.rm=T)
      knots_vl_2_boundary <- quantile(bootA0$lag2_log10vl[bootA0$y_obs==0], probs = c(0.05, 0.95), na.rm=T)

      knots_cd4 <- quantile(bootA0$cd4[bootA0$y_obs==0], probs = knots1, na.rm=T)
      knots_cd4_boundary <- quantile(bootA0$cd4[bootA0$y_obs==0], probs = c(0.05, 0.95), na.rm=T)
      knots_cd4_1 <- quantile(bootA0$lag1_cd4[bootA0$y_obs==0], probs = knots1, na.rm=T)
      knots_cd4_1_boundary <- quantile(bootA0$lag1_cd4[bootA0$y_obs==0], probs = c(0.05, 0.95), na.rm=T)
      knots_cd4_2 <- quantile(bootA0$lag2_cd4[bootA0$y_obs==0], probs = knots1, na.rm=T)
      knots_cd4_2_boundary <- quantile(bootA0$lag2_cd4[bootA0$y_obs==0], probs = c(0.05, 0.95), na.rm=T)

      assign("knots_age", knots_age, envir = .GlobalEnv)
      assign("knots_weight", knots_weight, envir = .GlobalEnv)
      assign("knots_height", knots_height, envir = .GlobalEnv)
      assign("knots_vl_1", knots_vl_1, envir = .GlobalEnv)
      assign("knots_vl_2", knots_vl_2, envir = .GlobalEnv)
      assign("knots_cd4", knots_cd4, envir = .GlobalEnv)
      assign("knots_cd4_1", knots_cd4_1, envir = .GlobalEnv)
      assign("knots_cd4_2", knots_cd4_2, envir = .GlobalEnv)

      assign("knots_age_boundary", knots_age_boundary, envir = .GlobalEnv)
      assign("knots_weight_boundary", knots_weight_boundary, envir = .GlobalEnv)
      assign("knots_height_boundary", knots_height_boundary, envir = .GlobalEnv)
      assign("knots_vl_1_boundary", knots_vl_1_boundary, envir = .GlobalEnv)
      assign("knots_vl_2_boundary", knots_vl_2_boundary, envir = .GlobalEnv)
      assign("knots_cd4_boundary", knots_cd4_boundary, envir = .GlobalEnv)
      assign("knots_cd4_1_boundary", knots_cd4_1_boundary, envir = .GlobalEnv)
      assign("knots_cd4_2_boundary", knots_cd4_2_boundary, envir = .GlobalEnv)

      logit_den_a0 <- glm(g_mod,
                          data=bootA0,
                          family='binomial'(link = "logit"))
      logit_num_a0 <- glm(n_mod,
                          data=bootA0,
                          family='binomial'(link = "logit"))
      bootA0$pr_c_den <- predict(logit_den_a0, bootA0, type='response') # predicted probabilities
      bootA0$pr_c_num <- predict(logit_num_a0, bootA0, type='response') # predicted probabilities

      bootData <- rbind(bootA1, bootA0)

      bootData <- bootData %>%
        group_by(ID) %>%
        mutate(ipcw = cumprod(1 / pr_c_den),           # Inverting to obtain weights
               s_ipcw = cumprod(pr_c_num / pr_c_den))

      # truncate weights to tmax if truncate=TRUE
      if (truncate==TRUE) {
        max <- truncate_max
        min <- 1/truncate_max
        bootData <- bootData %>% mutate(ipcw_t = case_when(s_ipcw > max ~ max,
                                                           s_ipcw < min ~ min,
                                                           .default = s_ipcw))
        bootData$s_ipcw <- bootData$ipcw_t
      }

      wt_dist[[1]][b,] <- summary(bootData$ipcw)[1:6] #IPCW
      wt_dist[[2]][b,] <- summary(bootData$s_ipcw)[1:6] # stablilized IPCW

      # Estimating a IPCW Kaplan-Meier
      km_ipcw <- survfit(Surv(t_start, t_end, y_obs) ~ tot_rx, data=bootData %>% filter(nltfu_t==1), weights=bootData$s_ipcw[bootData$nltfu_t==1],
                         robust = TRUE, id = ID)
      km_ipcw_dat <- fortify(km_ipcw)
      km_ipcw_dat$risk <- 1 - km_ipcw_dat$surv

      estimates_ipw[b,1] <- max(km_ipcw_dat$risk[km_ipcw_dat$strata=="DTG"])
      estimates_ipw[b,2] <- km_ipcw_dat$std.err[km_ipcw_dat$strata=="DTG" & km_ipcw_dat$time==96]
      estimates_ipw[b,3] <- max(km_ipcw_dat$risk[km_ipcw_dat$strata=="SOC"])
      estimates_ipw[b,4] <- km_ipcw_dat$std.err[km_ipcw_dat$strata=="SOC" & km_ipcw_dat$time==96]
      estimates_ipw[b,5] <- max(km_ipcw_dat$risk[km_ipcw_dat$strata=="DTG"])-max(km_ipcw_dat$risk[km_ipcw_dat$strata=="SOC"])
    }

    if (estimator=="gcomp"|estimator=="all") {
      #first save naive
      estimates_gcomp[b,6] <- r_itt_DTG
      estimates_gcomp[b,7] <- r_itt_SOC
      estimates_gcomp[b,8] <- itt

      # Set up call to lmtp package
      A <- vars_select(names(wideBoot), starts_with("A"))
      Lnodes <- list(vars_select(names(wideBoot), starts_with("Z") & contains("12")),
                     vars_select(names(wideBoot), starts_with("Z") & contains("24")),
                     vars_select(names(wideBoot), starts_with("Z") & contains("36")),
                     vars_select(names(wideBoot), starts_with("Z") & contains("48")),
                     vars_select(names(wideBoot), starts_with("Z") & contains("60")),
                     vars_select(names(wideBoot), starts_with("Z") & contains("72")),
                     vars_select(names(wideBoot), starts_with("Z") & contains("84")),
                     vars_select(names(wideBoot), starts_with("Z") & contains("96")))

      Ynodes <- vars_select(names(wideBoot), starts_with("Y"))
      Cnodes <- vars_select(names(wideBoot), starts_with("C_"))
      # Use lmtp to implement ICE g-comp
      wide_dtg <- wideBoot %>% filter(A12==1) %>% select(-ID)
      wide_soc <- wideBoot %>% filter(A12==0) %>% select(-ID)

      res_dtg <- tryCatch({lmtp_sub(data=wide_dtg,
                                    trt = A, outcome = Ynodes,
                                    time_vary = Lnodes, cens = Cnodes,
                                    shift = NULL, k=2, folds = 1,
                                    outcome_type = "survival")},
                          error=function(error_message) {
                            message("Error in DTG arm, return NA.")
                            message("Error message from R:")
                            message(error_message)
                            return(list(theta = NA,
                                        standard_error = NA))})

      res_soc <- tryCatch({lmtp_sub(data=wide_soc,
                                    trt = A, outcome = Ynodes,
                                    time_vary = Lnodes, cens = Cnodes,
                                    shift = NULL, k=2, folds = 1,
                                    outcome_type = "survival")},
                          error=function(error_message) {
                            message("Error in SOC arm, return NA.")
                            message("Error message from R:")
                            message(error_message)
                            return(list(theta = NA,
                                        standard_error = NA))})

      estimates_gcomp[b,1] <- 1-res_dtg$theta
      estimates_gcomp[b,2] <- res_dtg$standard_error
      estimates_gcomp[b,3] <- 1-res_soc$theta
      estimates_gcomp[b,4] <- res_soc$standard_error
      estimates_gcomp[b,5] <- res_soc$theta-res_dtg$theta
    }

    if (estimator=="TMLE"|estimator=="all") {
      #first save naive
      estimates_tmle[b,6] <- r_itt_DTG
      estimates_tmle[b,7] <- r_itt_SOC
      estimates_tmle[b,8] <- itt

      # Set up call to lmtp package
      A <- vars_select(names(wideBoot), starts_with("A"))
      Lnodes <- list(vars_select(names(wideBoot), starts_with("Z") & contains("12")),
                     vars_select(names(wideBoot), starts_with("Z") & contains("24")),
                     vars_select(names(wideBoot), starts_with("Z") & contains("36")),
                     vars_select(names(wideBoot), starts_with("Z") & contains("48")),
                     vars_select(names(wideBoot), starts_with("Z") & contains("60")),
                     vars_select(names(wideBoot), starts_with("Z") & contains("72")),
                     vars_select(names(wideBoot), starts_with("Z") & contains("84")),
                     vars_select(names(wideBoot), starts_with("Z") & contains("96")))

      Ynodes <- vars_select(names(wideBoot), starts_with("Y"))
      Cnodes <- vars_select(names(wideBoot), starts_with("C_"))
      # Use lmtp to implement ICE g-comp
      wide_dtg <- wideBoot %>% filter(A12==1) %>% select(-ID)
      wide_soc <- wideBoot %>% filter(A12==0) %>% select(-ID)

      tmle_dtg <- tryCatch({lmtp_tmle(data=wide_dtg,
                                    trt = A, outcome = Ynodes,
                                    time_vary = Lnodes, cens = Cnodes,
                                    shift = NULL, k=2, folds = 1,
                                    outcome_type = "survival")},
                          error=function(error_message) {
                            message("Error in DTG arm, return NA.")
                            message("Error message from R:")
                            message(error_message)
                            return(list(theta = NA,
                                        standard_error = NA))})

      tmle_soc <- tryCatch({lmtp_tmle(data=wide_soc,
                                    trt = A, outcome = Ynodes,
                                    time_vary = Lnodes, cens = Cnodes,
                                    shift = NULL, k=2, folds = 1,
                                    outcome_type = "survival")},
                          error=function(error_message) {
                            message("Error in SOC arm, return NA.")
                            message("Error message from R:")
                            message(error_message)
                            return(list(theta = NA,
                                        standard_error = NA))})

      estimates_tmle[b,1] <- 1-tmle_dtg$theta
      estimates_tmle[b,2] <- tmle_dtg$standard_error
      estimates_tmle[b,3] <- 1-tmle_soc$theta
      estimates_tmle[b,4] <- tmle_soc$standard_error
      estimates_tmle[b,5] <- tmle_soc$theta-tmle_dtg$theta
    }
  }

  #Output dataframes with results (as lists) according to specified estimators

  if (estimator=="gcomp") {
    gcomp <- list(estimates_gcomp, wt_dist, b_rand_seed)
    names(gcomp) <- c("G-computation", "Weight Distributions", "Seeds")
    return(gcomp)
  }
  if (estimator=="IPW") {
    if(boot==0) {ipw <- list(estimates_ipw, wt_dist, b_rand_seed, km_ipcw_dat)
    names(ipw) <- c("IPW", "Weight Distributions", "Seeds", "km_dat")}
    else {ipw <- list(estimates_ipw, wt_dist, b_rand_seed)
    names(ipw) <- c("IPW", "Weight Distributions", "Seeds")}
    return(ipw)
  }
  if (estimator=="AIPW") {
    aipw <- list(estimates_aipw, wt_dist, b_rand_seed)
    names(aipw) <- c("AIPW", "Weight Distributions", "Seeds")
    return(aipw)
  }
  if (estimator=="TMLE") {
    tmle <- list(estimates_tmle, wt_dist, b_rand_seed)
    names(tmle) <- c("TMLE", "Weight Distributions", "Seeds")
    return(tmle)
  }
  if (estimator=="all") {
    all <- list(estimates_gcomp, estimates_ipw, estimates_tmle, wt_dist, b_rand_seed)
    names(all) <- c("gcomp","IPW", "TMLE", "Weight Distributions", "Seeds")
    return(all)
  }
}
