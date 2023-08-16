import_cv <- function(path_to_files) {
  files <- list.files(path_to_files, full.names = T, pattern = ".tsv")
  dat <- read_tsv(files, id = "path") 
  return(dat)
}


any_convergence_warning <- function(path) {
  patterns <-
    sapply(
      list.files(path, pattern = ".err", full.names = TRUE),
      FUN = function(x) {
        any(grepl("ConvergenceWarning", readLines(x)))
      }
    )
  return(patterns)
}

median_performance <- function(df) {
  
  df %>%
    group_by(path) %>% 
    summarise(
      auc = median(auc),
      binaryaccuracy = median(binaryaccuracy),
      precision = median(precision),
      recall = median(recall)
    )
}


sd_performance <- function(df){
  
  df %>%
    group_by(path) %>% 
    summarise(
      auc = IQR(auc),
      binaryaccuracy = IQR(binaryaccuracy),
      precision = IQR(precision),
      recall = IQR(recall)
    )
  
}

transform_performance <- function(list_df) {
  
  list_median <- median_performance(list_df)
  
  performance <- list_median %>%
    mutate(parameter = gsub(".*miles_","", path)) %>%
    mutate(parameter = gsub(".tsv","", parameter)) %>%
    separate(parameter, into = c("sigma2", "c"), sep = "_")
  
  df_sd <- sd_performance(list_df) 
  colnames(df_sd) <-
    c("path_iqr",
      "auc_iqr",
      "binaryaccuracy_iqr",
      "precision_iqr",
      "recall_iqr")

  performance <- bind_cols(performance, df_sd)
  
  return(performance)
}

get_best_model <- function(df) {
  df1 <- df %>%
    filter(auc == max(auc))%>% 
    mutate(sigma2 = as.numeric(sigma2),
           c = as.numeric(c))
  if (nrow(df1) > 1) {
    df1 <- df1 %>%
      filter(binaryaccuracy == max(binaryaccuracy))
    if (nrow(df1) > 1) {
      df1 <- df1 %>%
        filter(precision  == max(precision))
      if (nrow(df1) > 1) {
        df1 <- df1 %>%
          filter(sigma2  == min(sigma2)) %>%
          filter(c == min(c))
      }
    }
  }
  return(df1)
}


summarise_results <- function(df){
  
  res_raw <- lapply(df$path, import_cv)
  
  names(res_raw) <- df$path
  
  res_transformed <- res_raw %>% 
    map(., transform_performance)
  
  best_models <- res_transformed %>%
    map(., get_best_model)
  
  result_summary <- bind_rows(best_models, .id = "id")
  
  res <- df %>% 
    left_join(result_summary, by = c("path" = "id"))
  
  return(res)
  
}


summarise_warning <- function(df){
  
  warn_raw <- lapply(df$path, any_convergence_warning) 
  
  warn_transformed <- lapply(warn_raw, function(x) tibble(warning = x))
  names(warn_transformed) <- df$path
  
  warn_transformed1 <- warn_transformed %>% map(., function(x) count(x, warning))
  warn_transformed2 <- bind_rows(warn_transformed1, .id = "path")
  warn_transformed3 <- warn_transformed2 %>% 
    group_by(path) %>%
    mutate(frac = n/sum(n)) %>% 
    filter(warning == TRUE)
  
  df1 <- df %>% 
    left_join(warn_transformed3, by = "path") 
  
  return(df1)
  
}


# FEATURE IMPORTANCE ------------------------------------------------------


import_feature_importance <- function(path_to_files){
  files <- list.files(path_to_files, pattern = ".csv", full.names = T)
  nams <- gsub(".csv", "", gsub( "/", "" , gsub(path_to_files, "", files)))
  dat <- files %>%
    map(read_delim) 
  names(dat) <- nams
  return(dat)
}

transform_performance_importance <- function(list_df){
  list_mean <- list_df %>%
    map(., median_performance)
  list_sd <- list_df %>%
    map(., sd_performance)
  performance <- list_mean %>%
    bind_rows(.id = "feature") 
  df_sd <- list_sd %>%
    bind_rows()
  colnames(df_sd) <- c("auc_iqr", "binaryaccuracy_iqr", "precision_iqr", "recall_iqr")
  performance <- bind_cols(performance, df_sd)
  return(performance)
}




error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

transform_mean_auc <- function(list_performance){
  mean_auc <- sapply(list_performance, function(x) median(x$auc))
  params <- str_split(names(mean_auc),"_")
  names(params) <- names(mean_auc)
  params <- lapply(params, function(x) {
    names(x) <- c("sigma2", "c")
    return(x)
  })
  params <- bind_rows(params)
  params$auc <- mean_auc
  return(params)
}

transform_mean_auc_importance <- function(list_performance){
  mean_auc <- sapply(list_performance, function(x) median(x$auc))
  params <- tibble(feature= names(mean_auc), auc = mean_auc)
  return(params)
}



return_feature_importance <- function(path){
  
  dat_imp <- import_feature_importance(path)
  
  sorted_names <- mixedsort(names(dat_imp))
  dat_imp <- dat_imp[sorted_names]
  mean_auc <- transform_mean_auc_importance(dat_imp)
  return(mean_auc)
  
}


plot_miles_load <- function(dat, entity_ = "MEL", col_method, main = NULL){
  
  
  dat <- dat %>%
    mutate(mutation_type =ifelse(mutation_type == "combined", "all", mutation_type)) %>%
    mutate(mutation_type = factor(
      mutation_type,
      levels = c("all", "SNV", "INDEL", "Fusion gene")
    ))
  
  pvals <- dat %>%
    group_by(mutation_type, entity) %>%
    wilcox_test(AUROC ~ type) 
  
  pvals <- pvals %>% 
    mutate(y.position = c(1.03)) %>%
    mutate(p_ = signif(p, digits = 2))
  
  dat %>%
    filter(entity == entity_) %>%
    ggplot(aes(type, AUROC))+
    geom_boxplot(aes(fill = type))+
    scale_fill_manual(values = col_method)+
    facet_wrap(~ mutation_type, ncol = 4)+
    stat_pvalue_manual(pvals %>% filter(entity == entity_), label = "p_", size = 2) +
    xlab("")+
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45,  hjust=1),
          panel.spacing = unit(0, "lines"),
          strip.text.x = element_text(size = 6))+
    scale_x_discrete(labels=c("MILES" = "MILES", "Neoantigen candidate load" = "Neoantigen\ncandidate load"))+
    ggtitle(main)+
    ylim(c(0,1.05))
}


plot_importance <- function(importance, type) {
  
  col_mut <- c("SNV" = "#88CCEE","INDEL" = "#CC6677", "Fusion gene" = "#DDCC77","combined" = "#117733")
  
  col_mut <- col_mut[which(names(col_mut) %in% importance$class)]
  
  importance %>%
    mutate(class = factor(class, levels = c("SNV", "INDEL", "Fusion gene", "combined"))) %>%
    filter(entity == type) %>%
    group_by(feature) %>%
    filter(any(`delta AUC` > 0.05)) %>%
    ungroup() %>%
    ggplot(aes(`delta AUC`, reorder(feature, `delta AUC`), fill = class)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = col_mut)+
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(type) + 
    ylab("") +
    xlab("delta AUC")
}


# compare best and random apporach 
compare_random_best <- function(path_random, path_best, class = ""){
  
  dat_random <- read_delim(path_random)
  dat_best <- read_delim(path_best)
  
  dat_comp <- bind_rows(list(Best = dat_best, Random = dat_random), .id = "ds")
  
  dat_comp %>%
    ggplot(aes(ds, auc, fill = ds)) +
    geom_boxplot() +
    ylab("AUC")+ theme(legend.position="none")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    stat_compare_means(size = 2) +
    ggtitle(class) + 
    xlab("")+
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, NA))
}

