import_cv <- function(path_to_files) {
  files <- list.files(path_to_files, full.names = T, pattern = ".tsv")
  nams <- gsub(path_to_files, "", files)
  nams <- gsub("/miles_", "", gsub(".tsv", "", nams))
  dat <- files %>%
    map(., read_delim, delim = "\t")
  names(dat) <- nams
  return(dat)
}

import_feature_importance <- function(path_to_files){
  files <- list.files(path_to_files, pattern = ".csv", full.names = T)
  nams <- gsub(".csv", "", gsub( "/", "" , gsub(path_to_files, "", files)))
  dat <- files %>%
    map(read_delim) 
  names(dat) <- nams
  return(dat)
}


median_performance <- function(df){
  apply(df[,2:5], 2, median, na.rm = T)
}


mean_performance <- function(df){
  apply(df[,2:5], 2, mean, na.rm = T)
}

sd_performance <- function(df){
  apply(df[,2:5], 2, IQR, na.rm = T)
}

transform_performance <- function(list_df){
  list_mean <- list_df %>%
    map(., median_performance)
  list_sd <- list_df %>%
    map(., sd_performance)
  performance <- list_mean %>%
    bind_rows(.id = "parameter") %>%
    separate(parameter, into = c("sigma2", "c"), sep = "_")
  df_sd <- list_sd %>%
    bind_rows()
  colnames(df_sd) <- c("auc_iqr", "binaryaccuracy_iqr", "precision_iqr", "recall_iqr")
  performance <- bind_cols(performance, df_sd)
  return(performance)
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


get_best_model <- function(df){
  df1 <- df %>% 
    filter(auc == max(auc))
  if(nrow(df1) > 1){
    df1 <- df1 %>% 
      filter(binaryaccuracy == max(binaryaccuracy))
    if(nrow(df1) > 1){
      df1 <- df1 %>% 
        filter(precision  == max(precision ))
      if(nrow(df1) > 1){
        df1 <- df1 %>% 
          mutate(sigma2 = as.numeric(sigma2))%>%
          mutate(c = as.numeric(c))%>%
          filter(sigma2  == min(sigma2)) %>%
          filter(c == min(c))
      }
    }
  }
  return(df1)
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


compare_random_best <- function(path_random, path_best, class = ""){
  
  dat_random <- read_delim(path_random)
  dat_best <- read_delim(path_best)
  
  dat_comp <- bind_rows(list(Best = dat_best, Random = dat_random), .id = "ds")
  
  dat_comp %>%
    ggplot(aes(ds, auc, fill = ds)) +
    geom_boxplot() +
    ylab("AUC")+ theme(legend.position="none")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    stat_compare_means() +
    ggtitle(class) + 
    xlab("")+
    ylim(c(0,1))
}


return_feature_importance <- function(path){
  
  dat_imp <- import_feature_importance(path)
  
  sorted_names <- mixedsort(names(dat_imp))
  dat_imp <- dat_imp[sorted_names]
  mean_auc <- transform_mean_auc_importance(dat_imp)
  return(mean_auc)
  
}