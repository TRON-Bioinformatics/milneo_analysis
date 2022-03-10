

# QUANTITATIVE OVERVIEW ---------------------------------------------------


# count candidates per patients 
get_numbers <- function(list, nam ="count"){
  lapply(list, function(x) {
    x%>% 
      count(patientIdentifier)%>%
      dplyr::rename("{nam}" := n)
  })
}

merge_counts <- function(number_snvs, number_indels, number_fusion_genes){
  # combined number candidates 
  number_summary <- mapply(function(snv, indel){
    if(snv == indel){
      number_snvs[[snv]] %>%
        full_join(number_indels[[indel]], by = "patientIdentifier")
    }
  }, snv = names(number_snvs), indel = names(number_indels), SIMPLIFY=F)
  number_summary <-lapply(names(number_summary), function(x){
    if(x %in% names(number_fusion_genes)){
      print(x)
      ind <- which(x == names(number_fusion_genes))
      print(ind)
      dat <- number_fusion_genes[[ind]]
      number_summary[[x]] %>%
        full_join(dat, by = "patientIdentifier")
    }else{
      number_summary[[x]] %>%
        mutate("fusion gene" = NA)
    }
  })
  number_summary <- number_summary %>%
    map(., function(x) x%>%mutate(all=rowSums(select(.,snv, INDEL,`fusion gene`), na.rm = T)))
  names(number_summary) <- names(number_snvs)
  
  return(number_summary)
}

# get numbers of expressed neoantigen candidates
filter_expressed <- function(x, nam = "snv/INDEL"){
  if(nam == "fusion gene"){
    x1 <- x%>%
      filter(sum_junc_span > 0)
  }else{
    x1 <- x%>%
      filter(Expression_mutated_transcript > 0)
  }
}

# get numbers of  neoantigen candidates with at least one predicted binder 
filter_mhc_binding <- function(x, thresh = 500){
  x1 <- x%>%
    filter(Best_affinity_MHCI_score < thresh | Best_affinity_MHCII_score < thresh)
}

# transform count list into tall skinny tibble 
transfrom_count2tall <- function(number_list){
  nice_names <- c("Van Allen","Hugo", "Riaz", "Miao", "McDermott")
  names(nice_names) <- c( "vanallen" ,"hugo"  ,"riaz_pre" ,"miao" ,"mcdermott_Atezo")
  technical_problems <- c("hugo" = "Pt8", "miao" = "PtRCC_106","mcdermott_Atezo" =  "Pt164370145747331", "riaz_pre" = "no", "vanallen" = "no")
  l <- lapply(number_list, function(x){
    gather(x, class, count, snv:all, factor_key=TRUE)
  })
  res <- bind_rows(l, .id = "cohort")
  res %>%
    mutate(class = as.character(class)) %>%
    mutate(class = ifelse(class == "snv", "SNV", ifelse(class == "fusion gene", "Fusion gene", class))) %>%
    mutate(class =factor(class, levels = c("SNV", "INDEL", "Fusion gene", "all"))) %>%
    mutate(count = ifelse(is.na(count) & patientIdentifier != technical_problems[cohort], 0, count )) %>%
    mutate(cohort = nice_names[cohort]) %>%
    
    mutate(cohort = factor(cohort, levels = c(
      "Van Allen",
      "Hugo", 
      "Riaz", 
      "Miao",
      "McDermott"
    )))
      
}


# PLOTTING ----------------------------------------------------------------

## proportion barplot
make_proportion_barplot <- function(x, main = ""){
  x1 <- x %>%
    arrange(percentage_snv, percentage_fusion, percentage_indel)
  x1df <- t(as.matrix(x1[,6:8]))
  barplot(x1df, names.arg =x1$id, las = 2, col= col_mutation, main = main, cex.names = 0.5)
  legend(x = "bottomright", legend = c("SNV", "INDEL", "Fusion gene"), bg = "white", pch  =15, col = col_mutation)
}

# compare numbers of all predicted and selected neoantigen candidates
boxplot_all_selected <- function(x, main = "") {
  x %>%
    ggplot(aes(y = count + 1, x = class, fill = category)) +
    geom_boxplot()  +
    geom_boxplot(outlier.fill = NULL, outlier.shape = 21) +
    scale_y_continuous(trans =                         'log10') +
    scale_fill_manual(values = col_expressed) +
    scale_color_manual(values =                         col_expressed) +
    theme(panel.grid.major = element_blank()) +
    stat_compare_means(aes(group = category), label = "p.format") +
    ggtitle(main)
  
}

boxplot_all_selected_facet <- function(x, main = "") {
  x %>%
    ggplot(aes(y = count + 1, x = class, fill = category)) +
    geom_boxplot()  +
    geom_boxplot(outlier.fill = NULL, outlier.shape = 21) +
    scale_y_continuous(trans = 'log10') +
    scale_fill_manual(values = col_expressed) +
    scale_color_manual(values = col_expressed) +
    theme(panel.grid.major = element_blank()) +
    stat_compare_means(aes(group = category), label = "p.format") +
    facet_grid(cohort ~ . ) +
    ggtitle(main)
  
}

# boxplot response analysis 
make_boxplot_response <- function(x, main = ""){
  col_response <- c(no=rgb(12,44,132, 255, maxColorValue = 255),
                    yes = rgb(178,34,34, 250 , maxColorValue = 255))
  # add wilcoxon p-value + number of patients per group 
  n_responder = paste0("responder\n(n=", length(which(x$binary_response =="yes")), ")")
  n_non_responder = paste0("non-responder\n(n=", length(which(x$binary_response =="no")), ")")
  # ylab = paste0("# neontigen candidates + 1")
  # 
  x %>%
    filter(!is.na(binary_response)) %>%
    ggplot(aes(binary_response, count + 1, fill=binary_response))+
    geom_boxplot() +
    scale_y_log10()+
    scale_fill_manual(values = col_response) +
    stat_compare_means(aes(group = binary_response), label = "p.format") +
    #ylab(ylab)+
    scale_x_discrete(labels=c(n_non_responder, n_responder))+
    theme(plot.title = element_text(size = 8))+ 
    ggtitle(main)
}

make_boxplot_response_facet <- function(x, class="snv", main = ""){
  col_response <- c(no=rgb(12,44,132, 255, maxColorValue = 255),
                    yes = rgb(178,34,34, 250 , maxColorValue = 255))
  # add wilcoxon p-value + number of patients per group 
  x1 <- x%>% distinct(cohort_pat, .keep_all = TRUE)
  n_responder = paste0("responder\n(n=", length(which(x1$binary_response =="yes")), ")")
  n_non_responder = paste0("non-responder\n(n=", length(which(x1$binary_response =="no")), ")")
  ylab = paste0("# neontigen candidates + 1")
  
  x %>%
    filter(!is.na(binary_response)) %>%
    ggplot(aes(binary_response, count + 1, fill=binary_response))+
    geom_boxplot() +
    scale_y_log10()+
    scale_fill_manual(values = col_response) +
    stat_compare_means(aes(group = binary_response), label = "p.format") +
    ylab(ylab)+
    facet_wrap(.~ class )+
    #facet_wrap(.~ class + cohort.x)
    scale_x_discrete(labels=c(n_non_responder, n_responder))
}



# response prediction  -------------------------------------------------------------

# count only candidates that meet defined thresholds 
count_neoantigens_quality2 <- function(df, feature1, feature2, threshold, all_patients){
  
  technical_problems <- c("hugo" = "Pt8", "miao" = "no","mcdermott_Atezo" =  "Pt164370145747331", "riaz_pre" = "no", "vanallen" = "no")
  cohort_patients <- unique(paste0(df$cohort, "_", df$patientIdentifier))
  
  df1 <- df %>%
    filter(!!sym(feature1) < threshold | !!sym(feature2) < threshold) %>%
    count(cohort, patientIdentifier, binary_response) %>%
    dplyr::rename(count=n) %>%
    mutate(cohort_patient = paste0(cohort, "_", patientIdentifier))
  
  df_no_counts <- all_patients%>%
    filter(!cohort_patient %in% df1$cohort_patient) %>%
    mutate(count = NA )%>%
    mutate(count = ifelse(is.na(count) & patientIdentifier != technical_problems[cohort], 0, count )) %>%
    select(-cohort_patient)
  df1 <- df1 %>%
    select(-cohort_patient) %>%
    mutate(binary_response = as.character(binary_response))
  df_res <- bind_rows(df1, df_no_counts)
  
  return(df_res)
}

neoag_qual_class_count <- function(snv, indel, fg, feature1, feature2, thres){
  nice_names <- c("Van Allen","Hugo", "Riaz", "Miao", "McDermott", "all")
  names(nice_names) <- c( "vanallen" ,"hugo"  ,"riaz_pre" ,"miao" ,"mcdermott_Atezo", "all")
  dat_transformed_df <- bind_rows(snv, .id = "cohort")
  dat_indel_transformed_df <- bind_rows(indel, .id = "cohort")
  dat_fusions_transformed_df <- bind_rows(fg, .id = "cohort")
  
  all_patients <- bind_rows(dat_transformed_df[,c("cohort", "patientIdentifier", "binary_response")], 
                            dat_indel_transformed_df[,c("cohort", "patientIdentifier", "binary_response")], 
                            dat_fusions_transformed_df[,c("cohort", "patientIdentifier", "binary_response")]) %>%
    select(cohort, patientIdentifier, binary_response) %>%
    distinct()%>%
    mutate(cohort_patient = paste0(cohort, "_", patientIdentifier))
  
  df_snv <- count_neoantigens_quality2(dat_transformed_df, feature1 = feature1, feature2 = feature2, threshold = thres, all_patients)
  df_indel <- count_neoantigens_quality2(dat_indel_transformed_df, feature1 = feature1, feature2 = feature2, threshold = thres, all_patients)
  df_fusion <- count_neoantigens_quality2(dat_fusions_transformed_df, feature1 = feature1, feature2 = feature2, threshold = thres, all_patients)

  
  df_all <- bind_rows(list(SNV = df_snv, INDEL = df_indel, "Fusion gene" = df_fusion), .id = "class")
  df_sum <- df_all %>%
    mutate(cohort_pat = paste0(cohort, "_", patientIdentifier))%>%
    group_by(cohort_pat) %>%
    summarise(cohort_pat, cohort,patientIdentifier, binary_response, count = sum(count, na.rm = TRUE)) %>%
    ungroup()%>%
    mutate(class = "all") %>% 
    distinct()%>%select(-cohort_pat)

  df_all <- bind_rows(df_all, df_sum)
  df_all$class <- factor(df_all$class, levels= c("SNV", "INDEL", "Fusion gene", "all"))

  return(df_all)
  
}

neoag_qual_class_pval <- function(df_all){
  nice_names <- c("Van Allen","Hugo", "Riaz", "Miao", "McDermott", "all")
  names(nice_names) <- c( "vanallen" ,"hugo"  ,"riaz_pre" ,"miao" ,"mcdermott_Atezo", "all")
  pval_mhc <- compare_means(count ~ binary_response, data = df_all, group.by = c("class", "cohort"))
  
  pval_mhc_all_cohorts <- compare_means(count ~ binary_response, data = df_all, group.by = c("class"))
  pval_mhc_all_cohorts$cohort <- "all"
  
  pval_mhc <- bind_rows(list(pval_mhc, pval_mhc_all_cohorts))
  pval_mhc <- pval_mhc %>%
    mutate(cohort = nice_names[cohort])
  return(pval_mhc)
}



# ROC-Analysis ------------------------------------------------------------

split_cv <- function(seed, df){
  set.seed(seed)
  indx <- caret::createFolds(df$binary_response , k = 10, list = TRUE, returnTrain = FALSE)
  indx_data <- lapply(indx, function(x){
    df[x,]
  })
  return(indx_data)
}

# call this to create folds for 10x10 nested cross validation 
nested_CV <- function(df){
  sd <- c(42, 23, 1991, 3, 5, 21, 7, 2022, 31, 8 )
  list_cv_splits <- lapply(sd, split_cv, df = df)
  list_cv_splits <- c(list_cv_splits[[1]], list_cv_splits[[2]], list_cv_splits[[3]], list_cv_splits[[4]],
                      list_cv_splits[[5]], list_cv_splits[[6]], list_cv_splits[[7]], list_cv_splits[[8]],
                      list_cv_splits[[9]], list_cv_splits[[10]])
  return(list_cv_splits)
}

roc_with_cv <- function(predictions_full, labels, col.background = "grey", col.average = "black", main = "") {
  
  pred <- prediction(predictions_full, labels)
  perf <- performance(pred,'tpr','fpr')
  auc <- unlist(performance(pred,'auc')@y.values)
  auc.mean <- signif(median(auc, na.rm = T), digits = 2)
  auc.sd <- signif(IQR(auc, na.rm = T), digits = 2)
  
  print(perf)
  
  plot(perf,
       col = col.background,
       las = 2,
       lwd=1, main = main)
  plot(perf,
       avg='vertical',
       lwd=3, add = T,
       col = col.average,
       las = 1)
  legend(x = "bottomright", legend = paste0("AUC = ", auc.mean, " +/- ", auc.sd))
  
  tib_auc = tibble(mean_auc = auc.mean, sd_auc = auc.sd, num = "1")
  return(tib_auc)
}
