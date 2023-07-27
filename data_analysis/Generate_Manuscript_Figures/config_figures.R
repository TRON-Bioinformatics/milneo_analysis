

# THEME -------------------------------------------------------------------

theme_adju <- function(){
  theme_bw()+ 
    theme(axis.text=element_text(colour="black"))+ 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "top",
          legend.justification="right",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-10),
          legend.title = element_blank(),
          legend.text = element_text(size = 7),
          title = element_text(size = 8),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 8),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
}

theme_set(theme_adju())



# COLOR DEFINITIONS -------------------------------------------------------

col_mutation <-
  c(SNV = "#6477B9",
    "INDEL" = "#D969A5",
    "Fusion gene" = "#D9AF6B",
    "combined" = "#85BF63")

col_response <- c(
  no = rgb(12, 44, 132, 255, maxColorValue = 255),
  yes = rgb(178, 34, 34, 250 , maxColorValue = 255)
)

col_method = c("#A1D7DF","#295960")


col_entity <- c("MEL+RCC" = "#117733", "MEL" = "#999933","RCC" =  "#6699CC")


# NEOANTIGEN FEATURES -----------------------------------------------------

features_full <-
  c(
    "Best_rank_MHCI_score",
    "Best_rank_MHCII_score",
    "rnaExpression",
    "Pathogensimiliarity_MHCI_9mer",
    "Pathogensimiliarity_MHCII",
    "Selfsimilarity_MHCI",
    "Selfsimilarity_MHCII",
    "rnaVariantAlleleFrequency",
    "Amplitude_MHCII_rank",
    "Amplitude_MHCI_affinity",
    "DAI_MHCI_affinity",
    "Dissimilarity_MHCI",
    "Generator_rate_MHCI",
    "Generator_rate_MHCII",
    "IEDB_Immunogenicity_MHCI",
    "IEDB_Immunogenicity_MHCII",
    "MixMHC2pred_best_rank",
    "MixMHCpred_best_rank",
    "PHBR_I",
    "PHBR_II",
    "vaxrank_binding_score",
    "Hex_alignment_score_MHCI",
    "Hex_alignment_score_MHCII",
    "Neoag_immunogenicity" ,
    "Priority_score",
    "Recognition_Potential_MHCI_9mer",
    "Tcell_predictor_score",
    "vaxrank_total_score",
    "PRIME_best_rank"
  )
