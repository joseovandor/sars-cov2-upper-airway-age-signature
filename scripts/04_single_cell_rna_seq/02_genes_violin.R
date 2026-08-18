# ============================================================
# PART 2
# VIOLIN PLOTS OF SIGNIFICANT RESULTS
#
# FDR < 0.05
#
# Pediatric = yellow
# Adult     = blue
#
# Each point = DONOR
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

library(dplyr)
library(ggplot2)
library(tibble)


# ============================================================
# 1. COLORS
# ============================================================

group_colors <- c(
  
  "Pediatric" =
    "#E9C46A",
  
  "Adult" =
    "#457B9D"
)


# ============================================================
# 2. SELECT SIGNIFICANT ONLY
# ============================================================

significant_results <- test_results %>%
  
  dplyr::filter(
    
    !is.na(FDR),
    
    FDR < 0.05
  ) %>%
  
  dplyr::arrange(
    FDR
  )


cat(
  "\nSignificant results FDR < 0.05:\n"
)


print(
  tibble::as_tibble(
    significant_results
  ),
  n = Inf
)


# ============================================================
# 3. STOP IF THERE ARE NO SIGNIFICANT RESULTS
# ============================================================

if (nrow(significant_results) == 0) {
  
  stop(
    "\nNo associations with FDR < 0.05 exist."
  )
}


# ============================================================
# 4. PREPARE DATA
# ============================================================

violin_data <- expr_primary %>%
  
  dplyr::mutate(
    
    Group_plot =
      dplyr::case_when(
        
        Group == "Ped" ~
          "Pediatric",
        
        Group == "Adult" ~
          "Adult",
        
        TRUE ~
          as.character(Group)
      )
  )


violin_data$Group_plot <- factor(
  
  violin_data$Group_plot,
  
  levels = c(
    "Pediatric",
    "Adult"
  )
)


# ============================================================
# 5. KEEP SIGNIFICANT GENE x CELL TYPE PAIRS
# ============================================================

violin_data <- violin_data %>%
  
  dplyr::inner_join(
    
    significant_results %>%
      
      dplyr::select(
        SYMBOL,
        cell_type,
        FDR,
        p_value
      ),
    
    by = c(
      "SYMBOL",
      "cell_type"
    )
  )


# ============================================================
# 6. PANEL
# ============================================================

violin_data <- violin_data %>%
  
  dplyr::mutate(
    
    panel =
      paste0(
        SYMBOL,
        "\n",
        cell_type
      )
  )


# ============================================================
# 7. SIGNIFICANCE LABEL
# ============================================================

violin_data <- violin_data %>%
  
  dplyr::mutate(
    
    significance =
      dplyr::case_when(
        
        FDR < 0.001 ~
          "***",
        
        FDR < 0.01 ~
          "**",
        
        FDR < 0.05 ~
          "*",
        
        TRUE ~
          ""
      ),
    
    
    annotation =
      paste0(
        
        significance,
        
        "   FDR = ",
        
        formatC(
          FDR,
          format = "f",
          digits = 3
        )
      )
  )


# ============================================================
# 8. ANNOTATION POSITION
# ============================================================

violin_annotation <- violin_data %>%
  
  dplyr::group_by(
    panel
  ) %>%
  
  dplyr::summarise(
    
    y_max =
      max(
        mean_expression,
        na.rm = TRUE
      ),
    
    
    y_min =
      min(
        mean_expression,
        na.rm = TRUE
      ),
    
    
    annotation =
      dplyr::first(
        annotation
      ),
    
    
    .groups =
      "drop"
  ) %>%
  
  dplyr::mutate(
    
    y_range =
      y_max - y_min,
    
    
    y_range =
      ifelse(
        
        y_range <= 0,
        
        pmax(
          abs(y_max),
          0.1
        ),
        
        y_range
      ),
    
    
    bracket_y =
      y_max +
      0.20 * y_range,
    
    
    annotation_y =
      y_max +
      0.37 * y_range
  )


# ============================================================
# 9. COMBINED VIOLIN PLOTS
# ============================================================

p_violin <- ggplot(
  
  violin_data,
  
  aes(
    x = Group_plot,
    y = mean_expression,
    fill = Group_plot
  )
  
) +
  
  # ----------------------------------------------------------
# VIOLIN
# ----------------------------------------------------------

geom_violin(
  trim = FALSE,
  width = 0.78,
  alpha = 0.72,
  linewidth = 0.5
) +
  
  
  # ----------------------------------------------------------
# BOXPLOT
# ----------------------------------------------------------

geom_boxplot(
  width = 0.13,
  outlier.shape = NA,
  fill = "white",
  alpha = 0.88,
  linewidth = 0.5
) +
  
  
  # ----------------------------------------------------------
# DONORS
# ----------------------------------------------------------

geom_jitter(
  
  aes(
    color = Group_plot
  ),
  
  width = 0.065,
  height = 0,
  
  size = 3,
  
  alpha = 0.88
) +
  
  
  # ----------------------------------------------------------
# BRACKET
# ----------------------------------------------------------

geom_segment(
  
  data = violin_annotation,
  
  aes(
    x = 1,
    xend = 2,
    y = bracket_y,
    yend = bracket_y
  ),
  
  inherit.aes = FALSE,
  
  linewidth = 0.55
) +
  
  
  geom_segment(
    
    data = violin_annotation,
    
    aes(
      x = 1,
      xend = 1,
      
      y = bracket_y,
      
      yend =
        bracket_y -
        0.04 * y_range
    ),
    
    inherit.aes = FALSE,
    
    linewidth = 0.55
  ) +
  
  
  geom_segment(
    
    data = violin_annotation,
    
    aes(
      x = 2,
      xend = 2,
      
      y = bracket_y,
      
      yend =
        bracket_y -
        0.04 * y_range
    ),
    
    inherit.aes = FALSE,
    
    linewidth = 0.55
  ) +
  
  
  # ----------------------------------------------------------
# SIGNIFICANCE + FDR
# ----------------------------------------------------------

geom_text(
  
  data = violin_annotation,
  
  aes(
    x = 1.5,
    y = annotation_y,
    label = annotation
  ),
  
  inherit.aes = FALSE,
  
  size = 5
) +
  
  
  # ----------------------------------------------------------
# FACETS
# ----------------------------------------------------------

facet_wrap(
  ~ panel,
  scales = "free_y",
  nrow = 1
) +
  
  
  # ----------------------------------------------------------
# COLORS
# ----------------------------------------------------------

scale_fill_manual(
  values = group_colors
) +
  
  
  scale_color_manual(
    values = group_colors
  ) +
  
  
  # ----------------------------------------------------------
# SPACE FOR FDR
# ----------------------------------------------------------

scale_y_continuous(
  
  expand = expansion(
    mult = c(
      0.05,
      0.43
    )
  )
) +
  
  
  labs(
    x = NULL,
    y = "Mean normalized expression"
  ) +
  
  
  theme_classic(
    base_size = 15
  ) +
  
  
  theme(
    
    legend.position =
      "none",
    
    
    axis.text.x =
      element_text(
        size = 14,
        face = "plain",
        color = "black"
      ),
    
    
    axis.text.y =
      element_text(
        size = 12,
        face = "plain",
        color = "black"
      ),
    
    
    axis.title.y =
      element_text(
        size = 14,
        face = "plain",
        color = "black"
      ),
    
    
    strip.background =
      element_blank(),
    
    
    strip.text =
      element_text(
        size = 15,
        face = "plain",
        color = "black"
      ),
    
    
    panel.spacing =
      unit(
        1.5,
        "lines"
      )
  )


print(
  p_violin
)


# ============================================================
# 10. DIRECTORY
# ============================================================

violin_dir <- file.path(
  figure_dir,
  "significant_violin_plots"
)


dir.create(
  violin_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 11. SAVE COMBINED TIFF FIGURE
# ============================================================

ggsave(
  
  filename =
    file.path(
      violin_dir,
      "Significant_genes_Pediatric_vs_Adult_violin.tiff"
    ),
  
  plot =
    p_violin,
  
  width =
    10,
  
  height =
    6,
  
  dpi =
    600,
  
  compression =
    "lzw"
)