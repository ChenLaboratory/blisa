plotLRIsum <- function(LR_out, top = NULL, pt.size = 4) {

  if(!is.null(top)){
    LR_out <- LR_out[1:top,]
  }

  LR_out$LR_pair <- rownames(LR_out)

  LR_out <- LR_out %>%
    arrange(desc(`sig_numbers`)) %>%
    mutate(LR_pair = factor(LR_pair, levels = rev(LR_pair))) # reverse so high at top

  p <- ggplot(LR_out, aes(x = `sig_numbers`, y = LR_pair, color = annotation)) +
    geom_point(size = pt.size)+
    scale_color_manual(
      values = c(
        "Secreted Signaling" = "#90a955",   # blue
        "ECM-Receptor"       = "#219ebc",   # green
        "Cell-Cell Contact"  = "#f7b801",   # red
        "Non-protein Signaling" = "#9f86c0" # purple
      )
    ) +
    scale_x_continuous(expand = expansion(add = 100)) +
    labs(x = "Sig Spot Numbers", y = "Ligand–Receptor Pair", color = "Annotation") +
    theme_minimal()+
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.text.y = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )+
    coord_cartesian(clip = "off")

  return(p)

}
