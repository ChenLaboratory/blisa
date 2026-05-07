# Suppress R CMD check NOTEs for non-standard evaluation column references.
# ggplot2 aes() and data.table := / .[] expressions reference column names as
# bare symbols, which the static checker cannot resolve.
utils::globalVariables(c(
  # ggplot2 aes() column references — plotLRIsum
  "sig_numbers", "LR_pair", "annotation",
  # ggplot2 aes() column references — plotLRI.sf
  "bin_status", "plot_val",
  # data.table column references — CCIspatial
  "hex_id", "receptor_expr", "ct", "ligand_expr",
  "hh_hex", "r_sum", "l_sum", "product",
  "cell_pair", "ct_l", "ct_r", ".I", "cell_pair_plot",
  # external dataset default argument — filterLRpairs / getLRpairs
  "CellChatDB.human"
))
