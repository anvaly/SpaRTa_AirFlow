a = "XXXX"
b = "XXXX"
c = "XXXX"


a_df <- readr::read_delim(a, delim="\t")
b_df <- readr::read_delim(b, delim="\t")
c_df <- readr::read_delim(c, delim="\t")

data_df <- a_df %>%
  dplyr::full_join(b_df, by="Symbol") %>%
  dplyr::full_join(c_df, by="Symbol_mouse") %>%
  dplyr::select(Symbol, Symbol_mouse, ClLin_DepMap_ExpressionMean, KP2SyngBulk_Mean, KP2SyngBulk_Ductal_2, KP2CLparentalExprNorm) %>%
  dplyr::mutate(HumanMouseDiff = ClLin_DepMap_ExpressionMean - KP2CLparentalExprNorm)

readr::write_delim(data_df, "XXXX")
