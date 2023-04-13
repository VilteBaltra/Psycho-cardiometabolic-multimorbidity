## tidies up MR results and adds intercept and heterogeneity output ##

my_MR_tests <- function(res, dat){
  # obtain MR Egger intercept 
  intercept = mr_pleiotropy_test(dat)
  intercept$method <- "MR Egger"
  names(intercept) <- c('id.exposure', 'id.outcome', 'outcome', 'exposure', 'egger_intercept', 'intercept_se', 'intercept_pval', 'method')
  order <- c('Inverse variance weighted', 'MR Egger', 'Simple mode', 'Weighted median', 'Weighted mode')
  res_ordered <- res %>% slice(match(order, method))
  res_ordered <- merge(res_ordered,intercept[, c('method', 'egger_intercept', 'intercept_se', 'intercept_pval')], by = "method", all = T)
  # obtain heterogeneity values 
  het = mr_heterogeneity(dat)
  merged <- merge(res_ordered,het[, c('method', 'Q', 'Q_df', 'Q_pval')], by = "method", all = T)
  
  # generate 95% OR
  generate_odds_ratios <- function(res) 
  {
    lo_ci <- res$b - 1.96 * res$se
    up_ci <- res$b + 1.96 * res$se
    or <- exp(res$b)
    or_lci95 <- exp(lo_ci)
    or_uci95 <- exp(up_ci)
    res$`OR (95% CI)` <- paste0(round(or, 2), " (", round(or_lci95, 2), "-", round(or_uci95, 2), ")")
    return(res)
  }
  
  merged <- generate_odds_ratios(merged)
  # return combined data frame
  return(merged)
}
