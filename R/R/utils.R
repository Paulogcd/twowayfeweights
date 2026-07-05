#' @importFrom rlang sym

printf <- function(...) print(sprintf(...))

fn_ctrl_rename <- function(x) paste("ctrl", x, sep="_")
get_controls_rename <- function(controls) unlist(lapply(controls, fn_ctrl_rename))
fn_treatment_rename <- function(x) paste("OT", x, sep="_")

get_treatments_rename <- function(treatments) {
  unlist(lapply(treatments, fn_treatment_rename))
}

fn_treatment_weight_rename <- function(x) paste("weight_", x, sep = "")
fn_random_weight_rename <- function(x) paste("RW", x, sep="_")
get_random_weight_rename <- function(ws) unlist(lapply(ws, fn_random_weight_rename))


##
# twowayfeweights_rename_var
twowayfeweights_rename_var <- function(df, Y, G, T, D, D0, controls, treatments, random_weights) {

  random_weight_df <- NULL
  D0_df <- NULL
  controls_df <- NULL
  treatments_df <- NULL

  # Generate a random character value that is not contained in the column names
  # of the dataframe to act as the identification column of the dataset:
  id.col = paste(sample(LETTERS, 10, TRUE), collapse = '')
  while(id.col %in% colnames(df)){
    id.col = paste(sample(LETTERS, 10, TRUE), collapse = '')
  }
  df <- df %>%
    dplyr::mutate(!!id.col := 1:nrow(df))

  # Y, G, T, D
  original_names = c(Y, G, T, D)
  df.base <- df %>%
    dplyr::select(dplyr::all_of(c(!!id.col, original_names)))
  
  names(original_names) = c("Y", "G", "T", "D")
  
  df.base <- df.base %>%
    dplyr::rename(original_names)
  
  controls_rename   <- get_controls_rename(controls)
  treatments_rename <- get_treatments_rename(treatments)

  if(!is.null(controls)){
    controls_df = df %>% dplyr::select(dplyr::all_of(c(controls, !!id.col)))
    
    old.names.controls         = c(controls, id.col)
    names(old.names.controls)  = c(controls_rename, id.col)
    controls_df <- controls_df %>% dplyr::rename(old.names.controls)

  }
  
  if(!is.null(treatments)){
    treatments_df = df %>% dplyr::select(dplyr::all_of(c(treatments, !!id.col)))
    
    old.names.treatments         = c(treatments, id.col)
    names(old.names.treatments)  = c(treatments_rename, id.col)
    treatments_df <- treatments_df %>% dplyr::rename(old.names.controls)

  }

  if (length(random_weights) > 0) {
    random_weight_rename <- get_random_weight_rename(random_weights)
    # random_weight_df <- df[, random_weights, drop = FALSE]
    random_weight_df <- df %>% dplyr::select(dplyr::all_of(c(random_weights, !!id.col)))
    
    old.names.weights         = c(random_weights, id.col)
    names(old.names.weights)  = c(random_weight_rename, id.col)
    random_weight_df <- random_weight_df %>% dplyr::rename(old.names.weights)
  }
  
  if (!is.null(D0)) {
    D0_df = df %>% dplyr::select(dplyr::all_of(c(D0, !!id.col)))
    
    old.names         = c(D0, id.col)
    names(old.names)  = c("D0", id.col)
    D0_df <- D0_df %>% dplyr::rename(old.names)
  }
  
  # This formulation does not select a column twice if there is a duplicate in the original_names vector.
  # df <- data.frame(df) %>% dplyr::select_at(dplyr::vars(original_names))
  
  df <- df.base

  # For each non null variable, we add the corresponding column:
  if(!is.null(D0_df)){
    df <- df %>%
      dplyr::left_join(D0_df, by = dplyr::join_by(!!id.col))
  }
  if(!is.null(random_weight_df)){
    df <- df %>%
      dplyr::left_join(random_weight_df, by = dplyr::join_by(!!id.col))
  }
  if(!is.null(treatments_df)){
    df <- df %>%
      dplyr::left_join(treatments_df, by = dplyr::join_by(!!id.col))
  }
  if(!is.null(controls_df)){
    df <- df %>%
      dplyr::left_join(controls_df, by = dplyr::join_by(!!id.col))
  }

  df <- df %>% dplyr::select(-c(!!id.col))
  
  return(df)
}


##
# twowayfeweights_transform
twowayfeweights_transform <- function(df, controls, weights, treatments) {
  
  .data = NULL
  
  ret = twowayfeweights_normalize_var(df, "D")
  if (ret$retcode) {
    df <- ret$df
    printf("The treatment variable in the regression varies within some group * period cells.")
    printf("The results in de Chaisemartin, C. and D'Haultfoeuille, X. (2020) apply to two-way fixed effects regressions")
    printf("with a group * period level treatment.")
    printf("The command will replace the treatment by its average value in each group * period.")
    printf("The results below apply to the two-way fixed effects regression with that treatment variable.")
  }
  
  for (control in controls) {
    ret = twowayfeweights_normalize_var(df, control)
    if (ret$retcode) {
      df <- ret$df
      printf("The control variable %s in the regression varies within some group * period cells.", control)
      printf("The results in de Chaisemartin, C. and D'Haultfoeuille, X. (2020) apply to two-way fixed effects regressions")
      printf("with controls apply to group * period level controls.")
      printf("The command will replace replace control variable %s by its average value in each group * period.", control)
      printf("The results below apply to the regression with control variable %s averaged at the group * period level.", control)
    }
  }
  
  for (treatment in treatments) {
    ret = twowayfeweights_normalize_var(df, treatment)
    if (ret$retcode) {
      df <- ret$df
      printf("The other treatment variable %s in the regression varies within some group * period cells.", treatment)
      printf("The results in de Chaisemartin, C. and D'Haultfoeuille, X. (2020) apply to two-way fixed effects regressions")
      printf("with several treatments apply to group * period level controls.")
      printf("The command will replace replace other treatment variable %s by its average value in each group * period.", treatment)
      printf("The results below apply to the regression with other treatment variable %s averaged at the group * period level.", treatment)
    }
  }
  
  if (is.null(weights)) {
    df$weights <- 1
  } else {
    df$weights <- weights
  }
  
  df$Tfactor <- factor(df$T)
  TfactorLevels <- length(levels(df$Tfactor))
  df <- df %>% dplyr::mutate(TFactorNum = as.numeric(factor(.data$Tfactor, labels = seq(1:TfactorLevels))))
  
  return(df)
}


##
# twowayfeweights_filter
twowayfeweights_filter <- function(df, Y, G, T, D, D0, cmd_type, controls, treatments) {
  .data = NULL
  # Remove rows with NA values
  if (cmd_type != "fdTR") {
    df <- df %>%
      dplyr::mutate(tag = rowSums(dplyr::across(.cols = c(Y, G, T, D, controls, treatments), .fns = is.na))) %>%
      dplyr::filter(.data$tag == 0) %>%
      dplyr::select(-.data$tag)
  } else {
    df <- df %>%
      dplyr::mutate(tag1 = rowSums(dplyr::across(.cols = c(D, T, Y), .fns = is.na))) %>%
      dplyr::mutate(tag2 = rowSums(dplyr::across(.cols = c(D0), .fns = is.na))) %>%
      dplyr::filter(.data$tag1 == 0 | .data$tag2 == 0)
    
    if (length(controls) > 0) {
      df <- df %>%
        dplyr::mutate(tag3 = rowSums(dplyr::across(.cols = controls, .fns = is.na))) %>%
        dplyr::filter(.data$tag1 == 1 | .data$tag3 == 0) %>%
        dplyr::select(-.data$tag3)
    }
    df <- df %>% dplyr::select(-.data$tag1, -.data$tag2)
  }
  return(df)
}



##
# twowayfeweights_summarize_weights
twowayfeweights_summarize_weights <- function(df, var_weight) {

  weight_plus <- df[[var_weight]][df[[var_weight]] > 0 & !is.na(df[[var_weight]])]
  nr_plus <- length(weight_plus)
  sum_plus <- sum(weight_plus, na.rm = TRUE)
  
  weight_minus <- df[[var_weight]][df[[var_weight]] < 0 & !is.na(df[[var_weight]])]
  nr_minus <- length(weight_minus)
  sum_minus <- sum(weight_minus, na.rm = TRUE)
  
  nr_weights <- nr_plus + nr_minus
  
  return(
    list(
      nr_plus    = nr_plus,
      nr_minus   = nr_minus,
      nr_weights = nr_weights,
      sum_plus   = sum_plus,
      sum_minus  = sum_minus
    )
  )

}

##
# twowayfeweights_test_random_weights
twowayfeweights_test_random_weights <- function(df, random_weights) {
  
  .data = NULL
  
  mat <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(mat) <- c("Coef", "SE", "t-stat", "Correlation")
  df_filtered <- df %>% dplyr::filter(is.finite(.data$W))
  df_filtered_sub <- subset(df_filtered, df_filtered$nat_weight != 0) #Modif. Diego: added extra line to solve note in R CMD Check
 
  for (v in random_weights) {
    formula <- sprintf("%s ~ W", v)
    # rw.lm = estimatr::lm_robust(formula = as.formula(formula), data = df_filtered_sub, weights = df_filtered_sub$nat_weight, clusters = df_filtered_sub$G, se_type = "stata")
    # beta <- rw.lm$coefficients[["W"]]
    # se <- rw.lm$std.error[["W"]]
    # r2 <- rw.lm$r.squared
    rw_lm = fixest::feols(fml = as.formula(formula), data = df_filtered_sub, weights = ~nat_weight, vcov = ~G)
    beta = stats::coef(rw_lm)[["W"]]
    se = sqrt(diag(stats::vcov(rw_lm)))[["W"]]
    r2 = fixest::r2(rw_lm)["r2"]
    
    mat[v, ] <- c(beta, se, beta/se, if (beta > 0) { sqrt(r2) } else { -sqrt(r2) })
  }
  
  return(mat)
}
