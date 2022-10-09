fit_all_gams <- function(data_set, outcome_variable, model_family, terms) {
  
  fit_gam = function(frm) {
    try(gam(formula=as.formula(frm), model_family, data_set, select=TRUE), silent=FALSE)
  }
  
  terms_grid = do.call(expand.grid, terms)
  
  #Create model forumulas from the grid
  formulas = apply(as.matrix(terms_grid), 1, function(row) paste(row, collapse = " + ")) %>%
    stri_replace_all_regex("\\s[\\+\\s]+\\s", " + ") %>%
    {paste(outcome_variable, "~", .)} %>%
    rearrange_formula %>%
    unique
  
  models = tibble(formula = formulas)
  
  n_cores = detectCores()
  n_cores_use = round(nrow(models) / ceiling(nrow(models) / (n_cores - 1)))
  options(mc.cores = n_cores_use)
  message("Using ", n_cores_use, " cores to fit ", nrow(models), " models")
  
  models_vec = mclapply(models$formula, fit_gam)
  
  models = models %>%
    mutate(model = models_vec)
  
  # Calculate models
  models = models %>%
    filter(map_lgl(model, ~ !("try-error" %in% class(.) | is.null(.)))) %>%
    mutate(aic = map_dbl(model, MuMIn::AICc),
           daic = aic - min(aic),
           weight = exp(-daic/2)) %>%
    arrange(aic)
  
  # Remove unused terms from models and reduce to unique ones
  models_reduced = models %>%
    dplyr::select(model) %>%
    mutate(formula = map_chr(model, ~rearrange_formula(rm_low_edf(.)))) %>%
    distinct(formula, .keep_all = TRUE)
  
  n_cores_use = round(nrow(models_reduced) / ceiling(nrow(models_reduced) / (n_cores - 1)))
  options(mc.cores = n_cores_use)
  message("Using ", n_cores_use, " cores to fit ", nrow(models_reduced), " reduced models")
  
  # Reduce the remaining models
  models_reduced = models_reduced %>%
    mutate(model = mclapply(model, reduce_model))
  
  
  models_reduced = models_reduced %>%
    filter(map_lgl(model, ~ !("try-error" %in% class(.) | is.null(.)))) %>%
    mutate(aic = map_dbl(model, MuMIn::AICc)) %>%
    mutate(formula = map_chr(model, ~rearrange_formula(rm_low_edf(.)))) %>%
    distinct(formula, .keep_all=TRUE) %>%
    arrange(aic) %>%
    mutate(daic = aic - min(aic),
           weight = exp(-daic/2),
           terms = shortform(map(model, ~ rearrange_formula(.$formula))),
           relweight = ifelse(daic > 2, 0, weight/sum(weight[daic < 2])),
           relweight_all = weight/sum(weight),
           cumweight = cumsum(relweight_all))
  
  return(models_reduced)
  
}

# Returns a model formula from a GAM with the low_edf terms removed
rm_low_edf <- function(mod, edf_cutoff = 0.0000000000001) {
  fr = as.character(formula(mod))
  lhs = fr[2]
  rhs = fr[3]
  edfs = pen.edf(mod)
  low_edfs = edfs[edfs < edf_cutoff]
  vars_to_remove = stri_replace_all_fixed(unique(stri_extract_first_regex(names(low_edfs), "(?<=s\\()[^\\)]+(?=\\))")), ",",", ")
  vars_regex = paste0("(", paste(vars_to_remove, collapse="|"), ")")
  new_rhs = stri_replace_all_regex(rhs, paste0("\\s*s\\(", vars_regex, "\\,[^\\)]+\\)\\s*\\+?"), "")
  new_rhs= stri_replace_all_fixed(new_rhs, "+, k = 7) ", "")
  new_rhs= stri_replace_all_fixed(new_rhs, "+ +s", "+ s")
  new_formula = paste(lhs, "~", new_rhs)
  new_formula = stri_replace_all_regex(new_formula, "[\\s\\n]+", " ")
  new_formula = stri_replace_all_regex(new_formula, "[+\\s]*$", "")
  return(new_formula)
}

#' Alphabetizes the right-hand side of a formula so as to compare formulas across models
rearrange_formula = function(formula) {
  if(class(formula) == "formula") {
    formula = as.character(formula)
    formula = paste(formula[2], "~", formula[3], collapse=" ")
  }
  lhs = stri_extract_first_regex(formula, "^[^\\s~]+")
  rhs = stri_replace_first_regex(formula, "[^~]+~\\s+", "")
  terms = stri_split_regex(rhs, "[\\s\\n]+\\+[\\s\\n]+")
  terms = lapply(terms, sort)
  new_formula = mapply(function(lhs, terms) {paste(lhs, "~", paste(terms, collapse = " + "))}, lhs, terms)
  new_formula = stri_replace_all_regex(new_formula, "[\\s\\n]+", " ")
  new_formula = stri_replace_all_fixed(new_formula, "+ +s", "+ s")
  new_formula = stri_replace_all_fixed(new_formula, "+ +", "+")
  names(new_formula) <- NULL
  return(stri_trim(new_formula))
}

# Re-fits a gam model, dropping terms that have been selected out
reduce_model <- function(mod, edf_cutoff = 0.001, recursive=TRUE) {
  low_edf_vars = any(pen.edf(mod) < edf_cutoff)
  if(recursive) {
    while(low_edf_vars) {
      mod = update(mod, formula = as.formula(rm_low_edf(mod, edf_cutoff)))
      low_edf_vars = any(pen.edf(mod) < edf_cutoff)
    }
  } else {
    if(low_edf_vars) {
      mod = update(mod, formula = as.formula(rm_low_edf(mod, edf_cutoff)))
    }
  }
  return(mod)
}

# Makes a reduced version of the RHS of a formula
shortform = function(formula) {
  rhs = stri_replace_first_regex(formula, "[^~]+~\\s+", "")
  rhs = stri_replace_all_regex(rhs, "s\\(([^\\,]+)\\,[^\\)]+\\)", "s($1)")
  rhs= stri_replace_all_fixed(rhs, "+ +s", "+ s")
  stri_replace_all_fixed(rhs, "(1 + | + 1)", "")
}

addSomeVirus = function(path_raw_virus, path_s2v, path_virus, virus_type){
  doAdd <- function(path_s2v_, path_virus_, virus_type_) {
    s2v = read.csv(path_s2v_)
    df = data.frame(read.table(path_virus_, encoding = "UTF-8"), virus=virus_type_)
    names(df) = c("Cluster_name", "virus")
    step1 = s2v %>%
      left_join(df, by="Cluster_name")
    step1_decast = dcast(step1, Sample_ID~virus)
    if ("NA" %in% names(step1_decast)) {step1_decast = step1_decast %>% dplyr::select(!c("NA"))}
    step2 = s2v %>%
      left_join(step1_decast, by="Sample_ID") %>%
      distinct(Sample_ID, .keep_all=TRUE) %>%
      dplyr::select(!Cluster_name)
    names(step2) = c("Sample_ID", paste("virus_", virus_type, sep=""))
    return(step2)
  }
  newV = doAdd(path_s2v_ = path_s2v, path_virus_ = path_virus, virus_type_ = virus_type)
  addedV = read.csv(path_raw_virus) %>%
    left_join(newV, by="Sample_ID")
  addedV[is.na(addedV)] = 0
  write.csv(addedV, path_raw_virus, row.names = FALSE)
}