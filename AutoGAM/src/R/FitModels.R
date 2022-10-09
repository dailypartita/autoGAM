fitGAM <- function(fit_term) {
  path_data_set = P('temp/GAM_INPUT.csv')
  set.seed(0)
  gam_dir = P(paste('output', today(), sep = ''))
  gam_file = P(paste('output', today(), '/GAM_', fit_term, '.rds', sep = ''))
  dir.create(gam_dir)
  data_raw = read.csv(path_data_set)
  names(data_raw)[names(data_raw) == "genus"] <- "Genus_"
  removeColsAllNa  = function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
  data_raw = removeColsAllNa(data_raw)
  data_raw = DropNA(data_raw)
  model_family <<- poisson
  # smooths
  smooth_terms = list(
    host_Traits1 = c("s(Mass, bs='tp', k=3)",
                     "s(RWL, bs='tp', k=3)",
                     "s(AR, bs='tp', k=3)"),
    
    host_Traits2 = c("s(Lifespan, bs='tp', k=3)"),
    
    host_Traits3 = c("s(DietBreadth, bs='tp', k=3)",
                     "s(TrophicLevel, bs='tp', k=3)"),
    
    host_Traits4 = c("s(LitterSize, bs='tp', k=3)",
                     "s(LittersPerYear, bs='tp', k=3)"),
    
    host_population = c("s(GroupSize, bs='tp', k=3)",
                        "s(PopulationGrpSize, bs='tp', k=3)"),
    
    ecology_human_pop = c("s(TotalPopulationSize, bs='tp', k=3)",
                          "s(DistanceToClosestTown, bs='tp', k=3)"),
    
    ecology_human_land = c("s(CroplandSize, bs='tp', k=3)",
                           "s(PastureSize, bs='tp', k=3)"),
    
    ecology_abio = c("s(Elevation, bs='tp', k=3)"),
    
    ecology_bio = c("s(MammalSympatry, bs='tp', k=3)"),
    
    seq_data = c("s(CleanDataSize, bs='tp', k=3)",
                 "s(CleanQ30, bs='tp', k=3)",
                 "s(TotalContigLength, bs='tp', k=3)",
                 "s(ContigCount, bs='tp', k=3)"),
    
    neighbor_data1 = c("s(neighborSampleNum, bs='tp', k=3)",
                       "s(neighborBatGenus, bs='tp', k=3)"),
    
    neighbor_data2 = c("s(neighborDNAvirusMean, bs='tp', k=3)",
                       "s(neighborRNAvirusMean, bs='tp', k=3)")
  )
  # dummys
  dummys = as.data.frame(with(data_raw, model.matrix(~Genus_))[,-1])
  dummy_terms = paste0("s(", names(dummys), ", bs = 're')")
  names(dummy_terms) = names(dummys)
  data_set <<- cbind(data_raw, dummys)
  
  terms = c(smooth_terms, dummy_terms)
  gam_item  = fit_all_gams(data_set,
                           fit_term,
                           poisson,
                           terms)
  print(paste("GAM output write to:", gam_file))
  saveRDS(gam_item, gam_file)
  return(gam_file)
}