inputTrimer = function(raw_sample_path, raw_virus_path) {
  trimed_virus_file = P('tmp/1_virus.csv')
  sample_list = read_excel(raw_sample_path) %>% dplyr::select(sample)
  virus_list  = read_excel(raw_virus_path) %>%
    filter(sample %in% sample_list$sample) %>%
    dplyr::select(sample, vfamily, vtype, vsname, abundance)
  trim_positive = dcast(virus_list, sample ~ vtype)
  trim_positive = cbind(trim_positive, ALL = rowSums(trim_positive %>% dplyr::select(!sample)))
  trim_all = sample_list %>%
    left_join(trim_positive, by='sample')
  trim_all[is.na(trim_all)] = 0
  write.csv(trim_all, trimed_virus_file, row.names=FALSE)
}

getGAMinput <- function() {
  #in
  clean_sample = P('data/0_bat_samples_959.csv')
  path_gps = P('data/1_GPS_2_ecology_38.csv')
  path_host= P('data/2_Genus_2_traits_8.csv')
  path_seq = P('data/3_id_2_sequence_959.csv')
  path_neighbor = P('data/4_nighbors959^2.csv')
  #out
  path_virus = P('tmp/1_virus.csv')
  data_temp = P('tmp/GAM_INPUT.csv')
  raw_samples = read.csv(clean_sample)
  gps = read.csv(path_gps)
  host = read.csv(path_host)
  seqence = read.csv(path_seq)
  neighbor = read.csv(path_neighbor)
  virus = read.csv(path_virus)
  samples = virus$sample
  
  neighbor_stat = neighbor %>%
    filter(root_ID %in% samples | neighbor %in% samples) %>%
    left_join(raw_samples, by=c('neighbor'='sample_id')) %>%
    group_by(root_ID) %>%
    mutate(neighborSampleNum = n()) %>%
    left_join(neighbor %>% 
                left_join(raw_samples, by=c('neighbor'='sample_id')) %>%
                distinct(root_ID, genus) %>%
                group_by(root_ID) %>%
                mutate(neighborBatGenus = n()) %>% 
                distinct(root_ID, neighborBatGenus), 
              by = 'root_ID') %>%
    left_join(neighbor %>% 
                left_join(virus, by=c('neighbor'='sample')) %>%
                mutate(ifDNA = ifelse(DNA!=0, 1, 0),
                       ifRNA = ifelse(RNA!=0, 1, 0)) %>%
                distinct(root_ID, neighbor, ifDNA, ifRNA) %>%
                group_by(root_ID) %>%
                mutate(neighborDNAvirusMean = sum(ifDNA)/n(),
                       neighborRNAvirusMean = sum(ifRNA)/n()) %>%
                distinct(root_ID, neighborDNAvirusMean, neighborRNAvirusMean), 
              by = 'root_ID') %>%
    distinct(root_ID,
             neighborSampleNum,
             neighborBatGenus,
             neighborDNAvirusMean,
             neighborRNAvirusMean)
  neighbor_stat[is.na(neighbor_stat)] = 0
  
  dataset = raw_samples %>%
    filter(sample_id %in% samples) %>%
    mutate(GPS_adjusted = paste(lat, lon, sep = ',')) %>%
    dplyr::select(!c(lat, lon)) %>%
    left_join(virus, by=c('sample_id' = 'sample')) %>%
    left_join(gps, by='GPS_adjusted') %>%
    left_join(host, by='genus') %>%
    left_join(seqence, by='sample_id') %>%
    left_join(neighbor_stat, by=c('sample_id' = 'root_ID')) %>%
    dplyr::select(!c('sample_id', 'GPS_adjusted'))
  write.csv(dataset, data_temp, row.names=FALSE)
}