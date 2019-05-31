get_labware_profile <- function(compare, profile, gene, delimiter){
  if(compare){
    molec_profile <- gene
  } else {
    molec_profile <- paste(profile, delimiter, gene, sep = "")
  }
  molec_profile
}
