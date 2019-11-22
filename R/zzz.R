# Fixing code shortcomings without complete rewrite


# Support functions
# Maybe put read_vcf and blastcmd in here some day


`%notin%` <- Negate(`%in%`)


# .onAttach <- function(libname, pkgname) {
#   packageStartupMessage("Welcome to WADE")
# }



# On load create global variables:
.onLoad <- function(libname, pkgname){
  utils::globalVariables(c("out_location", 
                           "Locus_id",
                           "Type",
                           "Ident",
                           "Mismatches",
                           "Gaps",
                           "Align",
                           "Allele2",
                           "gastests",
                           "gonotests",
                           "pneumotests"))
}



