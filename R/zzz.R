# Fixing code shortcomings without complete rewrite


# Support functions
# Maybe put read_vcf and blastcmd in here some day


`%notin%` <- Negate(`%in%`)



# Create filecache environment to store filename outside of all functions
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to WADE")
}

# On load fix the notes:

.onLoad <- function(libname, pkgname){
  utils::globalVariables(c("out_location", 
                           "Locus_id",
                           "Type",
                           "Ident",
                           "Mismatches",
                           "Gaps",
                           "Align",
                           "Allele2"))
}



