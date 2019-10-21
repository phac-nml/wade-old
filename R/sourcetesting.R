# Temporary Script to Source Necessary WADE Files During Testing
# 2019-10-15

# Run this in terminal
# for o in GAS GONO PNEUMO; do echo Rscript R/wade_cmd.R -o $o -t TEST -s ~/testswade/SC10-1841-A.fasta; don


# Call library statements
library(stringr)

# Fake CL args
# opt <- list(organism = "GAS",
#             test = "rRNA23S",
#             samples = "~/testswade/test.vcf",
#             locus="list",
#             outdir=NA)
# out_location <- here("output")
# path <- "./testdata/WADE/test.vcf"

# Shared functions, 'utilities.R' eventually
`%notin%` <- Negate(`%in%`)

# Scripts added as initial testing done (run with wade_cmd.R )
source(normalizePath('R/database_pipeline.R'))
source(normalizePath('R/emm.R'))
source(normalizePath('R/general_mlst.R'))
source(normalizePath('R/rRNA.R'))
# source(normalizePath('R/PneumoCaT.R')) # Doesn't accept a list of 


# 23S RNA test needs read_vcf from this file
source(normalizePath('R/read_vcf.R'))

# This blast() function is only used in emm and generalmlst. DBpipeline uses a system command
source(normalizePath('R/blast.R'))