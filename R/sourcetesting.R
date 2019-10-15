# Temporary Script to Source Necessary WADE Files During Testing
# 2019-10-15

# Scripts added as initial testing done (run with wade_cmd.R )
source(normalizePath('R/database_pipeline.R'))
source(normalizePath('R/emm.R'))
source(normalizePath('R/general_mlst.R'))

# This blast() function is only used in emm and generalmlst. DBpipeline uses a system command
source(normalizePath('R/blast.R'))