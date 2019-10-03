sink(stdout(), type = "message")

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(purrr)
  library(here)
})

source(normalizePath("R/database_pipeline.R"))

option_list = list(
  make_option(c("--organism"),
              type="character",
              default=NULL,
              help='organism',
              metavar="character"),
  make_option(c("-t", "--test"),
              type="character",
              default=NULL,
              help='test',
              metavar="characater"),
  make_option(c("-l", "--locus"),
              type="character",
              default=NULL,
              help='locus',
              metavar="character"),
  make_option(c("-s","--samples"),
              type="character",
              default=c(),
              help="samples to analyse (comma separated list of files)",
              metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$samples)) {
  print_help(opt_parser)
  stop("Samples is a required argument", call.=FALSE)
}

samples <- unlist(strsplit(opt$samples, ","))

samples <- samples %>% map(~ normalizePath(.x))

samples <- samples %>% map_df(~ data.frame(fullpath = .x,
                             parent_dir = base::dirname(base::dirname(.x)),
                             subdir_id = base::basename(base::dirname(.x)),
                             filename = base::basename(.x),
                             type = "file")
)

test <- 'AMR_DB'

org <- 'GAS'

locus <- 'list'

switch(test,
       AMR_DB = { database_pipeline(org, samples, FALSE) },
       AMR_LW = { labware_gono_amr() },
       EMM = { emm(org, samples, locus) },
       MASTER = { master_blastr(org, test, samples, locus) },
       MLST = { general_mlst_pipeline(org, samples, locus, test) },
       NGSTAR = { general_mlst_pipeline(org, samples, locus, test) },
       NGMAST = { general_mlst_pipeline(org, samples, locus, test) },
       rRNA23S = { rna_23s(org, samples) },
       SERO = { PneumoCaT_pipeline(samples) },
       VFDB = { database_pipeline(org, samples, TRUE) },
       { master_blastr(org, test, samples, locus) }
)