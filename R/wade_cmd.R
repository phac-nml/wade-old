sink(stdout(), type = 'message')

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(purrr)
  library(here)
})

source(normalizePath('R/sourcetesting.R'))


#' Check if value matches a pre-defined list of values
check_choose_from <- function(choices) {
  if (length(choices) == 0) {
    stop(paste0('Choice is empty for "', opt_flag, '".'))
  }
  function(opt, opt_flag, opt_value, parser, choices. = choices) {
    if (is.na(opt_value)) {
      opt_value <- choices.[1]
    } else if (length(opt_value) != 1) {
      stop(
        paste0('Expects a single value for "', opt_flag, '".'),
        call. = FALSE
      )
    } else if (! opt_value %in% choices.) {
      stop(
        paste0(
          'Supplied value "', opt_value, '" for "', opt_flag,
          '" is not one of [', paste(choices., collapse = ', '), '].'
        ),
        call. = FALSE
      )
    }
    opt_value
  }
}

option_list = list(
  make_option(c('-o', '--organism'),
              action='callback',
              type='character',
              default=NULL,
              metavar='character',
              callback=check_choose_from(choices = c('GAS','GONO','PNEUMO')),
              help='organism'),
  make_option(c('-t', '--test'),
              action='callback',
              type='character',
              default=NULL,
              metavar='character',
              callback=check_choose_from(choices = c('MLST','VIRULENCE','EMM','AMR_DB','VFDB','NGSTAR','rRNA23S','NGMAST')), # c('TOXINS','MLST','VIRULENCE','EMM','AMR_DB','VFDB','NGSTAR','rRNA23S','NGMAST')),
              help='test'),
  make_option(c('-l', '--locus'),
              type='character',
              default='list',
              help='locus',
              metavar='character'),
  make_option(c('-s','--samples'),
              type='character',
              default=c(),
              help='samples to analyse (comma separated list of files)',
              metavar='character'),
  make_option(c('-d','--outdir'),
              type='character',
              default=NA,
              help='Output directory',
              metavar='character')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
print(opt) # TESTING ONLY
if (is.null(opt$samples)) {
  print_help(opt_parser)
  stop('Samples is a required argument', call.=FALSE)
}

samples <- unlist(strsplit(opt$samples, ','))

samples <- samples %>% map(~ normalizePath(.x))

samples <- samples %>% map_df(~ data.frame(fullpath = .x,
                             parent_dir = base::dirname(base::dirname(.x)),
                             subdir_id = base::basename(base::dirname(.x)),
                             filename = base::basename(.x),
                             type = 'file')
)

test <- opt$test
org <- opt$organism
locus <- opt$locus
out_location <- if_else(is.na(opt$outdir), 
                        here("output"),
                        suppressWarnings(normalizePath(opt$outdir))) # Throws warnings for some reason. Seems like if_else() evaluates all code within it even if the logical dictates otherwise

switch(test,
       AMR_DB = { database_pipeline(org, samples, FALSE) }, # 
       VFDB = { database_pipeline(org, samples, TRUE) },
       EMM = { emm(org, samples, locus) },
       MLST = { general_mlst_pipeline(org, samples, locus, test) },
       NGSTAR = { general_mlst_pipeline(org, samples, locus, test) },
       NGMAST = { general_mlst_pipeline(org, samples, locus, test) },
       rRNA23S = { rna_23s(org, samples) } # Need proper testing data
       # Tested to here
       # SERO = { PneumoCaT_pipeline(samples) }, # Needs Samplenum issue (not sample list) addressed
       # MASTER = { master_blastr(org, test, samples, locus) },
       # AMR_LW = { labware_gono_amr() }, # Not to include for first run through
       #{ master_blastr(org, test, samples, locus) }
)