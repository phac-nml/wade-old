sink(stdout(), type = 'output')

suppressPackageStartupMessages({
  library(optparse)
  library(wade)
  library(dplyr)
  library(purrr)
  library(here)
  library(stringr)
})

#source(normalizePath('R/sourcetesting.R'))


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
              callback=check_choose_from(choices = c('MLST','VIRULENCE','EMM','AMR','VFDB','NGSTAR','rRNA23S','NGMAST')), # c('TOXINS','MLST','VIRULENCE','EMM','AMR_DB','VFDB','NGSTAR','rRNA23S','NGMAST')),
              help='test'),
  make_option(c('-l', '--locus'),
              type='character',
              default='list',
              help='locus name or list, default is list',
              metavar='character'),
  make_option(c('-s','--samples'),
              type='character',
              default=c(),
              help='samples to analyse (comma separated list of files)',
              metavar='character'),
  make_option(c('-s2','--samples2'),
              type='character',
              default=c(),
              help='VCF samples to analyse (comma separated list of files, only used for GONO/Labware)',
              metavar='character'),
  make_option(c('-d','--outdir'),
              type='character',
              default=NA,
              help='Output directory: default is ./output/',
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
                             type = 'file')) %>% 
  mutate_if(is.factor, as.character)

samples2 <- unlist(strsplit(opt$samples2, ','))
samples2 <- samples2 %>% map(~ normalizePath(.x))
samples2 <- samples2 %>% map_df(~ data.frame(fullpath = .x,
                                           parent_dir = base::dirname(base::dirname(.x)),
                                           subdir_id = base::basename(base::dirname(.x)),
                                           filename = base::basename(.x),
                                           type = 'file')) %>% 
  mutate_if(is.factor, as.character)


test <- opt$test
org <- opt$organism
locus <- opt$locus
out_location <- paste(if(is.na(opt$outdir)){ here("output")} else{normalizePath(opt$outdir)}, "/", sep="")

try(system(paste("mkdir", out_location)))

switch(test,
       AMR = { database_pipeline(org, samples, FALSE) }, # 
       VFDB = { database_pipeline(org, samples, TRUE) },
       EMM = { emm(org, samples, locus) },
       MLST = { general_mlst_pipeline(org, samples, locus, test) },
       NGSTAR = { general_mlst_pipeline(org, samples, locus, test) },
       NGMAST = { general_mlst_pipeline(org, samples, locus, test) },
       rRNA23S = { rna_23s(org, samples)}, 
       AMR_LW = { 
         if (org == "GONO") {
           labware_gono_amr(amrDF = database_pipeline(org, samples, is_vfdb = FALSE, stdout = TRUE), 
                            ngstarDF = general_mlst_pipeline(org, samples, locus, seq_type = "NGSTAR"), 
                            rnaDF = rna_23s(org, samples2)) # RNA USES VCF FILES NOT FASTA!
         } else{ stop('Organism must be GONO for AMR_LW', call.=FALSE) }
        } 
       # Tested to here
       # SERO = { PneumoCaT_pipeline(samples) }, # Needs Samplenum issue (not sample list) addressed
       # MASTER = { master_blastr(org, test, samples, locus) },
       #{ master_blastr(org, test, samples, locus) }
)