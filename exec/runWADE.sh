#!/bin/sh

Rscript --vanilla -e 'library(wade); source(file = system.file("exec/wade_cmd.R", package = "wade"))' $@
