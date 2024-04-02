#!/usr/bin/env Rscript

library(shiny)
library(readxl)
library(ggplot2)
library(plyr)

runGitHub( "ConlonLab-MMExpression", "drlaurenwasson", launch.browser = getOption("shiny.launch.browser", interactive()))


