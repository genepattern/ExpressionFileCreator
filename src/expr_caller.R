# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2003-2013) by the
# Broad Institute/Massachusetts Institute of Technology. All rights are
# reserved.

# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for its
# use, misuse, or functionality.

# Bridge script to call the ExpressionFileCreator code
args <- commandArgs(trailingOnly=TRUE)

vers <- "2.15"            # R version
libdir <- args[1]
server.dir <- args[2]
patch.dir <- args[3]

source(file.path(libdir, "loadRLibrary.R"))
load.packages(libdir, patch.dir, server.dir, vers)

# Call the actual code with only the expected arguments
source(file.path(libdir, "expr.R"))
parseCmdLine(args[4:NROW(args)])
