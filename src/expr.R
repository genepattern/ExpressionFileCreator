# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2003-2007) by the
# Broad Institute/Massachusetts Institute of Technology. All rights are
# reserved.

# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for its
# use, misuse, or functionality.

# this file contains a collection of functions to go from probe level data (Cel files) to expression measures 

string.to.boolean <- function(s) {
	if(s=="yes") {
		return(TRUE)	
	}
	return(FALSE)
}

message <- function (..., domain = NULL, appendLF = TRUE) {

}
do.install.packages <<- T

zip.file.name <- ''
output.data.file.name <- ''
clm.input.file <- ''
output.cls.file.name <- ''
probe.descriptions.file.name <- ''
exec.info <- 'gp_module_execution_log.txt'
mycdfenv <<- NULL
cdf.file <<- NULL

cleanup <- function() {
	files <- dir()
	for(file in files) {
		if(file != zip.file.name && file!=cdf.file && file!=output.data.file.name && file!=output.cls.file.name && file!=clm.input.file && file!=exec.info && file!=probe.descriptions.file.name && file!="stdout.txt" && file!="stderr.txt" && file != "cmd.out") {
         unlink(file, recursive=T)
      }
	}		
}

parseCmdLine <- function(...) {
	suppressMessages(.parseCmdLine(...))
}

.parseCmdLine <- function(...) {
	args <- list(...)
	input.file.name <- ''
	output.file.name <- ''
	method <- ''
	quantile.normalization <- ''
	background <- ''
	compute.calls <- ''
	normalization.method <- ''
	libdir <- ''
	clm.input.file <- ''
	annotate.probes <- ''
	value.to.scale.to <- NULL

	for(i in 1:length(args)) {
		flag <- substring(args[[i]], 0, 2)
		value <- substring(args[[i]], 3, nchar(args[[i]]))
		if(flag=='-i') {
			input.file.name <- value
			zip.file.name <<- input.file.name # for cleanup
		} else if(flag=='-o') {
			output.file.name <- value
		} else if(flag=='-m') {
			method <- value # RMA, mas5, dchip
		} else if(flag=='-q') {
			quantile.normalization <- value
		} else if(flag=='-b') { # whether to background correct when using RMA
			background <- value
		} else if(flag=='-c') {
			compute.calls <- value
		} else if(flag=='-n') {
			normalization.method <- value # none, mean, median, linear fit
			if(normalization.method=='') {
				normalization.method <- 'none'
			}
		} else if(flag=='-f') {
			clm.input.file <- value
		} else if(flag=='-l') {
			libdir <- value
      } else if(flag=='-v') {
         if(value!='') {
            value.to.scale.to <- as.integer(value)
         }  	
      } else if(flag=='-a') {
         annotate.probes <- value
      } else if(flag=='-e') {
         cdf.file <<- value
		} else {
			stop(paste("unknown flag ", flag, " value ", value, sep=""), .call=FALSE)
		} 
		
	}
	
	
	create.expression.file(input.file.name=input.file.name, output.file.name=output.file.name, method=method, quantile.normalization=quantile.normalization, background=background, compute.calls=compute.calls, normalization.method=normalization.method, clm.input.file=clm.input.file, libdir=libdir, value.to.scale.to=value.to.scale.to, annotate.probes=annotate.probes, cdf.file=cdf.file)

}

create.expression.file <- function(input.file.name, output.file.name, method, quantile.normalization, background, compute.calls, normalization.method, clm.input.file, libdir, value.to.scale.to, annotate.probes, cdf.file)  {
	source(paste(libdir, "common.R", sep=''))

	DEBUG <<- F
	zip.file.name <<- input.file.name # for cleanup
	quantile.normalization <- string.to.boolean(quantile.normalization)
	background <- string.to.boolean(background)
	compute.calls <- string.to.boolean(compute.calls)
	annotate.probes <- string.to.boolean(annotate.probes)
	clm.input.file <<- clm.input.file

	if(libdir!='') {
		setLibPath(libdir)
		on.exit(cleanup())
		install.required.packages(libdir, method)
	}
	
	suppressMessages(library(tools))
	suppressMessages(library(Biobase))
	suppressMessages(library(affy))
	suppressMessages(library(affyio))
	
	dataset <- NULL # list containing data and calls if isRes is true
	isRes <- compute.calls
	
	clm <- NULL
	if(clm.input.file!='') {
		clm <- read.clm(clm.input.file)
	}

	cel.file.names <- NULL
	if(input.file.name!='') {
		cel.file.names <- get.celfilenames(input.file.name)
		if(!is.null(clm)) {
			scan.names <- clm$scan.names
			new.cel.file.names <- vector("character")
			i <- 1
			scanIdx <- 1
			remove.scan.index <- c()
			for(scan in scan.names) {
				if(length(grep('cel$', scan, ignore.case=T, extended=F)) == 0) { # check if scan name ends with .cel
					s1 <- paste('^', scan, '.cel', "$",sep='')
					s2 <- paste('^', scan, '.cel.gz', "$",sep='')
					s <- paste(s1, "|", s2, sep="")
				} else {
					s <- paste('^', scan, "$",sep='')
				}
				index <- grep(s, cel.file.names, ignore.case=T, extended=T)

				if(length(index) == 0) 
				{
					cat(paste("Scan", scan, "in clm file was not found. \n"))
			        remove.scan.index <- c(remove.scan.index, scanIdx)
				} 
				else if(length(index)>1) 
				{
					cat(paste("Scan", scan, "in clm file matches more than one CEL file. \n"))
			        remove.scan.index <- c(remove.scan.index, scanIdx)
				} 
				else 
				{
				    new.cel.file.names[i] <- cel.file.names[index[1]]
					i <- i + 1
				}	
				scanIdx <- scanIdx + 1
			} # for
			cel.file.names <- new.cel.file.names
			
			#remove duplicate or missing scan names
			if(length(remove.scan.index) !=0)
			{
                clm$scan.names <- clm$scan.names[-remove.scan.index]
                clm$sample.names <- clm$sample.names[-remove.scan.index]
                
                if(!is.null(clm$factor))
                {
			        clm$factor <- clm$factor[-remove.scan.index, drop=TRUE]
			    }
			}
		} # clm
	} else {
		exit("Either a zip of CEL files or a clm file is required.")
	}

	if(length(cel.file.names) == 0) {
	   exit("No CEL files listed in clm file found.")
	}

	chip <- read.celfile.header(cel.file.names[[1]])$cdfName
	for(c in cel.file.names) {
	   if(chip != read.celfile.header(c)$cdfName) {
	      exit("CEL files are from different chips.")
	   } 
	}
	
	if(!is.null(cdf.file) && cdf.file!='') {
		if(do.install.packages && !is.package.installed(libdir, "makecdfenv")) {
			install.package(libdir, "makecdfenv_1.14.0.zip", "makecdfenv_1.14.0.tgz", "makecdfenv_1.14.0.tar.gz")
		}
		library(makecdfenv)
		mycdfenv <<- make.cdf.env(filename=basename(cdf.file), cdf.path=dirname(cdf.file), verbose=T)
		
	  # c <- read.cdffile.list(filename=basename(cdf.file),     cdf.path=dirname(cdf.file))
	  # custom.chip.name <- c$Chip$Name
	   
		#index <- grep(chip, custom.chip.name, ignore.case=T)
		#if(length(index) == 0) {
		 #  exit("The custom cdf file provided does not appear to be from the ", chip, " chip.")	
		#}
	}

	compressed <- is.compressed(cel.file.names)
	if(method=='dChip' || method=='RMA' || method=='GCRMA') {
		if(method=='dChip') {
			dataset <- gp.dchip(cel.file.names, compressed, compute.calls)
		} else if(method=='RMA'){
			dataset <- gp.rma(cel.file.names, compressed, quantile.normalization, background, compute.calls)
		} else if(method=='GCRMA') {
			library(gcrma, verbose=FALSE)
			dataset <- gp.gcrma(cel.file.names, compressed, quantile.normalization, compute.calls)
		}
		
	} else if(method=='MAS5'){
		dataset <- gp.mas5(cel.file.names, compressed, compute.calls)
	} else if(method=="FARMS") {
		dataset <- gp.farms(cel.file.names, compressed, compute.calls, libdir)
	} else {
		exit('Unknown method')	
	}
	
	if(!is.null(clm)) { 
		if(!is.null(clm$sample.names)) {
			colnames(dataset$data) <- clm$sample.names
		}
		
		factor <- clm$factor
		if(!is.null(factor)) {
			output.cls.file.name <<- get.cls.file.name(output.file.name)		
			write.factor.to.cls (factor, output.cls.file.name)
		}

	} else {  # remove .cel extension from names
		col.names <- colnames(dataset$data)		
		col.names <- sub(".[cC][eE][lL].gz$|.[cC][eE][lL]$", "", col.names)
		colnames(dataset$data) <- col.names
	}

	if(annotate.probes) {
		cdf <- whatcdf(filename=cel.file.names[1], compress=compressed)
		row.descriptions <- try(get.row.descriptions.csv(libdir, dataset$data, cdf))
		if(!is.null(row.descriptions) && class(row.descriptions)!="try-error") {
			dataset$row.descriptions <- row.descriptions
		}
	}

	dataset$column.descriptions <- vector("character", length=length(ncol(dataset$data)))

	if(method=='MAS5') {
		dataset <- gp.normalize(dataset=dataset, method=normalization.method, value.to.scale.to=value.to.scale.to)
	}
	if(isRes) {
		output.data.file.name <<- write.res(dataset, output.file.name)
	} else {
		output.data.file.name <<- write.gct(dataset, output.file.name)
	}
	if(clm.input.file!='') {
		return(list(output.data.file.name, output.cls.file.name))
	} else {
		return(list(output.data.file.name))
	}
}


####################################################
# EXPRESSION MEASURES

#normalize: logical value. If 'TRUE' normalize data using quantile
 #         normalization

#background: logical value. If 'TRUE' background correct using RMA
 #         background correction

#bgversion: integer value indicating which RMA background to use 1: use
 #         background similar to pure R rma background given in affy
  #        version 1.0 - 1.0.2 2: use background similar to pure R rma
   #       background given in affy version 1.1 and above
	
gp.rma <- function(cel.files, compressed, normalize, background, compute.calls=FALSE) {
    # samplenames <- gsub("^/?([^/]*/)*", "", unlist(cel.files), 
    # extended = TRUE)
    # n <- length(cel.files)
    # pdata <- data.frame(sample = 1:n, row.names = samplenames)
    # phenoData <- new("phenoData", pData = pdata, varLabels = list(sample = "arbitrary numbering"))

    if(is.null(mycdfenv)) {
        eset <- just.rma(filenames=cel.files, compress=compressed, normalize=normalize, background=background, verbose=FALSE)
    } else {
        eset <- just.rma(filenames=cel.files, compress=compressed, normalize=normalize, background=background, verbose=FALSE, cdfname='mycdfenv')   
    }
    
    data <- exprs(eset)
    data <- 2^data # rma produces values that are log scaled
    if(!compute.calls) {
        return(list(data=data))
    } else {
        r <- ReadAffy(filenames=cel.files, compress=compressed) 
        if(!is.null(mycdfenv)) {
            r@cdfName <- 'mycdfenv'
        }
        calls <- get.calls(r)
        return(list(data=data, calls=calls)) 
    }
}



gp.gcrma <- function(cel.files, compressed, normalize, compute.calls=FALSE) { 
	# samplenames <- gsub("^/?([^/]*/)*", "", unlist(cel.files), 
    #        extended = TRUE)
   # n <- length(cel.files)
	# pdata <- data.frame(sample = 1:n, row.names = samplenames)
	# phenoData <- new("phenoData", pData = pdata, varLabels = list(sample = "arbitrary numbering"))
	
	r <- read.affybatch(filenames=cel.files) 
	if(!is.null(mycdfenv)) {
		r@cdfName <- 'mycdfenv'
	}
	
	eset <- gcrma(object=r, normalize=normalize, verbose=F)
   #eset <- just.gcrma(filenames=cel.files, normalize=normalize, verbose=F)   
  
   data <- exprs(eset)
   data <- 2^data # produces values that are log scaled
   if(!compute.calls) {
   	return(list(data=data))
   } else {
   	r <- ReadAffy(filenames=cel.files, compress=compressed) 
   	if(!is.null(mycdfenv)) {
   		r@cdfName <- 'mycdfenv'
   	}
   	calls <- get.calls(r)
		return(list(data=data, calls=calls)) 
   }
}

gp.dchip <- function(cel.file.names, compressed, compute.calls=FALSE) {
	afbatch <- ReadAffy(filenames=cel.file.names, compress=compressed)
	if(!is.null(mycdfenv)) {
   		afbatch@cdfName <- 'mycdfenv'
   	}
   	
	eset <- expresso(afbatch, normalize.method="invariantset", bg.correct=FALSE, pmcorrect.method="pmonly",summary.method="liwong", verbose=FALSE) 
	data <- exprs(eset)
	if(!compute.calls) {
		return(list(data=data))
	} else {
		calls <- get.calls(afbatch)
		res <- list(data=data, calls=calls) 
		return(res)
	}
}

# r AffyBatch object
get.calls <- function(r) {
	calls.eset <- try(mas5calls.AffyBatch(r, verbose = FALSE))
	if(class(calls.eset)=="try-error") {
	   return(NULL)
	}
	calls <- exprs(calls.eset)
	return(calls)
}

gp.mas5 <- function(cel.file.names, compressed, compute.calls) {
	r <- read.affybatch(filenames=cel.file.names) 
	if(!is.null(mycdfenv)) {
		r@cdfName <- 'mycdfenv'
	}
	eset <- mas5(r, normalize=FALSE)
	
	data <- exprs(eset)
	if(!compute.calls) {
		return(list(data=data))
	} else {
		calls <- get.calls(r)
		res <- list(data=data, calls=calls) 
		return(res)
	}
}

gp.farms <- function(cel.file.names, compressed, compute.calls, libdir)  {
	if(do.install.packages && !is.package.installed(libdir, "farms")) {
		install.package(libdir, "farms.zip", "farms_1.0.0.tar.gz", "farms_1.0.0.tar.gz")
	}
	library(farms)
	r <- ReadAffy(filenames=cel.file.names, compress=compressed) 
	if(!is.null(mycdfenv)) {
		r@cdfName <- 'mycdfenv'
	}
	eset <- exp.farms(r, bgcorrect.method = "none", pmcorrect.method = "pmonly", normalize.method = "quantiles")
	data <- exprs(eset)
	if(!compute.calls) {
		return(list(data=data))
	} else {
		calls <- get.calls(r)
		res <- list(data=data, calls=calls) 
		return(res)
	}
}

########################################################
# MISC FUNCTIONS

install.required.packages <- function(libdir, method) {
	#if(!do.install.packages) {
	#   return
	#}
	if(!is.package.installed(libdir, "Biobase"))
	{
		install.package(libdir, "Biobase_2.2.2.zip", "Biobase_2.2.2.tgz", "Biobase_2.2.2.tar.gz")
	}
	
	if(!is.package.installed(libdir, "affyio"))
	{
		install.package(libdir, "affyio_1.14.0.zip", "affyio_1.14.0.tgz", "affyio_1.14.0.tar.gz")
	}

	if(!is.package.installed(libdir, "preprocessCore"))
	{
		install.package(libdir, "preprocessCore_1.4.0.zip", "preprocessCore_1.4.0.tgz", "preprocessCore_1.4.0.tar.gz")
	}

	if(!is.package.installed(libdir, "affy"))
	{
		install.package(libdir, "affy_1.20.2.zip", "affy_1.20.2.tgz","affy_1.20.2.tar.gz")
	}	

    if(!is.package.installed(libdir, "makecdfenv"))
    {
		install.package(libdir, "makecdfenv_1.20.0.zip", "makecdfenv_1.20.0.tgz", "makecdfenv_1.20.0.tar.gz")
	}

	if(method=='GCRMA' && !is.package.installed(libdir, "matchprobes"))
	{
		#if(isMac()) {
		#	Sys.putenv(MAKEFLAGS="LIBR= SHLIB_LIBADD= LIBS=")
		#}

		install.package(libdir, "matchprobes_1.14.1.zip", "matchprobes_1.14.1.tgz", "matchprobes_1.14.1.tar.gz")
	}

    #commented out to use already installed gcrma_2.15.0.tar.gz
	#if(method=='GCRMA' && !is.package.installed(libdir, "gcrma"))
	#{
		#if(isMac()) {
		#	Sys.putenv(MAKEFLAGS="LIBR= SHLIB_LIBADD= LIBS=")
		#}

	#    if(length(grep(R.version$os, "darwin9")) != 0)
	#    {
	#        mac_package <- "gcrma_2.16.0_leopard.tgz"
	#    }
	#    else
	#    {
	#        mac_package <- "gcrma_2.16.0_tiger.tgz"
	#    }

	#	install.package(libdir, "gcrma_2.16.0.zip", mac_package, "gcrma_2.16.0.tar.gz")
	#}
}


get.cls.file.name <- function(data.output.file.name) {
	gct.ext <- regexpr(paste(".gct","$",sep=""), tolower(data.output.file.name))
	if(gct.ext[[1]] != -1) {
		cls <- substring(data.output.file.name, 0, gct.ext[[1]]-1)
		return(paste(cls, ".cls", sep=""))
	}
	
	res.ext <- regexpr(paste(".res","$",sep=""), tolower(data.output.file.name))
	if(res.ext[[1]] != -1) {
		cls <- substring(data.output.file.name, 0, res.ext[[1]]-1)
		return(paste(cls, ".cls", sep=""))
	}
	return(paste(data.output.file.name, ".cls", sep=""))
}

get.celfilenames <- function(input.file.name) {
	if(!file.info(input.file.name)[['isdir']]) {
		isWindows <- Sys.info()[["sysname"]]=="Windows"
		if(isWindows) {
			zip.unpack(input.file.name, dest=getwd())
		} else {
			 zip <- getOption("unzip")
			 system(paste(zip, " -q '", input.file.name, "'", sep=''))
		}
		
		files <- list.files()
		
		for(file in files) {
		   if(file.info(file)[['isdir']]) {  # move all CEL files to working directory
		      subfiles <- list.files(file, full.names=T)
		      subfiles.names.only <- list.files(file)

              # if the directory is not empty then move files to working directory
              if(length(subfiles) != 0)
              {
		        for(j in 1:length(subfiles))
		        {
		            file.rename(subfiles[[j]], subfiles.names.only[[j]])
		        }
		      }
		   }
		}
	
	   return(list.celfiles(path = ".", recursive=F, full.names=F))
	   
	}
	files <- list.celfiles(path = input.file.name, recursive=F, full.names=TRUE)
	return(files)
	
}

is.compressed <- function(cel.files) {
	for(f in cel.files) {
		r <- grep(".[cC][eE][lL].gz$", f, extended=F)
		if(length(r) > 0) {
			return(TRUE)
		}
	}
	return(FALSE)
}




############################################################
# NORMALIZATION METHODS

# data a data frame
rank.normalize <- function(data) {
	for(i in 1:NCOL(data)) {
		data[,i] <- rank(data[,i], ties.method="mean")	
	}
	return(data)
}

# reference.column - reference column to use for scaling. Use median column if -1
# value.to.scale.to) - if mean/median scaling, scale all column means/medians to this value
gp.normalize <- function(dataset, method, reference.column=-1, value.to.scale.to=NULL) {
	if(reference.column!=-1 && reference.column > NCOL(dataset)) {
		stop("Invalid reference column")
	}
	use.p.p.genes <- F # FIXME
	if(!(method %in% c('none', 'target signal', 'quantile normalization', 'linear fit', 'mean scaling', 'median scaling'))) {
		stop("Unknown scaling method")
	}
	data <- dataset$data
	calls <- dataset$calls
	if(method=="none") {
		return(dataset)
	}	
	
	if(method=='target signal') {
		nf <- 1
		for (i in 1:ncol(data)) {
        	slg <- data[, i]
        	sf <- sc/mean(slg, trim = 0.02)
        	reported.value <- nf * sf * slg
        	data[, i] <- reported.value
    	}
    	dataset$data <- data
    	return(dataset) 
	} 
    
	if(method=="quantile normalization") {
	    row.names <- row.names(data)
	    col.names <- colnames(data)
		data <- normalize.quantiles(data)
		row.names(data) <- row.names
		colnames(data) <- col.names
		dataset$data <- data
		return(dataset)
	}
	
	if(reference.column == -1) {
		means <- apply(data, 2, mean)
		#if(method %in% c('mean scaling', 'median scaling')) {
		#	value.to.scale.to <- median(means)
		#} else {
			median.index <- as.integer(length(means)/2) + 1
			reference.column <- which(rank(means, ties="first")==median.index)
		#}
	}
	if(!is.null(value.to.scale.to) && method %in% c('mean scaling', 'median scaling')) {
		reference.column <- -1
	}
	
	scalingFactors <- vector(mode="numeric", length=NCOL(data));
	scalingFactors[reference.column] <- 1;
	row.indices <- matrix(TRUE, nrow=nrow(data), ncol=ncol(data))
	if(!is.null(calls) && use.p.p.genes) { # FIXME
		for(j in 1:NCOL(calls)) {
			if(j==reference.column) {
				next
			}
			row.indices[,j] <- calls[,j]=='P' && calls[,reference.column]=='P'
		}
	}

	if(method=="linear fit") {
		for(j in 1:NCOL(data)) {
			if(j==reference.column) {
				next
			}
			linear.fit <- linear.fit(data[,reference.column], data[,j])
   			scalingFactors[j] <- 1.0/linear.fit$m
		}
	} else if(method=='mean scaling' || method=='median scaling') {
		if(method=='mean scaling') {
			FUN <- 'mean'
		}
		else if(method=='median scaling') {
			FUN <- 'median'
		}
		function.to.apply <- get(FUN)
		if(is.null(value.to.scale.to)) {
			ref.mean.or.median <- function.to.apply(data[,reference.column])
		} else {
			ref.mean.or.median <- value.to.scale.to
		}
		
		for(j in 1:ncol(data)) {
			if(j!=reference.column) {
				scaling.factor <- function.to.apply(data[,j])/ref.mean.or.median
				scalingFactors[j] <- 1.0/scaling.factor
			}
		}
	}

	for(j in 1:NCOL(data)) {
		if(j==reference.column) {
			next
		}
		data[,j] <- data[,j]*scalingFactors[j]
		if(!is.null(dataset$column.descriptions) && length(dataset$column.descriptions) >= j) {
			if(length(dataset$column.descriptions[j]) <= 0 || dataset$column.descriptions[j]=='') {
				prev <- ''
			} else {
				prev <- paste(dataset$column.descriptions[j], ", ", sep='')
			}
			dataset$column.descriptions[j] <- paste(prev, "scaling factor=", scalingFactors[j], sep='')
		}
	}
	dataset$data <- data
	return(dataset)
}

# data - matrix 
# cdf - the cdf file for data, used to construct the URL to download
# t - result of read.table(csv), for testing
get.row.descriptions.csv <- function(libdir, data, cdf, t=NULL) {
    file.name <- paste(cdf, ".zip", sep='')
    absolute.file.name <- paste(libdir, file.name, sep='')
	if(is.null(t) && !file.exists(file.name)) {
		url <- paste("ftp://ftp.broadinstitute.org/pub/genepattern/csv/Affymetrix/2012-annotations", file.name, sep='')
		sink(stdout(), type = "message")
		try(suppressMessages(download.file(url, quiet=T, destfile=absolute.file.name, mode="wb")))
		sink(stderr(), type = "message")
		if(!file.exists(absolute.file.name) || file.info(absolute.file.name)[['size']] == 0) {
			cat(paste("No annotations found for chip ", cdf, "\n", sep=''))
			return(NULL)
		}
		isWindows <- Sys.info()[["sysname"]]=="Windows"
		if(isWindows) {
			rc <- zip.unpack(absolute.file.name, dest=getwd())
        }
        else {
		   rc <- .Internal(int.unzip(absolute.file.name, NULL, getwd()))
        }

		csv.file <- attr(rc, "extracted")
	    on.exit(unlink(csv.file))

		desc <- as.matrix(read.table(row.names=1, file=csv.file, header=T, quote='"', comment.char='#', fill=T, sep=","))
	}

    if(file.info(absolute.file.name)[['size']] == 0)
    {
        return(NULL)
    }

	gene.title.idx <-  match('Gene.Title', colnames(desc))
	gene.symbol.idx <-  match('Gene.Symbol', colnames(desc))
	probeids <- row.names(data)

	get.gene.info <- function(probe) {
	    return(paste(desc[[probe, gene.title.idx]], ", ", desc[[probe, gene.symbol.idx ]], sep=''))
	}
	row.descriptions <- vector(mode = "character", length = length(probeids))
	for(i in 1:length(probeids)) {
	   probe <- probeids[[i]]
       ann <- try(get.gene.info(probe), silent=T)
       if(class(ann)=="try-error" || is.na(ann) || ann == ', ' || ann=="---, ---") {
        next
      
        }
        row.descriptions[i] <- ann
	}

	return(row.descriptions)
}

get.row.descriptions.annaffy <- function(cdf, data) {
	library(annaffy)
	cdf <- tolower(sub("_", "", cdf))
	if(!require(cdf, quietly=TRUE, character.only=TRUE)) {
			info(paste("installing", cdf, "..."))
			source("http://www.bioconductor.org/getBioC.R")		
			getBioC(pkgs=c(cdf))
			info(paste(cdf, "installed"))
	}	
	affy.descriptions <- aafDescription(row.names(data),chip=cdf)
	affy.descriptions <- sapply(affy.descriptions, getText)	
	gene.symbols <- aafSymbol(row.names(data), chip=cdf)
	gene.symbols <- sapply(gene.symbols, getText)
	row.descriptions <- paste(affy.descriptions, gene.symbols, sep=", ")	
	return(row.descriptions)
}

linear.fit <- function(xpoints, ypoints) {
	n <- length(xpoints)
	xBar_yBar <- sum(xpoints*ypoints)
	xBar <- sum(xpoints)
	yBar <- sum(ypoints)
	x2Bar <- sum(xpoints*xpoints)
	
	xBar_yBar <- xBar_yBar / n
	xBar <- xBar / n
	yBar <- yBar / n
	x2Bar <- x2Bar / n
	deltaX2 <- x2Bar - xBar * xBar
	m <- (xBar_yBar - xBar * yBar) / deltaX2
	b <- yBar - m * xBar
	return(list("m"=m, "b"=b))
}





