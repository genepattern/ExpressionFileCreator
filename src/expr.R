# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2003-2006) by the
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

zip.file.name <- ''
output.data.file.name <- ''
clm.input.file <- ''
output.cls.file.name <- ''
probe.descriptions.file.name <- ''
exec.info <- 'gp_task_execution_log.txt'

cleanup <- function() {
	files <- dir()
	for(file in files) {
		if(file != zip.file.name && file!=output.data.file.name && file!=output.cls.file.name && file!=clm.input.file && file!=exec.info && file!=probe.descriptions.file.name && file!="stdout.txt" && file!="stderr.txt") {
			info(paste("removing", file))
         unlink(paste("'", file, "'", sep=''), recursive=TRUE) # XXX directories fail to delete on windows
      }
	}		
}

get.row.descriptions <- function(data, file) {
	row.descriptions <- vector("character", length=NROW(data))
	table <- read.table(file, colClasses="character", header=FALSE, sep="\t")
	for(i in 1:NROW(table)) {
		probe <- table[i, 1]
		index <- which(row.names(data)==probe)
		if(length(index) > 0 && index > 0) {
			row.descriptions[index] <- table[i, 2]
		}
	}
	return(row.descriptions)	
}

parseCmdLine <- function(...) {
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
	row.descriptions.file <- ''
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
		} else if(flag=='-r') {
			row.descriptions.file <- value
			probe.descriptions.file.name <<- value # for cleanup
		} else if(flag=='-a') {
			if (value !='' && .Platform$OS.type == "windows") {
      		memory.limit(size=as.numeric(value))
   		}
		}  else  {
			cat(args)
			cat("\n")
			stop(paste("unknown option", flag, sep=": "), .call=FALSE)
		} 
		
	}
	create.expression.file(input.file.name=input.file.name, output.file.name=output.file.name, method=method, quantile.normalization=quantile.normalization, background=background, compute.calls=compute.calls, normalization.method=normalization.method, clm.input.file=clm.input.file, libdir=libdir, row.descriptions.file=row.descriptions.file)

}

create.expression.file <- function(input.file.name, output.file.name, method, quantile.normalization, background, compute.calls, normalization.method, clm.input.file, libdir, row.descriptions.file)  {
	source(paste(libdir, "common.R", sep=''))
	DEBUG <<- F
	info(paste("normalization.method", normalization.method))
	options("warn"=-1)
	zip.file.name <<- input.file.name # for cleanup
	quantile.normalization <- string.to.boolean(quantile.normalization)
	background <- string.to.boolean(background)
	compute.calls <- string.to.boolean(compute.calls)
	clm.input.file <<- clm.input.file

	if(libdir!='') {
		setLibPath(libdir)
		on.exit(cleanup())
		install.required.packages(libdir, method)
	}

	library(affy, verbose=FALSE)
	result <- NULL # data.frame if calls are not computed, else list containing calls and data
	isRes <- compute.calls
	
	zipFileGiven <- TRUE
	clm <- NULL
	if(clm.input.file!='') {
		clm <- read.clm(clm.input.file)
	}
	if(input.file.name!='') {
		cel.file.names <- get.cel.file.names(input.file.name)
		if(clm.input.file!='') {
			scan.names <- clm$scan.names
			new.cel.file.names <- vector("character")
			i <- 1
			for(scan in scan.names) {
				s1 <- paste(scan, '.cel', "$",sep='')
				s2 <- paste(scan, '.cel.gz', "$",sep='')
				s <- paste(s1, "|", s2, "|", scan, sep="")
				index <- grep(s, cel.file.names, ignore.case=TRUE)
				if(length(index) == 0) {
					cat(paste("Scan ", scan, "in clm file not found.\n"))
				} else if(length(index)>1) {
					cat(paste("Scan ", scan, "in clm file matches more than one CEL file.\n"))
				} else {
					new.cel.file.names[i] <- cel.file.names[index[1]]
					i <- i + 1
				}
			}
			cel.file.names <- new.cel.file.names
			info(paste("cel.file.names", cel.file.names))
		}
	} else if(clm.input.file!='') {
		if(output.file.name=='') {
			output.file.name <- sub(".clm$", '', clm.input.file, ignore.case=TRUE)
		}
			
		zipFileGiven <- FALSE
		scan.names <- clm$scan.names
		cel.file.names <- paste(scan.names, ".CEL", sep='')
		bzip.names <- paste("/xchip/data/Affy/Genechip/RawCelFiles/", cel.file.names, ".bz2", sep='')
		for(i in 1:length(bzip.names)) {
			bz <- bzip.names[i]
			cel <- cel.file.names[i]
			if(!file.exists(bz)) {
				exit(paste("Unable to find scan", cel))
			}
			cmd <- paste("bzcat ", bz, " > ", cel, sep='')
			info(cmd)
			system(cmd)
		}
	} else {
		exit("Either a zip of CEL files or a clm file is required.")
	}
	
	is.compressed <- is.compressed(cel.file.names)
	
	if(method=='dChip' || method=='RMA' || method=='GCRMA') {
		if(method=='dChip') {
			info("running dChip")
			result <- gp.dchip(cel.file.names, is.compressed, compute.calls)
		} else if(method=='RMA'){
			info("running rma")
			info(paste("cel.file.names", cel.file.names))
			info(paste("is.compressed", is.compressed))
			info(paste("is.compressed", is.compressed))
			info(paste("quantile.normalization", quantile.normalization))
			info(paste("background", background))
			result <- gp.rma(cel.file.names, is.compressed, quantile.normalization, background, compute.calls)
		} else if(method=='GCRMA') {
			library(gcrma, verbose=FALSE)
			result <- gp.gcrma(cel.file.names, is.compressed, quantile.normalization, compute.calls)
		}
		info(paste("Finished running", method))
	} else if(method=='MAS5'){
		result <- gp.mas5(cel.file.names, is.compressed, compute.calls)
		info("Finished running mas5")
	} else if(method=="FARMS") {
		result <- gp.farms(cel.file.names, is.compressed, compute.calls, libdir)
		info("Finished running FARMS")
	} else {
		exit('Unknown method')	
	}
	
	if(!is.null(clm)) { 
		if(!is.null(clm$sample.names)) {
			info("setting names")
			info(clm$sample.names)
			if(isRes) {
				colnames(result$data) <- clm$sample.names
			} else {
				colnames(result) <- clm$sample.names	
			}
		}
		
		factor <- clm$factor
		if(!is.null(factor)) {
			cls <- list(labels=factor,names=levels(factor))
			info(paste("cls: ", cls))
			output.cls.file.name <<- get.cls.file.name(output.file.name)
			info(paste("saving cls file to", output.cls.file.name))
			write.cls(cls, output.cls.file.name)
			info("cls file saved")
		}
			
	} else {  # remove .cel extension from names
		if(isRes) {
			col.names <- colnames(result$data)
		} else {
			col.names <- colnames(result)
		}
		info(paste("removing .cel extension from sample names", colnames(result)))
		col.names <- sub(".[cC][eE][lL].gz$|.[cC][eE][lL]$", "", col.names)
		info(paste("sample names", col.names))
		info(paste("isRes", isRes))
		if(isRes) {
			colnames(result$data) <- col.names
		} else {
			colnames(result) <- col.names
		}
	}
	
	
	includeDescriptionsFromR <- FALSE
	row.descriptions <- NULL
	if(includeDescriptionsFromR) {
		if(!require("annaffy", quietly=TRUE, character.only=TRUE)) {
			source("http://www.bioconductor.org/getBioC.R")
			getBioC(pkgs=c("annaffy"), lib=libdir)
		}
		library(annaffy)
		cdf <- cleancdfname(afbatch@cdfName) # e.g "hgu133acdf"
		cdf <- substring(cdf, 0, nchar(cdf)-3)# remove cdf from end
		row.descriptions <- aafDescription(row.names(data),chip=cdf)
		row.descriptions <- sapply(row.descriptions, getText)
	} else {
		if(row.descriptions.file!='') {
			if(isRes) {
				row.descriptions <- get.row.descriptions(result$data, 	row.descriptions.file)
			} else {
				row.descriptions <- get.row.descriptions(result, 	row.descriptions.file)
			}
		}
	}
	
	if(isRes) {
		result$row.descriptions <- row.descriptions
		struc <- result
		struc$column.descriptions <- vector("character", length=length(ncol(result$data)))
	} else {
		struc <- list(data=result, row.descriptions=row.descriptions)		
	}

	if(method=='MAS5') {
		result <- gp.normalize(struc, normalization.method)
	}
	if(isRes) {
		output.data.file.name <<- write.res(struc, output.file.name)
	} else {
		output.data.file.name <<- write.gct(struc, output.file.name)
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
	info("creating samplenames")
	samplenames <- gsub("^/?([^/]*/)*", "", unlist(cel.files), 
            extended = TRUE)
   info(paste("samplenames", samplenames))
   n <- length(cel.files)
	pdata <- data.frame(sample = 1:n, row.names = samplenames)
	phenoData <- new("phenoData", pData = pdata, varLabels = list(sample = "arbitrary numbering"))
   
   info(paste("normalize", normalize))
   info(paste("background", background))
   eset <- just.rma(filenames=cel.files, compress=compressed, normalize=normalize, background=background, verbose=FALSE, phenoData=phenoData)   
   eset@exprs <- 2^eset@exprs # rma produces values that are info scaled
   data <- exprs(eset)
   if(!compute.calls) {
   	return(data)
   } else {
   	r <- ReadAffy(filenames=cel.files, compress=compressed) 
   	calls <- get.calls(r)
		res <- list(data=data, calls=calls) 
		return(res)
   }
}



gp.gcrma <- function(cel.files, compressed, normalize, compute.calls=FALSE) { 
	info("creating samplenames")
	samplenames <- gsub("^/?([^/]*/)*", "", unlist(cel.files), 
            extended = TRUE)
   info(paste("samplenames", samplenames))
   n <- length(cel.files)
	pdata <- data.frame(sample = 1:n, row.names = samplenames)
	phenoData <- new("phenoData", pData = pdata, varLabels = list(sample = "arbitrary numbering"))
	
 #  cdf <- cleancdfname(afbatch@cdfName) # e.g "hgu133acdf"
 #  cdf <- substring(cdf, 0, nchar(cdf)-3)# remove cdf from end
 #  pkg <- paste(cdf, "probe", sep='') 
  # if(!is.package.installed(libdir, pkg)) {
  # if(!library(package=pkg, lib.loc=libdir, logical.return=TRUE, version=NULL)) {
   #   repList <- getOptReposList()
    #  n <- repNames(repList)
      
     # get.index <- function(n) {
      #   for(index in 1:length(n)) {
       #     if(n[index] =='Bioconductor Probe Data Packages') {
        #       return(index)
         #   }
       #  }
       #  exit("Probe Data repository not found.")
      # }
      # isWindows <- Sys.info()[["sysname"]]=="Windows"
       
      # if(isWindows) {
      #   type <- "Win32"  
      # } else {
      #    type <- "Source"
      # }
      # r <- getRepEntry(repList, get.index(n))
    #  info(paste("installing package", pkg))
      # install.packages2(repEntry=r, pkgs=c(pkg), type=type)
   #}
   eset <- just.gcrma(filenames=cel.files, compress=compressed, normalize=normalize, verbose=TRUE, phenoData=phenoData)   
   eset@exprs <- 2^eset@exprs # produces values that are log scaled
   data <- exprs(eset)
   
   if(!compute.calls) {
   	return(data)
   } else {
   	r <- ReadAffy(filenames=cel.files, compress=compressed) 
   	calls <- get.calls(r)
		res <- list(data=data, calls=calls) 
		return(res)
   }
}

gp.dchip <- function(cel.file.names, is.compressed, compute.calls=FALSE) {
	info("running dchip")
	afbatch <- ReadAffy(filenames=cel.file.names, compress=is.compressed)
	eset <- expresso(afbatch, normalize.method="invariantset", bg.correct=FALSE, pmcorrect.method="pmonly",summary.method="liwong", verbose=FALSE) 
	data <- exprs(eset)
	if(!compute.calls) {
		return(data)
	} else {
		calls <- get.calls(afbatch)
		res <- list(data=data, calls=calls) 
		return(res)
	}
}

# r AffyBatch object
get.calls <- function(r) {
	calls.eset <- mas5calls.AffyBatch(r, verbose = FALSE) 
	calls <- exprs(calls.eset)
	return(calls)
}

gp.mas5 <- function(cel.file.names, is.compressed, compute.calls) {
	r <- ReadAffy(filenames=cel.file.names, compress=is.compressed) 
	info("running mas5...")
	eset <- mas5(r, normalize=FALSE)
	data <- exprs(eset)
	if(!compute.calls) {
		return(data)
	} else {
		calls <- get.calls(r)
		res <- list(data=data, calls=calls) 
		return(res)
	}
}

gp.farms <- function(cel.file.names, is.compressed, compute.calls, libdir)  {
	if(!is.package.installed(libdir, "farms")) {
		info("installing farms")
		install.package(libdir, "farms.zip", "farms_1.0.0.tar.gz", "farms_1.0.0.tar.gz")
	}
	library(farms)
	r <- ReadAffy(filenames=cel.file.names, compress=is.compressed) 
	info("running farms...")
	eset <- exp.farms(r, bgcorrect.method = "none", pmcorrect.method = "pmonly", normalize.method = "quantiles")
	data <- exprs(eset)
	if(!compute.calls) {
		return(data)
	} else {
		calls <- get.calls(r)
		res <- list(data=data, calls=calls) 
		return(res)
	}
}

########################################################
# MISC FUNCTIONS

install.required.packages <- function(libdir, method) {
	info(libdir)
	if(!is.package.installed(libdir, "reposTools")) {
		info("installing reposTools")
		install.package(libdir, "reposTools_1.8.0.zip", "reposTools_1.8.0.tar.gz", "reposTools_1.8.0.tar.gz")
	}
	
	if(!is.package.installed(libdir, "Biobase")) {
		info("installing Biobase")
		install.package(libdir, "Biobase_1.5.0.zip", "Biobase_1.5.0.tgz", "Biobase_1.5.0.tar.gz")
	}
	if(!is.package.installed(libdir, "affy")) {
		info("installing affy")
		install.package(libdir, "affy_1.5.8-1.zip", "affy_1.5.8.tgz","affy_1.5.8.tar.gz")
	}	
	
	if(method=='GCRMA' && !is.package.installed(libdir, "matchprobes")) {
		info("installing matchprobes")
		if(isMac()) {
			Sys.putenv(MAKEFLAGS="LIBR= SHLIB_LIBADD= LIBS=")
		}
		install.package(libdir, "matchprobes_1.0.22.zip", "matchprobes_1.0.22.tar.gz","matchprobes_1.0.22.tar.gz")
	}
	
#	if(!is.package.installed(libdir, "gcrma")) {
#		info("installing gcrma")
#		install.package(libdir, "gcrma_1.1.4.zip", "gcrma_1.1.4.tar.gz","gcrma_1.1.4.tar.gz")
#	}
	
	if(method=='GCRMA' && !is.package.installed(libdir, "gcrma")) {
		info("installing gcrma")
		if(isMac()) {
			Sys.putenv(MAKEFLAGS="LIBR= SHLIB_LIBADD= LIBS=")
		}
		install.package(libdir, "gcrma_2.2.1.zip", "gcrma_2.2.1.tar.gz","gcrma_2.2.1.tar.gz")
	}
	
	
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

get.cel.file.names <- function(input.file.name) {
	isWindows <- Sys.info()[["sysname"]]=="Windows"
	if(isWindows) {
		zip.unpack(input.file.name, dest=getwd())
	} else {
		 zip <- getOption("unzip")
		 system(paste(zip, "-q", input.file.name))
	}

	files <- list.files(recursive=TRUE)
   cel.files <- files[grep(".[cC][eE][lL].gz$|.[cC][eE][lL]$", files)]
	return(cel.files)
}

is.compressed <- function(cel.files) {
	for(f in cel.files) {
		r <- grep(".[cC][eE][lL].gz$", f)
		if(length(r) > 0) {
			return(TRUE)
		}
	}
	return(FALSE)
}


read.clm <- function(input.file.name) {
	s <- read.table(input.file.name, colClasses = "character", sep="\t", comment.char="") 
	columns <- ncol(s)
	scan.names <- s[,1]
	sample.names <- NULL
	if(columns > 1) {
		sample.names <- s[,2]
	}
	
	class.names <- NULL
	if(columns>2) {
		class.names <- s[, 3]
	}
	
	for(i in 1:length(scan.names)) { # check for duplicate scans
		scan.names[i] <- trim(scan.names[i])
		scan <- scan.names[i]
		l <- which(scan.names==scan)
		if(length(l) >= 2) {
			exit(paste("Duplicate scan:", scan))
		}
	}
	

	#for(name in names) {
	#	info(paste("checking if clm file contains:", name))
	#	if(!(name %in% scans)) {	
	#		exit(paste("Clm file missing scan", name))	
	#	}
	#}
	info("Finished reading clm file")
	f <- NULL
	if(columns > 2) {
		f <- factor(class.names)
	}
	list("factor"=f, "scan.names"=scan.names , "sample.names"=sample.names)
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

gp.normalize <- function(struc, method, reference.column=-1, use.p.p.genes=FALSE, sc=500) {
	if(!(method %in% c('none', 'target signal', 'quantile normalization', 'linear fit', 'mean scaling', 'median scaling'))) {
		stop("Unknown scaling method")
	}
	data <- struc$data
	calls <- struc$calls
	if(method=="none") {
		return(struc)
	}	
	
	if(method=='target signal') {
		nf <- 1
		for (i in 1:ncol(data)) {
        slg <- data[, i]
        sf <- sc/mean(slg, trim = 0.02)
        reported.value <- nf * sf * slg
        data[, i] <- reported.value
    	}
    	struc$data <- data
    	return(struc) 
	} 
    
	if(method=="quantile normalization") {
		data <- normalize.quantiles(data)
		struc$data <- data
		return(struc)
	}
	
	if(reference.column == -1) {
		index <- sort(apply(data, 2, median), index=TRUE)$ix
		reference.column <- index[(length(index)/2) + 1] 
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
		ref.mean.or.median <- function.to.apply(data[,reference.column])
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
		if(!is.null(struc$column.descriptions)) {
			if(struc$column.descriptions[j]=='') {
				prev <- ''
			} else {
				prev <- paste(struc$column.descriptions[j], ", ", sep='')
			}
			struc$column.descriptions[j] <- paste(prev, "scale factor=", scalingFactors[j], sep='')
		}
	}
	struc$data <- data
	return(struc)
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


      

