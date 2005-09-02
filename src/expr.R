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

cleanup <- function() {
	files <- dir()
	for(file in files) {
		if(file != zip.file.name && file!=output.data.file.name && file!=output.cls.file.name && file!=clm.input.file) {
			log(paste("removing", file))
         unlink(file, recursive=TRUE)
      }
	}		
}
# Scales all arrays so they have the same mean or median value
# new.value.i <- sample.i/mean.of.sample.i * ref.sample
# constant <- mean.of.sample.i * ref.sample

normalize <- function(data, method, reference.sample.name) {
	if(method=='') {
		return(data)
	}
	
	# get index of reference.sample.name in data
	colnames <- colnames(data)
	index <- 1
	refindex <- NULL
	for(name in colnames) {
		if(name==reference.sample.name) {
			refindex <- index
			break
		}
		index <- index + 1
	}
	
	log(paste("refindex", refindex))
	if(is.null(refindex)) {
		warning("Could not find reference sample name. Skipping scaling.")
		return(data)
	}
	if(method=='mean scaling' || method=='median scaling') {
		if(method=='mean scaling') {
			FUN <- 'mean'
		}
		else if(method=='median scaling') {
			FUN <- 'median'
		}
		function.to.apply <- get(FUN)
		ref.mean.or.median <- function.to.apply(data[,refindex])
		for(i in 1:ncol(data)) {
			if(i!=refindex) {
				scaling.factor <- function.to.apply(data[,i])/ref.mean.or.median
				data[, i] <- data[, i] / scaling.factor
			}
		}
		return(data)
	} else if(method=='rank normalization') {
		return(rank.normalize(data))
	} else {
		warning(paste("Unknown normalization method:", method, sep=''))	
	}
}

parseCmdLine <- function(...) {
	args <- list(...)
	input.file.name <- ''
	output.file.name <- ''
	method <- ''
	quantile.normalization <- ''
	background <- ''
	scale <- ''
	compute.calls <- ''
	normalization.method <- ''
	reference.sample.name <- ''
	libdir <- ''
	clm.input.file <- ''
	
	for(i in 1:length(args)) {
		flag <- substring(args[[i]], 0, 2)
		value <- substring(args[[i]], 3, nchar(args[[i]]))
		if(flag=='-i') {
			input.file.name <- value
			zip.file.name <<- input.file.name # for cleanup
		} else if(flag=='-o') {
			output.file.name <- value
		} else if(flag=='-m') {
			method <- value
		} else if(flag=='-q') {
			quantile.normalization <- value
		} else if(flag=='-b') { # whether to background correct when using RMA
			background <- value
		} else if(flag=='-s') {
			scale <- value
		} else if(flag=='-c') {
			compute.calls <- value
		} else if(flag=='-n') {
			normalization.method <- value
		} else if(flag=='-x') {
			reference.sample.name <- value
		} else if(flag=='-f') {
			clm.input.file <- value
		} else if(flag=='-l') {
			libdir <- value
		}  else  {
			exit(paste("unknown option", flag, sep=": "))
		} 
	}
	create.expression.file(input.file.name, output.file.name, method, 	quantile.normalization, background, scale, compute.calls, normalization.method, reference.sample.name, clm.input.file, libdir)

}

create.expression.file <- function(input.file.name, output.file.name, method, quantile.normalization, background, scale, compute.calls, normalization.method, reference.sample.name, clm.input.file, libdir)  {
	
	ret <- try(.create.expression.file(input.file.name, output.file.name, method, quantile.normalization, background, scale, compute.calls, normalization.method, reference.sample.name, clm.input.file, libdir)
		)
	if(class(ret)=="try-error") {
      exit("An error occurred while running the module.")
   }
}

.create.expression.file <- function(input.file.name, output.file.name, method, quantile.normalization, background, scale, compute.calls, normalization.method, reference.sample.name, clm.input.file, libdir)  {
	source(paste(libdir, "common.R", sep=''))
	options("warn"=-1)
	zip.file.name <<- input.file.name # for cleanup
	quantile.normalization <- string.to.boolean(quantile.normalization)
    background <- string.to.boolean(background)
    scale <- as.integer(scale)
	compute.calls <- string.to.boolean(compute.calls)
	clm.input.file <<- clm.input.file

	if(libdir!='') {
		setLibPath(libdir)
		on.exit(cleanup())
		install.required.packages(libdir)
	}

	library(affy, verbose=FALSE)
	
	result <- NULL
	isRes <- FALSE
	if(method=='dChip' || method=='RMA' || method=='GCRMA') {
		log("reading zip file")
		#afbatch <- gp.readAffyBatch(input.file.name)
      f <- get.cel.file.names(input.file.name)
		if(method=='dChip') {
			eset <- gp.dchip(afbatch)
		} else if(method=='RMA'){
			log(paste("running with", f$cel.files))
			eset <- gp.rma(f$cel.files, f$compressed, quantile.normalization, background)
		} else {
         eset <- gp.gcrma(afbatch)
      }
		data <- as.data.frame(exprs(eset))
		result <- data
		log(paste("Finished running", method))
	} else if(method=='MAS5'){
		isRes <- compute.calls
		result <- gp.mas5(input.file.name, compute.calls, scale)
		log("Finished running mas5")
	} else {
		exit('Unknown method')	
	}
	
	
	if(clm.input.file!='') { 
		if(isRes) {
			# reorder data
			clm <- apply.clm(result$data, clm.input.file)
			result$data <- clm$data
			
			# reorder calls TODO only need to read the clm file once
			clm <- apply.clm(result$calls, clm.input.file)
			result$calls <- clm$data
			
			factor <- clm$factor		
			if(!is.null(factor)) {
				cls <- list(labels=factor,names=levels(factor))
				class(cls) <- "cls"
				output.cls.file.name <<- get.cls.file.name(output.file.name)
				log(paste("saving cls file to", output.cls.file.name))
				write.cls(cls, output.cls.file.name)
			}
		} else {
			clm <- apply.clm(result, clm.input.file)
			result <- clm$data
			factor <- clm$factor	
			if(!is.null(factor)) {
				cls <- list(labels=factor,names=levels(factor))
				output.cls.file.name <<- get.cls.file.name(output.file.name)
				log(paste("saving cls file to", output.cls.file.name))
				write.cls(cls, output.cls.file.name)
			}
		}
	}
	if(isRes) {
		log("normalizing res file...")
		result$data <- normalize(result$data, normalization.method, reference.sample.name)
		log("finished normalizing res file")
	} else {
		result <- normalize(result, normalization.method, reference.sample.name)
	}
	if(isRes) {
		log("writing res file...")
		output.data.file.name <<- write.res(result, output.file.name)
		log("finished writing res file")
	} else {
		gct <- list(data=result, row.descriptions='')
		output.data.file.name <<- write.gct(gct, output.file.name)
		log(paste("wrote gct file to",output.data.file.name)) 
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
	
gp.rma <- function(cel.files, compressed, normalize, background) {
	samplenames <- gsub("^/?([^/]*/)*", "", unlist(cel.files), 
            extended = TRUE)
   n <- length(cel.files)
	pdata <- data.frame(sample = 1:n, row.names = samplenames)
	phenoData <- new("phenoData", pData = pdata, varLabels = list(sample = "arbitrary numbering"))
        
   eset <- just.rma(filenames=cel.files, compress=compressed, normalize=normalize, background=background, verbose=FALSE, phenoData=phenoData)
  # r <- ReadAffy(filenames=cel.files, compress=compressed) 
   #eset <- rma(r, normalize=normalize, background=background, verbose=FALSE, compress=compressed)
   
   eset@exprs <- 2^eset@exprs # rma produces values that are log scaled
   return(eset)
}


gp.dchip <- function(afbatch) {
	eset <- expresso(afbatch, normalize.method="invariantset", bg.correct=FALSE, pmcorrect.method="pmonly",summary.method="liwong", verbose=FALSE) 
	return(eset)
}

#normalize: logical. If 'TRUE' scale normalization is used after we obtain an instance of 'exprSet-class'

#sc: Value at which all arrays will be scaled to.

#analysis: should we do absolute or comparison analysis, although "comparison" is still not implemented.

gp.mas5 <- function(input.file.name, compute.calls, scale) {
	r <- gp.readAffyBatch(input.file.name)
	if(!is.na(scale)) {
		log("mas5: normalizing and scaling data")
		eset <- mas5(r, normalize=TRUE, sc=scale)
	} else {
		log("mas5: normalizing and scaling data")
		eset <- mas5(r, normalize=FALSE)
	}
	
	if(!compute.calls) {
		data <- as.data.frame(exprs(eset))
		names(data) <- colnames((exprs(eset)))
		return(data)
	} else {
		calls.eset <- mas5calls.AffyBatch(r, verbose = FALSE) 
		data <- as.data.frame(exprs(eset))
		calls <- as.data.frame(exprs(calls.eset))
		names(data) <- colnames((exprs(eset)))
		names(calls) <- colnames((exprs(eset)))
		res <- list(row.descriptions='', column.descriptions='', data=data, calls=calls) 
		return(res)
	}
}

gp.gcrma <- function(afbatch) {
   if(!require("gcrma", quietly=TRUE)) {
      install.package(libdir, "gcrma_1.0.0.zip", "gcrma_1.0.0.zip", "gcrma_1.0.0.zip")
   }
   
   if(!require("matchprobes", quietly=TRUE)) {
       install.package(libdir, "matchprobes_1.0.0.zip", "matchprobes_1.0.0.zip", "matchprobes_1.0.0.zip")
   }
   
   cdf <- cleancdfname(afbatch@cdfName) # e.g "hgu133acdf"
   cdf <- substring(cdf, 0, nchar(cdf)-3)# remove cdf from end
   pkg <- paste(cdf, "probe", sep='') 
   if(!library(package=pkg, lib.loc=libdir, logical.return=TRUE, version=NULL)) {
      repList <- getOptReposList()
      n <- repNames(repList)
      
      get.index <- function(n) {
         for(index in 1:length(n)) {
            if(n[index] =='Bioconductor Probe Data Packages') {
               return(index)
            }
         }
         exit("Probe Data repository not found.")
       }
       isWindows <- Sys.info()[["sysname"]]=="Windows"
       
       if(isWindows) {
         type <- "Win32"  
       } else {
          type <- "Source"
       }
       r <- getRepEntry(repList, get.index(n))
       install.packages2(repEntry=r, pkgs=c(pkg), type='Win32')
   }
 #  affinity.info <- compute.affinities(afbatch@cdfName, verbose=FALSE)
   eset <- gcrma(afbatch, verbose=TRUE)
   eset@exprs <- 2^eset@exprs # produces values that are log scaled
   return(eset)
   
}

########################################################
# MISC FUNCTIONS

install.required.packages <- function(libdir) {
	log(libdir)
	if(!is.package.installed(libdir, "reposTools")) {
		log("installing reposTools")
		install.package(libdir, "reposTools_1.5.2.zip", "reposTools_1.5.19.tgz", "reposTools_1.5.19.tar.gz")
	}
	
	if(!is.package.installed(libdir, "Biobase")) {
		log("installing Biobase")
		install.package(libdir, "Biobase_1.5.0.zip", "Biobase_1.5.0.tgz", "Biobase_1.5.0.tar.gz")
	}
	if(!is.package.installed(libdir, "affy")) {
		log("installing affy")
		install.package(libdir, "affy_1.5.8-1.zip", "affy_1.5.8.tgz","affy_1.5.8.tar.gz")
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

gp.readAffyBatch <- function(input.file.name) {
	f <- get.cel.file.names(input.file.name)
	r <- ReadAffy(filenames=f$cel.files, compress=f$compressed) 
	return(r)
}

get.cel.file.names <- function(input.file.name) {
	isWindows <- Sys.info()[["sysname"]]=="Windows"
	if(isWindows) {
		zip.unpack(input.file.name, dest=getwd())
	} else {
		 zip <- getOption("unzip")
		 system(paste(zip, "-q", input.file.name))
	}

	d <- dir(recursive=TRUE)
	cel.files <- vector()
	sampleNames <- vector()
	
	index <- 1
	compressed <- FALSE
	cel.files.only <- FALSE

	for(i in 1:length(d)) {
      if(d[i]==input.file.name) {
         next  
      }
		gzipped <- regexpr(paste(".cel.gz","$",sep=""), tolower(d[i]))
		zipped <- regexpr(paste(".cel.zip","$",sep=""), tolower(d[i]))
		if(gzipped[[1]] != -1 || zipped[[1]] != -1) {
			cel.files[index] <- d[i]
			sampleNames[index] <- basename(d[i])
			if(gzipped[[1]] != -1) {
				dot.index <- regexpr(paste(".cel.gz","$",sep=""), tolower(sampleNames[index]))[[1]]
			} else {
				dot.index <- regexpr(paste(".cel.zip","$",sep=""), tolower(sampleNames[index]))[[1]]
			}
			# sampleNames[index] <- substring(sampleNames[index], 0, dot.index-1)
			index <- index + 1
			compressed <- TRUE
		} else {
			result <- regexpr(paste(".cel","$",sep=""), tolower(d[i]))
			if(result[[1]] != -1) {
				cel.files[index] <- d[i]
				sampleNames[index] <- basename(d[i])
				dot.index <- regexpr(paste(".cel","$",sep=""), tolower(sampleNames[index]))[[1]]
				# sampleNames[index] <- substring(sampleNames[index], 0, dot.index-1)
				index <- index + 1
				cel.files.only <- TRUE
			}	
		}
	}

	if(compressed && cel.files.only) {
		exit("Both compressed and uncompressed .cel files found.")	
	}
	class(sampleNames) <- 'character'
	list("cel.files"=cel.files, "compressed"=compressed)
}


# reorders the given data, substitutes sample for scan names, and creates a factor based on class names
apply.clm <- function(data, clm.file.name) {
	clm <- read.clm(clm.file.name, names(data))
   
	reordered.scan.names <- clm$scan.names
   	reordered.sample.names <- clm$sample.names
   	
   	i <- 1
   	order <- list()
   	for(reordered.scan in reordered.scan.names) {
   		order[reordered.scan] <- i
      	i <- i+1
   	}
	log("reordering...")
	new.data <- reorder.data.frame(order, data)
	
	# substitute sample names for scan names
	names(new.data) <- reordered.sample.names
	
	# note we don't need to reorder the factor because reordered data is now in same order as factor
	list("factor"=clm$factor, "data"=new.data)
}

# order - a list containing new order of data e.g. 3, 4, 1, 2
reorder.data.frame <- function(order, data) {
	new.data <- data.frame(data)
	
	for(j in 1:NCOL(data)) {
		
		name <- names(data)[j]
	
		i <- order[name]
		i <- i[[1]]
		new.data[, i] <- data[,j]
		names(new.data)[i] <- name
	}
	new.data
}

# names - the column names in the expression dataset
read.clm <- function(input.file.name, names) {
	s <- read.table(input.file.name, colClasses = "character", sep="\t") 
	columns <- ncol(s)
	class.names <- vector() 
	sample.names <- vector()
	scans <- vector()	
	found.scans <- list() # keeps track of scan names in clm file
	for(i in 1:NROW(s)) {
		scans[i] <- s[i, 1]
		if(is.null(found.scans[[scans[i]]])==FALSE) {
			exit("Class file contains scan ", scans[i], " more than once")
		}
		found.scans[scans[i]] <- 1
		sample.names[i] <- s[i, 2]
		if(columns>2) {
			class.names[i] <- s[i, 3]
		}
	}
	log(paste("names", names ))
	for(name in names) {
		log(paste("checking if clm file contains:", name))
		if(is.null(found.scans[[scans[i]]])) {	
			exit("Class file missing scan ", name)	
		}
	}
	log("Finished reading clm file")
	f <- NULL
	if(columns > 2) {
		f <- factor(class.names)
	}
	list("factor"=f, "scan.names"=scans, "sample.names"=sample.names)
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


