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
exec.log <- 'gp_task_execution_log.txt'

cleanup <- function() {
	files <- dir()
	for(file in files) {
		if(file != zip.file.name && file!=output.data.file.name && file!=output.cls.file.name && file!=clm.input.file && file!=exec.log && file!=probe.descriptions.file.name) {
			log(paste("removing", file))
         unlink(paste("'", file, "'", sep=''), recursive=TRUE)
      }
	}		
}

get.median.index <- function(data) {
	medians <- vector(mode="numeric") # find median scan in data
	for(col in 1:ncol(data)) {
		medians[col] <- median(data[,col])
	}

	index <- sort(medians, index=TRUE)$ix
	refindex <- index[length(index)/2]
	return(refindex)
}

# returns a matrix with the same number of rows as calls and one column containing the sum of the number of P calls for that row
get.p.indices <- function(calls) {
	 temp <- apply(calls, 1, function(val) {
     	 return(sum(val=="P"))
    })
    expected.number <- ncol(calls)
    temp <- as.matrix(temp)
    ok <- (temp >= expected.number)
    return(ok)
}

# returns a new data matrix containing only those genes which have all P calls
get.p.only.data <- function(data, calls) {
	return(data[get.p.indices(calls), ])
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


# Scales all arrays so they have the same mean or median value
# new.value.i <- sample.i/mean.of.sample.i * ref.sample
# constant <- mean.of.sample.i * ref.sample
my.normalize <- function(data, method, reference.sample.name='', sc=NULL, calls=NULL) {
	log(paste("my normalize function...", method))
	if(method=='none' || method=='') {
		log("not normalizing, returning data")
		return(data)
	}

	
	if(method=='target signal') {
		nf <- 1
		for (i in 1:ncol(data)) {
        slg <- data[, i]
        sf <- sc/mean(slg, trim = 0.02)
        reported.value <- nf * sf * slg
        data[, i] <- reported.value
    	}
    	return(data) 
    } 
    
	refindex <- NULL
	

	if(reference.sample.name=='') {
		refindex <- get.median.index(data)
	} else {
		# get index of reference.sample.name in data
		colnames <- colnames(data)
		refindex <- which(colnames==name)
	}
	
	log(paste("refindex", refindex))
	if(!isTRUE(refindex > 0)) {
		refindex <- get.median.index(data)
		warning("Could not find reference scan name. Using median scan.")
	}
	
	if(method=='mean scaling' || method=='median scaling') {
    	scale.factor <- vector()
    	scale.factor[refindex] <- "scale factor=1"
		if(method=='mean scaling') {
			FUN <- 'mean'
		}
		else if(method=='median scaling') {
			FUN <- 'median'
		}
		function.to.apply <- get(FUN)
		ref.mean.or.median <- function.to.apply(data[,refindex])
		if(!is.null(calls)) {
			for(i in 1:ncol(data)) {
				if(i!=refindex) {
					p.data <- get.p.only.data(cbind(data[,refindex], data[, i]), cbind(calls[,refindex], calls[, i]))
					ref.mean.or.median <- function.to.apply(p.data[,refindex])
					scaling.factor <- function.to.apply(p.data[,i])/ref.mean.or.median 
					data[, i] <- data[, i] / scaling.factor
					scale.factor[i] <- paste("scale factor=", scaling.factor, sep='')
				} 
			}
			
		} else {
			for(i in 1:ncol(data)) {
				if(i!=refindex) {
					scaling.factor <- function.to.apply(data[,i])/ref.mean.or.median
					data[, i] <- data[, i] / scaling.factor
					scale.factor[i] <- paste("scale factor=", scaling.factor, sep='')
				} 
			}
			attr(data, "scale.factor") <-  scale.factor
		}
		return(data)
	} else if(method=='rank normalization') {
		return(rank.normalize(data))
	} else if(method=='linear fit') {
		scale.factor <- vector()
		scale.factor[i] <- "scale factor=1"	
				
 		if(!is.null(calls)) {
			for(i in 1:ncol(data)) {
				if(i!=refindex) {
					p.data <- get.p.only.data(cbind(data[,refindex], data[, i]), cbind(calls[,refindex], calls[, i]))
					scaling.factor <- linear.fit(p.data[,refindex], p.data[,i])$m
					data[, i] <- data[, i] * scaling.factor	
					scale.factor[i] <- paste("scale factor=", scaling.factor, sep='')
				} 
			}
			
		} else {
			for(i in 1:ncol(data)) {
				if(i!=refindex) {
					scaling.factor <- linear.fit(data[,refindex], data[,i])$m
					data[, i] <- data[, i] * scaling.factor
					scale.factor[i] <- paste("scale factor=", scaling.factor, sep='')
				} 
			}
		}
		attr(data, "scale.factor") <-  scale.factor
		return(data)
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
	use.p.p.genes <- FALSE
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
		} else if(flag=='-s') {
			scale <- value
		} else if(flag=='-c') {
			compute.calls <- value
		} else if(flag=='-n') {
			normalization.method <- value # none, mean, median
		} else if(flag=='-x') {
			reference.sample.name <- value
		} else if(flag=='-f') {
			clm.input.file <- value
		} else if(flag=='-l') {
			libdir <- value
		} else if(flag=='-p') {
			use.p.p.genes <- value
		} else if(flag=='-r') {
			row.descriptions.file <- value
			probe.descriptions.file.name <<- value # for cleanup
		} else if(flag=='-a') {
			if (value !='' && .Platform$OS.type == "windows") {
      		memory.limit(size=as.numeric(value))
   		}
		}  else  {
			stop(paste("unknown option", flag, sep=": "), call.=FALSE)
		} 
		
	}
	
	
	create.expression.file(input.file.name, output.file.name, method, 	quantile.normalization, background, scale, compute.calls, normalization.method, reference.sample.name, clm.input.file, libdir, use.p.p.genes, row.descriptions.file)

}

create.expression.file <- function(input.file.name, output.file.name, method, quantile.normalization, background, scale, compute.calls, normalization.method, reference.sample.name, clm.input.file, libdir, use.p.p.genes, row.descriptions.file)  {
	source(paste(libdir, "common.R", sep=''))
	DEBUG <<- FALSE
	log(paste("normalization.method", normalization.method))
	options("warn"=-1)
	zip.file.name <<- input.file.name # for cleanup
	quantile.normalization <- string.to.boolean(quantile.normalization)
	background <- string.to.boolean(background)
	use.p.p.genes <- string.to.boolean(use.p.p.genes)
   scale <- as.integer(scale)
   if(normalization.method=='target signal' && is.na(scale)) {
		stop("Target signal value required when scaling method is target signal", call.=FALSE)
	}
	compute.calls <- string.to.boolean(compute.calls)
	if(method=='MAS5' && use.p.p.genes && !compute.calls) {
		stop("Must create res file if using only P-P genes for scaling.", call.=FALSE)
	}
	if(use.p.p.genes && !compute.calls) {
		use.p.p.genes <- FALSE
	}
	
	clm.input.file <<- clm.input.file

	if(libdir!='') {
		setLibPath(libdir)
		on.exit(cleanup())
		install.required.packages(libdir)
	}

	library(affy, verbose=FALSE)
	result <- NULL # data.frame if calls are not computed, else list containing calls and data
	isRes <- compute.calls
	
	zipFileGiven <- TRUE
	clm <- NULL
	if(input.file.name!='') {
		cel.file.names <- get.cel.file.names(input.file.name)
	} else if(clm.input.file!='') {
		if(output.file.name=='') {
			output.file.name <- sub(".clm$", '', clm.input.file, ignore.case=TRUE)
		}
			
		zipFileGiven <- FALSE
		clm <- read.clm(clm.input.file)
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
			log(cmd)
			system(cmd)
		}
	} else {
		exit("Either a zip of CEL files or a clm file is required.")
	}
	
	
	
	if(zipFileGiven && clm.input.file!='') { # reorder scan names
		clm <- read.clm(clm.input.file)
		scan.names <- clm$scan.names
		log(paste("scan.names", scan.names))
		new.cel.file.names <- vector("character")
		index <- 1
		for(scan in scan.names) {
			s1 <- paste(scan, '.cel', "$",sep='')
			s2 <- paste(scan, '.cel.gz', "$",sep='')
			s <- paste(s1, "|", s2, sep="")
			index <- grep(s, cel.file.names, ignore.case=TRUE)
			if(length(index) == 0) {
				cat(paste("Scan ", scan, "in clm file not found.\n"))
			} else if(length(index)>1) {
				cat(paste("Scan ", scan, "in clm file matches more than one CEL file.\n"))
			} else {
				new.cel.file.names[index] <- cel.file.names[index[1]]
				index <- index + 1
			}
		}
		cel.file.names <- new.cel.file.names
	}
	
	is.compressed <- is.compressed(cel.file.names)
	
	
	if(method=='dChip' || method=='RMA' || method=='GCRMA') {
		log("reading zip file")
		if(method=='dChip') {
			log("running dChip")
			result <- gp.dchip(cel.file.names, is.compressed, compute.calls)
		} else if(method=='RMA'){
			result <- gp.rma(cel.file.names, is.compressed, quantile.normalization, background, compute.calls)
		} 
		log(paste("Finished running", method))
	} else if(method=='MAS5'){
		result <- gp.mas5(cel.file.names, is.compressed, compute.calls)
		log("Finished running mas5")
	} else {
		exit('Unknown method')	
	}
	
	if(!is.null(clm)) { 
		if(!is.null(clm$sample.names)) {
			log("setting names")
			log(clm$sample.names)
			if(isRes) {
				colnames(result$data) <- clm$sample.names
			} else {
				colnames(result) <- clm$sample.names	
			}
		}
		
		factor <- clm$factor
		if(!is.null(factor)) {
			cls <- list(labels=factor,names=levels(factor))
			log(paste("cls: ", cls))
			output.cls.file.name <<- get.cls.file.name(output.file.name)
			log(paste("saving cls file to", output.cls.file.name))
			write.cls(cls, output.cls.file.name)
			log("cls file saved")
		}
			
	} else {  # remove .cel extension from names
		if(isRes) {
			col.names <- colnames(result$data)
		} else {
			col.names <- colnames(result)
		}
		log("removing .cel extension")
		col.names <- sub(".[cC][eE][lL].gz$|.[cC][eE][lL]$", "", col.names)
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
		log("writing res file...")
		#result$column.descriptions <- attr(result$data, "scale.factor")
		result$column.descriptions <- NULL
		result$row.descriptions <- row.descriptions
		output.data.file.name <<- write.res(result, output.file.name)
		log("finished writing res file")
	} else {
		log("writing gct file...")
		log(paste("dim of result:", dim(result)))
		gct <- list(data=result, row.descriptions=row.descriptions)
		output.data.file.name <<- write.gct(gct, output.file.name)
		log(paste("wrote gct file to", output.data.file.name)) 
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
	samplenames <- gsub("^/?([^/]*/)*", "", unlist(cel.files), 
            extended = TRUE)
   n <- length(cel.files)
	pdata <- data.frame(sample = 1:n, row.names = samplenames)
	phenoData <- new("phenoData", pData = pdata, varLabels = list(sample = "arbitrary numbering"))
   
   log(paste("normalize", normalize))
   log(paste("background", background))
   eset <- just.rma(filenames=cel.files, compress=compressed, normalize=normalize, background=background, verbose=FALSE, phenoData=phenoData)   
   eset@exprs <- 2^eset@exprs # rma produces values that are log scaled
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
	log("running dchip")
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
	log("running mas5...")
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
		install.package(libdir, "reposTools_1.5.19.zip", "reposTools_1.5.19.tgz", "reposTools_1.5.19.tar.gz")
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
	#	log(paste("checking if clm file contains:", name))
	#	if(!(name %in% scans)) {	
	#		exit(paste("Clm file missing scan", name))	
	#	}
	#}
	log("Finished reading clm file")
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



