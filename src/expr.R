# this file contains a collection of functions to go from probe level data (Cel files) to expression measures 

#input.file.name a zip file containing cel or gzipped or zipped cel files
#output.file.name name of the gct output file


string.to.boolean <- function(s) {
	if(s=="yes") {
		return(TRUE)	
	}
	return(FALSE)
}

# Scales all arrays so they have the same mean or median value
# new.value.i <- sample.i/mean.of.sample.i * ref.sample
# constant <- mean.of.sample.i * ref.sample

normalize <- function(data, method, refindex) {
	if(method=='') {
		return(data)
	} else if(method=='mean scaling' || method=='median scaling') {
		if(refindex=='') {
			warning("No sample index specified.")	
			return(data)
		} 
		if(refindex < 1 || refindex > NCOL(refindex)) {
			warning("Sample index out of range.")	
			return(data)
		}
		b <- new ('AffyBatch', exprs=as.matrix(data))
		if(method=='mean scaling') {
			FUN <- 'mean'
		}
		else if(method=='median scaling') {
			FUN <- 'median'
		}
		norm <- normalize.AffyBatch.constant(b, refindex=refindex, get(FUN))
	} else if(method=='rank normalization') {
		return(rank.normalize(data))
	} else {
		warning(paste("Unknown normalization method:", method, sep=''))	
	}
}

create.expression.file <- function(...) {
	args <- list(...)
	classes.file <- ''
	
	#optional args
	normalize <- FALSE
	background <- FALSE
	quantile.normalization <- FALSE
	compute.calls <- FALSE
	scale <- NULL
	normalization.method <- NULL
	classes.file <- NULL
	refindex <- NULL
	options(warn=-1)
	for(i in 1:length(args)) {
		flag <- substring(args[[i]], 0, 2)
		if(flag=='-i') {
			input.file.name <- substring(args[[i]], 3, nchar(args[[i]]))
		} else if(flag=='-o') {
			output.file.name <- substring(args[[i]], 3, nchar(args[[i]]))
		} else if(flag=='-m') {
			method <- substring(args[[i]], 3, nchar(args[[i]]))
		} else if(flag=='-q') {
			quantile.normalization <- substring(args[[i]], 3, nchar(args[[i]]))
			quantile.normalization <- string.to.boolean(quantile.normalization)
		} else if(flag=='-b') { # whether to background correct when using RMA
			background <- substring(args[[i]], 3, nchar(args[[i]]))
			background <- string.to.boolean(background)
		} else if(flag=='-s') {
			scale <- substring(args[[i]], 3, nchar(args[[i]]))
			scale <- as.integer(scale)
		} else if(flag=='-c') {
			compute.calls <- substring(args[[i]], 3, nchar(args[[i]]))
			compute.calls <- string.to.boolean(compute.calls)
		} else if(flag=='-n') {
			normalization.method <- substring(args[[i]], 3, nchar(args[[i]]))
		} else if(flag=='-x') {
			refindex <- substring(args[[i]], 3, nchar(args[[i]]))
			refindex <- as.integer(refindex)
		} else if(flag=='-f') {
			classes.file <- substring(args[[i]], 3, nchar(args[[i]]))
		} else if(flag=='-l') {
			libdir <- substring(args[[i]], 3, nchar(args[[i]]))
		}  else  {
			stop(paste("unknown option", flag, sep=": "))
		} 
	}
	
	on.exit(unlink(substring(input.file.name, 0, nchar(input.file.name)-4), recursive=TRUE))
	if(exists("libdir")) {
		install.affy.packages(libdir)
	}
	library(affy, verbose=FALSE)
	library(GenePattern, verbose=FALSE)
	
	if(method=='dChip' || method=='RMA') {
		my.list <- gp.readAffyBatch(input.file.name)
		afbatch <- my.list[[1]]
		sampleNames <- my.list[[2]]
		if(method=='dChip') {
			eset <- dchip(afbatch)
		} else {
			eset <- gp.rma(afbatch, quantile.normalization, background)
		}
		data <- as.data.frame(exprs(eset))
		names(data) <- sampleNames

		if(classes.file!='') { 
			reorder <- reorder(data, classes.file)
			data <- reorder$data
			factor <- reorder$factor		
			cls <- list(labels=factor,names=levels(factor))
			class(cls) <- "cls"
			cls.file.name <- get.cls.file.name(output.file.name)
			save.cls(cls, cls.file.name)
		}
		
		data <- normalize(data, normalization.method, refindex)
	
		gct.file.name <- save.data.as.gct(data, output.file.name)
		if(classes.file!='') { 
			return(list(gct.file.name, cls.file.name))
		}
		return(list(gct.file.name))
	} else if(method=='MAS5'){
		return(gp.mas5(input.file.name, output.file.name, compute.calls, scale, classes.file, post.normalization=normalization.method, refindex=refindex))
	} else {
		stop('Unknown method')	
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
	
gp.rma <- function(afbatch, normalize, background) {
	eset <- rma(afbatch, normalize=normalize, background=background, bgversion="2", verbose=FALSE)
	return(eset)
}


dchip <- function(afbatch) {
	eset <- expresso(afbatch, normalize.method="invariantset", bg.correct=FALSE, pmcorrect.method="pmonly",summary.method="liwong", verbose=FALSE) 
	return(eset)
}


#normalize: logical. If 'TRUE' scale normalization is used after we obtain an instance of 'exprSet-class'

#sc: Value at which all arrays will be scaled to.

#analysis: should we do absolute or comparison analysis, although "comparison" is still not implemented.

gp.mas5 <- function(input.file.name, output.file.name, compute.calls, scale, classes.file, post.normalization, refindex=refindex) {
	
	ab <- gp.readAffyBatch(input.file.name)
	r <- ab[[1]]
	sampleNames <- ab[[2]]
	if(scale!='') {
		eset <- mas5(r, normalize=TRUE, sc=scale)
	} else {
		eset <- mas5(r, normalize=FALSE)
	}
	
	if(!compute.calls) {
		data <- as.data.frame(exprs(eset))
		names(data) <- sampleNames

		if(classes.file!='') {
			reorder <- reorder(data, classes.file)
			data <- reorder$data
			factor <- reorder$factor		
			cls <- list(labels=factor,names=levels(factor))
			class(cls) <- "cls"
			cls.file.name <- get.cls.file.name(output.file.name)
			save.cls(cls, cls.file.name)
		}
		normalize(data, post.normalization, refindex)
		gct.file.name <- save.data.as.gct(data, output.file.name)
		if(classes.file!='') {
			return(list(gct.file.name, cls.file.name))
		}
		return(list(gct.file.name))
	} else {
		calls.eset <- mas5calls.AffyBatch(r, verbose = FALSE) 
		# write.res in GenePattern library is broken
		res.ext <- regexpr(paste(".res","$",sep=""), tolower(output.file.name))
		if(res.ext[[1]] == -1) {
			output.file.name <- paste(output.file.name, ".res", sep="") # ensure correct file extension
		}
	
		data <- as.data.frame(exprs(eset))
		calls <- as.data.frame(exprs(calls.eset))
		names(data) <- sampleNames
		names(calls) <- sampleNames
		
		if(classes.file!='') {
			reorder <- reorder(data, classes.file)
			data <- reorder$data
			factor <- reorder$factor
			
			calls <- reorder(calls, classes.file)$data
			
			cls <- list(labels=factor,names=levels(factor))
			class(cls) <- "cls"
			cls.file.name <- get.cls.file.name(output.file.name)
			save.cls(cls, cls.file.name)
		}
		normalize(data, post.normalization, refindex)
		res <- new ("res", gene.descriptions='', sample.descriptions='', data=data, calls=calls) 
		my.write.res(res, output.file.name)
		if(classes.file!='') {
			return(list(output.file.name, cls.file.name))
		}
		return(list(output.file.name))
	}
}


########################################################
# MISC FUNCTIONS


install.affy.packages <- function(libdir) {
	if(!require("reposTools", quietly=TRUE)) {	
	  isWindows <- Sys.info()[["sysname"]]=="Windows"
	  if(isWindows) {
		  f <- paste(libdir, "reposTools_1.3.6.zip", sep="")
		  .install.windows(f)
	  } else { # install from source
		  f <- paste(libdir, "reposTools_1.3.6.tar.gz", sep="")
		  install(f)
	  }
	}
	library(reposTools)
	if(!require("Biobase", quietly=TRUE)) {	
	  isWindows <- Sys.info()[["sysname"]]=="Windows"
	  if(isWindows) {
		  f <- paste(libdir, "Biobase_1.4.0.zip", sep="")
		  .install.windows(f)
	  } else { # install from source
		  f <- paste(libdir, "Biobase_1.4.0.tar.gz", sep="")
		  install(f)
	  }
	}
	library(Biobase)
	if(!require("affy", quietly=TRUE)) {	
	  isWindows <- Sys.info()[["sysname"]]=="Windows"
	  if(isWindows) {
		  f <- paste(libdir, "affy_1.3.28.zip", sep="")
		  .install.windows(f)
	  } else { # install from source
		  f <- paste(libdir, "affy_1.3.28.tar.gz", sep="")
		  install(f)
	  }
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

save.cls <- function(cls, output.file.name) {
	cls.ext <- regexpr(paste(".cls","$",sep=""), tolower(output.file.name))
	if(cls.ext[[1]] == -1) {
		output.file.name <- paste(output.file.name, ".cls", sep="") # ensure correct file extension
	}
	write.cls(file=output.file.name, cls)
}

install <- function(pkg) {
    lib <- .libPaths()[1]
    cmd <- paste(file.path(R.home(), "bin", "R"), "CMD INSTALL")
    cmd <- paste(cmd, "-l", lib)
    cmd <- paste(cmd, " '", pkg, "'", sep = "")
    status <- system(cmd)
    if (status != 0) 
    	cat("\tnpackage installation failed\n")
}

.install.windows <- function(file) {
	install.packages(file, .libPaths()[1], CRAN=NULL)
}

gp.readAffyBatch <- function(input.file.name) {
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
			sampleNames[index] <- substring(sampleNames[index], 0, dot.index-1)
			index <- index + 1
			compressed <- TRUE
		} else {
			result <- regexpr(paste(".cel","$",sep=""), tolower(d[i]))
			if(result[[1]] != -1) {
				cel.files[index] <- d[i]
				sampleNames[index] <- basename(d[i])
				dot.index <- regexpr(paste(".cel","$",sep=""), tolower(sampleNames[index]))[[1]]
				sampleNames[index] <- substring(sampleNames[index], 0, dot.index-1)
				index <- index + 1
				cel.files.only <- TRUE
			}	
		}
	}

	if(compressed && cel.files.only) {
		stop("Both compressed and uncompressed .cel files found.")	
	}
	class(sampleNames) <- 'character'
	r <- ReadAffy(filenames=cel.files, sampleNames=sampleNames, compress=compressed) 
	return(list(r, sampleNames))
}



save.data.as.gct <- function(data, output.file.name) {
	gct.ext <- regexpr(paste(".gct","$",sep=""), tolower(output.file.name))
	if(gct.ext[[1]] == -1) {
		output.file.name <- paste(output.file.name, ".gct", sep="") # ensure correct file extension
	}
	
	gct <- new ("gct", row.descriptions='', data=data) 
	
	write.gct(gct, output.file.name)
	
	return(output.file.name)	
}

my.write.res <-
#
# write a res structure as a file
#
function(res, filename)
{
	if(!inherits(res,"res")) {
		stop("argument `res' must be a res structure.")
	}
	f <- file(filename, "w")
	
	# write the labels
	cat("Description\tAccession\t", file=f, append=TRUE)
	cat(names(res@data), sep="\t\t", file=f, append=TRUE)
	cat("\n", file=f, append=TRUE)
	
	# write the descriptions
	cat("\t", file=f, append=TRUE)
	cat(res@sample.descriptions, sep="\t\t", file=f, append=TRUE)
	cat("\n", file=f, append=TRUE)
	
	# write the size
	cat(NROW(res@data), "\n", sep="", file=f, append=TRUE)
	
	# write the data
	# 1st combine matrices
	dim <- dim(res@data)
	dim[2] <- dim[2]*2
	
	m <- matrix(nrow=dim[1], ncol=dim[2]+2)
	m[,1] <- res@gene.descriptions
	m[,2] <- row.names(res@data)
	
	index <- 3
	for(i in 1:dim(res@data)[2]) {
		m[,index] <- res@data[,i]
		index <- index + 2
	}
	index <- 4
	
	for(i in 1:dim(res@calls)[2]) {
		m[,index] <- as.character(res@calls[,i])
		index <- index + 2
	}
	write.table(m, file=f, col.names=FALSE, row.names=FALSE, append=TRUE, quote=FALSE, sep="\t", eol="\n")
	close(f)
	return(filename)
}

reorder <- function(data, cls.file.name) {
	l <- read.class(cls.file.name, names(data))
	order <- order.by.class(l)
	new.data <- reorder.by.class(order, data)
	list("factor"=l$factor, "data"=new.data)
}

reorder.by.class <- function(order, data) {
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

order.by.class <- function(list) {
	scans <- list$scans
	factor <- list$factor
	index <- 1
	order <- list(length=length(factor))
	for(l in levels(factor)) {
		# find indices in factor with value l
		for(i in 1:length(factor)) {
			if(factor[i]==l) {
				order[scans[i]] <- index
				index <- index + 1
			}
		}	
	}
	order
}

read.class <- function(input.file.name, names) {
	s <- read.table(input.file.name, colClasses=c('character', 'character'), sep="\t")
	x <- vector()
	scans <- vector()	
	found.scans <- list()
	for(i in 1:NROW(s)) {
		scans[i] <- s[i, 1]
		if(is.null(found.scans[[scans[i]]])==FALSE) {
			stop("Class file contains scan ", scans[i], " more than once")
		}
		found.scans[scans[i]] <- 1
		x[i] <- s[i, 2]
	}
	for(name in names) {
		if(is.null(found.scans[[scans[i]]])) {	
			stop("Class file missing scan ", name)	
		}
	}
	
	f <- factor(x)
	list("factor"=f, "scans"=scans)
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


