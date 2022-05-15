library("KEGGREST")

# FUNCTIONS FOR PARSING INPUT
promptInput <- function(){
	cat("Options:\n
		Make network from input with KO ids: c [filename] [outputdir]\n
		Make network from input with only gene names: cg [filename] [outputdir]\n")
	cat("Please enter command: ")
	userinp <- readLines("stdin",n=1)
	parseInput(strsplit(userinp," "))
}
parseInput <- function(s){
	# check for each option
	if(s[1] == "cg"){ # get KO ids and make network
		dir <- makeOutDir(s[3])
		getKO(s[2],dir)
		makeNetwork(paste(dir,"KO_genelist.txt",sep=""),dir)
	}else if(s[1] == "c"){ # already have KO ids
		dir <- makeOutDir(s[2])
		makeNetwork(s[1],dir)
	}else{ # exit error
		stop(paste("Invalid argument:",s[1]))
	}
	summary(s)
}
makeOutDir <- function(outdir){ # makes output and formats dir name
	if(is.na(outdir)) # exit error
		stop("Output directory missing")
	dir.create(outdir)
	if(substring(outdir,nchar(outdir))=="/")
		return(outdir)
	else
		return(paste(outdir,"/",sep="")) # format outdir with backslash
}

# FUNCTIONS FOR GETTING KO IDS
getKO <- function(fname,outdir){
	fin <- file(fname,open="r") # read file with one gene name per line
	cat(paste(outdir,"KO_genelist.txt\n",sep=""))
	fout <- file(paste(outdir,"KO_genelist.txt",sep=""),open="w")
	genes <- readLines(fin)
	for (i in 1:length(genes)){
		ko <- names(keggFind("ko",genes[i]))
		if(length(ko)==0){
			cat(paste("KO for",genes[i],"not found, line skipped\n"))
			next
		}
		writeLines(paste(genes[i],ko[1],sep="\t"),fout)
	}
	close(fin)
	close(fout)
}

# FUNCTIONS FOR MAKING CYTOSCAPE OUTPUT
makeNetwork <- function(fname,outdir){
	cat(paste(outdir,"KO_",fname,sep = ""))
}

# MAIN
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) promptInput() else parseInput(args)