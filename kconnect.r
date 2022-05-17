library("KEGGREST")
library("KEGGgraph")

# FUNCTIONS FOR PARSING INPUT
promptInput <- function(){
	cat("Options:\n
		Make network from input with KO ids: c [organism] [filename] [outputdir]\n
		Make network from input with only gene names: cg [organism] [filename] [outputdir]\n")
	cat("Please enter command: ")
	userinp <- readLines("stdin",n=1)
	parseInput(strsplit(userinp," "))
}
parseInput <- function(s){
	dir <- makeOutDir(s[4])
	# check for each option
	if(s[1] == "cg"){ # get ids and make network
		getIDs(s[3],s[2],dir)
		makeNetwork(paste(dir,s[2],"_genelist.txt",sep=""),s[2],dir)
	}else if(s[1] == "c"){ # already have ids
		makeNetwork(s[3],s[2],dir)
	}else{ # exit error
		stop(paste("Invalid argument:",s[1]))
	}
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

# FUNCTIONS FOR GETTING IDS
getIDs <- function(fname,org,outdir){
	fin <- file(fname,open="r") # read file with one gene name per line
	fout <- file(paste(outdir,org,"_genelist.txt",sep=""),open="w")
	genes <- readLines(fin)
	for(i in 1:length(genes)){
		found = FALSE
		kf <- keggFind(org,genes[i])
		for(hit in names(kf)){ # check for exact match across each entry of keggfind
			gnames <- unlist(strsplit(kf[hit],";"))[1] # KEGG puts known aliases before a semicolon, then separates with commas
			if(genes[i] %in% unlist(strsplit(gnames,", "))){
				found = TRUE
				writeLines(paste(genes[i],hit,sep="\t"),fout)
				break
			}
		}
		if(!found)
			cat(paste("ID for",genes[i],"in organism",org,"not found, line skipped\n"))
	}
	close(fin)
	close(fout)
	cat(paste("Wrote new ids of genes to ",outdir,org,"_genelist.txt\n",sep=""))
}

# FUNCTIONS FOR MAKING CYTOSCAPE OUTPUT
makeNetwork <- function(fname,org,outdir){
	# read input file and save gene ids
	genes <- c()
	namekey <- new.env(hash = TRUE)
	fin <- file(fname,open="r")
	lines <- readLines(fin)
	for(i in 1:length(lines)){
		lsplit <- unlist(strsplit(lines[i],"\t"))
		genes <- append(genes,lsplit[2])
		namekey[[lsplit[2]]] <- lsplit[1] 
	}
	close(fin)
	# map between pathways and target genes within them
	pathsfound <- keggLink("pathway",genes)
	pathsmapped <- new.env(hash = TRUE)
	for(i in seq_along(pathsfound)){
		p <- pathsfound[[i]]
		if(! p %in% names(pathsmapped)) # if pathway isn't initialized yet, initialize it
			pathsmapped[[p]] <- c()
		pathsmapped[[p]] <- append(pathsmapped[[p]],names(pathsfound)[i])
	}
	# get interactions between target genes in each pathway found
	cytoedges <- file(paste(outdir,"cytograph.sif",sep=""),open="w")
	allfoundpaths <- c()
	for(targpath in names(pathsmapped)){
		if(length(pathsmapped[[targpath]])<2) # only care about pathways with two or more genes
			next
		tmp <- tempfile()
		kgml <- retrieveKGML(targpath, organism=org, destfile=tmp, method="wget", quiet=TRUE)
		parsed <- parseKGML2Graph(kgml,expandGenes=TRUE, genesOnly = TRUE)
		gedges <- edges(parsed)
		# traverse graph downstream from each gene and make directed connection if it exists
		cat("taking new path\n")
		genes <- pathsmapped[[targpath]]
		for(groot in genes){
			newg <- new.env()
			newg$edges <- gedges
			newg$root <- groot
			newg$targets <- genes[genes!=groot]
			newg$foundpaths <- allfoundpaths
			traverseGene(newg,groot,parsed,namekey,cytoedges) # find if path exists from gene to any target
			allfoundpaths <- newg$foundpaths
		}
	}
	close(cytoedges)
}
traverseGene <- function(myenv,gene,mapkG,namekey,fout){
	if(length(myenv$edges[[gene]]) == 0){ # reached a dead end
		return(myenv)
	}
	children <- myenv$edges[[gene]]
	myenv$edges[[gene]] <- c()  # remove node after traversed
	for(c in children){
		if(c %in% myenv$targets){
			edat <- getKEGGedgeData(mapkG,paste(myenv$root,"~",c,sep="")) # get edge info so we can mark type
			if(is.null(edat))
				etype <- "indirect"
			else
				etype <- gsub(" ","",getName(getSubtype(edat)[[1]]))
			newconnect <- paste(namekey[[myenv$root]],etype,namekey[[c]],sep=" ")
			if(! newconnect %in% myenv$foundpaths){
				myenv$foundpaths <- append(myenv$foundpaths,newconnect)
				writeLines(newconnect,fout)
			}
		}
		# go as far down as possible, updating graph to avoid retracing
		traverseGene(myenv,c,mapkG,namekey,fout)
	}
}


# MAIN
args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0) promptInput() else parseInput(args)