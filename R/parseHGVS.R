
#' HGVS Parser 
#'
#' Parses HGVS strings
#' @param strings A character vector containing the HGVS strings
#' @param aacode allowed values: 1, 3, or NA. Determines whether 1-letter codes or 3-letter codes should be forced. NA uses input format.
#' @return A \code{data.frame} with the following columns: 
#' @keywords HGVS parsing
#' @export
#' @examples
#' parseHGVS("g.1318G>T")

parseHGVS <- function(strings,aacode=c(NA,1,3)) {

	#Check that parameters are valid
	if (!is.character(strings)) {
		stop("Input for 'parse' function must be a character vector! Found '",class(strings),"' instead.")
	}

	aacode <- aacode[[1]]
	if (!is.na(aacode) && !(aacode %in% c(1,3))) {
		warning("Invalid aacode parameter, defaulting to NA!")
		aacode <- NA
	}
	
	#Helper function: turns a list of lists (lol) in to a dataframe
	to.df <- function(lol) {
		colnames <- unique(do.call(c,lapply(lol,names)))
		columns <- lapply(colnames,function(cn) sapply(lol,function(row) {
			if (cn %in% names(row)) row[[cn]] else NA
		}))
		names(columns) <- colnames
		empty <- which(sapply(columns,function(xs)all(is.na(xs))))
		columns[empty] <- NULL
		do.call(data.frame,columns)
	}

	# ###
	# # Binds matrices of same size together to a 3D matrix, analogously
	# # to cbind and rbind.
	# #
	# zbind <- function(...) {
	# 	x <- list(...)
	# 	y <- array(0,dim=c(nrow(x[[1]]),ncol(x[[1]]),length(x)),dimnames=dimnames(x[[1]]))
	# 	for (i in 1:length(x)) y[,,i] <- x[[i]]
	# 	y
	# }


	###
	# Function to *locally* excise regex groups from string vectors.
	# I.e. only extract the first occurrence of each group within each string.
	# x = string vector
	# re = regular expression with groups
	#
	extract.groups <- function(x, re) {
		matches <- regexpr(re,x,perl=TRUE)
		start <- attr(matches,"capture.start")
		end <- start + attr(matches,"capture.length") - 1
		do.call(cbind,lapply(1:ncol(start), function(i) {
			sapply(1:nrow(start),function(j){
				if (start[j,i] > -1) substr(x[[j]],start[j,i],end[j,i]) else NA
			})
		}))
	}

	# ###
	# # Function to *globally* excise regex groups from string vectors.
	# # x = string vector
	# # re = regular expression with groups
	# #
	# global.extract.groups <- function(x,re) {
	#     all.matches <- gregexpr(re,x,perl=TRUE)
	#     mapply(function(matches,x) {
	#         start <- attr(matches,"capture.start")
	#         end <- start + attr(matches,"capture.length") - 1
	#         apply(zbind(start,end),c(1,2),function(pos) substr(x,pos[[1]],pos[[2]]) )
	#     },matches=all.matches,x=x,SIMPLIFY=FALSE)
	# }


	###
	# Helper function:
	# Given an HGVS body and a list of regexes corresponding to types,
	# find the (first) matching type.
	findType <- function(body,types) {
		i <- 0
		found <- done <- FALSE
		while (!found && !done) {
			found <- regexpr(types[[i <- i+1]],body) > 0
			done <- i >= length(types)
		}
		if (found) {
			return(names(types)[[i]])
		} else {
			return("invalid")
		}
	}

	out <- lapply(strings,function(s) {

		if (regexpr("^[gcnmrp]\\.",s) < 1) {
			return(list(hgvs=s,type="invalid"))
		}

		body <- substr(s,3,nchar(s))

		subjects <- c(
			g="genomic",c="coding",n="noncoding",
			m="mitochondrial",r="rna",p="protein"
		)
		subject <- subjects[[substr(s,1,1)]]

		if (subject=="genomic") {

			types <- c(
				substitution="\\d+[ACGT]>[ACGT]", singledeletion="^\\d+del$",
				deletion="\\d+_\\d+del$",inversion="\\d+_\\d+inv",
				duplication="\\d+_\\d+dup",insertion="\\d+_\\d+ins[ATCG]+",
				conversion="\\d+_\\d+con\\d+_\\d+",delins="\\d+_\\d+delins[ATCG]+",
				amplification="\\d+_\\d+\\[\\d+\\]"
			)
			
			type <- findType(body,types)

			if (type == "substitution") {
				groups <- extract.groups(body,"(\\d+)([ACGT])>([ACGT])")[1,]
				position <- as.integer(groups[[1]])
				ancestral <- groups[[2]]
				variant <- groups[[3]]
				return(list(hgvs=s,subject=subject,type=type,start=position,
					ancestral=ancestral,variant=variant))

			} else if (type == "singledeletion") {
				groups <- extract.groups(body,"(\\d+)del")[1,]
				position <- as.integer(groups[[1]])
				return(list(hgvs=s,subject=subject,type=type,start=position))

			} else if (type == "deletion") {
				groups <- extract.groups(body,"(\\d+)_(\\d+)del")
				start <- as.integer(groups[[1]])
				end <- as.integer(groups[[2]])
				return(list(hgvs=s,subject=subject,type=type,start=start,end=end))

			} else if (type == "inversion") {
				groups <- extract.groups(body,"(\\d+)_(\\d+)inv")
				start <- as.integer(groups[[1]])
				end <- as.integer(groups[[2]])
				return(list(hgvs=s,subject=subject,type=type,start=start,end=end))

			} else if (type == "duplication") {
				groups <- extract.groups(body,"(\\d+)_(\\d+)dup")
				start <- as.integer(groups[[1]])
				end <- as.integer(groups[[2]])
				return(list(hgvs=s,subject=subject,type=type,start=start,end=end))

			} else if (type == "insertion") {
				groups <- extract.groups(body,"(\\d+)_(\\d+)ins([ATCG]+)")
				start <- as.integer(groups[[1]])
				end <- as.integer(groups[[2]])
				if (abs(end-start)!=1) {
					warning("Invalid insertion definition: 
						Start and end positions must be adjacent!")
				}
				variant <- groups[[3]]
				return(list(hgvs=s,subject=subject,type=type,start=start,end=end,variant=variant))

			} else if (type == "conversion") {
				groups <- extract.groups(body,"(\\d+)_(\\d+)con(\\d+)_(\\d+)")
				start <- as.integer(groups[[1]])
				end <- as.integer(groups[[2]])
				tStart <- as.integer(groups[[3]])
				tEnd <- as.integer(groups[[4]])
				return(list(hgvs=s,subject=subject,type=type,start=start,end=end,
					templateStart=tStart,templateEnd=tEnd))

			} else if (type == "delins") {
				groups <- extract.groups(body,"(\\d+)_(\\d+)delins([ATCG]+)")
				start <- as.integer(groups[[1]])
				end <- as.integer(groups[[2]])
				variant <- groups[[3]]
				return(list(hgvs=s,subject=subject,type=type,start=start,end=end,variant=variant))

			} else if (type == "amplification") {
				groups <- extract.groups(body,"(\\d+)_(\\d+)\\[(\\d+)\\]")
				start <- as.integer(groups[[1]])
				end <- as.integer(groups[[2]])
				copies <- as.integer(groups[[3]])
				return(list(hgvs=s,subject=subject,type=type,start=start,end=end,copies=copies))

			} else {
				return(list(hgvs=s,type="invalid"))
			}

		} else if (subject=="coding") {

			types <- c(
				substitution="\\d+([+-]\\d+)?[ACGT]>[ACGT]", 
				singledeletion="^\\d+([+-]\\d+)?del$",
				deletion="\\d+([+-]\\d+)?_\\d+([+-]\\d+)?del$",
				inversion="\\d+([+-]\\d+)?_\\d+([+-]\\d+)?inv",
				duplication="\\d+([+-]\\d+)?_\\d+([+-]\\d+)?dup",
				insertion="\\d+([+-]\\d+)?_\\d+([+-]\\d+)?ins[ATCG]+",
				conversion="\\d+([+-]\\d+)?_\\d+([+-]\\d+)?con\\d+([+-]\\d+)?_\\d+([+-]\\d+)?",
				delins="\\d+([+-]\\d+)?_\\d+([+-]\\d+)?delins[ATCG]+",
				amplification="\\d+([+-]\\d+)?_\\d+([+-]\\d+)?\\[\\d+\\]"
			)
			
			type <- findType(body,types)

			if (type == "substitution") {
				groups <- extract.groups(body,"(\\d+)([+-]\\d+)?([ACGT])>([ACGT])")[1,]
				position <- as.integer(groups[[1]])
				intronOffset <- as.integer(groups[[2]])
				ancestral <- groups[[3]]
				variant <- groups[[4]]
				return(list(hgvs=s,subject=subject,type=type,start=position,
					startIntron=intronOffset,ancestral=ancestral,variant=variant))

			} else if (type == "singledeletion") {
				groups <- extract.groups(body,"(\\d+)([+-]\\d+)?del")[1,]
				position <- as.integer(groups[[1]])
				intronOffset <- as.integer(groups[[2]])
				return(list(hgvs=s,subject=subject,type=type,start=position,startIntron=intronOffset))

			} else if (type == "deletion") {
				groups <- extract.groups(body,"(\\d+)([+-]\\d+)?_(\\d+)([+-]\\d+)?del")
				start <- as.integer(groups[[1]])
				intronOffset <- as.integer(groups[[2]])
				end <- as.integer(groups[[3]])
				intronOffset2 <- as.integer(groups[[4]])
				return(list(hgvs=s,subject=subject,type=type,start=start,startIntron=intronOffset,
					end=end,endIntron=intronOffset2))

			} else if (type == "inversion") {
				groups <- extract.groups(body,"(\\d+)([+-]\\d+)?_(\\d+)([+-]\\d+)?inv")
				start <- as.integer(groups[[1]])
				intronOffset <- as.integer(groups[[2]])
				end <- as.integer(groups[[3]])
				intronOffset2 <- as.integer(groups[[4]])
				return(list(hgvs=s,subject=subject,type=type,start=start,startIntron=intronOffset,
					end=end,endIntron=intronOffset2))

			} else if (type == "duplication") {
				groups <- extract.groups(body,"(\\d+)([+-]\\d+)?_(\\d+)([+-]\\d+)?dup")
				start <- as.integer(groups[[1]])
				intronOffset <- as.integer(groups[[2]])
				end <- as.integer(groups[[3]])
				intronOffset2 <- as.integer(groups[[4]])
				return(list(hgvs=s,subject=subject,type=type,start=start,startIntron=intronOffset,
					end=end,endIntron=intronOffset2))

			} else if (type == "insertion") {
				groups <- extract.groups(body,"(\\d+)([+-]\\d+)?_(\\d+)([+-]\\d+)?ins([ATCG]+)")
				start <- as.integer(groups[[1]])
				intronOffset <- as.integer(groups[[2]])
				end <- as.integer(groups[[3]])
				intronOffset2 <- as.integer(groups[[4]])
				variant <- groups[[5]]
				return(list(hgvs=s,subject=subject,type=type,start=start,startIntron=intronOffset,
					end=end,endIntron=intronOffset2,variant=variant))

			} else if (type == "conversion") {
				groups <- extract.groups(body,"(\\d+)([+-]\\d+)?_(\\d+)([+-]\\d+)?con(\\d+)([+-]\\d+)?_(\\d+)([+-]\\d+)?")
				start <- as.integer(groups[[1]])
				intronOffset <- as.integer(groups[[2]])
				end <- as.integer(groups[[3]])
				intronOffset2 <- as.integer(groups[[4]])
				tStart <- as.integer(groups[[5]])
				intronOffset3 <- as.integer(groups[[6]])
				tEnd <- as.integer(groups[[7]])
				intronOffset4 <- as.integer(groups[[8]])
				return(list(hgvs=s,subject=subject,type=type,start=start,startIntron=intronOffset,
					end=end,endIntron=intronOffset2,templateStart=tStart,
					templateStartIntron=intronOffset3,templateEnd=tEnd,
					templateEndIntron=intronOffset4))

			} else if (type == "delins") {
				groups <- extract.groups(body,"(\\d+)([+-]\\d+)?_(\\d+)([+-]\\d+)?delins([ATCG]+)")
				start <- as.integer(groups[[1]])
				intronOffset <- as.integer(groups[[2]])
				end <- as.integer(groups[[3]])
				intronOffset2 <- as.integer(groups[[4]])
				variant <- groups[[5]]
				return(list(hgvs=s,subject=subject,type=type,start=start,startIntron=intronOffset,
					end=end,endIntron=intronOffset2,variant=variant))

			} else if (type == "amplification") {
				groups <- extract.groups(body,"(\\d+)([+-]\\d+)?_(\\d+)([+-]\\d+)?\\[(\\d+)\\]")
				start <- as.integer(groups[[1]])
				intronOffset <- as.integer(groups[[2]])
				end <- as.integer(groups[[3]])
				intronOffset2 <- as.integer(groups[[4]])
				copies <- as.integer(groups[[3]])
				return(list(hgvs=s,subject=subject,type=type,start=start,startIntron=intronOffset,
					end=end,endIntron=intronOffset2,copies=copies))

			} else {
				return(list(hgvs=s,type="invalid"))
			}

		} else if (subject=="protein") {

			one2three <- c(A="Ala",C="Cys",D="Asp",E="Glu",F="Phe",G="Gly",H="His",
				I="Ile",K="Lys",L="Leu",M="Met",N="Asn",P="Pro",Q="Gln",R="Arg",
				S="Ser",T="Thr",V="Val",W="Trp",Y="Tyr",`*`="Ter")
			three2one <- c(Ala="A",Arg="R",Asn="N",Asp="D",Cys="C",Gln="Q",Glu="E",
				Gly="G",His="H",Ile="I",Leu="L",Lys="K",Met="M",Phe="F",Pro="P",
				Ser="S",Thr="T",Trp="W",Tyr="Y",Val="V",Ter="*")
			codes <- paste(c(one2three,three2one[-21],"\\*"),collapse="|")

			types <- list(
				substitution=paste0("^(",codes,")(\\d+)(",codes,")$"),
				singledeletion=paste0("^(",codes,")(\\d+)del$"),
				deletion=paste0("^(",codes,")(\\d+)_(",codes,")(\\d+)del$"),
				duplication=paste0("^(",codes,")(\\d+)_(",codes,")(\\d+)dup$"),
				insertion=paste0("^(",codes,")(\\d+)_(",codes,")(\\d+)ins((",codes,")+)$"),
				delins=paste0("^(",codes,")(\\d+)_(",codes,")(\\d+)delins((",codes,")+)$"),
				frameshift1=paste0("^(",codes,")(\\d+)fs$"),
				frameshift2=paste0("^(",codes,")(\\d+)(",codes,")fs(Ter|\\*)(\\d+)$")
			)
			
			type <- findType(body,types)

			if (type == "substitution") {
				groups <- extract.groups(body,types$substitution)
				aa1 <- groups[[1]]
				pos <- as.integer(groups[[2]])
				aa2 <- groups[[3]]
				if (aa1 %in% c(one2three,three2one) && aa2 %in% c(one2three,three2one)) {
					if (is.na(aacode)) {
						#do nothing
					} else if (aacode == 1) {
						if (nchar(aa1) == 3) aa1 <- three2one[[aa1]]
						if (nchar(aa2) == 3) aa2 <- three2one[[aa2]]
					} else if (aacode ==3) {
						if (nchar(aa1) == 1) aa1 <- one2three[[aa1]]
						if (nchar(aa2) == 1) aa2 <- one2three[[aa2]]
					} else {
						#this should never happen, as it's supposed to be detected at start of function
						stop("Invalid aacode. If you see this, report this as a bug!")
					}
					return(list(hgvs=s,subject=subject,type=type,start=pos,
						ancestral=aa1,variant=aa2))
				} else {#not valid amino acid
					return(list(hgvs=s,type="invalid"))
				}

			} else if (type == "singledeletion") {
				groups <- extract.groups(body,types$singledeletion)
				aa1 <- groups[[1]]
				pos <- as.integer(groups[[2]])
				if (is.na(aacode)) {
					#do nothing
				} else if (aacode == 1) {
					if (nchar(aa1) == 3) aa1 <- three2one[[aa1]]
				} else if (aacode == 3) {
					if (nchar(aa1) == 1) aa1 <- one2three[[aa1]]
				} else {
					#this should never happen, as it's supposed to be detected at start of function
					stop("Invalid aacode. If you see this, report this as a bug!")
				}
				return(list(hgvs=s,subject=subject,type=type,start=pos,ancestral=aa1))
			} else if (type == "deletion") {
				groups <- extract.groups(body,types$deletion)
				aa1 <- groups[[1]]
				pos <- as.integer(groups[[2]])
				aa2 <- groups[[3]]
				pos2 <- as.integer(groups[[4]])
				if (is.na(aacode)) {
					#do nothing
				} else if (aacode == 1) {
					if (nchar(aa1) == 3) aa1 <- three2one[[aa1]]
					if (nchar(aa2) == 3) aa2 <- three2one[[aa2]]
				} else if (aacode == 3) {
					if (nchar(aa1) == 1) aa1 <- one2three[[aa1]]
					if (nchar(aa2) == 1) aa2 <- one2three[[aa2]]
				} else {
					#this should never happen, as it's supposed to be detected at start of function
					stop("Invalid aacode. If you see this, report this as a bug!")
				}
				return(list(hgvs=s,subject=subject,type=type,start=pos,
					ancestral=aa1,end=pos2,ancestral2=aa2))	

			} else if (type == "duplication") {
				groups <- extract.groups(body,types$duplication)
				aa1 <- groups[[1]]
				pos <- as.integer(groups[[2]])
				aa2 <- groups[[3]]
				pos2 <- as.integer(groups[[4]])
				if (is.na(aacode)) {
					#do nothing
				} else if (aacode == 1) {
					if (nchar(aa1) == 3) aa1 <- three2one[[aa1]]
					if (nchar(aa2) == 3) aa2 <- three2one[[aa2]]
				} else if (aacode == 3) {
					if (nchar(aa1) == 1) aa1 <- one2three[[aa1]]
					if (nchar(aa2) == 1) aa2 <- one2three[[aa2]]
				} else {
					#this should never happen, as it's supposed to be detected at start of function
					stop("Invalid aacode. If you see this, report this as a bug!")
				}
				return(list(hgvs=s,subject=subject,type=type,start=pos,
					ancestral=aa1,end=pos2,ancestral2=aa2))	

			} else if (type == "insertion") {
				groups <- extract.groups(body,types$insertion)
				aa1 <- groups[[1]]
				pos <- as.integer(groups[[2]])
				aa2 <- groups[[3]]
				pos2 <- as.integer(groups[[4]])
				insert <- groups[[5]]
				#TODO: Implement code conversion 
				return(list(hgvs=s,subject=subject,type=type,start=pos,
					ancestral=aa1,end=pos2,ancestral2=aa2,variant=insert))	

			} else if (type == "delins") {
				groups <- extract.groups(body,types$delins)
				aa1 <- groups[[1]]
				pos <- as.integer(groups[[2]])
				aa2 <- groups[[3]]
				pos2 <- as.integer(groups[[4]])
				insert <- groups[[5]]
				#TODO: Implement code conversion 
				return(list(hgvs=s,subject=subject,type=type,start=pos,
					ancestral=aa1,end=pos2,ancestral2=aa2,variant=insert))	

			} else if (type == "frameshift1") {
				groups <- extract.groups(body,types$frameshift1)
				aa1 <- groups[[1]]
				pos <- as.integer(groups[[2]])
				#TODO: Implement code conversion 
				return(list(hgvs=s,subject=subject,type="frameshift",start=pos,ancestral=aa1))

			} else if (type == "frameshift2") {
				groups <- extract.groups(body,types$frameshift2)
				aa1 <- groups[[1]]
				pos <- as.integer(groups[[2]])
				aa2 <- groups[[3]]
				term <- as.integer(groups[[5]])
				#TODO: Implement code conversion 
				return(list(hgvs=s,subject=subject,type="frameshift",start=pos,ancestral=aa1,variant=aa2,end=term))

			} else {#unmatched type
				return(list(hgvs=s,type="invalid"))
			}


		} else if (subject=="noncoding") {
			return("Not implemented! If you see this, report it as a bug!")
		} else if (subject=="mitochondrial") {
			return("Not implemented! If you see this, report it as a bug!")
		} else if (subject=="rna") {
			return("Not implemented! If you see this, report it as a bug!")
		} else {#unmatched subject, shouldn't happen
			stop("Unmatched subject! If you see this, report it as a bug!")
		}
	})
	to.df(out)
}