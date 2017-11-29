
#' Genomic HGVS Builder 
#'
#' A constructor for a genomic HGVS builder object. The object contains a collection of functions
#' for building genomic HGVS strings.
#' @return A \code{hgvs.builder.g} object with functions for building genomic HGVS strings.
#' @keywords HGVS builder
#' @export
#' @examples
#' builder <- new.hgvs.builder.g()
#' string <- builder$substitution(123,"A","G")

new.hgvs.builder.g <- function() {

	substitution <- function(pos,ancestral,variant) {
		if (!is.numeric(pos) || pos < 1) stop("position must be a positive integer")
		if (!is.character(ancestral) || !(ancestral %in% c("A","C","G","T"))) stop("ancestral must be single nucleotide")
		if (!is.character(variant) || !(variant %in% c("A","C","G","T"))) stop("variant must be single nucleotide")
		paste0("g.",pos,ancestral,">",variant)
	}

	deletion <- function(start,stop) {
		if (!is.numeric(start)) stop("start must be an integer")
		if (!is.numeric(stop)) stop("stop must be an integer")
		paste0("g.",start,"_",stop,"del")
	}

	inversion <- function(start,stop) {
		if (!is.numeric(start)) stop("start must be an integer")
		if (!is.numeric(stop)) stop("stop must be an integer")
		paste0("g.",start,"_",stop,"inv")
	}

	duplication <- function(start,stop) {
		if (!is.numeric(start)) stop("start must be an integer")
		if (!is.numeric(stop)) stop("stop must be an integer")
		paste0("g.",start,"_",stop,"dup")
	}

	insertion <- function(start,seq) {
		if (!is.numeric(start)) stop("start must be an integer")
		if (!is.character(seq) || regexpr("^[ACGT]+$",seq) < 1) stop("variant must be nucleotide sequence")
		paste0("g.",start,"_",start+1,"ins",seq)
	}

	delins <- function(start,stop,seq) {
		if (!is.numeric(start)) stop("start must be an integer")
		if (!is.numeric(stop)) stop("stop must be an integer")
		if (!is.character(seq) || regexpr("^[ACGT]+$",seq) < 1) stop("variant must be nucleotide sequence")
		paste0("g.",start,"_",stop,"delins",seq)
	}

	return(structure(list(
		substitution=substitution,
		deletion=deletion,
		inversion=inversion,
		duplication=duplication,
		insertion=insertion,
		delins=delins
	),class="hgvs.builder.g"))
}

new.hgvs.builder.c <- function() {

	substitution <- function(pos,ancestral,variant,posOffset=0) {
		if (!is.numeric(pos) || pos < 1) stop("position must be a positive integer")
		if (!is.character(ancestral) || !(ancestral %in% c("A","C","G","T"))) stop("ancestral must be single nucleotide")
		if (!is.character(variant) || !(variant %in% c("A","C","G","T"))) stop("variant must be single nucleotide")
		offsetStr <- if (posOffset==0) {
			""
		} else if (posOffset > 0) {
			paste0("+",posOffset)
		} else if (posOffset < 0) {
			as.character(posOffset)
		}
		paste0("c.",pos,offsetStr,ancestral,">",variant)
	}

	# deletion <- function(start,stop,startOffset=0,stopOffset=0) {
	# 	if (!is.numeric(start)) stop("start must be an integer")
	# 	if (!is.numeric(stop)) stop("stop must be an integer")
	# 	paste0("c.",start,"_",stop,"del")
	# }

	# inversion <- function(start,stop,startOffset=0,stopOffset=0) {
	# 	if (!is.numeric(start)) stop("start must be an integer")
	# 	if (!is.numeric(stop)) stop("stop must be an integer")
	# 	paste0("c.",start,"_",stop,"inv")
	# }

	# duplication <- function(start,stop,startOffset=0,stopOffset=0) {
	# 	if (!is.numeric(start)) stop("start must be an integer")
	# 	if (!is.numeric(stop)) stop("stop must be an integer")
	# 	paste0("c.",start,"_",stop,"dup")
	# }

	# insertion <- function(start,seq,startOffset=0,stopOffset=0) {
	# 	if (!is.numeric(start)) stop("start must be an integer")
	# 	if (!is.character(seq) || regexpr("^[ACGT]+$",seq) < 1) stop("variant must be nucleotide sequence")
	# 	paste0("c.",start,"_",start+1,"ins",seq)
	# }

	# delins <- function(start,stop,seq,startOffset=0,stopOffset=0) {
	# 	if (!is.numeric(start)) stop("start must be an integer")
	# 	if (!is.numeric(stop)) stop("stop must be an integer")
	# 	if (!is.character(seq) || regexpr("^[ACGT]+$",seq) < 1) stop("variant must be nucleotide sequence")
	# 	paste0("c.",start,"_",stop,"delins",seq)
	# }

	return(structure(list(
		substitution=substitution#,
		# deletion=deletion,
		# inversion=inversion,
		# duplication=duplication,
		# insertion=insertion,
		# delins=delins
	),class="hgvs.builder.c"))
}
