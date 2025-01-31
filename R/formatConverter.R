#' Convert GInteractions to linkSet
#' @aliases Convert
#' @description Convert other data formats to linkSet. Currently supported: GInteractions, data.frame.
#' 
#' @param x A GInteractions object
#' @param ... Additional arguments (not used)
#' 
#' @rdname Convert
#' 
#' @return A linkSet object
#' @export
#' @importFrom S4Vectors mcols mcols<-
#' @examples
#' library(InteractionSet)
#' gi <- GInteractions(anchor1 = c(1, 2), anchor2 = c(3, 4), regions = GRanges(seqnames = c("chr1", "chr1", "chr2", "chr2"),
#' ranges = IRanges(start = c(100, 200, 300, 400), width = 50)))
#' ls <- Convert(gi)
#' ls
#' 
setMethod("Convert", signature(x = "GInteractions"), function(x, baitCol = NULL, ...) {
  anchor1 <- x@anchor1
  anchor2 <- x@anchor2
  regions <- x@regions
  metadata <- as.data.frame(x@elementMetadata)

  if (!is.null(baitCol) && !baitCol %in% colnames(metadata)) {
    warning("baitCol not found in metadata, using first regions as bait")
    baitCol <- NULL
  }

  if (is.null(baitCol)) {
    anchor1Region <- regions[anchor1]
    nameBait <- paste0(anchor1Region)
  } else {
    nameBait <- metadata[[baitCol]]
  }
  ls <- linkSet(
    anchor1 = regions[anchor1],
    anchor2 = regions[anchor2],
    specificCol =  nameBait
  )
  mcols(ls) <- metadata
  return(ls)
})

#' Convert string intervals to GRanges
#' @importFrom S4Vectors Rle
#' @keywords internal
.convert_to_grange <- function(intervals) {
  # convert "chr1.816066.816566" or "chr1:816066-816566" to grange format
  parts <- strsplit(intervals, "[.:\\-]")
  
  # Check if all parts have exactly 3 elements
  if (!all(vapply(parts, length, FUN.VALUE = numeric(1)) == 3)) {
    stop("Please input peak format like chr1.816066.816566 or chr1:816066-816566")
  }
  
  # Extract components with proper type specification
  chromosomes <- vapply(parts, `[`, FUN.VALUE = character(1), 1)  # First element is character (chromosome)
  starts <- as.numeric(vapply(parts, `[`, FUN.VALUE = character(1), 2), USE.NAMES = FALSE)  # Convert to numeric after extraction
  ends <- as.numeric(vapply(parts, `[`, FUN.VALUE = character(1), 3), USE.NAMES = FALSE)    # Convert to numeric after extraction
  # Create a GRanges object
  granges_obj <- GRanges(
    seqnames = Rle(chromosomes),
    ranges = IRanges(start = starts, end = ends)
  )
  return(granges_obj)
}
#' Convert data.frame to linkSet
#'
#' @param x A data.frame object
#' @param source The source of the data frame, either "data.frame" or "chicane"
#' @param baitCol The column name in the data frame that contains the bait information
#' @param oeCol The column name in the data frame that contains the other end information
#' @param ... Additional arguments (not used)
#' @rdname Convert
#' @return A linkSet object
#' @export
#' @examples
#' df <- data.frame(
#'   gene = c("gene1", "gene2"),
#'   peak = c("chr1:1000-2000", "chr2:1500-2500"),
#'   score = c(5.5, 6.0)
#' )
#' ls <- Convert(df, source = "data.frame", baitCol = "gene", oeCol = "peak")
#' ls
setMethod("Convert", signature(x = "data.frame"), function(x, source = "data.frame", baitCol = "gene", oeCol = "peak", ...) {
  if (source == "chicane") {
    required_cols <- c("target.id", "bait.id", "bait.chr", "bait.start", "bait.end", "target.chr", "target.start", "target.end")
    if (!all(required_cols %in% colnames(x))) {
      stop("chicane data.frame must contain the columns:\n",
           paste(required_cols, collapse = ' '))
    }
    # Create GRanges for bait and other end
    bait_gr <- GRanges(
      seqnames = x$bait.chr,
      ranges = IRanges(start = x$bait.start, end = x$bait.end)
    )
    oe_gr <- GRanges(
      seqnames = x$target.chr,
      ranges = IRanges(start = x$target.start, end = x$target.end)
    )
        # Create GRanges for bait and other end
    bait_gr <- GRanges(
      seqnames = x$bait.chr,
      ranges = IRanges(start = x$bait.start, end = x$bait.end)
    )
    oe_gr <- GRanges(
      seqnames = x$target.chr,
      ranges = IRanges(start = x$target.start, end = x$target.end)
    )

    # Create linkSet object
    ls <- linkSet(
      anchor1 = bait_gr,
      anchor2 = oe_gr,
      specificCol = x[["bait.id"]]
    )

    # Add metadata
    metadata_cols <- setdiff(colnames(x), c("bait.chr", "bait.start", "bait.end", "target.chr", "target.start", "target.end"))
    if (length(metadata_cols) > 0) {
      mcols(ls) <- x[, metadata_cols, drop = FALSE]
    }
    return(ls)
  }
  # Implementation for data.frame conversion
  # This is a placeholder and needs to be implemented based on your data.frame structure
  else if (source == "data.frame") {

    if (!baitCol %in% colnames(x) | !oeCol %in% colnames(x)) {
      stop("baitCol and oeCol must be columns in the data.frame")
    }
    bait <- x[[baitCol]]
    oe <- x[[oeCol]]
    oeGrange <- .convert_to_grange(oe)
    metadata <- x[ , !(colnames(x) %in% c(baitCol, oeCol)), drop = FALSE]
    ls <- linkSet(
      anchor1 = bait,
      anchor2 = oeGrange
    )
    mcols(ls) <- metadata
    return(ls)
  }
})

#' Convert Pairs to linkSet
#' @param x A Pairs object
#' @param baitCol A character string specifying the column to use for bait naming
#' 
#' @rdname Convert
#'
#' @return A linkSet object
#' @export
setMethod("Convert", signature(x = "Pairs"), function(x,baitCol = NULL, ...) {
  # Extract first and second GRanges
  anchor1 <- x@first
  anchor2 <- x@second
  
  # Extract metadata
  metadata <- as.data.frame(x@elementMetadata)
  
  # Create specificCol (using 'symbol' if available, otherwise use ranges)
  if ("symbol" %in% colnames(metadata)) {
    specificCol <- metadata$symbol
  } else {
    specificCol <- paste0(anchor1)
  }

  if (is.null(baitCol)) {
    nameBait <- paste0(anchor1)
  } else {
    nameBait <- metadata[[baitCol]]
  }
  # Create linkSet object
  ls <- linkSet(
    anchor1 = anchor1,
    anchor2 = anchor2,
    specificCol = nameBait
  )
  
  # Add metadata to linkSet
  mcols(ls) <- metadata
  
  return(ls)
})

#' Default conversion method
#'
#' @param x An object of unsupported class
#' @param ... Additional arguments (not used)
#' @param baitCol A character string specifying the column to use for bait naming
#'
#' @rdname Convert
#' @return Nothing, throws an error
#' @export
setMethod("Convert", signature(x = "ANY"), function(x, baitCol = NULL, ...) {
  if (inherits(x, "GenomicInteractions")) {
    if (!requireNamespace("GenomicInteractions", quietly = TRUE)) {
      stop("Package 'GenomicInteractions' is needed for this function to work. Please install it.",
         call. = FALSE)
    }
    # Extract anchor1 and anchor2 GRanges
    anchor1 <- x@regions[x@anchor1]
    anchor2 <- x@regions[x@anchor2]
  
    # Extract metadata
    metadata <- as.data.frame(x@elementMetadata)
    
    # Create specificCol (using 'InteractionID' if available, otherwise use ranges)
    if ("InteractionID" %in% colnames(metadata)) {
      specificCol <- metadata$InteractionID
    } else {
      specificCol <- paste0(anchor1)
    }

    if (is.null(baitCol)) {
      nameBait <- paste0(anchor1)
    } else if (baitCol %in% colnames(metadata)) {
      nameBait <- metadata[[baitCol]]
    } else {
      warning("baitCol not found in metadata, using first regions as bait")
      nameBait <- paste0(anchor1)
    }

    # Create linkSet object
    ls <- linkSet(
      anchor1 = anchor1,
      anchor2 = anchor2,
      specificCol = nameBait
    )
    
    # Add metadata to linkSet
    mcols(ls) <- metadata
    
    return(ls)
  } else if (inherits(x, "chicagoData")) {
    .exportToLinkSet(x, ...)
  }
   else {
    stop(paste("Conversion from", class(x), "to linkSet is not supported"))
  }
})



###############################################################
#' @rdname baitGInteractions
#' @param geneSymbol A character vector with same length as geneGr or column name in mcols(geneGr) for gene symbols
#' @examples
#' # Example usage:
#' library(GenomicRanges)
#' library(InteractionSet)
#' 
#' # Create example GRanges objects for genes and peaks
#' geneGr <- GRanges(seqnames = "chr1", 
#'                   ranges = IRanges(start = c(100, 200), end = c(150, 250)), 
#'                   geneSymbol = c("Gene1", "Gene2"))
#' peakGr <- GRanges(seqnames = "chr1", 
#'                   ranges = IRanges(start = c(300, 400), end = c(350, 450)))
#' 
#' # Create example GInteractions object
#' gi <- GInteractions(anchor1 = geneGr, anchor2 = peakGr)
#' 
#' # Convert to linkSet
#' linkSetObj <- baitGInteractions(gi, geneGr, peakGr, geneSymbol = "geneSymbol")
#' 
#' # Print the linkSet object
#' print(linkSetObj)
#' @export
setMethod("baitGInteractions", signature(x = "GInteractions", geneGr = "GRanges", peakGr = "GRanges"), function(x, geneGr, peakGr, geneSymbol=NULL) {
  gi <- x
  # validate geneSymbol in in mcols(geneGr)
  # validate length of geneSymbol is 1 or same length as geneGr
  if (length(geneSymbol) != 1 & length(geneSymbol) != length(geneGr)) {
    warning("geneSymbol must be a column name in mcols(geneGr) or a vector of the same length as geneGr, 
    using range to represent gene instead.")
    geneSymbol = NULL
  } else if (length(geneSymbol) == length(geneGr)) {
    geneSymbol <- geneSymbol
  } else if (!geneSymbol %in% colnames(mcols(geneGr))){
    warning("geneSymbol not found in geneGr, using range to represent gene instead.")
    geneSymbol = NULL
  } else {
    geneSymbol <- mcols(geneGr)[[geneSymbol]]
  }

  if (is.null(geneSymbol)) {
    geneSymbol <- paste0(geneGr)
  }
  out <- InteractionSet::linkOverlaps(gi, geneGr, peakGr)
  linkGene <- geneSymbol[out$subject1]
  linkGeneRange <- geneGr[out$subject1]
  linkPeak <- peakGr[out$subject2]
  linkSetObj <- linkSet(linkGeneRange, linkPeak, linkGene)
  
  # Add mcols from gi to linkSetObj
  mcols(linkSetObj) <- cbind(mcols(linkSetObj), mcols(gi)[out$query, ])
  
  return(linkSetObj)
})



###############################################################
#' Read validPairs file to GInteractions
#' @param file A character string specifying the path to the validPairs file or 4DN pairs file
#' @param format A character string specifying the format of the file, either "validPairs" or "pair". Pair format should be "readID chr1 pos1 chr2 pos2 strand1 strand2". And validPairs should be "readID chr1 pos1 strand1 chr2 pos2 strand2".
#' @param njobs An integer specifying the number of threads to use for reading the file
#' @return A GInteractions object
#' @rdname Convert
#' @export
#' 
#' 
readvalidPairs <- function(file, njobs = 1, format = "validPairs") {
  # Validate format parameter
  if (!format %in% c("validPairs", "pair")) {
    stop("format must be either 'validPairs' or 'pair'")
  }

  # Read file with data.table
  message("Reading file...")
  
  # Check if file is gzipped
  is_gzipped <- grepl("\\.gz$", file)
  
  # First check the header structure
  if(is_gzipped) {
    header_lines <- system(paste("zcat", file, "| grep '^#'"), intern = TRUE)
  } else {
    header_lines <- system(paste("grep '^#'", file), intern = TRUE)
  }
  skip_lines <- length(header_lines)
  
  # Check column structure from header
  if(format == "pair" && length(header_lines) > 0) {
    col_line <- header_lines[grep("#columns:", header_lines)]
    if(length(col_line) > 0) {
      cols <- strsplit(gsub("#columns:\\s*", "", col_line), "\\s+")[[1]]
    } else {
      cols <- c("readID", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2")
    }
  }
  
  if(format == "pair") {
    # Define base columns and classes
    base_cols <- c("readID", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2")
    base_classes <- c(
      "character", # readID
      "character", # chr1
      "numeric",   # pos1
      "character", # chr2
      "numeric",   # pos2
      "character", # strand1
      "character"  # strand2
    )
    
    dt <- data.table::fread(
      cmd = if(is_gzipped) paste("zcat", file) else file,
      skip = skip_lines,
      check.names = TRUE,
      nThread = njobs
    )
    
    # Rename columns based on position
    if(ncol(dt) >= length(base_cols)) {
      setnames(dt, 1:length(base_cols), base_cols)
    } else {
      stop("File has fewer columns than expected minimum")
    }
    
  } else {
    dt <- data.table::fread(
      cmd = if(is_gzipped) paste("zcat", file) else file,
      skip = skip_lines,
      col.names = c("readID", "chr1", "pos1", "strand1", "chr2", "pos2", "strand2"),
      colClasses = c(
        "character", # readID
        "character", # chr1
        "numeric",   # pos1
        "character", # strand1
        "character", # chr2
        "numeric",   # pos6
        "character"  # strand2
      ),
      check.names = TRUE,
      nThread = njobs
    )
  }
  
  # Filter for valid rows
  valid_rows <- dt
  if(nrow(dt) == 0) {
    stop("No valid interactions found after filtering.")
  }
  
  message("Converting to GInteractions...")
  
  # Create GRanges for anchor1
  gr1 <- GRanges(
    seqnames = valid_rows$chr1,
    ranges = IRanges(start = valid_rows$pos1, width = 1),
    strand = valid_rows$strand1
  )
  
  # Create GRanges for anchor2
  gr2 <- GRanges(
    seqnames = valid_rows$chr2,
    ranges = IRanges(start = valid_rows$pos2, width = 1),
    strand = valid_rows$strand2
  )
  
  # Create GInteractions object
  gi <- InteractionSet::GInteractions(anchor1 = gr1, anchor2 = gr2)
  
  # Add additional metadata if columns exist
  if(format == "pair") {
    if("pair_type" %in% names(dt)) mcols(gi)$pair_type <- dt$pair_type
    if("mapq1" %in% names(dt)) mcols(gi)$mapq1 <- dt$mapq1
    if("mapq2" %in% names(dt)) mcols(gi)$mapq2 <- dt$mapq2
  }
  
  # Clean up
  rm(dt, valid_rows, gr1, gr2)
  #gc()
  
  return(gi)
}


#' coerce linkSet to DataFrame
#' @param x A linkSet object
#' @return A DataFrame object
#' @export
#' @examples
#' # Create a linkSet object
#' data(linkExample)
#' # Convert linkSet to DataFrame
#' df <- as.data.frame(linkExample)
#' print(df)
setMethod("as.data.frame", "linkSet", function(x) {
  df <- data.frame(
    bait = bait(x),
    bait_region = paste0(regionsBait(x)),
    oe_region = paste0(oe(x)),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(mcols(x))) {
    df <- cbind(df, as.data.frame(mcols(x)))
  }
  
  return(df)
})


###############################################################
## Hidden "read" functions from chicago package ----------------------

.readRmap = function(s){
  fread(s$rmapfile, colClasses = list(character=1))
}

.readBaitmap = function(s){
  fread(s$baitmapfile, colClasses = list(character=1))
}

.readNPBfile = function(s){
  
  # Reads a pre-made text file containing the numbers of fragments per bait per distance bin 
  # within the interval maxl, given binsize.
  # The file can be generated by countNperBin.py and its first line should start with # and 
  # contain the parameter settings used. In addition to maxl and binsize, it defines the 
  # filtering parameters for restriction fragment minsize, maxsize as well as a boolean
  # variable removeb2b specifying whether bait2bait interactions should also not be counted.
  # (In fact, they probably should never be).
  
  # s is the current chicagoData object's settings list
  
  message("Reading NPerBin file...")
  header = readLines(s$nperbinfile, n=1)
  params = sapply(sapply(strsplit(header, "\t")[[1]],function(x)strsplit(x,"=")[[1]]), function(x)x[2])
  params = params[2:length(params)]
  names(params) = gsub("(\\S+)=.+", "\\1", names(params))
  minsize = as.numeric(params[["minFragLen"]])
  if (minsize != s$minFragLen){
    stop("The minFragLen in the NPerBin file header is not equal to minFragLen defined in experiment settings. Amend either setting (and if needed, generate a new NPerBin file) before running the analysis\n")
  }
  maxsize = as.numeric(params[["maxFragLen"]])
  if (maxsize != s$maxFragLen){
    stop("The maxFragLen in the NPerBin file header is not equal to maxFragLen defined in experiment settings. Amend either setting (and if needed, generate a new NPerBin file) before running the analysis\n")
  }
  maxl = as.numeric(params[["maxLBrownEst"]])
  if (maxl != s$maxLBrownEst){
    stop("The maxLBrownEst in the NPerBin file header is not equal to maxLBrownEst defined in experiment settings. Amend either setting (and if needed, generate a new NPerBin file) before running the analysis\n")
  }
  binsz = as.numeric(params[["binsize"]]) 
  if (binsz != s$binsize){
    stop("The binsize in the NPerBin file header is not equal to binsize defined in experiment settings. Amend either setting (and if needed, generate a new NPerBin file) before running the analysis\n")
  }  
  if (params[["removeb2b"]]!="True"){
    stop("The NPerBin file must be generated with removeb2b==True. Please generate a new file.\n")
  }
  if ( (params[["removeAdjacent"]]=="True" & !s$removeAdjacent) | (params[["removeAdjacent"]]!="True" & s$removeAdjacent)  ){
    stop("The removeAdjacent parameter settings used for generating NPerBin file (according to its header) and defined in experiment settings do not match. Amend either setting (and if needed, generate a new NPerBin file) before running the analysis\n")
  }  
  if(basename(params[["rmapfile"]]) != basename(s$rmapfile)){
    stop("The .rmap file used for generating the NPerBin file (according to the NPerBin header) and the one defined in experiment settings do not match. Amend either setting (and if needed, generate a new NPerBin file) before running the analysis\n")
  }
  
## Not checking this for now as we have a mixup of _baits and _baits_ID files used at different times...
#   if(basename(params[["baitmapfile"]]) != basename(baitmapfile)){
#     stop("Bait files used for generating the NfragPerBin file and defined here do not match. 
#          Amend either setting before running the analysis\n")
#   }

  npb = data.table::fread(s$nperbinfile, skip=1L)
  setnames(npb, names(npb)[1], "baitID")
  for(i in 2:ncol(npb)){
    setnames(npb, names(npb)[i], paste0("bin", i-1))    
  }
  npb
}

.readNbaitsPBfile = function(s){
  
  # Reads a pre-made text file containing the numbers of baits per other end per distance bin 
  # within the interval maxl, given binsize.
  # The file can be generated by countNBaitsPerBin.py and its first line should start with # and 
  # contain the parameter settings used. 
  
  # s is the current chicagoData object's settings list
  
  message("Reading NBaitsPerBin file...")
  header = readLines(s$nbaitsperbinfile, n=1)
  params = sapply(sapply(strsplit(header, "\t")[[1]],function(x)strsplit(x,"=")[[1]]), function(x)x[2])
  params = params[2:length(params)]
  names(params) = gsub("(\\S+)=.+", "\\1", names(params))

  maxl = as.numeric(params[["maxLBrownEst"]])
  if (maxl != s$maxLBrownEst){
    stop("The maxLBrownEst in the NBaitsPerBin file header is not equal to maxLBrownEst defined in experiment settings. Amend either setting (and if needed, generate a new NBaitsPerBin file) before running the analysis\n")
  }

  # Currently binsize is called bin, but should correct this
  binsz = as.numeric(params[["binsize"]]) 
  if (binsz != s$binsize){
    stop("The binsize in the NBaitsPerBin file header is not equal to binsize defined in experiment settings. Amend either setting (and if needed, generate a new NBaitsPerBin file) before running the analysis\n")
  }  
  
  # Currently not in the file
#   if (params[["removeb2b"]]!="True"){
#     stop("The NfragPerBin file must be generated with removeb2b==True\n")
#   }
#   if ( (params[["removeAdjacent"]]=="True" & !removeAdjacent) | (params[["removeAdjacent"]]!="True" & removeAdjacent)  ){
#     stop("The removeAdjacent parameter settings used for generating NfragPerBin file and defined here do not match. 
#          Amend either setting before running the analysis\n")
#   }  
  
  if(basename(params[["rmapfile"]]) != basename(s$rmapfile)){
    stop("The .rmap file used for generating the NBaitsPerBin file (according to the NBaitsPerBin header) and the one defined in experiment settings do not match. Amend either setting (and if needed, generate a new NBaitsPerBin file) before running the analysis\n")
  }
  
  ## Not checking this for now as we have a mixup of _baits and _baits_ID files used at different times...
  #   if(basename(params[["baitmapfile"]]) != basename(baitmapfile)){
  #     stop("Bait files used for generating the NfragPerBin file and defined here do not match. 
  #          Amend either setting before running the analysis\n")
  #   }
  
  nbpb = data.table::fread(s$nbaitsperbinfile, skip=1L)
  setnames(nbpb, names(nbpb)[1], "otherEndID")
  for(i in 2:ncol(nbpb)){
    setnames(nbpb, names(nbpb)[i], paste0("bin", i-1))    
  }
  nbpb
}


.readProxOEfile <- function(s){
  
  # Reads a pre-computed text file that denotes which other ends are in the proximal 
  # range relative to each bait, and gives that distance. 
  # Note that fragments that are too small/too large have already been removed.
  # s is the current chicagoData object's settings list
  
  message("Reading ProxOE file...")
  header = readLines(s$proxOEfile, n=1)
  params = sapply(sapply(strsplit(header, "\t")[[1]],function(x)strsplit(x,"=")[[1]]), function(x)x[2])
  params = params[2:length(params)]
  names(params) = gsub("(\\S+)=.+", "\\1", names(params))
  minsize = as.numeric(params[["minFragLen"]])
  if (minsize != s$minFragLen){
    stop("The minFragLen specified in the ProxOE file header is not equal to minFragLen defined in experiment settings. Amend either parameter setting (and if needed, generate a new ProxOE file) before running the analysis\n")
  }
  maxsize = as.numeric(params[["maxFragLen"]])
  if (maxsize != s$maxFragLen){
    stop("The maxFragLen specified in the ProxOE file header is not equal to maxFragLen defined in experiment settings. Amend either parameter setting (and if needed, generate a new ProxOE file) before running the analysis\n")
  }
  maxl = as.numeric(params[["maxLBrownEst"]])
  if (maxl != s$maxLBrownEst){
    stop("The maxLBrownEst specified in the ProxOE file header is not equal to maxLBrownEst defined in experiment settings. Amend either parameter setting (and if needed, generate a new ProxOE file) before running the analysis\n")
  }
  binsz = as.numeric(params[["binsize"]]) 
  if (binsz != s$binsize){
    stop("The binsize specified in the ProxOE file header is not equal to binsize defined in experiment settigs. Amend either parameter setting (and if needed, generate a new ProxOE file) before running the analysis\n")
  }  
  if (params[["removeb2b"]]!="True"){
    stop("The ProxOE file must be generated with removeb2b==True. Please generate a new file.\n")
  }
  if ( (params[["removeAdjacent"]]=="True" & !s$removeAdjacent) | (params[["removeAdjacent"]]!="True" & s$removeAdjacent)  ){
    stop("The removeAdjacent parameter settings used for generating ProxOE file (according to its header) and defined in experiment settings do not match. Amend either setting (and if needed, generate a new ProxOE file) before running the analysis\n")
  }  
  if(basename(params[["rmapfile"]]) != basename(s$rmapfile)){
    stop("The .rmap files used for generating the ProxOE file (according to the ProxOE header) and the one defined in experiment settings do not match. Amend either setting (and if needed, generate a new ProxOE file) before running the analysis\n")
  }
  ## Not checking this for now as we have a mixup of _baits and _baits_ID files used at different times...
  #   if(basename(params[["baitmapfile"]]) != basename(baitmapfile)){
  #     stop("Bait files used for generating the ProxOE file and defined here do not match. 
  #          Amend either setting before running the analysis\n")
  #   }
  proxOE = data.table::fread(s$proxOEfile, skip=1L)
  data.table::setnames(proxOE, 1:3, c("baitID", "otherEndID", "dist"))
  proxOE
  }

#' Export to linkSet format
#' @keywords internal
#' @importFrom rlang .data
.exportToLinkSet <- function(cd, scoreCol="score", cutoff=0, b2bcutoff=NULL,
                       order=c("position", "score")[1], removeMT=TRUE)
{
  if (any(c("rChr", "rStart", "rEnd", "rID", "bChr", "bStart", "bEnd", "bID") %in% colnames(cd@x))){
    stop ("Colnames x shouldn't contain rChr, rStart, rEnd, rID, bChr, bStart, bEnd, bSign, bID\n") 
  }

  if (! order %in% c("position","score")){
    stop ("Order must be either position (default) or score\n")
  }
  if (! removeMT %in% c(TRUE,FALSE)){
    stop("removeMT must be TRUE or FALSE")
  }
  
  message("Reading the restriction map file...")
  rmap = .readRmap(cd@settings)
  data.table::setnames(rmap, "V1", "rChr")
  data.table::setnames(rmap, "V2", "rStart")
  data.table::setnames(rmap, "V3", "rEnd")
  data.table::setnames(rmap, "V4", "otherEndID")
  
  message("Reading the bait map file...")
  baitmap = .readBaitmap(cd@settings)
  
  data.table::setnames(baitmap, "V1", "baitChr")
  data.table::setnames(baitmap, "V2", "baitStart")
  data.table::setnames(baitmap, "V3", "baitEnd")
  data.table::setnames(baitmap, cd@settings$baitmapFragIDcol, "baitID")
  data.table::setnames(baitmap, cd@settings$baitmapGeneIDcol, "promID")
  
  message("Preparing the output table...")
  if (is.null(b2bcutoff)){
    x = cd@x[ cd@x[[scoreCol]]>=cutoff, ]
  }
  else{
    x = cd@x[ (cd@x$isBait2bait==TRUE & cd@x[[scoreCol]]>=b2bcutoff ) | 
                ( cd@x$isBait2bait==FALSE & cd@x[[scoreCol]]>=cutoff ), ]
  }
  
  x = x[, c("baitID", "otherEndID", "N", scoreCol,"distSign","isBait2bait"), with=FALSE]
  
  data.table::setkey(x, "otherEndID")
  data.table::setkey(rmap, "otherEndID")
  
  x = merge(x, rmap, by="otherEndID", allow.cartesian = TRUE)
  data.table::setkey(x, "baitID")
  
  data.table::setkey(baitmap, "baitID")  
  x = merge(x, baitmap, by="baitID", allow.cartesian = TRUE)
  
  # note that baitmapGeneIDcol has been renamed into "promID" above 
  bm2 = baitmap[,c ("baitID", "promID"), with=FALSE]
  
  data.table::setDF(x)
  data.table::setDF(bm2)
  
  # this way we can be sure that the new column will be called promID.y  
  out = merge(x, bm2, by.x="otherEndID", by.y="baitID", all.x=TRUE, all.y=FALSE, sort=FALSE)
  out[is.na(out$promID.y), "promID.y"] = "."
  
  out = out[,c("baitChr", "baitStart", "baitEnd", "baitID", "promID.x", "rChr", "rStart", "rEnd", "otherEndID", scoreCol, "N", "promID.y","isBait2bait","distSign")]
  
  names(out) = c("bait_chr", "bait_start", "bait_end", "bait_ID", "bait_name", "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_ID", "score", "N_reads", "otherEnd_name","isBait2bait","distSign")
  
  out$N_reads [ is.na(out$N_reads) ] = 0
  out$score = round(out$score,2)
  
  if (order=="position"){
    out = out[order(out$bait_chr, out$bait_start, out$bait_end, out$otherEnd_chr, out$otherEnd_start, out$otherEnd_end), ]
  }
  if (order=="score"){
    out = out[order(out$score, decreasing=TRUE), ]
  }
  
  if(removeMT)
  {
    ##Remove mitochondrial DNA
    selMT <- tolower(out$bait_chr) == c("chrmt")
    if(any(selMT))
    {
      out <- out[!selMT,]
    }
  }
  
  #out
  
  ##convert out to a GI
  anchor.one = with(out, GenomicRanges::GRanges(as.character(bait_chr), IRanges::IRanges(start=bait_start, end=bait_end)))
  anchor.two = with(out, GenomicRanges::GRanges(as.character(otherEnd_chr), IRanges::IRanges(start=otherEnd_start, end=otherEnd_end)))
  linkSet(anchor.one, anchor.two, specificCol = out$bait_name,
                      counts=out$N_reads, baitName=out$bait_name, otherEndName=out$otherEnd_name,
                      baitID = out$bait_ID,oeID = out$otherEnd_ID,
                      score= out$score,distSign = out$distSign)
}


#######################################################
#' @rdname as.GInteractions
#' @export
setMethod("as.GInteractions", "linkSet", function(x) {
  anchor.one = regionsBait(x)
  anchor.two = oe(x)
  metadata = mcols(x)
  gi = InteractionSet::GInteractions(anchor.one,anchor.two)
  mcols(gi) = metadata
  return(gi)
})


#' @rdname exportInterBed
#' @importFrom utils write.table
#' @export
setMethod("exportInterBed", "linkSet", function(x, outfile) {
  gr1 <- as.data.frame(regionsBait(x))
  colnames(gr1) <- paste0("bait_",colnames(gr1))
  gr2 <- as.data.frame(oe(x))
  colnames(gr2) <- paste0("otherEnd",colnames(gr2))
  gr1 = gr1[,1:3]
  gr2 = gr2[,1:3]
  gr1Name = bait(x)
  gr2Name = paste0(oe(x))

  metaDf <- as.data.frame(mcols(x))

  out <- cbind(gr1, bait_name = gr1Name, gr2, otherEnd_name = gr2Name, metaDf)

  utils::write.table(out, outfile, sep="\t", quote=FALSE, row.names=FALSE)
})

#' @rdname exportWashU
#' @importFrom utils write.table
#' @export
setMethod("exportWashU", "linkSet", function(x, outfile) {
  gr1 <- as.data.frame(regionsBait(x))
  colnames(gr1) <- paste0("bait_",colnames(gr1))
  gr2 <- as.data.frame(oe(x))
  colnames(gr2) <- paste0("otherEnd",colnames(gr2))
  gr1 = gr1[,1:3]
  gr2 = gr2[,1:3]
  gr1Name = bait(x)
  gr2Name = paste0(oe(x))

  metaDf <- as.data.frame(mcols(x))
  if ("score" %in% colnames(metaDf)) {
    out = cbind(gr1, bait_name = gr1Name, gr2, otherEnd_name = gr2Name, metaDf$score)
  } else {
    out = cbind(gr1, bait_name = gr1Name, gr2, otherEnd_name = gr2Name)
  }

  utils::write.table(out, outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
})