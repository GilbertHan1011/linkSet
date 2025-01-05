# Annotation--------------
# Create a new environment to store database connections
.dbCache <- new.env(parent = emptyenv())

# Helper function to get or create database connection
.getDBConnection <- function(genome) {
  # Create a key for the connection
  conn_key <- paste0("conn_", genome)

  # Check if connection exists and is valid
  if (exists(conn_key, envir = .dbCache)) {
    src <- get(conn_key, envir = .dbCache)
    if (DBI::dbIsValid(src$con)) {
      return(src)
    } else {
      # Remove invalid connection
      rm(list = conn_key, envir = .dbCache)
    }
  }

  # Create new connection
  src <- if (genome == "hg38") {
    if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE)) {
    stop("Package TxDb.Hsapiens.UCSC.hg38.knownGene needed for this function")
  }
    Organism.dplyr::src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
  } else if (genome == "mm10") {
    if (!requireNamespace("TxDb.Mmusculus.UCSC.mm10.knownGene", quietly = TRUE)) {
    stop("Package TxDb.Mmusculus.UCSC.mm10.knownGene needed for this function")
  }
    Organism.dplyr::src_organism("TxDb.Mmusculus.UCSC.mm10.knownGene")
  } else {
    stop("Unsupported genome. Please use 'hg38' or 'mm10'.")
  }

  # Store the connection
  assign(conn_key, src, envir = .dbCache)

  return(src)
}

# Function to cleanup connections
.cleanupConnections <- function() {
  for (conn_key in ls(.dbCache)) {
    src <- get(conn_key, envir = .dbCache)
    if (inherits(src, "src_dbi") && DBI::dbIsValid(src$con)) {
      tryCatch({
        DBI::dbDisconnect(src$con)
      }, error = function(e) {
        warning(paste("Failed to disconnect from database:", e$message))
      })
    }
  }
  rm(list = ls(.dbCache), envir = .dbCache)
}


#' Annotate the link set with txDb. Give a gene list, and return a
#'
#' @aliases annotatePromoter
#' @param x linkSet
#' @param genome the genome you want to annotate
#' @param keyType the key type. You can check with AnnotationDbi::keytypes
#' @param upstream The upstream base from the gene
#' @param overwrite Whether to overwrite the regionsBait if it already exists
#'
#' @return linkSet object
#' @export
#'
#' @examples
#'   gr1 <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
#'                 ranges = IRanges(start = c(1000, 2000, 3000), width = 100),
#'                 strand = "+", symbol = c("BRCA1", "TP53", "NONEXISTENT"))
#'   gr2 <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
#'                 ranges = IRanges(start = c(5000, 6000, 7000), width = 100),
#'                 strand = "+")
#'   ls <- linkSet(gr1, gr2, specificCol = "symbol")
#'
#'   # Test annotatePromoter
#'   annotated_ls <- suppressWarnings(annotatePromoter(ls, genome = "hg38", upstream = 500,overwrite = TRUE))
#'
#'
# Modified annotatePromoter method
setMethod("annotatePromoter", "linkSet", function(x, genome = "hg38",
                                                  keyType = "symbol", upstream = 5000,
                                                  overwrite = FALSE) {
  if (!is.null(regionsBait(x)) && !overwrite) {
    warning("regionsBait already exists, set overwrite = TRUE to overwrite it")
    return(x)
  }
  if (!is.null(regionsBait(x)) && overwrite) {
    warning("regionsBait is being overwritten")
  }

  tryCatch({
    # Get cached or new connection
    src <- .getDBConnection(genome)

    genes <- bait(x)
    geneGr <- Organism.dplyr::genes(src, filter = ~(symbol %in% genes))
    promoterGr <- IRanges::promoters(geneGr, upstream = upstream)

    index <- match(genes, geneGr$symbol)
    if (any(is.na(index))) {
      warning("Some genes are not found in the txDb, they will be set to chrNULL")
      newIndex <- which(!is.na(index))
      grMatch <- promoterGr[index[newIndex]]
      gr <- GRanges(seqnames = rep("chrNULL", length(genes)),
                    ranges = IRanges(rep(0, length(genes)), rep(0, length(genes))))
      mcols(grMatch) <- NULL
      gr[newIndex] <- grMatch
    } else {
      gr <- promoterGr[index]
    }

    regionsBait(x) <- gr
    return(x)
  }, error = function(e) {
    warning(paste("An error occurred:", e$message))
    return(x)
  })
})

# Register cleanup function to be called when R exits
reg.finalizer(.dbCache, function(e) {
  .cleanupConnections()
}, onexit = TRUE)

#' Execute Database Operation with Automatic Connection Management
#'
#' @param x The genome name or object to operate on
#' @param expr Expression to evaluate with database connection
#' @param ... Additional arguments
#' @return Result of the database operation
#' @export
#' @examples
#' requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene")
#' requireNamespace("org.Hs.eg.db")
#' result <- withTxDb("hg38", function(src) {
#'   genes <- Organism.dplyr::genes(src)
#'   return(genes)
#' })
#' print(result)
setMethod("withTxDb", signature(x = "character", expr = "function"),
  function(x, expr, ...) {
    if (!x %in% c("mm10", "hg38")) {
      stop("Unsupported genome. Please use 'hg38' or 'mm10'.")
    }
    
    src <- .getDBConnection(x)
    
    tryCatch({
      result <- expr(src, ...)
      return(result)
    }, error = function(e) {
      stop(paste("Database operation failed:", e$message))
    })
  }
)
