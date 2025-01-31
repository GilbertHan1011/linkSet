#' @rdname run_chicane
#' @param replicate.merging.method Method for merging replicates (default: 'sum')
#' @param distribution Distribution to use for modeling (default: 'negative-binomial')
#' @param include.zeros How to handle zero counts (default: 'none')
#' @param bait.filters Vector of length 2 for bait filtering thresholds (default: c(0,1))
#' @param target.filters Vector of length 2 for target filtering thresholds (default: c(0,1))
#' @param distance.bins Number of distance bins (default: NULL)
#' @param multiple.testing.correction Method for multiple testing correction (default: 'bait-level')
#' @param adjustment.terms Additional terms for model adjustment (default: NULL)
#' @param remove.adjacent Whether to remove adjacent fragments (default: FALSE)
#' @param temp.directory Directory for temporary files (default: NULL)
#' @param keep.files Whether to keep temporary files (default: FALSE)
#' @param maxit Maximum iterations for model fitting (default: 100)
#' @param epsilon Convergence threshold (default: 1e-8)
#' @param cores Number of CPU cores to use (default: 1)
#' @param trace Whether to print trace information (default: FALSE)
#' @param verbose Whether to print progress information (default: FALSE)
#' 
#' @examples
#' # Create example data
#' gr1 <- GRanges(seqnames = c("chr1", "chr3", "chr3"),
#'                ranges = IRanges(start = c(1000, 2000, 3000), width = 100),
#'                strand = "+", symbol = c("BRCA1", "TP53", "NONEXISTENT"))
#' gr2 <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
#'                ranges = IRanges(start = c(5000, 6000, 7000), width = 100),
#'                strand = "+")
#' ls <- linkSet(gr1, gr2, specificCol = "symbol")
#' 
#' # Annotate and prepare data
#' annotated_ls <- suppressWarnings(
#'   annotatePromoter(ls, genome = "hg38", upstream = 500, overwrite = TRUE)
#' )
#' annotated_ls <- countInteractibility(annotated_ls)
#' annotated_ls <- linkSet::pairdist(annotated_ls)
#' 
#' # Run analysis
#' result_ls <- run_chicane(
#'   annotated_ls, 
#'   replicate.merging.method = 'sum',
#'   bait.filters = c(0, 1),
#'   target.filters = c(0, 1),
#'   distance.bins = NULL,
#'   multiple.testing.correction = 'bait-level',
#'   verbose = TRUE
#' )
#' 
#' @import data.table
#' @export
setMethod("run_chicane", "linkSet", function(linkSet, 
	replicate.merging.method = 'sum',
	distribution = 'negative-binomial',
	include.zeros = 'none',
	bait.filters = c(0, 1),
	target.filters = c(0, 1),
	distance.bins = NULL,
	multiple.testing.correction = c('bait-level', 'global'),
	adjustment.terms = NULL,
	remove.adjacent = FALSE,
	temp.directory = NULL,
	keep.files = FALSE,
	maxit = 100,
	epsilon = 1e-8,
	cores = 1,
	trace = FALSE,
	verbose = FALSE
	) {
	# TO DO:
	#	- check format of linkSet object if passed directly
	linkSet = .verify.linkSet(linkSet)
	.check.packages(distribution)
	### MAIN #############################################################
	replicate.merging.method <- match.arg(replicate.merging.method);
	multiple.testing.correction <- match.arg(multiple.testing.correction);

	# Convert linkSet to data.table for processing
	#interaction.data <- as.data.table(linkSet)
	linkSet <- .filter.fragments(
		linkSet,
		bait.filters = bait.filters,
		target.filters = target.filters,
		verbose = verbose
		);
	
	if( verbose ) {
		cat('FITTING MODEL\n');
	}
	
	chicane.results <- fit.model(
		linkSet, 
		distance.bins = distance.bins, 
		distribution = distribution,
		adjustment.terms = adjustment.terms,
		verbose = verbose,
		cores = cores,
		maxit = maxit,
		epsilon = epsilon,
		trace = trace
		);

	chicane.results <- multiple.testing.correct(
		chicane.results,
		bait.level = 'bait-level' == multiple.testing.correction
		);
	#print(dim(chicane.results))
	# sort by q-value and p-value
	chicane.results <- chicane.results[order(chicane.results$q.value, chicane.results$p.value), ];
	#print(dim(chicane.results))
	# Convert results back to linkSet
	result_linkSet <- linkSet
  #print(head(chicane.results))
	mcols(result_linkSet)$expected <- chicane.results$expected
	mcols(result_linkSet)$p.value <- chicane.results$p.value
	mcols(result_linkSet)$q.value <- chicane.results$q.value
	return(result_linkSet);
})


#' verify.linkSet
#'
#' @description 
#'	Verify that linkSet object is in expected format. Throws an error if object does not fit requirements.
#'
#' @param linkSet Object to be verified.
#' @importFrom methods is
#' @importFrom S4Vectors mcols mcols<-
#' @keywords internal
#' @rdname chicane
#' @return None
#'
.verify.linkSet <- function(linkSet) {
	# LinkSet object
	if( !is(linkSet, "linkSet") ) {
		stop('linkSet must be a LinkSet object');
	}

	required.columns <- c('distance', 'count', 'bait.trans.count', 'target.trans.count');

	if( !all(required.columns %in% colnames(mcols(linkSet))) ) {
		error.message <- paste(
			'linkSet must contain the columns:\n',
			paste(required.columns, collapse = ' '),
			"\nPlease check wether you have run the function countInteractibility, pairdist and countInteractions before."
			);
		stop(error.message);
	}
	if (!"bait.id" %in% colnames(mcols(linkSet))) {
		mcols(linkSet)$bait.id = bait(linkSet)
		print("Not found column 'bait.id', adding bait name as default.")
	}
	if (!"bait.to.bait" %in% colnames(mcols(linkSet))) {
		mcols(linkSet)$bait.to.bait = FALSE
		print("Not found column 'bait.to.bait' Set 'FALSE' as default.")
	}
	if (any(is.na(linkSet$bait.id))){
		print("Found NA in 'bait.id',Filtering...")
		linkSet = linkSet[!is.na(linkSet$bait.id)]
	}
	return(linkSet)
}

#' filter.fragments
#'
#' @description 
#'	Filter low and high-interacting restriction fragments based on the total number of trans counts
#'
#' @param linkSet LinkSet object containing interactions
#'
#' @return LinkSet object containing fragments that passed all filters
#' @rdname chicane
#' @keywords internal
#' @noRd
#' @import linkSet
#' @export

.filter.fragments <- function(
	linkSet,
	bait.filters = c(0, 1),
	target.filters = c(0, 1),
	verbose = FALSE
	) {

	### MAIN ##################################################################	
	#print("filtering fragments....")
	apply.target.filters <- !identical( c(0, 1), target.filters);
	apply.bait.filters <- !identical( c(0, 1), target.filters);

	if( verbose && (apply.bait.filters || apply.target.filters) ) {
		cat('Applying filters\n');
		cat('\tbaits:', paste0(round(100*bait.filters[1], 2), '% - ', round(100*bait.filters[2], 2), '%'), '\n' );
		cat('\ttargets:', paste0(round(100*bait.filters[1], 2), '% - ', round(100*bait.filters[2], 2), '%'), '\n' );
	}

	filtered.ls <- linkSet;

	# filter targets
	if( apply.target.filters ) {

		# get counts per fragment
		target.counts <- unique( mcols(linkSet)[, c("target.id", "target.trans.count")] );
	
		# get absolute thresholds corresponding to percentiles
		target.lower.cutoff <- stats::quantile(
			target.counts$target.trans.count,
			prob = target.filters[1]
			);

		target.upper.cutoff <- stats::quantile(
			target.counts$target.trans.count,
			prob = target.filters[2]
			);

		filtered.ls <- filtered.ls[ filtered.ls$target.trans.count >= target.lower.cutoff & filtered.ls$target.trans.count <= target.upper.cutoff ];

	}
	# filter baits
	if( !apply.bait.filters ) {

		bait.counts <- unique( mcols(linkSet)[, c("bait.id", "bait.trans.count")] );

		# get absolute thresholds corresponding to percentiles
		bait.lower.cutoff <- stats::quantile(
			bait.counts$bait.trans.count,
			prob = bait.filters[1]
			);

		bait.upper.cutoff <- stats::quantile(
			bait.counts$bait.trans.count,
			prob = bait.filters[2]
			);

		filtered.ls <- filtered.ls[ (filtered.ls$bait.trans.count >= bait.lower.cutoff & filtered.ls$bait.trans.count) <= bait.upper.cutoff ,];
	}

	if( verbose && (apply.bait.filters || apply.target.filters) ) {
		cat('\ttargets:', nrow(filtered.ls), 'interactions remain after filtering');
	}

	return(filtered.ls);
}

#' .model.try.catch
#'
#' @description
#' 	Internal function for fitting model within a tryCatch loop, handling numerical errors gracefully.
#'
#' @param model.formula formula
#' @param data model data
#' @param init.theta Initial value of theta in negative binomial model
#' @param start starting values of coefficients in linear predictor
#' 
#' @return List with elements
#'  \item{model}{model object. Set to NULL if no model could be fit.}
#' 	\item{expected.values}{vector of expected values for each element in original data, or vector of NAs if no model could be fit}
#' 	\item{p.values}{vector of p-values for test of significantly higher response than expected, or vector of NAs if no model could be fit}
#' @keywords internal
#' @noRd
#' @rdname chicane
#' @export
.model.try.catch <- function(
	model.formula, 
	data,
	distribution = 'negative-binomial',
	maxit = 100,
	epsilon = 1e-8,
	init.theta = NULL, 
	start = NULL,
	trace = FALSE,
	verbose = FALSE
	) {

	# Negative binomial GLM only supports overdispersion. 
	# When the variance does not exceed the mean, this causes problems with model fitting/
	# These problems manifest in one of two ways:
	#	1) An error in while ((it <- it + 1) < limit && abs(del) > eps)
	#	2) A warning about NaNs produced in sqrt(1/i)
	#
	# Handle these cases by trying to fit a Poisson distribution instead.
	model <- tryCatch({
		# try fitting distribution as requested by user
		fit.glm( 
			model.formula, 
			data,
			distribution = distribution,
			maxit = maxit,
			epsilon = epsilon,
			trace = trace,
			init.theta = init.theta,
			start = start
 			); 
		}, error = function(e) {

			# if problem was negative binomial, try Poisson
			if( 'negative-binomial' == distribution && .is.glm.nb.theta.error(e) ) {

				if(verbose) cat('\nDispersion error - running Poisson\n');
				
				temp.model <- fit.glm( 
					model.formula, 
					data,
					distribution = 'poisson',
					maxit = maxit,
					epsilon = epsilon,
					trace = trace
 					); 
				return(temp.model);
				
			} else if ( grepl('no valid set of coefficients has been found: please supply starting values', e$message, fixed = TRUE) ) {
				# enter error as couldn't fit model	
				if(verbose) cat('\nNo valid coefficients error - skipping ahead\n');
				model.data <- list(
					model = NULL,
					expected.values = rep(NA, nrow(data)),
					p.values = rep(NA, nrow(data))
					);
				return( model.data );

			} else if( grepl("NA/NaN/Inf in 'x'", e$message, fixed = TRUE ) ) {
				if(verbose) cat('\nNA/NaN/Inf in x error - skipping ahead\n');
				model.data <- list(
					model = NULL,
					expected.values = rep(NA, nrow(data)),
					p.values = rep(NA, nrow(data))
					);
				return( model.data );
			} else {
				if(verbose) {
					cat('\nUnknown error - skipping ahead\n');
					cat(e$message, '\n');
				}

				model.data <- list(
					model = NULL,
					expected.values = rep(NA, nrow(data)),
					p.values = rep(NA, nrow(data))
					);
			}
		}, warning = function(w) {
			dispersion.problem <- FALSE;
			if( 'negative-binomial' == distribution && ( .is.glm.nb.maxiter.warning(w) || .is.glm.nb.theta.warning(w) ) ) {
				if(verbose) cat('Caught a warning - checking for dispersion problems\n');

				# See if problem is lack of overdispersion
				# Get estimate of theta after a low number of iterations
				# 	=> if no evidence for overdispersion, fit Poisson

				negbin.fit <- suppressWarnings(
					fit.glm( 
						model.formula, 
						data,
						distribution = 'negative-binomial',
						maxit = 25,
						epsilon = epsilon,
						trace = trace,
						init.theta = init.theta,
						start = start
 						)
					); 
				max.mu <- max(negbin.fit$model$fitted.values);
				var.at.max.mu <- max.mu + (max.mu^2)/negbin.fit$model$theta;

				# assess difference between NB variance and Poisson one
				if( (var.at.max.mu - max.mu)/max.mu < 0.001 ) {
					dispersion.problem <- TRUE;
				}
			} 

			# if there is a dispersion problem, try Poisson
			# if not, raise original warning and fit original distribution
			if( dispersion.problem ) {

				if(verbose) cat('Dispersion problem detected - fitting Poisson\n');
				distribution <- 'poisson';
			} else {
				warning(w);
			}

			temp.model <- tryCatch({
				fit.glm( 
					model.formula, 
					data,
					distribution = distribution,
					maxit = maxit,
					epsilon = epsilon,
					trace = trace
 					); 
				}, error = function(e) {
					if(verbose) cat('\nUnknown error - skipping ahead\n');
					list(
						model = NULL,
						expected.values = rep(NA, nrow(data)),
						p.values = rep(NA, nrow(data))
						);
				});

			return(temp.model);
		});
	return(model);
}


#' Fit Statistical Model for Interaction Analysis
#' 
#' @title Fit Model for Interaction Analysis
#' @description Fit negative binomial model to obtain p-values for interactions.
#' 
#' @param linkSet A linkSet object
#' @param distance.bins Distance bins for model fitting (default: NULL)
#' @param distribution Type of distribution to use (default: 'negative-binomial')
#' @param adjustment.terms Character vector of extra terms to adjust for in the model fit
#' @param maxit Maximum number of iterations (default: 100)
#' @param epsilon Convergence threshold (default: 1e-8)
#' @param cores Number of cores to use (default: 1)
#' @param trace Whether to print trace information (default: FALSE)
#' @param verbose Whether to print progress reports (default: FALSE)
#' 
#' @return A data.table with fitted model results including expected counts and p-values
#' @keywords internal
#' @noRd
fit.model <- function(
	linkSet, 
	distance.bins = NULL,
	distribution = 'negative-binomial', 
	adjustment.terms = NULL,
	maxit = 100,
	epsilon = 1e-8,
	cores = 1,
	trace = FALSE,
	verbose = FALSE
	) {
	### MAIN ##################################################################
  print("fitting model....")
	# filter out low and high-interacting fragments


	# split into bait to bait and other
	b2b.data <- linkSet[ linkSet$bait.to.bait ,];
	non.b2b.data <- linkSet[ !linkSet$bait.to.bait ,];
	b2b.data <- as.data.table(mcols(b2b.data))
	non.b2b.data <- as.data.table(mcols(non.b2b.data))
	# free up memory
	rm(linkSet);
	#gc();
	#browser()
	# fit models separately
	b2b.results <- NULL;
	non.b2b.results <- NULL;
	#before_fitting_time <- Sys.time()
	#print(paste0("before_fitting_time: ", before_fitting_time))
	if( nrow(b2b.data) > 0 ) {
		b2b.results <- run.model.fitting(
			b2b.data, 
			distance.bins = distance.bins, 
			distribution = distribution, 
			verbose = verbose,
			bait.to.bait = TRUE,
			adjustment.terms = adjustment.terms,
			cores = cores,
			maxit = maxit,
			epsilon = epsilon,
			trace = trace
			);
	}

	if( nrow(non.b2b.data) > 0 ) {
		non.b2b.results <- run.model.fitting(
			non.b2b.data, 
			distance.bins = distance.bins, 
			distribution = distribution,
			verbose = verbose,
			bait.to.bait = FALSE,
			adjustment.terms = adjustment.terms,
			cores = cores,
			maxit = maxit, 
			epsilon = epsilon,
			trace = trace
			);
	}

	# combine results
	# hold off on sorting to optimize for speed
	#  - everything will get jumbled for the multiple testing correction anyways
	combined.data <- rbind(b2b.results, non.b2b.results);

	return(combined.data);
}




#' run.model.fitting
#'
#' @description
#' 	Run model fitting procedure for either bait-to-bait or other interactions.
#'	Meant for internal use only. 
#'
#' @param bait.to.bait Logical indicating if model should be fit as bait-to-bait
#' @param adjustment.terms Characted vector of extra terms to adjust for in the model fit
#'
#' @rdname chicane
#' @noRd
#' @return Interactions data with expeceted number of interactions and p-values added.
#' @keywords internal
#' @importFrom foreach %dopar%
#' @importFrom iterators icount
#' @importFrom stats logLik
#' @importFrom rlang .data
#' @importFrom data.table :=
run.model.fitting <- function(
	interaction.data,
	distance.bins = NULL, 
	distribution = 'negative-binomial',
	bait.to.bait = FALSE,
	adjustment.terms = NULL,
	maxit = 100,
	epsilon = 1e-8,
	cores = 1,
	trace = FALSE,
	verbose = FALSE
	) {

	# TO DO:
	# 	- see if you can avoid cis/ trans repetitiveness
	
	### INPUT TESTS ###########################################################
	#interaction.data <- mcols(linkSet)
	# be extremely strict about these things to avoid bugs
	if( bait.to.bait && !all( interaction.data$bait.to.bait ) ) {
		stop('Cannot fit bait-to-bait model when not all interactions are bait-to-bait');
	}

	if( !bait.to.bait && any(interaction.data$bait.to.bait) ) {
		stop('Cannot fit non-bait-to-bait model on bait-to-bait interactions');
	}

	if( !is.null(adjustment.terms) && !is.character(adjustment.terms) ) {
		stop('adjustment.terms must be a character vector');
	}

	### MAIN ##################################################################

	# figure out formula to use based on whether it's bait-to-bait or not
	# need separate cis and trans formulas because of distance adjustment
	if( bait.to.bait ) {
		cis.formula <- stats::as.formula(count ~ log(distance) + log(bait.trans.count + 1)*log(target.trans.count + 1));
		trans.formula <- stats::as.formula(count ~ log(bait.trans.count + 1)*log(target.trans.count + 1));
	} else {
		cis.formula <- stats::as.formula(count ~ log(distance) + log(bait.trans.count + 1) );
		trans.formula <- stats::as.formula(count ~ log(bait.trans.count + 1) )
	}

	# if requested, update model with user-requested terms
	if( !is.null(adjustment.terms) ) {
		adjustment.string <- paste(adjustment.terms, collapse = ' + ');

		cis.formula <- stats::update.formula(cis.formula, paste0('~ . + ', adjustment.string) );
		trans.formula <- stats::update.formula(trans.formula, paste0('~ . + ', adjustment.string) );

		# graceful error handling – make sure all variables are in the input data
		# do this here in case user specifies something like log(x) in adjustment.terms
		formula.vars <- unique( c(all.vars(cis.formula), all.vars(trans.formula)) );
		if( !all(formula.vars %in% names(interaction.data)) ) {
			error.message <- paste(
				'The following variables were not found in the data:', 
				paste(formula.vars[ !(formula.vars %in% names(interaction.data) ) ], collapse = ' ')
				);
			stop(error.message);
		}
	}
	trans.data <- interaction.data[ is.na(interaction.data$distance) ,];

	# Fit models separately in each quantile of distance
	cis.data <- interaction.data[ !is.na(interaction.data$distance) ,];
	cis.data <- cis.data[ order(cis.data$distance), ];

	# free up memory
	rm(interaction.data);
	#gc();

	# list of data.tables, where each element corresponds to 
	# a specific distance
	distance.binned.data <- .distance.split(
		cis.data, 
		distance.bins = distance.bins, 
		verbose = verbose
		);

	# store interaction data after fitting models
	p.value.data <- list();

	# speed up model fitting by passing starting values and theta between iterations 
	# hopefully this will also help with stability ?
	init.theta <- NULL;
	start <- NULL;
	#before_parallel_time <- Sys.time()
	#print(paste0("before_parallel_time: ", before_parallel_time))
	if( cores > 1 ) {
		computing.cluster <- parallel::makeCluster( cores );
		doParallel::registerDoParallel( computing.cluster );
	} else {
		foreach::registerDoSEQ();
	}
	iter.i <- NULL;
	if (distribution == 'poisson'|| distribution == 'negative-binomial') {
		packages <- c('MASS', 'data.table')
	} else {
		packages <- c('MASS', 'data.table', 'gamlss', 'gamlss.tr')
	}
		p.value.data <- foreach::foreach(
		temp.data = distance.binned.data,
		iter.i = icount(),
		.packages = packages,
    .export = c(".model.try.catch", "create.modelfit.plot")
		) %dopar% {
		
		# progress meter
		if(verbose) cat('*');

		# fit model through helper function that gracefully handles numerical errors
		model <- .model.try.catch(
			cis.formula, 
			get("temp.data"),
			distribution = distribution,
			maxit = maxit,
			epsilon = epsilon,
			trace = trace,
			init.theta = init.theta,
			start = start
			);
    #browser()
		data.table::set(get("temp.data"), j = "expected", value = model$expected.values)
		data.table::set(get("temp.data"), j = "p.value", value = model$p.values)
		# clear memory
		#for (gc.i in 1:5) { gc(); }
		return(get("temp.data"));
	}

	# fit trans-interactions
	# (same model, but no distance correction)
	if( nrow(trans.data) > 0  ) {

		if(verbose) {
			cat('\n\ttrans interactions\n');
		}

		trans.model <- .model.try.catch(
			trans.formula, 
			trans.data,
			distribution = distribution,
			maxit = maxit,
			epsilon = epsilon,
			trace = trace,
			init.theta = init.theta,
			start = start
			);
		data.table::set(trans.data, j = "expected", value = trans.model$expected.values)
		data.table::set(trans.data, j = "p.value", value = trans.model$p.values)

		# add to p-value data frame
		p.value.data[[ length(p.value.data) + 1 ]] <- trans.data; 
	}

	# clean up parallel computing
	if ( cores > 1 ) {
		foreach::registerDoSEQ();
		parallel::stopCluster(computing.cluster);
		remove(computing.cluster);
	}

	p.value.data <- do.call(rbind, p.value.data);

	return(p.value.data);

}


#' .distance.split
#'
#' @description
#'	Split interaction data into subsets that are large enough for the chicane model to be fit (see Details), 
#'	based on distance. This step allows the distance term in the model to be fit in a piecewise linear fashion. 
#'
#' @details
#'  Fitting \code{glm.nb} fails when there is a lack of overdispersion in the data. The chicane method
#' 	contains logic to catch these errors and instead fit a Poisson model. However, to avoid this happening
#' 	more than necessary, an attempt is made to avoid distance splits that will clearly result in numerical errors.
#' 	This includes bins of data where the count is the same for all rows, 
#'	or a covariate is a perfect predictor of count.  
#' 	
#' @param interaction.data 
#'	Data table of interaction data, typically from \code{prepare.data}
#' @param distance.bins 
#'	Number of distance bins desired. If NULL, a number is chosen to ensure that the negative binomial can be fit in all bins.
#' @param min.rows.bin
#' 	The minimum number of expected rows in a distance bin. Ignored if distance.bins is set
#' @param verbose 
#'	Logical indicating whether to print progress reports
#' @rdname chicane
#' @return 
#'	List where each element corresponds to a specified distance bin, and the final one corresponding to trans-interactions (if present)
#' @keywords internal
#' @noRd
#' @export
.distance.split <- function(
	interaction.data, 
	distance.bins = NULL, 
	min.rows.bin = 50,
	verbose = FALSE
	) {

	### INPUT TESTS ###########################################################

	if( !is.data.table(interaction.data) ) {
		stop('interaction.data must be a data.table object');
	}

	if( !is.null(distance.bins) && !is.numeric(distance.bins) ) {
		stop('distance.bins must be a positive integer');
	} 

	if( !is.null(distance.bins) && ( 0 != distance.bins %% 1 || distance.bins < 1 ) ) {
		stop('distance.bins must be a positive integer');
	}

	### MAIN ##################################################################

	cis.data <- interaction.data[!is.na(interaction.data$distance),];
	trans.data <- interaction.data[is.na(interaction.data$distance),];

	cis.data <- cis.data[order(cis.data$distance),];

	if( nrow(cis.data) < 50 || !.check.model.numerical.fit(cis.data) ) {
		# if fewer than 50 rows in cis-data, or already too small for a decent numerical fit, keep as one item
		split.data <- list( cis.data );

	} else {
		if( is.null(distance.bins) ) {

			distance.bins <- min( round( nrow(cis.data)/min.rows.bin, 1), 100 );
			numerical.fit <- FALSE;

			if( verbose ) {
				cat('Splitting data for model fitting\n');
				cat('\tchecking distance.bins =', distance.bins, '\n');
			}
			

			# as long as model cannot be fit, try more distance.bins
			while( distance.bins >= 1 && !numerical.fit ) {

				distance.bins <- round(distance.bins/2);
				if( verbose ) cat('\tchecking distance.bins =', distance.bins, '\n');
				split.data <- .smart.split(cis.data, bins = distance.bins);
				# get indicator of whether the model can be fit in each of the split data parts
				numerical.fit <- .check.split.data.numerical.fit(split.data);
			}

		} else {
			# user has requested a specified number of distance bins
			# split into this number of groups, and throw an error if the model cannot be fit

			split.data <- .smart.split(cis.data, bins = distance.bins);
			numerical.fit <- .check.split.data.numerical.fit(split.data);
			if( !numerical.fit ) {
				stop('Model cannot be fit with the specified number of distance bins. Try using fewer bins.');
			}
		} # end distance bins if/else
	
	}

	# add data on trans interactions if it exists
	if( nrow(trans.data) > 0 ) {
		split.data[[ length(split.data) + 1 ]] <- trans.data;
	}
	
	return(split.data);

}


#' check.glm.nb.theta.error
#'
#' @description
#' 	Check if an error matches the error raised by \code{glm.nb} due to an inflated theta estimate.
#'	This happens when the variance of the negative binomial does not exceed the mean (i.e. there is no overdispersion).
#'  In such cases, the Poisson distribution may be a suitable alternative.
#' 
#' @param e Error object
#' @rdname chicane
#' @return Boolean indicating if error matches
#'
.is.glm.nb.theta.error <- function(e) {
	
	error.code <- deparse(e$call)[1]; # will be a character vector if it really is theta error
	error.message <- e$message;

	# run checks for non-compatibility with glm.nb error
	# testing for equality seems to be less robust than running 
	if( grepl('missing value where TRUE/FALSE needed', error.message, fixed = TRUE) ) {

		if( grepl('while ((it <- it + 1) < limit && abs(del) > eps)', error.code, fixed = TRUE) ) return(TRUE);

		if( grepl('while ((iter <- iter + 1) <= control$maxit && (abs(Lm0 - Lm)/d1', error.code, fixed = TRUE) ) return(TRUE);

		if( grepl('if (t0 < 0)', error.code, fixed = TRUE) ) return(TRUE);
	}

	return(FALSE);
}


#' check.model.identifiability
#'
#' @description
#' 	Check if chicane model can be fit on a given dataset. 
#'	\code{glm.nb} does not work when all responses are constant, or there are only two unique values and a covariate is a perfect predictor.
#'
#' @param interaction.data Data table of interaction data on which model is to be fit
#' @rdname chicane
#' @return boolean indicating if model can be fit
#' 
.check.model.numerical.fit <- function(interaction.data) {

	### INPUT TESTS ###########################################################

	if( !is.data.table(interaction.data) ) {
		stop('interaction.data must be a data.table object');
	}

	### MAIN ##################################################################

	is.trans <- all( is.na(interaction.data$distance ));
	is.b2b <- all( interaction.data$bait.to.bait );

	# figure out terms to include in model
	model.terms <- c('bait.trans.count');
	
	if( !is.trans ) {
		model.terms <- c( model.terms, 'distance' );
	}

	if( is.b2b ) {
		model.terms <- c( model.terms, 'target.trans.count' );
	}

	unique.counts <- length(unique( interaction.data$count) );

	if( 1 ==  unique.counts) {
		return(FALSE);
	}

	if( unique.counts < 3 ) {
		count.levels <- as.numeric( as.factor( interaction.data$count ) );
		
		for( model.term in model.terms ) {
			model.levels <- as.numeric( as.factor( interaction.data[[ model.term ]] ) );
			
			if( identical(count.levels, model.levels) ) {
				return(FALSE);
			}
		}
	}

	return(TRUE);
}

#' multiple.testing.correct
#'
#' @description
#'	Perform multiple testing correction on p-values from interaction test.
#' 	By default, multiple testing correction is applied per bait. To change this
#' 	to a global multiple testing correction, set \code{bait.level = FALSE}.
#'
#' @param interaction.data 
#' 	Data table of interaction calls. Must contain columns p.value and bait.id.
#' @param bait.level 
#'	Logical indicating whether multiple testing correction should be performed per bait.
#' @rdname chicane
#' @return Original data table with new column
#' 	\item{q.value}{FDR-corrected p-value}

multiple.testing.correct <- function(
	interaction.data,
	bait.level = TRUE
	) {

	### INPUT TESTS ###########################################################

	.verify.interaction.data(interaction.data);

	if( !all(c('bait.id', 'p.value') %in% names(interaction.data)) ) {
		stop('interaction.data must contain columns bait.id and p.value');
	}	

	### MAIN ##################################################################


	if( bait.level ) {
		# split input data by locus, and perform multiple testing correction for 
		# each locus separately
		locus.data <- split(
			interaction.data, 
			interaction.data$bait.id
			);
		
		q.value.data <- list();

		for(i in seq_along(locus.data) ) {

			temp.data <- locus.data[[i]];

			temp.data$q.value <- stats::p.adjust(
				temp.data$p.value,
				method = 'fdr'
				);
			temp.data$q.value[is.na(temp.data$q.value)] <- 1
			q.value.data[[ i ]] <- temp.data;
		}
		q.value.data <- do.call(rbind, q.value.data);
	} else {
		# perform global FDR correction
		q.value.data <- interaction.data;
		q.value.data$q.value <- stats::p.adjust(
			q.value.data$p.value,
			method = 'fdr'
			);
	}
	return(q.value.data);
}

#' smart.split
#'
#' @description
#'	Split a data frame into a prespecified number of bins, using 
#'	\code{split} and \code{cut}. Unlike the default R functions, this does not
#'	fail when asked to split the data into a single bin.
#'
#' @param dat Data frame or data table to be split
#' @param bins Number of bins to split data into
#' @rdname chicane
#' @return 
#'	List with \code{bins} elements. Each element corresponds to one portion 
#'	of the data  
#'
.smart.split <- function(dat, bins) {

	### INPUT TESTS ###########################################################

	if( !is.data.frame(dat) ) {
		stop('dat must be a data frame');
	}

	if( !is.numeric(bins) || bins %% 1 != 0 || bins < 1 ) {
		stop('bins must be a positive integer');
	}

	### MAIN ##################################################################

	if( bins > 1 ) {
		split.data <- split(
			dat,
			cut( seq_len( nrow(dat) ), breaks = bins)
			);
	} else {
		# only a single bin - cut function will fail
		# simply return a list with one element containing all of the data
		split.data <- list( dat );
	}

	return(split.data);

}

#' check.split.data.numerical.fit
#'
#' @description
#'	Helper function to check if the chicane model can be fit on each element of a split data list.
#'
#' @param split.data List of data.table objects with fragment interaction data
#' 
#' @return Logical indicating if the model can be fit
#'
.check.split.data.numerical.fit <- function(split.data) {

	### INPUT TESTS ###########################################################

	if( !is.list(split.data) ) {
		stop('split.data must be a list');
	}

	element.is.data.table <- vapply(split.data, is.data.table, FUN.VALUE = FALSE);
	if( !all(element.is.data.table) ) {
		stop('Each element of split.data should be a data.table object');
	}

	### MAIN ##################################################################

	element.numerical.fit <- vapply(
		split.data,
		.check.model.numerical.fit,
		FUN.VALUE = FALSE
		);

	return( all(element.numerical.fit) );
}


#' .verify.interaction.data
#'
#' @description 
#'	Verify that interaction.data object is in expected format. Throws an error if object does not fit requirements.
#'
#' @param interaction.data Object to be verified.
#' @rdname chicane
#' @return None
#'
.verify.interaction.data <- function(interaction.data) {

	# data.table object
	if( !is.data.table(interaction.data) ) {
		stop('interaction.data must be a data.table object');
	}

	required.columns <- c('distance', 'count', 'bait.trans.count', 'target.trans.count');

	if( !all(required.columns %in% names(interaction.data)) ) {
		error.message <- paste(
			'interaction.data must contain the columns:\n',
			paste(required.columns, collapse = ' ')
			);
		stop(error.message);
	}
}

#' fit.glm
#'
#' @description
#'	Fit GLM according to a specified distribution. This needs to be done separately from \code{glm}
#' 	in order to include negative binomial and truncated distributions as options.
#'
#' @param formula 
#' 	Formula specifying model of interest
#' @param data 
#'	Data frame containing variables specified in formula
#' @param distribution 
#' 	Name of distribution of the counts. Options are 'negative-binomial', 
#'	'poisson', 'truncated-poisson', and 'truncated-negative-binomial'
#' @param start 
#' 	Starting values for model coefficients
#' @param init.theta
#' 	Initial value of theta if fitting the negative binomial distribution
#' @param maxit 
#'	Maximum number of IWLS iterations for fitting the model (passed to \code{glm.control})
#' @param epsilon
#'	Positive convergence tolerance for Poisson and negative binomial models. Passed to \code{glm.control}
#' @param trace
#' 	Logical indicating if output should be produced for each of model fitting procedure. Passed to \code{glm.control} or \code{gamlss.control}
#' @keywords internal
#' @rdname chicane
#' @return List with elements
#'  \item{model}{model object}
#' 	\item{expected.values}{vector of expected values for each element in original data}
#' 	\item{p.values}{vector of p-values for test of significantly higher response than expected}
#' @importFrom utils getFromNamespace
fit.glm <- function(
	formula, 
	data, 
	distribution = c('negative-binomial', 'poisson', 'truncated-poisson', 'truncated-negative-binomial'),
	start = NULL,
	init.theta = NULL,
	maxit = 100,
	epsilon = 1e-8,
	trace = FALSE
	) {

	distribution <- match.arg(distribution);
	if (distribution %in% c('truncated-poisson', 'truncated-negative-binomial')) {
    if (!requireNamespace("gamlss.tr", quietly = TRUE)) {
            stop("Package 'gamlss.tr' needed for truncated distributions. Please install it.",
                 call. = FALSE)
        }
    }
	### MAIN ##################################################################

	# get observed counts
	# might want to account for case where data is missing? 
	observed.count <- data[[ all.vars(formula)[1] ]];

	# negative binomial can be tricky to fit
	# increase the max number of iterations to help
	glm.control <- stats::glm.control(
		maxit = maxit, 
		epsilon = epsilon, 
		trace = trace
		);

	if( 'negative-binomial' == distribution ) {

		# unfortunately init.theta uses a missing() construct
		# need a different function call if init.theta = NULL
		if( is.null(init.theta) ) {
			model <- MASS::glm.nb(
				formula, 
				data,
				control = glm.control,
				start = start
				);
		} else {
			model <- MASS::glm.nb(
				formula, 
				data,
				control = glm.control,
				init.theta = init.theta,
				start = start
				);
		}
		# sanity check that no rows were lost due to missing data
		.model.rows.sanity.check(data, model);

		expected.values <- model$fitted.values;
		p.values <- stats::pnbinom(
			observed.count - 1, # probability of an observation at least this large
			mu = expected.values,
			size = model$theta,
			lower = FALSE
			);
	} else if( 'poisson' == distribution ) {
		model <- stats::glm(
			formula, 
			data, 
			family = 'poisson',
			control = glm.control,
			start = start
			);

		# Make sure no rows have been lost 
		# (causes problems when adding to the data frame)
		.model.rows.sanity.check(data, model);

		expected.values <- model$fitted.values;
		p.values <- stats::ppois(
			observed.count - 1, # probability of an observation at least this large
			lambda = expected.values,
			lower = FALSE
			);

	} else {
		gamlss.installed <- requireNamespace('gamlss', quietly = TRUE)
		gamlss.tr.installed <- requireNamespace('gamlss.tr', quietly = TRUE)

		if (!gamlss.installed || !gamlss.tr.installed) {
			stop('Truncated distributions depend on the GAMLSS and GAMLSS.TR packages - please install them and try again')
		}

		gamlss.control <- gamlss::gamlss.control(
			trace = trace,
			c.crit = 0.1 # see if this speeds up model fitting
		)

		# make a temporary data frame containing only columns needed for model we want to specify
		# GAMLSS throws an error when there are ANY NAs in the data.. even if we don't use that column
		temp.data <- stats::get_all_vars(formula, data = data);

		if( 'truncated-poisson' == distribution ) {
			gamlss.tr::gen.trun( 0, family = 'PO' );

			# TO DO: 
			#	- extend sanity check to accommodate S4 objects

			model <- gamlss::gamlss(
				formula,
				data = temp.data,
				family = utils::getFromNamespace("POtr", "gamlss.tr")(), 
				control = gamlss.control
				);

			# mean of truncated Poisson is no longer equal to lambda
			# see: https://en.wikipedia.org/wiki/Zero-truncated_Poisson_distribution
			expected.values <- ( model$mu.fv*exp(model$mu.fv) )/(exp(model$mu.fv) - 1);
			p.values <- mapply(	
				FUN = function(observed, mu) {
					# GAMLSS distribution functions do not accept arguments out of range of the distribution
					# Solution: set to 1 if observed count is lowest it can be
					p.value <- ifelse(
						observed >= 2,
						utils::getFromNamespace("pPOtr", "gamlss.tr")(observed - 1, mu = mu, lower.tail = FALSE),
						1
						);

					return(p.value);
					},
				observed.count,
				model$mu.fv
				);

		} else if( 'truncated-negative-binomial' == distribution ) {
			

			gamlss.tr::gen.trun( 0, family = 'NBI' );

			model <- gamlss::gamlss(
				formula,
				data = temp.data,
				family = utils::getFromNamespace("NBItr", "gamlss.tr")(), 
				control = gamlss.control
				);

			# calculate mean of a zero-truncated negative binomial with parameters mu, sigma
			expected.values <- model$mu.fv/( 1 - (1 + model$sigma.fv*model$mu.fv)^(-1/model$sigma.fv) );
			p.values <- mapply(	
				FUN = function(observed, mu, sigma) {
					# GAMLSS distribution functions do not accept arguments out of range of the distribution
					# Solution: set to 1 if observed count is lowest it can be
					p.value <- ifelse(
						observed >= 2,
						utils::getFromNamespace("pNBItr", "gamlss.tr")(
							observed - 1, 
							mu = mu, 
							sigma = sigma,
							lower.tail = FALSE
							),
						1
						);

					return(p.value);
					},
				observed.count,
				model$mu.fv,
				model$sigma.fv
				);	
		} 
	}

	model.data <- list(
		model = model,
		expected.values = expected.values,
		p.values = p.values
		);

	return(model.data);
}

#' model.rows.sanity.check
#' 
#' @description
#' 	Check that the model fit contains the same number of rows as the data used to fit it, 
#'	and throw an error if not
#'
#' @param model.data Data used to fit model
#' @param model Resulting negative binomial model object
#' @rdname chicane
#' @return None
#' 
.model.rows.sanity.check <- function(model.data, model) {

	if( nrow(model.data) != length(model$fit) ) {

		error.message <- paste(
			'Internal bug - row number mismatch.\n', 
			'Data used to fit model contains', nrow(model.data), 'rows\n',
			'Model fit has length', length(model$fit)
			);

		stop(error.message);
	}
}

#' is.glm.nb.maxiter.warning
#'
#' @description
#'	Check if a warning object is an iteration limit reached warning from \code{glm.nb}
#' 
#' @param w Warning object
#' @keywords internal
#' @rdname chicane
#' @return Logical indicating if warning matches iteration limit reached warning
.is.glm.nb.maxiter.warning <- function(w) {

	if( 'iteration limit reached' != w$message ) {
		return(FALSE);
	} 
	if( 'theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace = control$trace > ' != deparse(w$call)[1] ) {
		return(FALSE);
	}
	
	return(TRUE);
}

#' is.glm.nb.theta.warning
#'
#' @description
#' 	Check if a warning matches the square root warning raised by \code{glm.nb} due to an inflated theta estimate.
#'	This happens when the variance of the negative binomial does not exceed the mean (i.e. there is no overdispersion).
#'  In such cases, the Poisson distribution may be a suitable alternative.
#'
#' @param w Warning object
#' @rdname chicane
#' @return Boolean indicating if warning matches
#' 
.is.glm.nb.theta.warning <- function(w) {

	warning.code <- deparse(w$call)[1]

	if( 'sqrt(1/i)' == warning.code && 'NaNs produced' == w$message ) {
		return(TRUE);
    } else if( 'theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace = control$trace > ' == warning.code && 'estimate truncated at zero' == w$message ) {
    	return( TRUE );
    } else {
    	return(FALSE);
    }
}

# check if packages are installed for truncated distributions
.check.packages <- function(distribution) {
	if (distribution == 'poisson' || distribution == 'negative-binomial') {
		return(TRUE)
	}
	gamlss.installed <- requireNamespace('gamlss', quietly = TRUE)
	gamlss.tr.installed <- requireNamespace('gamlss.tr', quietly = TRUE)

	if (!gamlss.installed || !gamlss.tr.installed) {
		stop('Truncated distributions depend on the GAMLSS and GAMLSS.TR packages - please install them and try again')
	}
}