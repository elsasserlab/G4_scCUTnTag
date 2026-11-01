# Parse command-line options. Existing BED files are used by default.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L || (length(args) == 1L && args[[1]] != "--run-pqsfinder")) {
  stop("Usage: Rscript pqs_mm10.R [--run-pqsfinder]")
}
run_pqsfinder <- "--run-pqsfinder" %in% args

# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (run_pqsfinder && !requireNamespace("pqsfinder", quietly = TRUE))
    BiocManager::install("pqsfinder")

if (!requireNamespace("rtracklayer", quietly = TRUE))
    BiocManager::install("rtracklayer")

if (!requireNamespace("BSgenome", quietly = TRUE))
    BiocManager::install("BSgenome")

if (run_pqsfinder)
    library(pqsfinder)
library(rtracklayer)
library(GenomicRanges)
library(BSgenome)

# Get chromosome lengths from BSgenome (no external FASTA needed)
library(BSgenome.Mmusculus.UCSC.mm10)
mm10_genome <- BSgenome.Mmusculus.UCSC.mm10
seq_lengths <- seqlengths(mm10_genome)

# Resolve paths relative to this script so it can be run from anywhere
script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)])
script_dir <- if (length(script_file) == 1L) dirname(normalizePath(script_file)) else getwd()

bed_file <- file.path(script_dir, "PQS_scores.mm10.bed")
if (run_pqsfinder) {
  # Explicitly requested: run the expensive genome-wide scan and overwrite BED.
  cat("Running pqsfinder because --run-pqsfinder was supplied.\n")

  # Initialize an empty GRanges object to store results
  combined_granges <- GRanges()

  # Process each chromosome
  for (chrom in names(mm10_genome)) {
    cat("Processing", chrom, "...\n")

    # Extract DNA sequence for the current chromosome
    chrom_seq <- mm10_genome[[chrom]]

    # Run pqsfinder on the chromosome sequence
    pqs_results <- pqsfinder(chrom_seq, min_score = 20)

    # Convert pqsfinder results to a GRanges object
    gr <- GRanges(seqnames = chrom,
                  ranges = IRanges(start = start(pqs_results),
                                   end = end(pqs_results)),
                  score = score(pqs_results))

    # Append to the combined GRanges object
    combined_granges <- c(combined_granges, gr)
  }

  # Write the non-reduced GRanges to a BED file
  cat("Writing to BED file:", bed_file, "\n")
  export(combined_granges, bed_file, format = "BED")
} else if (file.exists(bed_file)) {
  # Reuse existing PQS calls and skip the expensive genome-wide scan.
  cat("Loading existing BED file:", bed_file, "\n")
  combined_granges <- import(bed_file, format = "BED")
} else {
  stop("BED file not found: ", bed_file,
       ". Supply --run-pqsfinder to generate it.")
}

# Assign sequence lengths to the GRanges object
seqlengths(combined_granges) <- seq_lengths[seqlevels(combined_granges)]

# Reduce GRanges to ensure no overlapping ranges for BigWig
reduced_granges <- reduce(combined_granges, with.revmap = TRUE, min.gapwidth = 0)

# Include chromosomes with no PQS hits so gaps() can create genome-wide zero regions.
seqlevels(reduced_granges) <- names(seq_lengths)
seqlengths(reduced_granges) <- seq_lengths

# Calculate average scores for the reduced PQS ranges
revmap <- mcols(reduced_granges)$revmap
avg_scores <- numeric(length(reduced_granges))

mapped_indices <- which(lengths(revmap) > 0)
avg_scores[mapped_indices] <- sapply(revmap[mapped_indices], function(idx) {
  if (length(idx) > 0) {
    mean(mcols(combined_granges[idx])$score, na.rm = TRUE)
  } else {
    0
  }
})
mcols(reduced_granges)$score <- avg_scores

# Add all regions not covered by a PQS interval with score 0.  Assigning zero
# only to existing ranges does not create the inter-interval regions.
strand(reduced_granges) <- "*"
zero_ranges <- gaps(reduced_granges)
zero_ranges <- zero_ranges[strand(zero_ranges) == "*"]
mcols(zero_ranges)$score <- 0
reduced_granges <- sort(c(reduced_granges, zero_ranges), ignore.strand = TRUE)

cat("final check...\n")
# Replace any remaining NA scores with 0 (final check)
mcols(reduced_granges)$score[is.na(mcols(reduced_granges)$score)] <- 0

# Assign sequence lengths to the reduced GRanges object
seqlengths(reduced_granges) <- seq_lengths[seqlevels(reduced_granges)]

# Export as a BigWig track
bigwig_file <- file.path(script_dir, "PQS_scores.mm10.bw")
cat("Writing to BigWig file:", bigwig_file, "\n")
export(reduced_granges, bigwig_file, format = "BigWig")
