#!/usr/bin/env Rscript
# Export relevant ASCAT tables to TSVs for downstream Python processing
# Usage:
#   Rscript unpack_RDS_ascat_helper.R --rdata path/to/ASCAT_objects.Rdata --out_prefix LTX0001_Tumour1 --output_dir path/to/outputs

args <- commandArgs(trailingOnly = TRUE)
parseArgs <- function(args){
  out <- list()
  i <- 1
  while(i <= length(args)){
    key <- args[i]
    if(substr(key, 1, 2) == "--"){
      name <- substring(key, 3)
      val <- if(i+1 <= length(args)) args[i+1] else TRUE
      out[[name]] <- val
    }
    i <- i + 2
  }
  out
}

opt <- parseArgs(args)
if(is.null(opt$rdata)) stop("--rdata is required")
if(is.null(opt$out_prefix)) opt$out_prefix <- "ascat_export"

# optional output directory to place all generated files
if(is.null(opt$output_dir)) opt$output_dir <- "."
if(!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)

load(opt$rdata)
if(!exists("ascat.bc") || !exists("ascat.output")) stop("RData must contain ascat.bc and ascat.output")

ASCATobj <- ascat.bc
ascat_result <- ascat.output

samples <- ASCATobj$samples

segments_raw <- ascat_result$segments_raw

seg_list <- list()
snp_list <- list()
meta_list <- list()

SNPpos <- ASCATobj$SNPpos
snp_ids <- rownames(SNPpos)
chr <- SNPpos[, 1]
pos <- SNPpos[, 2]

for(i in seq_along(samples)){
  sample_name <- samples[i]
  sample_idx <- i

  # segments (unrounded/raw) for this sample
  if(is.data.frame(segments_raw)){
    seg_sample <- segments_raw[segments_raw$sample == sample_name, , drop = FALSE]
  } else {
    seg_sample <- as.data.frame(segments_raw)
    if("sample" %in% colnames(seg_sample)){
      seg_sample <- seg_sample[seg_sample$sample == sample_name, , drop = FALSE]
    }
  }
  if(nrow(seg_sample) == 0){
    # create empty stub with sample info if no segments
    seg_sample <- data.frame(sample = sample_name, stringsAsFactors = FALSE)
  } else {
    seg_sample$sample <- sample_name
  }
  seg_list[[length(seg_list) + 1]] <- seg_sample

  # purity / ploidy / psi / gender / sexchromosomes
  purity_val <- as.numeric(ascat_result$purity[sample_name])
  ploidy_val <- as.numeric(ascat_result$ploidy[sample_name])
  psi_val <- as.numeric(ascat_result$psi[sample_name])
  gender <- ASCATobj$gender[sample_idx]
  sexchromosomes <- ASCATobj$sexchromosomes
  X_nonPAR <- ASCATobj$X_nonPAR
  meta <- data.frame(sample = sample_name,
                     purity = purity_val,
                     ploidy = ploidy_val,
                     psi = psi_val,
                     gender = gender,
                     sexchromosomes = paste(sexchromosomes, collapse=","),
                     X_nonPAR = if(is.null(X_nonPAR)) NA else paste(X_nonPAR, collapse=","),
                     stringsAsFactors = FALSE)
  meta_list[[length(meta_list) + 1]] <- meta

  # per-SNP table for this sample
  LogR <- ASCATobj$Tumor_LogR[, sample_idx]
  BAF <- ASCATobj$Tumor_BAF[, sample_idx]
  LogR_seg <- if(!is.null(ASCATobj$Tumor_LogR_segmented)) ASCATobj$Tumor_LogR_segmented[, sample_idx] else rep(NA, length(LogR))

  # Tumor_BAF_segmented[[sample_idx]] contains only het probes — expand to full length
  BAF_seg_sample <- NULL
  if(!is.null(ASCATobj$Tumor_BAF_segmented)){
    BAF_seg_sample <- ASCATobj$Tumor_BAF_segmented[[sample_idx]]
  }
  BAF_seg_full <- rep(NA, length(LogR))
  names(BAF_seg_full) <- snp_ids
  if(!is.null(BAF_seg_sample)){
    if(is.matrix(BAF_seg_sample) || is.data.frame(BAF_seg_sample)){
      vals <- as.numeric(BAF_seg_sample[,1])
      ids <- rownames(BAF_seg_sample)
      if(!is.null(ids)) BAF_seg_full[ids] <- vals
    } else if(!is.null(names(BAF_seg_sample))){
      BAF_seg_full[names(BAF_seg_sample)] <- as.numeric(BAF_seg_sample)
    } else if(length(BAF_seg_sample) == length(BAF_seg_full)){
      BAF_seg_full[] <- as.numeric(BAF_seg_sample)
    }
  }

  snp_df <- data.frame(
    sample = sample_name,
    SNP = snp_ids,
    chr = chr,
    pos = pos,
    LogR = as.numeric(LogR),
    BAF = as.numeric(BAF),
    LogR_seg = as.numeric(LogR_seg),
    BAF_seg = as.numeric(BAF_seg_full),
    stringsAsFactors = FALSE
  )

  # assign segment ids by position (segment index within this sample)
  snp_df$segment_id <- NA_integer_
  if(nrow(seg_sample) > 0 && all(c('chr','startpos','endpos') %in% colnames(seg_sample))){
    for(j in seq_len(nrow(seg_sample))){
      seg <- seg_sample[j, ]
      in_seg <- which(snp_df$chr == seg$chr & snp_df$pos >= seg$startpos & snp_df$pos <= seg$endpos)
      if(length(in_seg)>0) snp_df$segment_id[in_seg] <- j
    }
  }
  snp_list[[length(snp_list) + 1]] <- snp_df
}

# combine and write
seg_all_df <- do.call(rbind, seg_list)
snp_all_df <- do.call(rbind, snp_list)
meta_all_df <- do.call(rbind, meta_list)

seg_file <- file.path(opt$output_dir, paste0(opt$out_prefix, "_segments.tsv"))
meta_file <- file.path(opt$output_dir, paste0(opt$out_prefix, "_purity_ploidy.tsv"))
snp_file <- file.path(opt$output_dir, paste0(opt$out_prefix, "_snps.tsv"))

write.table(seg_all_df, file = seg_file, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(meta_all_df, file = meta_file, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(snp_all_df, file = snp_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Wrote:", seg_file, "\n")
cat("Wrote:", meta_file, "\n")
cat("Wrote:", snp_file, "\n")
