###G_matrix creation 
##Creation of genomic similarity matrix for a list of animal and a race:
.libPaths("/bao/lib_R3.5.0")
# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
})

##To be able to launch qsub:
Sys.setenv(SGE_ROOT = "/opt/sge")

# -------------------------------
# Arguments
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2){
  stop("Usage: Rscript script.R <file> <Race>")
}

listeANI <- args[1]
RACE <- args[2]
#listeANI="/espace_projets/inrae_gabi/rumigen/DATA/Liste_anim_2026-01-14.csv"

# RACE=66
# -------------------------------
# Parameters
# -------------------------------
nbBTA <- 29  # Default number of chromosomes

# [TODO] Set your working directories
USER <- Sys.getenv("USER")
TRAVdata <- paste0("/travail/",USER,"/GRM/data")
TRAVlog <- paste0("/travail/",USER,"/GRM/log")
TRAVmat <- paste0("/travail/",USER,"/GRM/mat")
DATE_BIG_WORK <- system("ls /big_work/INFOGENO/IMPUTATIONS/50K+/ | sort -rn | head -1", intern = TRUE)
PHASES=paste0("/big_work/INFOGENO/IMPUTATIONS/50K+/",DATE_BIG_WORK,"/",RACE,"/phasesnumeq")
CARTE <- paste0("/big_work/INFOGENO/IMPUTATIONS/50K+/",DATE_BIG_WORK,"/CARTE/snpIndexation_chromo")
nomMAT <- paste0("GRM_r", RACE)

# -------------------------------
# Extract phases
# -------------------------------
message("Extracting phases in ", TRAVdata)
dir.create(TRAVdata,recursive = T,showWarnings = F)
dir.create(TRAVlog,recursive = T,showWarnings = F)
setwd(TRAVdata)

for(chr in 1:nbBTA){
  phase_file <- paste0(PHASES, chr)
  output_file <- paste0("phases", chr)
  
  # Build awk command to filter animals in listeANI
  awk_cmd <- paste0(
    "awk '{if(FILENAME==ARGV[1]){p[$1]=1};",
    "if(FILENAME==ARGV[2]){if(p[$1]==1){print $0}}}' ",
    listeANI, " ", phase_file, " > ", output_file
  )
  
  # Run awk system call
  system(awk_cmd)
}



# Preallocate as list (MUCH faster than repeated cbind)
all_chr_list <- vector("list", nbBTA)
ID <- NULL

for (chr in 1:nbBTA) {
  
  phases <- fread(paste0("/travail/fbesnard/GRM/data/phases", chr))
  
  # Precompute row indices
  r1 <- seq(1, nrow(phases), 2)
  r2 <- r1 + 1
  
  # Compute genotype matrix for chromosome (fully vectorized)
  geno <- as.matrix(phases[r1, 3:ncol(phases), with = FALSE] +
                      phases[r2, 3:ncol(phases), with = FALSE] - 2)
  
  # Store ID only once
  if (is.null(ID)) {
    ID <- phases[r1, 1]
  }
  
  all_chr_list[[chr]] <- geno
}

# Final big matrix assembly (fast, in one step)
all_chr_typages <- data.frame(do.call(cbind, all_chr_list))

# Assign rownames
rownames(all_chr_typages) <- ID$V1

ID_list=fread(listeANI,header=F)

# Keep only animals in ID_com
all_chr_typages <- all_chr_typages[ID_list$V1, , drop = FALSE]

#Remove if non existing obs (only NA )
all_chr_typages <- all_chr_typages[rowSums(!is.na(all_chr_typages)) > 0, , drop = FALSE]

# Calculate allele frequency p for each SNP:
# mean(genotype)/2   →   fast columnMeans

p=apply(all_chr_typages,2,mean)/2 

s2pq=sum(2*p*(1-p)) 

grm=tcrossprod(scale(all_chr_typages,center = T,scale = F))/s2pq




fwrite(as.data.table(grm),"/espace_projets/inrae_gabi/rumigen/DATA/bglr/grm",sep=" " )

quit()
