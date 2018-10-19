# Walk directories with a starting point. Concatenate bed files. Run annotation on human.

# Delete extended files

# Find human Genome from bioconductor

# Walk, find ".bed" file

# Add back in HepG2

# Set wd to the histogram files folder

library(tidyverse)
library(geneXtendeR)

basedir <- getwd()
path <- "./../results/sorted"

# human <- readGFF("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz")
human.basic <- readGFF("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.basic.annotation.gtf.gz")

# Function to walk and concat all bed files within the results folders (less of a mess)
walk_and_concat <- function(path, skipping = TRUE) {
  dirs <- list.dirs(path)
  for (dir in dirs) {
      files <- list.files(dir)
      files <- files[grepl("*.bed", files)]
      concat(dir, files, name = paste0(basename(dir), ".bed"), skip = skipping)
  }
}

# Concat files in a path and bundle into a bed file with a specific name and specific folder
concat <- function(path, files, name = "all_bed.bed", skip = TRUE, folder = "folderbeds") {
    # Make sure not empty
    if (length(files) == 0) {
        message("Skipping this empty dir: ", path)
        return()
    }
    # Check if they are already done:
    if (name %in% list.files(folder) & skip) {
        message("This directory is already completed. Skipping :", path)
        return()
    }
    
    filefull <- file.path(path, files)
    # Import all the data into a linked-list and bind_rows, then write
    data <- tibble()
    for (fp in filefull) {
        data <- bind_rows(data, read_tsv(fp, col_names = FALSE, col_types = list(col_character(), col_integer(), col_integer())))
    }
    write_tsv(data, file.path(folder, name), append = FALSE, col_names = TRUE)
    message("Finished this dir: ", path)
}

# Walk the result folders and combine all the data within those files, then move to a folder called "folderbeds"
walk_and_concat(path)

# Concat all folderbed files into a new folder called singularity (file called super.bed)
concat("folderbeds", list.files("folderbeds"), name = "super.bed", skip = TRUE, folder = "singularity")

# Import the file to be split up by chromosome
data <- read_tsv(file.path("singularity", "super.bed"), col_names = TRUE, col_types = list(col_character(), col_integer(), col_integer())) %>%
    mutate(chr = as.factor(chr))
lev <- levels(data$chr)

path = "singularity"

# peaksInput("singularity/super.bed")

# Split the file by chromosome and save as tsv
for (level in lev) {
    fp <- file.path("singularity", paste0("chrom", level))
    if (!fp %in% list.dirs(path)) {
        dir.create(fp)
    }
    write_tsv(data %>% filter(chr == level), file.path(fp, paste0("chrom", level, ".tsv")), col_names = TRUE, append = FALSE)
}

# Function to create peaks files for each chromosome and save, then for each chromosome annotate and save, then create a histogram with.
generate_files <- function(bdir, path = "singularity") {
    directories <- list.dirs(file.path(bdir, path))
    directories <- directories[2: length(directories)]
    all_cut <- as_tibble()
    for (dir in directories) {
        setwd(dir)
        file <- paste0(basename(dir), ".tsv")
        if (!file.exists("peaks.txt")) {
            peaksInput(file)
        }
        message(file, " peaks finished. Moving onto annotation_n")
        if (!file.exists("annotated_1500_3.txt")) {
            message(file, " Not annotated yet, annotating 1500 3")
            dump <- annotate_n(human.basic, 1500, 3)
        }
        message(file, " completed entirely, grab data later")
        
        if (file.exists("annotated_1500_3.txt")) {
            message("Now using ", basename(dir), " and binning data")
            data <- read_tsv("annotated_1500_3.txt", col_types = c(col_integer(), col_integer(), col_integer(),
                                                                   col_integer(), col_integer(), col_integer(),
                                                                   col_character(), col_character(),
                                                                   col_integer(), col_integer()))
            
            data <- data %>% dplyr::select(`Peak-Num`, rank, `Minimum-Distance-to-Gene`)
            
            # Bohdans' way -- Relative distance. per peak do (rank2 - rank1) and count.
            db1 <- data %>% filter(rank == 1) %>%
                dplyr::select(`Minimum-Distance-to-Gene`)
            db2 <- data %>% filter(rank == 2) %>%
                dplyr::select(`Minimum-Distance-to-Gene`)
            
            cut <- (db2 - db1) %>% filter(`Minimum-Distance-to-Gene` <= 50000)
            
            all_cut <- bind_rows(all_cut, cut)
            #My way
            # db3 <- data %>% filter(rank == 3) %>%
            #     dplyr::select(`Minimum-Distance-to-Gene`)
        }
    }
    
    grap <- ggplot(all_cut, aes(x=`Minimum-Distance-to-Gene`)) +
        geom_histogram(binwidth = 500, color="black", fill="white") +
        labs(title = "Relative distance between first and second closest genes to each peak", x="Relative distance", y="Count")
    grap
    grap2 <- ggplot(all_cut, aes(x=`Minimum-Distance-to-Gene`)) +
        geom_histogram(binwidth = 2000, color="black", fill="white") +
        labs(title = "Relative distance between first and second closest genes to each peak", x="Relative distance", y="Count")
    grap2
    
    return(list(grap, grap2))
}

#keep base dir
basedir <- getwd()
graph <- generate_files(bdir = basedir)
#reset base dir after it is finished (takes quite a lot of time total, easier to do in chunks because each step stores the data)
setwd(basedir)



