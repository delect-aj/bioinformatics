Untitled
================

# Step 1: Obtain Initial Metadata from Each Project

To obtain the initial metadata: 1. Search for the BioProject ID (e.g.,
“PRJEBXXXX”) on the [NCBI SRA
website](https://www.ncbi.nlm.nih.gov/Traces/study/). 2. In the “Select”
table, click on “Metadata” and download the file. 3. The downloaded file
is typically named “SraRunTable.csv” and contains metadata for the
project.

``` r
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(gridExtra)

# Set working directory (update this path to your local directory)
setwd("D:/oral_project/PRJNA782768")

# Read the initial metadata files
initial_meta <- read.csv("SraRunTable.csv", header = TRUE)
initial_meta1 <- read.csv("SraRunTable1.csv", header = TRUE)

# Combine the two metadata files
meta <- rbind(initial_meta, initial_meta1)

# Remove entries where Organism is "gut metagenome" (non-oral samples)
meta <- subset(meta, Organism != "gut metagenome")

# Check the number of samples containing "HT" in Sample.Name (e.g., healthy controls)
nrow(subset(meta, grepl("HT", Sample.Name)))
```

    ## [1] 34

``` r
# Display the first few rows of the metadata
head(meta)
```

    ##            Run Assay.Type AvgSpotLen    Bases  BioProject    BioSample
    ## 11 SRR17012673   AMPLICON        493 44859056 PRJNA782768 SAMN23394734
    ## 12 SRR17012674   AMPLICON        494 38840750 PRJNA782768 SAMN23394725
    ## 13 SRR17012675   AMPLICON        495 38367945 PRJNA782768 SAMN23394724
    ## 14 SRR17012676   AMPLICON        492 34512816 PRJNA782768 SAMN23394783
    ## 15 SRR17012677   AMPLICON        493 35535440 PRJNA782768 SAMN23394782
    ## 16 SRR17012678   AMPLICON        494 31619458 PRJNA782768 SAMN23394781
    ##                 BioSampleModel    Bytes              Center.Name
    ## 11 Metagenome or environmental 14879583 SCHOOL OF BASIC MEDICINE
    ## 12 Metagenome or environmental 13064180 SCHOOL OF BASIC MEDICINE
    ## 13 Metagenome or environmental 12895958 SCHOOL OF BASIC MEDICINE
    ## 14 Metagenome or environmental 11630135 SCHOOL OF BASIC MEDICINE
    ## 15 Metagenome or environmental 11872474 SCHOOL OF BASIC MEDICINE
    ## 16 Metagenome or environmental 10519269 SCHOOL OF BASIC MEDICINE
    ##    Collection_Date Consent DATASTORE.filetype DATASTORE.provider
    ## 11      2021-01-05  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 12      2021-01-05  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 13      2021-01-05  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 14      2021-01-06  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 15      2021-01-06  public   fastq,run.zq,sra         gs,ncbi,s3
    ## 16      2021-01-06  public   fastq,run.zq,sra         gs,ncbi,s3
    ##                        DATASTORE.region  Experiment geo_loc_name_country
    ## 11 gs.us-east1,ncbi.public,s3.us-east-1 SRX13202881                China
    ## 12 gs.us-east1,ncbi.public,s3.us-east-1 SRX13202880                China
    ## 13 gs.us-east1,ncbi.public,s3.us-east-1 SRX13202879                China
    ## 14 gs.us-east1,ncbi.public,s3.us-east-1 SRX13202877                China
    ## 15 gs.us-east1,ncbi.public,s3.us-east-1 SRX13202876                China
    ## 16 gs.us-east1,ncbi.public,s3.us-east-1 SRX13202875                China
    ##    geo_loc_name_country_continent   geo_loc_name    HOST     Instrument
    ## 11                           Asia China:Shanghai missing Illumina MiSeq
    ## 12                           Asia China:Shanghai missing Illumina MiSeq
    ## 13                           Asia China:Shanghai missing Illumina MiSeq
    ## 14                           Asia China:Shanghai missing Illumina MiSeq
    ## 15                           Asia China:Shanghai missing Illumina MiSeq
    ## 16                           Asia China:Shanghai missing Illumina MiSeq
    ##    isolation_source          lat_lon Library.Name LibraryLayout
    ## 11          missing 31.19 N 121.60 E         T11a        PAIRED
    ## 12          missing 31.19 N 121.60 E         T02a        PAIRED
    ## 13          missing 31.19 N 121.60 E         T01a        PAIRED
    ## 14          missing 31.19 N 121.60 E         Tb20        PAIRED
    ## 15          missing 31.19 N 121.60 E         Tb19        PAIRED
    ## 16          missing 31.19 N 121.60 E         Tb18        PAIRED
    ##    LibrarySelection LibrarySource        Organism Platform          ReleaseDate
    ## 11              PCR   METAGENOMIC oral metagenome ILLUMINA 2022-03-19T00:00:00Z
    ## 12              PCR   METAGENOMIC oral metagenome ILLUMINA 2022-03-19T00:00:00Z
    ## 13              PCR   METAGENOMIC oral metagenome ILLUMINA 2022-03-19T00:00:00Z
    ## 14              PCR   METAGENOMIC oral metagenome ILLUMINA 2022-03-19T00:00:00Z
    ## 15              PCR   METAGENOMIC oral metagenome ILLUMINA 2022-03-19T00:00:00Z
    ## 16              PCR   METAGENOMIC oral metagenome ILLUMINA 2022-03-19T00:00:00Z
    ##             create_date version Sample.Name SRA.Study
    ## 11 2021-11-23T02:11:00Z       1        T11a SRP347352
    ## 12 2021-11-23T02:11:00Z       1        T02a SRP347352
    ## 13 2021-11-23T02:11:00Z       1        T01a SRP347352
    ## 14 2021-11-23T02:11:00Z       1        Tb20 SRP347352
    ## 15 2021-11-23T02:11:00Z       1        Tb19 SRP347352
    ## 16 2021-11-23T02:11:00Z       1        Tb18 SRP347352

# Step 2: Create a Uniform Metadata Template

The uniform metadata should include the following columns:

1.  **BioProject**: Project ID (e.g., “PRJNA339212”)
2.  **Run**: Sample ID (e.g., “ERRXXXXXX” or “SRRXXXXXX”)
3.  **sample_name**: Original sample name submitted by the authors (not
    the database-assigned ID)
4.  **Age**: Age of the subject (if available)
5.  **BMI**: Body Mass Index (if available)
6.  **Sex**: Sex of the subject (“F” or “M”)
7.  **Region**: Targeted 16S rRNA region (e.g., “V3-V4”)
8.  **Platform**: Sequencing platform (e.g., “Illumina MiSeq”)
9.  **Country**: Country of sample origin (e.g., “China”)
10. **Site**: Sample site (e.g., “saliva”, “dental plaque”)
11. **public_time**: Data release date (format: “YYYY-MM-DD”)
12. **PCR_primers**: PCR primers used (e.g., “515F-806R”)
13. **title**: Title of the associated publication
14. **Disease_study**: Whether the study involves a disease (“YES” or
    “NO”)
15. **Status**: Disease status (e.g., “T2DM”, “HC”). Note: If the study
    includes both disease and healthy samples, this must be specified
    per sample.

``` r
# Create a uniform metadata dataframe
uniform_meta <- data.frame(
  BioProject = meta$BioProject,
  Run = meta$Run,
  sample_name = meta$Sample.Name,
  Age = NA,  # To be filled Update based on you metadata
  BMI = NA,  # To be filled based on you metadata
  Sex = NA,  # To be filled based on you metadata
  Region = "V3-V4",  # Update if different for your study
  Platform = meta$Instrument,
  Country = meta$geo_loc_name_country,
  Site = meta$isolation_source, # other names of this columns such as env_median, env_material
  public_time = meta$ReleaseDate,
  stringsAsFactors = FALSE
)
```

# Step 3: Manually Add Information from the Publication

Some information (e.g., Sex, PCR primers, study title, disease status)
must be manually extracted from the publication and added to the
metadata.

``` r
# Add information obtained from the publication
uniform_meta$PCR_primers <- "338F-806R"  # Update based on the publication
uniform_meta$title <- "Distribution characteristics of oral microbiota and its relationship with intestinal microbiota in patients with type 2 diabetes mellitus"
uniform_meta$Disease_study <- "YES"

# Assign disease status based on sample_name patterns
uniform_meta <- uniform_meta %>%
  mutate(Status = case_when(
    grepl("DT", sample_name) ~ "T2DM",
    grepl("HT|Tb", sample_name) ~ "HC"
  ))
```

# Step 4: Summarize and Validate the Metadata

``` r
# Check the number of unique runs and sample names
print(paste("Number of unique runs:", length(unique(uniform_meta$Run))))
```

    ## [1] "Number of unique runs: 298"

``` r
print(paste("Number of unique sample names:", length(unique(uniform_meta$sample_name))))
```

    ## [1] "Number of unique sample names: 298"

# Step 5: Visualize the Metadata

``` r
# Select columns for visualization
plot_data <- uniform_meta %>%
  select(-BioProject, -Run, -sample_name, -title, -PCR_primers, -Age, -BMI, -Sex)

# Function to create simple plots for each column
create_simple_plots <- function(data) {
  plots <- list()
  
  for (col in colnames(data)) {
    if (is.numeric(data[[col]])) {
      # Numeric variable: histogram
      p <- ggplot(data, aes(x = .data[[col]])) +
        geom_histogram(fill = "blue", alpha = 0.6, bins = 20) +
        labs(title = col) +
        theme_minimal()
    } else {
      # Categorical variable: bar plot
      count_data <- as.data.frame(table(data[[col]]))
      p <- ggplot(count_data, aes(x = reorder(Var1, -Freq), y = Freq)) +
        geom_col(fill = "red", alpha = 0.6) +
        labs(title = col, x = NULL, y = "Count") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
    plots[[col]] <- p
  }
  
  return(plots)
}

# Generate and arrange plots
plots <- create_simple_plots(plot_data)
grid.arrange(grobs = plots, ncol = 3, top = "Metadata Summary")
```

![](metadata_collection_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->
\# Step 6: Finalize and Export the Metadata

``` r
# Reorder columns to match the desired format
uniform_meta <- uniform_meta[, c(
  "BioProject", "Run", "sample_name", "Age", "BMI", "Sex", 
  "Region", "Platform", "Country", "Site", "public_time", 
  "PCR_primers", "title", "Disease_study", "Status"
)]

# Export the uniform metadata to a CSV file
write.csv(uniform_meta, file = "uniform_meta.csv", quote = FALSE)
```

## Notes for Students:

1.  **Metadata Accuracy**: Always cross-check the metadata with the
    original publication to ensure accuracy.
2.  **Manual Entries**: Columns like `Age`, `BMI`, `Sex`, `PCR_primers`,
    and `Status` often require manual entry based on the publication.
3.  **Disease Status**: Be careful when assigning `Status`—some studies
    include both disease and healthy samples.
4.  **Visualization**: The plots help quickly identify the distribution
    of values in each column, which is useful for quality control.
5.  **Reproducibility**: Update the working directory and file paths to
    match your local environment.
