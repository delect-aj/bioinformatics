### MaAsLin2参考流程，https://github.com/biobakery/biobakery/wiki/maaslin2
library(Maaslin2)

input_data = system.file("extdata", "HMP2_taxonomy.tsv", package="Maaslin2") # The abundance table file

input_metadata = system.file("extdata", "HMP2_metadata.tsv", package="Maaslin2") # The metadata table file

df_input_data = read.table(file = input_data,
                           header = TRUE,
                           sep = "\t", 
                           row.names = 1,
                           stringsAsFactors = FALSE, 
                           check.names = FALSE)

df_input_metadata = read.csv(file = input_metadata,
                           header = TRUE,
                           sep = "\t", 
                           row.names = 1,
                           stringsAsFactors = FALSE, 
                           check.names = FALSE)
### run Maaslin2
fit_data2 = Maaslin2(input_data     = df_input_data, 
                     input_metadata = df_input_metadata, 
                     min_prevalence = 0,
                     normalization  = "NONE",
                     output         = "difference_analysis/demo_output2", 
                     fixed_effects  = c("diagnosis", "dysbiosis"),
                     reference      = c("diagnosis,nonIBD"))


