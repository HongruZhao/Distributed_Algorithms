my_lib <- "/panfs/jay/groups/27/shenx/zhao1118/distributed_algorithm/R_libs"
dir.create(my_lib, showWarnings = FALSE, recursive = TRUE)


install.packages(c("MASS", "glmtlp"),
                 lib   = my_lib,
                 repos = "https://cloud.r-project.org")


#.libPaths(c("/panfs/jay/groups/27/shenx/zhao1118/distributed_algorithm/R_libs", .libPaths()))

#library(MASS)     # loads from your folder if present
#library(glmtlp)   # idem
