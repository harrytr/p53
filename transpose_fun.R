transpose_fun <- function(df)
  
library(data.table)

t_df <- data.table::transpose(df)

# get row and colnames in order
colnames(t_df) <- rownames(df)
rownames(t_df) <- colnames(df)


return(df)