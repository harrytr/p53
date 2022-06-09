shorten_names <- function() {
  temp_dir <- getwd()
  for (i in 1: length(list.files(temp_dir)) ) {
    
    
    
    
    temp <- list.files(temp_dir)[i]
    file.rename(temp,substr(temp, 8, length(temp)))
    
  }
}