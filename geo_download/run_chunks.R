## parallelizing functions

chunkRun <- function(list.items, my_function, out_dir, chunk_num, CHUNK_SIZE=100, GROUP_SIZE=10, log_prefix="chunk"){
  logfile <- sprintf("logs/%s_%s.log", log_prefix, chunk_num)
  cat(" \n", file=logfile)
  
  if (chunk_num*CHUNK_SIZE >length(list.items)){
    chunk_gses <- list.items[((chunk_num-1)*CHUNK_SIZE):length(list.items)]
  } else {
    chunk_gses <- list.items[((chunk_num-1)*CHUNK_SIZE):(chunk_num*CHUNK_SIZE)]
  }
  print(chunk_gses)
  
  NUM_GROUPS= ceiling((length(chunk_gses))/GROUP_SIZE)
  
  # run in chunks, then save R object
  for (i in 1:NUM_GROUPS){
    tryCatch({
      if (i==NUM_GROUPS){
        chunk_gses_i <- chunk_gses[((i-1)*GROUP_SIZE+1):length(chunk_gses)]
      } else {
        chunk_gses_i <- chunk_gses[((i-1)*GROUP_SIZE+1):(i*GROUP_SIZE)]
      }
      chunk_ds <- lapply(chunk_gses_i, my_function)
      names(chunk_ds) <- chunk_gses_i
      save(chunk_ds, file=sprintf("%s/chunk%s_%s.RData", out_dir, chunk_num, i))
    },error = function(err){
      print(err)
    })
  }
}




