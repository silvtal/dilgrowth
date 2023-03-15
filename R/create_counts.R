#' This function parses and puts together in a list of data.frames all the 
#' count data for time series abundance data. For each timestep / transfer, it
#' saves the mean value of all available replicates.
#' 
#' If there are no available replicates for a given transfer, it saves the data
#' for the previous transfer (in other words, assumes no abundance changes took
#' place). This assumption can be overridden with the option 
#' "allow_empty_transfers", that will set "NA" instead.
#' 
#' If select_transfers is not NULL, is has to be a vector. If you want to pick
#' only transfers 0, 2 and 4, you must input c(0, 2, 4)
#' 
#' @param exp Abundance table
#' @param map Metadata about m_inic, replicates, sample names...
#' @param m_inic Vector of samples to pick (by sample name)
#' @param allow_empty_transfers 
#' @param select_transfers Vector of transfers to pick from the abundance table
#' @param orig Original sample tag column name in map
#' @param transfer Transfer tag column name in map
#' @param sa Sample name tag column name in map
#'
#' @return
#' @export
create_counts <- function (exp,
                           map,
                           m_inic,
                           allow_empty_transfers = FALSE,
                           select_transfers = NULL,
                           orig = "ORIG",
                           transfer = "T",
                           sa = "SA") { 
  counts_ <- list()
  if (!is.null(select_transfers)){
    transfers=select_transfers
  } else {
    transfers=0:length(unique(map[[transfer]]))      # includes transfer "0" (original)
  }
  for (i in m_inic) {
    counts_[[i]] <- data.frame(matrix(ncol = 0, nrow = nrow(exp))) # prepare empty DF
    for (t in transfers){
      ### separo por muestra inicial los datos finales también
      temp <- exp[map[map[orig]==i & map[transfer]==t,][,sa]]
      if (ncol(temp)==0){  # if we don't have data for this m_inic in this transfer, 
                           # skip and repeat the data for the last transfer instead.
        if (!allow_empty_transfers){
          counts_[[i]]=cbind(counts_[[i]],
                             setNames(counts_[[i]][ncol(counts_[[i]])],nm=t)) 
        } else {
          counts_[[i]]=cbind(counts_[[i]],
                             setNames(as.data.frame(rep(NA,nrow(counts_[[i]]))),nm=t))
        }
      } else {
        ### calculate the mean for all replicates for that transfer
        r <- rowSums(temp)/dim(temp)[2] # might return float
        ### finally, I add the mean data for that transfer to the final counts_ list
        counts_[[i]]=cbind(counts_[[i]],
                           setNames(as.data.frame(r),nm=t)) #añado la nueva transfer
        counts_[[i]][is.na(counts_[[i]])]=0
      }
    }                                                                                   
  }  
  return(counts_) # FIX en el script usar create counts para parsear el resultado de la simulación, guardado ya en dataframe...
}