#' Scoring System
#'
#' @param scoredf The data frame of containing all assessment content
#' @param SimiRetaincutoff The cutoff of the structure similarity. The default value is 0.8.
#' @import tibble
#'
#' @return a dataframe with scores for each embedding.
#' @export

score_cal <- function(scoredf, SimiRetaincutoff=0, GOFcutoff=0.53, URATEcutoff=0.029) {
  if (all(is.na(scoredf))) {
    final_scoredf <- scoredf
    final_scoredf$score <- rep(NA, nrow(final_scoredf))
    final_scoredf$decision <- rep("No Trajectory", nrow(final_scoredf))
    final_scoredf <- tibble::rownames_to_column(final_scoredf, "Row.names")
    final_scoredf$note <- rep(NA, nrow(final_scoredf))
    return(final_scoredf)
  }
  if(sum(scoredf$DCcheck & scoredf$SimiRetain>=SimiRetaincutoff)>0) {
    df <- scoredf[scoredf$DCcheck & scoredf$SimiRetain>=SimiRetaincutoff,]
    df_norm <- data.frame(Scaled_GOF = (GOFcutoff-df$GOF)/GOFcutoff, Scaled_USHAPE=(URATEcutoff-df$USHAPE)/URATEcutoff)
    df_norm$score <- rowSums(df_norm)*df$SimiRetain
    df_norm <- round(df_norm, 3)
    rownames(df_norm) <- rownames(df)
    final_scoredf <- merge(scoredf, df_norm, by="row.names", all=T)
    final_scoredf$ranking[!is.na(final_scoredf$score)] <- rank(-final_scoredf$score, na.last = NA)
    final_scoredf$decision <- ifelse(final_scoredf$score>0, "Recommended Embeddings", "Non-recommended Embeddings")
    final_scoredf$decision[is.na(final_scoredf$decision)] <- "very bad"
    final_scoredf$note <- rep(NA, nrow(final_scoredf))
  }else {
    final_scoredf <- scoredf
    final_scoredf$score <- rep(NA, nrow(final_scoredf))
    final_scoredf$decision <- rep("very bad", nrow(final_scoredf))
    final_scoredf <- tibble::rownames_to_column(final_scoredf, "Row.names")
    final_scoredf$note <- rep(NA, nrow(final_scoredf))
  }
  if(sum(final_scoredf$decision=="very bad")>0) {
    dc <- final_scoredf$DCcheck[final_scoredf$decision=="very bad"]
    simi <- final_scoredf$SimiRetain[final_scoredf$decision=="very bad"]
    final_scoredf$note <- rep(NA, nrow(final_scoredf))
    final_scoredf$note[final_scoredf$decision=="very bad"] <- ifelse(!(dc|simi), "Disconnected Clusters & Poor Cell Relationship",
                                                                     ifelse(!dc, "Disconnected Clusters", "Poor Cell Relationship"))
    final_scoredf$decision[final_scoredf$decision=="very bad"] <- "Non-recommended Embeddings"
  }
  final_scoredf <- final_scoredf[order(final_scoredf$score, decreasing = T), ]
  return(final_scoredf)
}

