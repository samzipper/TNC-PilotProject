## MODFLOW_HTC_Navarro_SummarizeLeakage.R
#' This script is intended to calculate summarize output from the postprocess_thisRun.py
#' scripts into a single file. 
#' 
#' You need to delete all of the mf*-*.out files from your directly first
#' after running MODFLOW_HTC_CheckFailure.py

## choose stream boundary condition and modflow version
stream_BC <- "RIV"       # "RIV" or "SFR"
modflow_v <- "mfnwt"     # "mfnwt" or "mf2005"
timeType  <- "Transient" # "SteadyState" or "Transient"

## define which directory you are interested in
dir.runs <- file.path("modflow", "HTC", "Navarro", timeType, stream_BC, modflow_v)

## list of all postprocessed files
WellNums <- as.numeric(gsub("mf", "", list.files(dir.runs, "mf")))

## scroll through files and bind
start.flag <- T
for (w in WellNums){
  # load file
  df <- read.csv(file.path(dir.runs, paste0("mf", w), "Navarro-Transient_postprocess.csv"))
  
  # add WellNum
  df$WellNum <- w
  
  # create output data frame
  if (start.flag){
    df.all <- df
    start.flag <- F
  } else {
    df.all <- rbind(df.all, df)
  }
  
  print(paste0(w, " complete"))
}

## save output
write.csv(df.all, file.path(dir.runs, paste0(stream_BC, "-SummarizeLeakage.csv")))
