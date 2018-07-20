## Navarro_DepletionApportionment_CompareMethods.R

source(file.path("src", "paths+packages.R"))
require(streamDepletr)

## load depletion apportionment output
# AdjacentOnly
df_adj <- 
  file.path("results", "Navarro_DepletionApportionment_AdjacentOnly_AllMethods+Wells+Reaches.csv") %>% 
  read.csv(stringsAsFactors = F) %>% 
  dplyr::select(WellNum, SegNum, f.InvDistSq, f.WebSq, f.TPoly) %>% 
  transform(apportionment = "AdjacentOnly")

# LocalArea
df_loc <- 
  file.path("results", "Navarro_DepletionApportionment_LocalArea_AllMethods+Wells+Reaches.csv") %>% 
  read.csv(stringsAsFactors = F) %>% 
  dplyr::select(WellNum, SegNum, f.InvDistSq, f.WebSq, f.TPoly) %>% 
  transform(apportionment = "LocalArea")

# Dynamic - use 'Transient' glover model at end of simulation
df_dyn <- 
  file.path("results", "Depletion_Analytical_Transient_Dynamic_AllMethods+Wells+Reaches.csv") %>% 
  read.csv(stringsAsFactors = F) %>% 
  subset(Time==max(Time) & analytical=="glover") %>% 
  dplyr::select(WellNum, SegNum, f.InvDistSq, f.WebSq, f.TPoly) %>% 
  transform(apportionment = "Dynamic")

## first: just compare TPoly in all approaches
df_TPoly <- 
  full_join(subset(df_adj[,c("WellNum", "SegNum", "f.TPoly")], f.TPoly != 0),
            subset(df_loc[,c("WellNum", "SegNum", "f.TPoly")], f.TPoly != 0),
            by=c("WellNum", "SegNum")) %>% 
  full_join(subset(df_dyn[,c("WellNum", "SegNum", "f.TPoly")], f.TPoly != 0),
            by=c("WellNum", "SegNum")) %>% 
  set_colnames(c("WellNum", "SegNum", "AdjacentOnly", "LocalArea", "Dynamic")) %>% 
  subset(WellNum %in% df_dyn$WellNum)
df_TPoly[is.na(df_TPoly)] <- 0

ggplot(df_TPoly, aes(x=AdjacentOnly, y=LocalArea)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, color="red")

ggplot(df_TPoly, aes(x=AdjacentOnly, y=Dynamic)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, color="red")

## compare web squared
df_WebSq <- 
  full_join(subset(df_adj[,c("WellNum", "SegNum", "f.WebSq")], f.WebSq != 0),
            subset(df_loc[,c("WellNum", "SegNum", "f.WebSq")], f.WebSq != 0),
            by=c("WellNum", "SegNum")) %>% 
  full_join(subset(df_dyn[,c("WellNum", "SegNum", "f.WebSq")], f.WebSq != 0),
            by=c("WellNum", "SegNum")) %>% 
  set_colnames(c("WellNum", "SegNum", "AdjacentOnly", "LocalArea", "Dynamic")) %>% 
  subset(WellNum %in% df_dyn$WellNum)
df_WebSq[is.na(df_WebSq)] <- 0

ggplot(df_WebSq, aes(x=AdjacentOnly, y=LocalArea)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, color="red")

ggplot(df_WebSq, aes(x=AdjacentOnly, y=Dynamic)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, color="red")
