library(tidyverse)
library(Matrix)


## load data: population age structure, contact matrices, daily accumulated cases,
##            and interventions (start and end dates)

loadPopData <- TRUE
loadContactMatrices <- TRUE
loadCaseData <- TRUE
loadInterventions <- TRUE

County <- "Dublin"

# File paths for Dublin
pop_file_path <- "data/Dublin_pop_2019.csv"
contacts_file_path <- "data/contacts_IRL.Rdata"
cases_file_path <- "data/Covid19CountyStatisticsHPSCIreland.csv"
interventions_file_path <- "data/Dublin_Interventions.csv" ## Interventions adopted in Dublin

# 1) population data
if(loadPopData) 
{ 
  population <- read_csv(pop_file_path)
}

# 2) (projected) contact matrices 
if(loadContactMatrices)
{
  load(contacts_file_path)
  contacts <- contacts_IRL
  rm(contacts_IRL)
}


# 3) population data

if(loadCaseData) 
{ 
  cumulative_cases <- read_csv(cases_file_path) %>% 
      filter(CountyName == County) %>%
      select(TimeStamp,ConfirmedCovidCases) %>% 
      rename(date = TimeStamp , cases = ConfirmedCovidCases) %>%
      mutate(date = as.Date(date)) 
}

# 4) Interventions

if(loadInterventions) 
{ 
  interventions_info <- read_csv(interventions_file_path) %>% 
    mutate(start = as.Date(start, "%d/%m/%y"),end = as.Date(end, "%d/%m/%y")) 
}

rm(loadContactMatrices,loadPopData, loadCaseData, loadInterventions)
