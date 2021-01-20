## load data: population age structure, contact matrices, daily accumulated cases,
##            and interventions (start and end dates)

load_data <- function(data_path = 'data/', county = 'Dublin', 
  load_pop_data = TRUE, load_contact_matrices = TRUE, load_case_data = TRUE, 
  load_interventions = TRUE) {
  
  # Empty list for output - allows for variable output size
  out <- list()
  
  # Load Population Data
  if (isTRUE(load_pop_data)) {
    pop_file <- paste0(data_path, county, '_pop_2019.csv')
    out[[length(out) + 1]] <- readr::read_csv(pop_file, col_types = cols())
    names(out)[[length(out)]] <- 'population'
  }
  
  # Load (Projected) Contact Matrices  
  if (isTRUE(load_contact_matrices)) {
    contacts_file <- paste0(data_path, 'contacts_IRL.Rdata')
    out[[length(out) + 1]] <- get(load(contacts_file))
    names(out)[[length(out)]] <- 'contacts'
  }
  
  # Load Case Data
  if (isTRUE(load_case_data)) {
    cases_file <- paste0(data_path, 'Covid19CountyStatisticsHPSCIreland.csv')
    out[[length(out) + 1]] <- read_csv(cases_file, col_types = cols()) %>% 
      filter(CountyName == county) %>%
      select(TimeStamp, ConfirmedCovidCases) %>% 
      rename(date = TimeStamp , cases = ConfirmedCovidCases) %>%
      mutate(date = as.Date(date)) 
    names(out)[[length(out)]] <- 'cumulative_cases'
  }
  
  # Load Intervention Data
  if (isTRUE(load_interventions)) { 
    interventions_file <- paste0(data_path, county, '_Interventions.csv')
    out[[length(out) + 1]] <- read_csv(
        interventions_file, col_types = cols()
      ) %>% 
      mutate(
        start = as.Date(start, "%d/%m/%y"),
        end = as.Date(end, "%d/%m/%y")
      ) 
    names(out)[[length(out)]] <- 'interventions_info'
  }
  
  return(out)
}


