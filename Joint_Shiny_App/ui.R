sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Summary", tabName = "summary", icon = icon("dashboard"),
             selectInput(inputId = "I_type", label = "Type of infected to display:",
                         choices = list("All", "Asymptomatic", "Pre-symptomatic", "Immediate Isolation", 
                                        "Awaiting Test Results", "Isolating After Test", "Not Isolating"),
                         selected = "All"),
             numericInput(inputId = "R0", label = "Basic Reproduction Number:", min = 0.1, value = 3.7),
             numericInput(inputId = "L", label = "Latent Period:", min = 0.1, value = 3.7, step = 0.1),
             numericInput(inputId = "Cv", label = "Incubation Period:", min = 0.1, value = 5.8, step = 0.1),
             numericInput(inputId = "Dv", label = "Infectious Period:", min = 0.1, value = 11.6, step = 0.1),
             numericInput(inputId = "h", label = "Asymptomatic Transmission:", min = 0.01, max = 1, value = 0.55,
                          step = 0.01),
             numericInput(inputId = "f", label = "Proportion of Asymptomatic:", min = 0.01, max = 1, value = 0.2,
                          step = 0.01),
             numericInput(inputId = "tv", label = "Proportion of Tested:", min = 0.01, max = 1, value = 0.8,
                          step = 0.01), 
             numericInput(inputId = "q", label = "Proportion of Isolated:", min = 0.01, max = 1, value = 0.1,
                          step = 0.01),
             numericInput(inputId = "TT", label = "Average Test Result Time:", min = 0.1, value = 2, step = 0.1),
             numericInput(inputId = "Estart", label = "Exposed Starting Value:", min = 0, value = 16),
             numericInput(inputId = "Istart", label = "Infected Starting Value:", min = 0, value = 1),
             numericInput(inputId = "Rstart", label = "Recovered Starting Value:", min = 0, value = 0),
             checkboxGroupInput(inputId = "age_sel", label = "Select Ages to Display:", choices = as.character(1:16),
                                selected = as.character(1:16)),
             dateRangeInput(inputId = "dates", label = "Time Horizon:", start = as.Date('2020-02-28'), 
                            end = as.Date('2020-10-01'), language = "en-GB", format = "dd/mm/yyyy"),
             fileInput(inputId = "pop_file", label = "Population File:", multiple = FALSE, accept = ".csv",
                       placeholder = "csv file..."),
             fileInput(inputId = "House_CM", label = "Household Contact Matrix Input:", multiple = FALSE, accept = ".csv",
                       placeholder = "csv file..."),
             fileInput(inputId = "Work_CM", label = "Work Contact Matrix Input:", multiple = FALSE, accept = ".csv",
                       placeholder = "csv file..."),
             fileInput(inputId = "School_CM", label = "School Contact Matrix Input:", multiple = FALSE, accept = ".csv",
                       placeholder = "csv file..."),
             fileInput(inputId = "Other_CM", label = "Other Contact Matrix Input:", multiple = FALSE, accept = ".csv",
                       placeholder = "csv file..."),
             fileInput(inputId = "Lockdown_Info", label = "Lockdown Scenario Input:", multiple = FALSE, accept = ".csv",
                       placeholder = "csv file...")
             ),
    
    menuItem("Intervention Effect", tabName = "ie", icon = icon("dashboard"),
             fileInput(inputId = "PopData", label = "Population Structure", #"Choose CSV File",
                       multiple = FALSE,
                       accept = c("text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv")),
             fileInput(inputId = "CountData", label = "Cumulative Number of Cases per Day", #"Choose CSV File",
                       multiple = FALSE,
                       accept = c("text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv")),
             
             fileInput(inputId = "InterInfo", label = "Interventions Start and End dates", #"Choose CSV File",
                       multiple = FALSE,
                       accept = c("text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv")),
             
             fileInput(inputId = "House_CM", label = "Household Contact Matrix Input", 
                       multiple = FALSE, accept = ".csv"),
             # placeholder = "csv file..."),
             fileInput(inputId = "Work_CM", label = "Work Contact Matrix Input", 
                       multiple = FALSE, accept = ".csv"),
             #placeholder = "csv file..."),
             fileInput(inputId = "School_CM", label = "School Contact Matrix Input", 
                       multiple = FALSE, accept = ".csv"),
             #placeholder = "csv file..."),
             fileInput(inputId = "Other_CM", label = "Other Contact Matrix Input", 
                       multiple = FALSE, accept = ".csv"),
             # placeholder = "csv file..."),
             
             # Horizontal line ----
             tags$hr(),
             
             actionButton("go", "Retrain Model")
    )
  )
)

dashboardBody(
  #Tab 1 Overall View
  tabItems(
    tabItem(tabName = "summary", 
            h2("Overall Objectives"),
            fluidRow(
              box(
                title = "Overall", width = 6, solidHeader = TRUE,
                plotlyOutput("summary_plot")
              ),
              box(
                title = "Cost Category Breakdown", width = 6, solidHeader = TRUE,
                fluidRow(
                  column(6, 
                         plotlyOutput("S_age_plot")),
                  column(6, 
                         plotlyOutput("E_age_plot")),
                  column(6, 
                         plotlyOutput("I_age_plot")),
                  column(6, 
                         plotlyOutput("R_age_plot"))
                )
              )
            )
    ),
    #Tab 2 Wellpad Decision
    tabItem(tabName = "ie", 
            h2("Wellpad Operator"),
            fluidRow(
              column(6, 
                     plotlyOutput("PlotCc")%>% withSpinner(color = "#0dc5c1", size = 2)),
              column(6, 
                     plotlyOutput("PlotRt")%>% withSpinner(color = "#0dc5c1", size = 2))
            ) 
            
    )
  )
)

# Put them together into a dashboardPage
dashboardPage(
  dashboardHeader(title = "Simple tabs"),
  sidebar,
  body
)

