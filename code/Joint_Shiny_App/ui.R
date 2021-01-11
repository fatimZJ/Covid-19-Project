ui <- navbarPage(
  # Application title
  "Projecting the Effect of Lockdown Measures",
  tabPanel( "Simulation",
            fluidPage(
              titlePanel("Age Structured SEIR Model"),
              sidebarLayout(
                sidebarPanel(
                  
                  # Main Inputs ----
                  selectInput(inputId = "I_type", label = "Type of infected to display:",
                              choices = list("All", "Asymptomatic", "Pre-symptomatic", "Immediate Isolation", 
                                             "Awaiting Test Results", "Isolating After Test", "Not Isolating"),
                              selected = "All"),
                  numericInput(inputId = "R0", label = "Basic Reproduction Number:", min = 0.1, value = 3.4),
                  numericInput(inputId = "L", label = "Latent Period:", min = 0.1, value = 3.8, step = 0.1),
                  numericInput(inputId = "Cv", label = "Incubation Period:", min = 0.1, value = 5.8, step = 0.1),
                  numericInput(inputId = "Dv", label = "Infectious Period:", min = 0.1, value = 13.5, step = 0.1),
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
                            placeholder = "csv file..."),
                ),
                
                # Plot Output ----
                mainPanel(tabsetPanel(
                  tabPanel("Overall", plotlyOutput(outputId = "summary_plot")),
                  tabPanel("Age Breakdown", 
                           fluidRow(
                             column(6, 
                                    plotlyOutput("S_age_plot")),
                             column(6, 
                                    plotlyOutput("E_age_plot")),
                             column(6, 
                                    plotlyOutput("I_age_plot")),
                             column(6, 
                                    plotlyOutput("R_age_plot"))
                           ))
                )
                ))
            )
  ),
  
  tabPanel( "Estimating Intervention Effect",
            fluidPage(
              
              # App title ----
              #  titlePanel("The Effect of Interventions"),
              
              # Sidebar layout with input and output definitions ----
              sidebarLayout(
                
                # Sidebar panel for inputs ----
                sidebarPanel(
                  # Input: Select a file ----
                  
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
                  
                  actionButton("go", "Retrain Model")),
                
                # Main panel for displaying outputs ----
                mainPanel(
                  
                  fluidRow(
                    column(6, 
                           plotlyOutput("PlotCc")%>% withSpinner(color = "#0dc5c1", size = 2)),
                    column(6, 
                           plotlyOutput("PlotRt")%>% withSpinner(color = "#0dc5c1", size = 2))),
                  br(),#br(),#br(),#br(),#br(),br(),
                  DT::dataTableOutput("Estimates")#tableOutput("Estimates")#,
                  #)
                )
                
              )
            ))
  
)
