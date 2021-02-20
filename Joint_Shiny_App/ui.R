sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Info", tabName = "info", icon = icon("info")),
    menuItem("SEIR Model Dashboard", tabName = "seir", icon = icon("dashboard")),
    menuItem("Forecast Settings", tabName = "forecast", icon = icon("wrench"))
    #menuItem("Comparison", tabName = "compare", icon = icon("balance-scale-left"))
  )
)

body <- dashboardBody(
  tabItems(
    #Tab 1: Info
    tabItem(tabName = "info",
            box(title = "1. About This Application", solidHeader = TRUE,
                status = "primary", width = 12,
              textOutput("General_Intro"),
              uiOutput("pp_link")
            ),
            box(title = "2. Model Dashboard Tab", solidHeader = TRUE,
                status = "primary", width = 12,
              textOutput("Mod_Dash_Tab"),
            ),
            box(title = "3. Model Forecast Tab", solidHeader = TRUE,
                status = "primary", width = 12,
              textOutput("Mod_Fore_Tab"),
            ),
            
            box(title = "4. Acknowledgement", solidHeader = TRUE,
                status = "primary", width = 12,
              textOutput("Acknowledgement"),
              uiOutput("sfi_link")
            )
    ),
    #Tab 2: Model Dashboard
    tabItem(tabName = "seir", 
            fluidRow(
              box( title = "Inputs", width = 2, solidHeader = TRUE, status = "warning",
                   numericInput(inputId = "R0", label = "Basic Reproduction Number:", min = 0.1, value = 3.7),
                   numericInput(inputId = "L", label = "Latent Period:", min = 0.1, value = 3.7, step = 0.1),
                   numericInput(inputId = "Cv", label = "Incubation Period:", min = 0.1, value = 5.8, step = 0.1),
                   numericInput(inputId = "Dv", label = "Infectious Period:", min = 0.1, value = 11.6, step = 0.1),
                   numericInput(inputId = "h", label = "Asymptomatic Transmission:", min = 0.01, max = 1, value = 0.55,
                                step = 0.01),
                   numericInput(inputId = "k", label = "Isolated Transmission:", min = 0.01, max = 1, value = 0.05,
                                step = 0.01),
                   numericInput(inputId = "f", label = "Proportion of Asymptomatic:", min = 0.01, max = 1, value = 0.2,
                                step = 0.01),
                   numericInput(inputId = "tv", label = "Proportion of Tested:", min = 0.01, max = 1, value = 0.8,
                                step = 0.01), 
                   numericInput(inputId = "q", label = "Proportion of Isolated:", min = 0.01, max = 1, value = 0.1,
                                step = 0.01),
                   numericInput(inputId = "TT", label = "Average Test Result Time:", min = 0.1, value = 7, step = 0.1),
                   tags$hr(),
                   selectInput(inputId = "I_type", label = "Type of infected to display:",
                               choices = list("All", "Symptomatic", "Pre-symptomatic", "Asymptomatic"),
                               selected = "All"),
                   dateInput(inputId = "end_date", label = "End Date:", value = td,
                           min = as.Date('2020-02-29'), max = td),
                   checkboxInput(inputId = "shade", label = "Shade Intervention Dates", value = TRUE),
                   checkboxGroupInput(inputId = "age_sel", label = "Select Ages to Display:", 
                                      choices = def_age_groups, selected = def_age_groups)
              ),
              box( width = 10, title = "Dublin SEIR Model Outputs", solidHeader = TRUE,
                   status = "primary",
                                fluidRow(
                                  column(12,
                                         plotlyOutput("summary_plot_dub") %>% 
                                           withSpinner(color = "#0dc5c1", size = 2, hide.ui = FALSE)),
                                  column(6,
                                         plotlyOutput("S_age_plot_dub") %>% 
                                           withSpinner(color = "#0dc5c1", size = 1, hide.ui = FALSE)
                                  ),
                                  column(6, 
                                         plotlyOutput("E_age_plot_dub") %>% 
                                           withSpinner(color = "#0dc5c1", size = 1, hide.ui = FALSE)
                                  ),
                                  column(6, 
                                         plotlyOutput("I_age_plot_dub") %>% 
                                           withSpinner(color = "#0dc5c1", size = 1, hide.ui = FALSE)
                                  ),
                                  column(6, 
                                         plotlyOutput("R_age_plot_dub") %>% 
                                           withSpinner(color = "#0dc5c1", size = 1, hide.ui = FALSE)
                                  )
                                )
                      ),
                      
                      )
              
    ),
    #Tab 3: Forecast Settings
    tabItem(tabName = "forecast",
            fluidRow(
              box( title = "Dublin Forecast", width = 10, solidHeader = TRUE,
                   status = "primary",
                   plotlyOutput("forecast_plot_dub") %>%
                                withSpinner(color = "#0dc5c1", size = 2, hide.ui = FALSE)),
              box( title = "Display", width = 2, solidHeader = TRUE, status = "warning",
                 selectInput(inputId = "Disp_comp", label = "Compartment:",
                             choices = list("Susceptible", "Exposed", "All Infected", 
                                            "Symptomatic Infected", "Pre-symptomatic Infected", 
                                            "Asymptomatic Infected", "Removed"),
                             selected = "Infected:All"),
                 dateInput(inputId = "start_date", label = "Start Date:", 
                           value = td, min = as.Date('2020-02-29'), max = td),
                 checkboxGroupInput(inputId = "forecast_age_sel", label = "Age Groups:", 
                                      choices = forecast_age_groups, selected = forecast_age_groups)
                 )
            ),
            fluidRow(
              infoBox("", "Select lockdown levels:", 
                      icon = icon("chevron-right"), width = 2, fill = TRUE),
              box( title = "Restriction Weeks 1 & 2:", width = 2, solidHeader = TRUE, status = "info",
                   selectInput(inputId = "res1", label = "Restriction:",
                               choices = lockdown_measures),
                   textOutput("Cost_12")
                   ),
              box( title = "Restriction Weeks 3 & 4:", width = 2, solidHeader = TRUE, status = "info",
                   selectInput(inputId = "res2", label = "Restriction:",
                               choices = lockdown_measures),
                   textOutput("Cost_34")
              ),
              box( title = "Restriction Weeks 5 & 6:", width = 2, solidHeader = TRUE, status = "info",
                   selectInput(inputId = "res3", label = "Restriction:",
                               choices = lockdown_measures),
                   textOutput("Cost_56")
              ),
              box( title = "Restriction Weeks 7 & 8:", width = 2, solidHeader = TRUE, status = "info",
                   selectInput(inputId = "res4", label = "Restriction:",
                               choices = lockdown_measures),
                   textOutput("Cost_78")
              ),
              box( title = "Fit Choices", width = 2, solidHeader = TRUE, status = "success", 
                   actionButton("fit_forecast", label = "Create Forecast", icon = icon("play-circle"))
              )
            ),
            fluidRow(
              box( title = "Expected Deaths by Age Group", width = 10, solidHeader = TRUE,
                   status = "danger",
                   plotlyOutput("forecast_exp_deaths") %>%
                                withSpinner(color = "#0dc5c1", size = 2, hide.ui = FALSE)),
              box( title = "Expected Deaths", width = 2, solidHeader = TRUE, status = "danger",
                     tableOutput('deaths') ),
              infoBoxOutput('TotalCostBox', width = 2)
            )
    )
    
    #Tab 4: Comparison
    #tabItem(tabName = "compare")
  )
)

# Put them together into a dashboardPage
dashboardPage( skin = "blue",
  dashboardHeader(title = "Ireland COVID-19 Incidence Modelling",
                  titleWidth = 380),
  sidebar,
  body
)

