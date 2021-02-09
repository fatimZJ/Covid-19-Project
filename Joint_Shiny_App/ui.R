sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Info", tabName = "info", icon = icon("info")),
    menuItem("SEIR Model Dashboard", tabName = "seir", icon = icon("dashboard")),
    menuItem("Forecast Settings", tabName = "forecast", icon = icon("wrench")),
    menuItem("Comparison", tabName = "compare", icon = icon("balance-scale-left"))
  )
)

body <- dashboardBody(
  tabItems(
    #Tab 1: Info
    tabItem(tabName = "info"),
    #Tab 2: SEIR Model Output
    tabItem(tabName = "our_stats"),
    #Tab 3: Model Dashboard
    tabItem(tabName = "seir", 
            fluidRow(
              box( title = "Inputs", width = 2, solidHeader = TRUE, status = "primary",
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
                   numericInput(inputId = "Estart", label = "Exposed Starting Value:", min = 0, value = 16),
                   numericInput(inputId = "Istart", label = "Infected Starting Value:", min = 0, value = 1),
                   numericInput(inputId = "Rstart", label = "Removed Starting Value:", min = 0, value = 0),
                   selectInput(inputId = "I_type", label = "Type of infected to display:",
                               choices = list("All", "Symptomatic", "Pre-symptomatic", "Asymptomatic"),
                               selected = "All"),
                   dateRangeInput(inputId = "dates", label = "Dates to Display:", start = as.Date('2020-02-28'), 
                                  end = td, language = "en-GB", format = "dd/mm/yyyy"),
                   checkboxInput(inputId = "shade", label = "Shade Intervention Dates", value = TRUE),
                   checkboxGroupInput(inputId = "age_sel", label = "Select Ages to Display:", 
                                      choices = def_age_groups, selected = def_age_groups)
                   #sliderInput(inputId = "scalar_quantiles", label = "Lockdown Effect Quantiles",
                   #            min = 0, max = 100, value = c(2.5, 97.5), step = 0.5)
              ),
              tabBox( width = 10, title = "SEIR Model Outputs",
                      tabPanel( title = "Dublin",
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
                      tabPanel(title = "Ireland",
                               fluidRow(
                                 column(12,
                                        plotlyOutput("summary_plot_irl") %>% 
                                          withSpinner(color = "#0dc5c1", size = 2, hide.ui = FALSE)),
                                 column(6,
                                        plotlyOutput("S_age_plot_irl") %>% 
                                          withSpinner(color = "#0dc5c1", size = 1, hide.ui = FALSE)
                                 ),
                                 column(6, 
                                        plotlyOutput("E_age_plot_irl") %>% 
                                          withSpinner(color = "#0dc5c1", size = 1, hide.ui = FALSE)
                                 ),
                                 column(6, 
                                        plotlyOutput("I_age_plot_irl") %>% 
                                          withSpinner(color = "#0dc5c1", size = 1, hide.ui = FALSE)
                                 ),
                                 column(6, 
                                        plotlyOutput("R_age_plot_irl") %>% 
                                          withSpinner(color = "#0dc5c1", size = 1, hide.ui = FALSE)
                                 )
                               )
                      )
              
              )
              
    )),
    #Tab 4: Forecast Settings
    tabItem(tabName = "forecast",
            fluidRow(
              tabBox( title = "Forecast", width = 10,
                      tabPanel(title = "Dublin",
                              plotlyOutput("forecast_plot_dub")),
                      tabPanel(title = "Ireland",
                              plotlyOutput("forecast_plot_irl"))
              ),
              box( title = "Expected Outcome", width = 2, solidHeader = TRUE, status = "danger",
                  h2("Expected Deaths and Costs here"))
              ),
            box( title = "Display", width = 2, solidHeader = TRUE, status = "warning",
                 selectInput(inputId = "Disp_comp", label = "Compartment:",
                             choices = list("Susceptible", "Exposed", "Infected: All", 
                                            "Infected: Asymptomatic", "Infected: Pre-symptomatic", 
                                            "Infected: Symptomatic", "Removed"),
                             selected = "Infected:All"),
                 dateInput(inputId = "start_date", label = "Start Date:", value = td,
                           min = as.Date('2020-02-28'), max = td)
                 ),
            box( title = "Restriction Weeks 1 & 2:", width = 2, solidHeader = TRUE, status = "primary",
                 selectInput(inputId = "res1", label = "Restriction:",
                             choices = interventions_info$policy,
                             selected = "No Intervention")
                 ),
            box( title = "Restriction Weeks 3 & 4:", width = 2, solidHeader = TRUE, status = "primary",
                 selectInput(inputId = "res2", label = "Restriction:",
                             choices = interventions_info$policy,
                             selected = "No Intervention")
            ),
            box( title = "Restriction Weeks 5 & 6:", width = 2, solidHeader = TRUE, status = "primary",
                 selectInput(inputId = "res3", label = "Restriction:",
                             choices = interventions_info$policy,
                             selected = "No Intervention")
            ),
            box( title = "Restriction Weeks 7 & 8:", width = 2, solidHeader = TRUE, status = "primary",
                 selectInput(inputId = "res4", label = "Restriction:",
                             choices = interventions_info$policy,
                             selected = "No Intervention")
            )
            ),
    
    #Tab 5: Comparison
    tabItem(tabName = "compare")
  )
)

# Put them together into a dashboardPage
dashboardPage( skin = "blue",
  dashboardHeader(title = "Ireland COVID-19 Incidence Modelling",
                  titleWidth = 380),
  sidebar,
  body
)

