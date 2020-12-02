##### Shiny Application

library(shiny)
library(plotly)
library(readr)
source('code/Shiny_Function_Source.R')

# Define UI ----
ui <- fluidPage(
  titlePanel("Age Structured SEIR Model (WORK IN PROGRESS)"),
  sidebarLayout(
    sidebarPanel(
      
      # Main Inputs ----
      selectInput(inputId = "I_type", label = "Type of infected to display:",
                  choices = list("All", "Asymptomatic", "Pre-symptomatic", "Immediate Isolation", 
                                 "Awaiting Test Results", "Isolating After Test", "Not Isolating"),
                  selected = "All"),
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
      numericInput(inputId = "tmax", label = "How Many Days to Project Ahead:", min = 1, value = 225, step = 1),
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
      tabPanel("Age_Breakdown", 
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
    ))
  )
)

# Define server ----
server <- function(input, output) {
  
  # Process Contact Matrices Input ----
  input_h_CM <- reactive({
    if (is.null(input$House_CM)) { return(contacts$home) }
    as.matrix( read.csv(input$House_CM$datapath) )
  })
  
  input_w_CM <- reactive({
    if (is.null(input$Work_CM)) { return(contacts$work) }
    as.matrix( read.csv(input$Work_CM$datapath) )
  })
  
  input_s_CM <- reactive({
    if (is.null(input$School_CM)) { return(contacts$school) }
    as.matrix( read.csv(input$School_CM$datapath) )
  })
  
  input_o_CM <- reactive({
    if (is.null(input$Other_CM)) { return(contacts$others) }
    as.matrix( read.csv(input$Other_CM$datapath) )
  })
  
  input_CM <- reactive({
    all <- input_h_CM() + input_w_CM() + input_s_CM() + input_o_CM()
    list(home = input_h_CM(), work = input_w_CM(), school = input_s_CM(), 
         others = input_o_CM(), all = all)
  })
  
  # Process Lockdown Info ----
  linfo <- reactive({
    if (is.null(input$Lockdown_Info)) { return(input$Lockdown_Info) }
    read_csv(input$Lockdown_Info$datapath)
  })
  
  # Simulate SEIR Model ----
  sol <- reactive({
    simulation_SEIR_model(pars = c(input$L, input$Cv, input$Dv, input$h, input$i, 
                                   input$j, input$f, input$tv, input$q, input$TT),
                          contacts_ireland = input_CM(),
                          tmax = input$tmax, 
                          lockdown_information = linfo())$sol_out 
  })
  
  # Define Shaded Regions for Plots ----
  shaded_regions <- reactive({
    
    if( is.null(input$Lockdown_Info) ) { return(NULL) }
    
    max_val <- rowSums(sol()[1, ]) * 1.05
    scales <- linfo()[[3]]
    N <- nrow(linfo())
    shaded_list <- vector("list", length = N)
    
    for (i in 1:N) {
      
      shaded_list[[i]] <- list(type = "rect",
                               fillcolor = "gray", line = list(color = "gray"), 
                               opacity = 0.5 * (1 - scales[i]), x0 = linfo()[[1]][i], 
                               x1 = linfo()[[2]][i], xref = "x",
                               y0 = 0, y1 = max_val, yref = "y")
      
    }
    
    shaded_list
    
  })
  
  # Create Plots ----
  xval <- reactive({
    as.Date('2020-02-28') + sol()$time
  })
  
  output$summary_plot <- renderPlotly({
    
    # Extract Compartments
    S <- sol()[grepl('S_',names(sol()))]
    E <- sol()[grepl('Ev_',names(sol()))]
    I_Pr <- sol()[grepl('Ip_',names(sol()))]
    I_As <- sol()[grepl('IA_',names(sol()))]
    I_Im <- sol()[grepl('Ii_',names(sol()))]
    I_Aw <- sol()[grepl('It_',names(sol()))]
    I_Is <- sol()[grepl('Iti_',names(sol()))]
    I_No <- sol()[grepl('Iq_',names(sol()))]
    I_Al <- I_Pr + I_As + I_Im + I_Aw + I_Is + I_No
    I <- get( paste0("I_", substr(input$I_type, start = 1, stop = 2)) )
    R <- sol()[grepl('R_',names(sol()))]
    
    ### Draw the Plot
    plot_ly(x = ~xval(), y = ~rowSums(S), name = 'Susceptible', type = 'scatter', mode = 'lines') %>%
      add_trace(y = ~rowSums(E), name = 'Exposed', mode = 'lines') %>% 
      add_trace(y = ~rowSums(I), name = 'Infected', mode = 'lines') %>%
      add_trace(y = ~rowSums(R), name = 'Recovered', mode = 'lines') %>% 
      layout(shapes = shaded_regions()) %>%
      layout(xaxis = list(title = "Time"), 
             yaxis = list(title = "Compartment Size", 
                          range = c(0, rowSums(sol()[1, ]) * 1.01)))
    
  })
  
  output$S_age_plot <- renderPlotly({
    
    S <- sol()[grepl('S_',names(sol()))]
    
    ### Draw the Plot
    p_S <- plot_ly(x = ~xval(), y = ~S[[1]], name = 'Age 1', type = 'scatter', mode = 'lines')
    
    for ( i in 2:16 ) {
      p_S <- p_S %>% add_trace(y = S[[i]], name = paste0('Age ', i), mode = 'lines')
    }
    
    p_S %>% layout(shapes = shaded_regions()) %>%
      layout(xaxis = list(title = "Time"), 
             yaxis = list(title = "Susceptible", 
                          range = c(0, max(S) * 1.01)))
    
  })
  
  output$E_age_plot <- renderPlotly({
  
    E <- sol()[grepl('Ev_',names(sol()))]
    
    ### Draw the Plot
    p_E <- plot_ly(x = ~xval(), y = ~E[[1]], name = 'Age 1', type = 'scatter', mode = 'lines')
    
    for ( i in 2:16 ) {
      p_E <- p_E %>% add_trace(y = E[[i]], name = paste0('Age ', i), mode = 'lines')
    }
    
    p_E %>% layout(shapes = shaded_regions()) %>%
      layout(xaxis = list(title = "Time"), 
             yaxis = list(title = "Exposed", 
                          range = c(0, max(E) * 1.01)))
    
  })
  
  output$I_age_plot <- renderPlotly({
    
    I_Pr <- sol()[grepl('Ip_',names(sol()))]
    I_As <- sol()[grepl('IA_',names(sol()))]
    I_Im <- sol()[grepl('Ii_',names(sol()))]
    I_Aw <- sol()[grepl('It_',names(sol()))]
    I_Is <- sol()[grepl('Iti_',names(sol()))]
    I_No <- sol()[grepl('Iq_',names(sol()))]
    I_Al <- I_Pr + I_As + I_Im + I_Aw + I_Is + I_No
    I <- get( paste0("I_", substr(input$I_type, start = 1, stop = 2)) )
    
    ### Draw the Plot
    p_I <- plot_ly(x = ~xval(), y = ~I[[1]], name = 'Age 1', type = 'scatter', mode = 'lines')
    
    for ( i in 2:16 ) {
      p_I <- p_I %>% add_trace(y = I[[i]], name = paste0('Age ', i), mode = 'lines')
    }
    
    p_I %>% layout(shapes = shaded_regions()) %>%
      layout(xaxis = list(title = "Time"), 
             yaxis = list(title = paste0("Infected: ", input$I_type), 
                          range = c(0, max(I) * 1.01)))
    
  })
  
  output$R_age_plot <- renderPlotly({
    
    R <- sol()[grepl('R_',names(sol()))]
    
    ### Draw the Plot
    p_R <- plot_ly(x = ~xval(), y = ~R[[1]], name = 'Age 1', type = 'scatter', mode = 'lines')
    
    for ( i in 2:16 ) {
      p_R <- p_R %>% add_trace(y = R[[i]], name = paste0('Age ', i), mode = 'lines')
    }
    
    p_R %>% layout(shapes = shaded_regions()) %>%
      layout(xaxis = list(title = "Time"), 
             yaxis = list(title = "Removed", 
                          range = c(0, max(R) * 1.01)))
    
  })
  
}

# Run App ----
shinyApp(ui, server)
