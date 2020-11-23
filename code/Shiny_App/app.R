##### Shiny Application
# https://skinsella.shinyapps.io/euprocesses/?_ga=2.42999395.711233149.1605794510-294588088.1605794510

library(shiny)
library(plotly)
source('code/Shiny_Function_Source.R')

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  titlePanel("Age Structured SEIR Model"),
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId = "I_type", label = "Type of infected to display:",
                  choices = list("All", "Asymptomatic", "Pre-symptomatic", "Immediate Isolation", 
                                 "Awaiting Test Results", "Isolating After Test", "Not Isolating"),
                  selected = "All"),
      numericInput(inputId = "L", label = "Latent Period:", min = 0.1, value = 3.8, step = 0.1),
      numericInput(inputId = "Cv", label = "Incubation Period:", min = 0.1, value = 5.8, step = 0.1),
      numericInput(inputId = "Dv", label = "Infectious Period:", min = 0.1, value = 13.5, step = 0.1),
      numericInput(inputId = "h", label = "Asymptomatic Transmission:", min = 0.01, max = 1, value = 0.55,
                   step = 0.01),
      numericInput(inputId = "i", label = "Isolation 1:", min = 0.01, max = 1, value = 0.05, step = 0.01),
      numericInput(inputId = "j", label = "Isolation 2:", min = 0.01, max = 1, value = 0.05, step = 0.01),
      numericInput(inputId = "f", label = "Proportion of Asymptomatic:", min = 0.01, max = 1, value = 0.2,
                   step = 0.01),
      numericInput(inputId = "tv", label = "Proportion of Tested:", min = 0.01, max = 1, value = 0.8,
                   step = 0.01),
      numericInput(inputId = "q", label = "Proportion of Isolated:", min = 0.01, max = 1, value = 0.1,
                   step = 0.01),
      numericInput(inputId = "TT", label = "Average Test Result Time:", min = 0.1, value = 2, step = 0.1)
    ),
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

server <- function(input, output) {
  
  sol <- reactive({ 
    simulation_SEIR_model(pars = c(input$L, input$Cv, input$Dv, input$h, input$i, 
                                   input$j, input$f, input$tv, input$q, input$TT))$sol_out 
  })
  
  output$summary_plot <- renderPlotly({
    
    xval <- sol()$time
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
    plot_ly(x = ~xval, y = ~rowSums(S), name = 'Susceptible', type = 'scatter', mode = 'lines') %>% 
      add_trace(y = ~rowSums(E), name = 'Exposed', mode = 'lines') %>% 
      add_trace(y = ~rowSums(I), name = 'Infected', mode = 'lines') %>%
      add_trace(y = ~rowSums(R), name = 'Recovered', mode = 'lines') %>% 
      layout(xaxis = list(title = "Time"), yaxis = list(title = "Compartment Size"))
    
  })
  
  output$S_age_plot <- renderPlotly({
    
    xval <- sol()$time
    S <- sol()[grepl('S_',names(sol()))]
    
    ### Draw the Plot
    p_S <- plot_ly(x = ~xval, y = ~S[[1]], name = 'Age 1', type = 'scatter', mode = 'lines')
    
    for ( i in 2:16 ) {
      p_S <- p_S %>% add_trace(y = S[[i]], name = paste0('Age ', i), mode = 'lines')
    }
    
    p_S %>% layout(xaxis = list(title = "Time"), yaxis = list(title = "Susceptible"))
    
  })
  
  output$E_age_plot <- renderPlotly({
  
    xval <- sol()$time
    E <- sol()[grepl('Ev_',names(sol()))]
    
    ### Draw the Plot
    p_E <- plot_ly(x = ~xval, y = ~E[[1]], name = 'Age 1', type = 'scatter', mode = 'lines')
    
    for ( i in 2:16 ) {
      p_E <- p_E %>% add_trace(y = E[[i]], name = paste0('Age ', i), mode = 'lines')
    }
    
    p_E %>% layout(xaxis = list(title = "Time"), yaxis = list(title = "Exposed"))
    
  })
  
  output$I_age_plot <- renderPlotly({
    
    xval <- sol()$time
    I_Pr <- sol()[grepl('Ip_',names(sol()))]
    I_As <- sol()[grepl('IA_',names(sol()))]
    I_Im <- sol()[grepl('Ii_',names(sol()))]
    I_Aw <- sol()[grepl('It_',names(sol()))]
    I_Is <- sol()[grepl('Iti_',names(sol()))]
    I_No <- sol()[grepl('Iq_',names(sol()))]
    I_Al <- I_Pr + I_As + I_Im + I_Aw + I_Is + I_No
    I <- get( paste0("I_", substr(input$I_type, start = 1, stop = 2)) )
    
    ### Draw the Plot
    p_I <- plot_ly(x = ~xval, y = ~I[[1]], name = 'Age 1', type = 'scatter', mode = 'lines')
    
    for ( i in 2:16 ) {
      p_I <- p_I %>% add_trace(y = I[[i]], name = paste0('Age ', i), mode = 'lines')
    }
    
    p_I %>% layout(xaxis = list(title = "Time"), yaxis = list(title = paste0("Infected: ", input$I_type)))
    
  })
  
  output$R_age_plot <- renderPlotly({
    
    xval <- sol()$time
    R <- sol()[grepl('R_',names(sol()))]
    
    ### Draw the Plot
    p_R <- plot_ly(x = ~xval, y = ~R[[1]], name = 'Age 1', type = 'scatter', mode = 'lines')
    
    for ( i in 2:16 ) {
      p_R <- p_R %>% add_trace(y = R[[i]], name = paste0('Age ', i), mode = 'lines')
    }
    
    p_R %>% layout(xaxis = list(title = "Time"), yaxis = list(title = "Removed"))
    
  })
  
}

shinyApp(ui, server)
