##### Shiny Application

library(shiny)
library(gridExtra)
source('code/Shiny_Function_Source.R')

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  titlePanel("Age Structured SEIR Model"),
  
  sidebarLayout(
    
    sidebarPanel(
      
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
      
      numericInput(inputId = "TT", label = "Average Test Result Time:", min = 0.1, value = 2, step = 0.1),
      
      checkboxGroupInput(inputId = "age_groups", label = "Choose Age Groups:",
                         choices = as.list(1:16), selected = 1:16),
      
    ),
    
    mainPanel(
      
      plotOutput(outputId = "outplot")
      
    )
  )
)

server <- function(input, output) {
  
  output$outplot <- renderPlot({
    
    sol <- simulation_SEIR_model(pars = c(input$L, input$Cv, input$Dv, input$h, input$i, 
                                          input$j, input$f, input$tv, input$q, input$TT))$sol_out
    
    S <- sol[grepl('S_',names(sol))]
    Ev <- sol[grepl('Ev_',names(sol))]
    Ip <- sol[grepl('Ip_',names(sol))]
    IA <- sol[grepl('IA_',names(sol))]
    Ii <- sol[grepl('Ii_',names(sol))]
    It <- sol[grepl('It_',names(sol))]
    Iti <- sol[grepl('Iti_',names(sol))]
    Iq <- sol[grepl('Iq_',names(sol))]
    I <- cbind(Ip, IA, Ii,It,Iti,Iq)
    R <- sol[grepl('R_',names(sol))]
    xval <- sol$time
    
    S_plot <- ggplot(data = NULL, aes(x=xval, y=rowSums(S))) + 
      geom_line() +
      labs(x="Time (days)", y = "Daily no. of susceptibles") +
      theme_light()
    
    E_plot <- ggplot(data = NULL, aes(x=xval, y=rowSums(Ev))) + 
      geom_line() +
      labs(x="Time (days)", y = "Daily no. of exposed") +
      theme_light()
    
    I_plot <- ggplot(data = NULL, aes(x=xval, y=rowSums(I))) + 
      geom_line() +
      labs(x="Time (days)", y = "Daily no. of infections") +
      theme_light()
    
    R_plot <- ggplot(data = NULL, aes(x=xval, y=rowSums(R))) + 
      geom_line() +
      labs(x="Time (days)", y = "Daily no. of removed") +
      theme_light()
    
    grid.arrange(grobs=list(S_plot, E_plot, I_plot, R_plot), nrow=4)
    
  }, height = 1000)
  
}

shinyApp(ui, server)
