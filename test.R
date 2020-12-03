
#### Load required packages ------------------------------------------------####

list.of.packages <- c("shiny", "mlbench", "plotly", "shinythemes", "dplyr", 
                      "shinyWidgets", "ggplot2", "shinydashboard", "data.table",
                      "optimx", "shinycssloaders", "DT")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(new.packages) 

lapply(as.list(list.of.packages), require, character.only = TRUE)

#### Load required source files --------------------------------------------####
source("code/getbeta.R")
source("Code/function_IEMAGSEIR_lsoda_scalars.R")
source('code/1_loadData.r')

# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("The Effect of Lockdown Measures"),
  
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
      #fluidRow(
      # Output: Data file ----
      #tableOutput("contents")
      plotlyOutput("Plot", width = "100%")%>% withSpinner(color = "#0dc5c1", size = 2.5),
      br(),br(),br(),br(),br(),br(),
      
      DT::dataTableOutput("Estimates")#tableOutput("Estimates")#,
      #)
    )
    
  )
)

# Define server logic to read selected file ----
server <- function(input, output) {
  
  count_data <- reactive({
    req(input$CountData)
    
    read.csv(input$CountData$datapath, sep = ",")
    
  })
  
  pop_data <- reactive({
    req(input$PopData)
    
    read.csv(input$PopData$datapath, sep = ",")
    
  })
  
  
  Ests <- eventReactive(input$go, {
    
    req(count_data())
    req(pop_data())
    
    nlm_fun <- function(scalars, Acc_Cases, Pop_Struc) {
      
      scalars <- exp(scalars)  
      
      test <- simulation_SEIR_model(R0t = 3.4,
                                    dateStartSchoolClosure = as.Date('2020-03-12') , #schools closed before imtense lockdown
                                    dateStartIntenseIntervention = as.Date('2020-03-27') , #Intense intervention: starts at Wuhan Lockdown
                                    dateEndIntenseIntervention = as.Date('2020-05-18'), #date we begin relaxing intense intervention
                                    dateStart = as.Date('2020-02-28'), #start date for epidemic in Ireland
                                    POP = Pop_Struc,
                                    numWeekStagger = c(3,6,9,12,15),
                                    contacts_ireland = contacts,
                                    dt = 1,  
                                    tmax = 225,#length(Acc_Cases),
                                    scalars = scalars)
      
      test_Cc <- rowSums(test$sol_out[grepl('Cc_',names(test$sol_out))])
      
      sum((Acc_Cases - (test_Cc[-1]))^2) 
    }
    
   
    ests <- nlm(nlm_fun,log( c(2.21000000, 0.97171108, 0.20812880, 0.07611468, 0.22169280, 
                               0.40394184, 0.44685362, 0.36467352, 0.42729580)),#c(1.91365769, 1.16790530, 0.20215852, 0.06962944, 
                             # 0.18328977, 0.43466483, 0.60612880, 0.32053803, 
                              #0.45793943)), 
      Acc_Cases = count_data()$cases[1:225],
                Pop_Struc = pop_data(),
                stepmax = 0.5, iterlim = 1)$est
    
    ests
    
  })
  
  Sim <- reactive({
    req(Ests())
    
    sim <- simulation_SEIR_model(R0t = 3.4,
                                 dateStartSchoolClosure = as.Date('2020-03-12') , #schools closed before imtense lockdown
                                 dateStartIntenseIntervention = as.Date('2020-03-27') , #Intense intervention: starts at Wuhan Lockdown
                                 dateEndIntenseIntervention = as.Date('2020-05-18'), #date we begin relaxing intense intervention
                                 dateStart = as.Date('2020-02-28'), #start date for epidemic in Ireland
                                 POP = pop_data(),
                                 numWeekStagger = c(3,6,9,12,15),
                                 contacts_ireland = contacts,
                                 dt = 1,#0.1,  
                                 tmax = 225,
                                 scalars = exp(Ests()))
    
    S <- sim$sol_out[grepl('S_',names(sim$sol_out))]
    Ev <- sim$sol_out[grepl('Ev_',names(sim$sol_out))]
    Ip <- sim$sol_out[grepl('Ip_',names(sim$sol_out))]
    IA <- sim$sol_out[grepl('IA_',names(sim$sol_out))]
    Ii <- sim$sol_out[grepl('Ii_',names(sim$sol_out))]
    It <- sim$sol_out[grepl('It_',names(sim$sol_out))]
    Iti <- sim$sol_out[grepl('Iti_',names(sim$sol_out))]
    Iq <- sim$sol_out[grepl('Iq_',names(sim$sol_out))]
    R <- sim$sol_out[grepl('R_',names(sim$sol_out))]
    Cc <- sim$sol_out[grepl('Cc_',names(sim$sol_out))]
    
    as.data.frame(cbind("time" = sim$sol_out$time, "S" = rowSums(S), 
                        "Ev" = rowSums(Ev), "Ip" = rowSums(Ip), 
                        "IA" = rowSums(IA), "Ii" = rowSums(Ii), 
                        "It" = rowSums(It), "Iti" = rowSums(Iti), 
                        "Iq" = rowSums(Iq), "R" = rowSums(R),
                        "Cc" = rowSums (Cc)))
    
    
  })
  
  output$Plot <- renderPlotly({  
    req(Sim())
    req(count_data())
    p <- plot_ly(data = Sim(),  x = ~time, y = ~(Cc), name = 'Estimates',
                 type = "scatter", mode = 'lines',  height = 480,
                 line = list(color = 'tomato')) %>%
      layout(yaxis = list(title = "Daily cumulative no. of cases"),
             xaxis = list(title ="Time(days)")) %>% 
      add_trace(data = count_data()[c(1:225),], x = ~seq_along(date), y = ~cases, 
                name = 'Reported',
                mode = 'lines', 
                line = list(color = 'black') )
    
    p
  })
  
  output$Estimates <- DT::renderDataTable(
    DT::datatable({
      req(Ests())
      ests <- exp(c(Ests()[9], Ests()[7],Ests()[6],Ests()[5],Ests()[4],Ests()[3]))
      
      names(ests) <- c(paste("Lockdown lvl",seq(0,5)))
      
      t(as.data.frame(ests))}, 
      colnames = c(paste("Lockdown lvl",seq(0,5))),
      caption = 'The effect of lockdown levels on overall contacts.',
      rownames = FALSE,
      extensions = 'Buttons',
      
      options = list(
        paging = FALSE,#TRUE,
        searching = FALSE,#TRUE,
        fixedColumns = TRUE,
        autoWidth = TRUE,
        ordering = FALSE,#TRUE,
        dom = 'tB',
        buttons = c('copy', 'csv')
      ),
      
      class = "display"))
  
}

# Create Shiny app ----
shinyApp(ui, server)