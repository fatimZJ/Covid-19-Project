
#### Load required packages ------------------------------------------------####

list.of.packages <- c("shiny", "mlbench", "plotly", "shinythemes", "dplyr", 
                      "shinyWidgets", "ggplot2", "shinydashboard", "data.table",
                      "optimx")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(new.packages) 

lapply(as.list(list.of.packages), require, character.only = TRUE)

#### Load required source files --------------------------------------------####
source("code/getbeta.R")
source("Code/Test_fun.R")
source("Code/function_IEMAGSEIR_lsoda_scalars.R")
source('code/1_loadData.r')
# Define colours to be used
all_cols <- c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", 
              "#42d4f4", "#f032e6", "#bfef45", "#fabebe", "#469990", "#e6beff", 
              "#9A6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", 
              "#000075", "#a9a9a9", "#000000","deeppink3","turquoise4")

GR_COL = rgb(0,0.5,0.5)
YL_COL = rgb(0.75,0.75,0)

# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("The Effect of Lockdown Measures"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      # Input: Select a file ----
      fileInput(inputId = "CountData", label = "Cumulative Number of Cases per Day", #"Choose CSV File",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      tags$hr(),
      
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
      
      # Output: Data file ----
      #tableOutput("contents")
      plotOutput("Plot"),
      tableOutput("Estimates")#, height = "555px")
      
    )
    
  )
)

# Define server logic to read selected file ----
server <- function(input, output) {
  
  base <- reactive({
    req(input$CountData)
    df <- read.csv(input$CountData$datapath, sep = ",")
    
    Base <- simulation_SEIR_model_ContRt(Rt = df[,1],
                                         dateStartSchoolClosure = as.Date('2020-03-12') , #schools closed before imtense lockdown
                                         dateStartIntenseIntervention = as.Date('2020-03-27') , #Intense intervention: starts at Wuhan Lockdown
                                         dateEndIntenseIntervention = as.Date('2020-05-18'), #date we begin relaxing intense intervention
                                         dateStart = as.Date('2020-02-28'), #start date for epidemic in Ireland
                                         POP = Irlpop,
                                         numWeekStagger = c(3,6,9,12,15),
                                         contacts_ireland = contacts,
                                         dt = 0.1,  
                                         tmax = 225,
                                         scalars = c(1, 1, 1, 1, 1, 1, 1, 1, 1)) 
    
    ## Plotting the solution
    
    Base_S <- Base$sol_out[grepl('S_',names(Base$sol_out))]
    Base_Ev <- Base$sol_out[grepl('Ev_',names(Base$sol_out))]
    Base_Ip <- Base$sol_out[grepl('Ip_',names(Base$sol_out))]
    Base_IA <- Base$sol_out[grepl('IA_',names(Base$sol_out))]
    Base_Ii <- Base$sol_out[grepl('Ii_',names(Base$sol_out))]
    Base_It <- Base$sol_out[grepl('It_',names(Base$sol_out))]
    Base_Iti <- Base$sol_out[grepl('Iti_',names(Base$sol_out))]
    Base_Iq <- Base$sol_out[grepl('Iq_',names(Base$sol_out))]
    Base_R <- Base$sol_out[grepl('R_',names(Base$sol_out))]
    
    as.data.frame(cbind("time" = Base$sol_out$time, "Base_S" = rowSums(Base_S), 
                        "Base_Ev" = rowSums(Base_Ev), "Base_Ip" = rowSums(Base_Ip), 
                        "Base_IA" = rowSums(Base_IA), "Base_Ii" = rowSums(Base_Ii), 
                        "Base_It" = rowSums(Base_It), "Base_Iti" = rowSums(Base_Iti), 
                        "Base_Iq" = rowSums(Base_Iq), "Base_R" = rowSums(Base_R)))
  })
  
  Ests <- eventReactive(input$go, {
    req(base())
    #req(input$CountData)
    #df <- read.csv(input$CountData$datapath, sep = ",")
    
    nlm_fun <- function(scalars, Base) {
      scalars <- exp(scalars)
      
      Base_S <- Base$Base_S
      ##
      test <- simulation_SEIR_model(R0t = 3.4,#JG_Rt[1],
                                    dateStartSchoolClosure = as.Date('2020-03-12') , #schools closed before imtense lockdown
                                    dateStartIntenseIntervention = as.Date('2020-03-27') , #Intense intervention: starts at Wuhan Lockdown
                                    dateEndIntenseIntervention = as.Date('2020-05-18'), #date we begin relaxing intense intervention
                                    dateStart = as.Date('2020-02-28'), #start date for epidemic in Ireland
                                    POP = Irlpop,
                                    numWeekStagger = c(3,6,9,12,15),
                                    contacts_ireland = contacts,
                                    dt = 0.1,
                                    tmax = 225,
                                    scalars = scalars)
      
      test_S <- rowSums(test$sol_out[grepl('S_',names(test$sol_out))])
      
      norm((Base_S - test_S), type = "2")
      
    }
    
    ests <- nlm(nlm_fun,log(c(1.91365769, 1.16790530, 0.20215852, 0.06962944, 
                              0.18328977, 0.43466483, 0.60612880, 0.32053803, 
                              0.45793943)), Base = base(),
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
                                 POP = Irlpop,
                                 numWeekStagger = c(3,6,9,12,15),
                                 contacts_ireland = contacts,
                                 dt = 0.1,  
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
    
    as.data.frame(cbind("time" = sim$sol_out$time, "S" = rowSums(S), 
                        "Ev" = rowSums(Ev), "Ip" = rowSums(Ip), 
                        "IA" = rowSums(IA), "Ii" = rowSums(Ii), 
                        "It" = rowSums(It), "Iti" = rowSums(Iti), 
                        "Iq" = rowSums(Iq), "R" = rowSums(R)))
    
    
  })
  
  output$Plot <- renderPlot({  
    req(base())
    req(Sim())
    p <- ggplot() + geom_line(data = base(), aes(time, 
                                                 (rowSums(cbind(Base_Ip, Base_IA, Base_Ii, 
                                                                Base_It,Base_Iti,Base_Iq))))
                              , color = "black", size = 1.5) +
      geom_line(data = Sim(), aes(time, 
                                  (rowSums(cbind(Ip, IA, Ii, 
                                                 It, Iti,Iq)))), color = "tomato",
                size = 1.5) + theme_minimal() + labs(x = "Time(days)", 
                                                     y = "Daily no. of infections") 
    p
  })
  
  output$Estimates <- renderTable({
    req(Ests())
    ests <- exp(c(Ests()[9], Ests()[7],Ests()[6],Ests()[5],Ests()[4],Ests()[3]))
    
    names(ests) <- c(paste("Lockdown lvl",seq(0,5)))
    
    t(as.data.frame(ests))
    
    xtable(t(as.data.frame(ests)), caption = NULL, label = NULL, align = NULL, digits = NULL,
           display = NULL, auto = FALSE, ...)
  

  })
  
}

# Create Shiny app ----
shinyApp(ui, server)