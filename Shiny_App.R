
#### Load required packages ------------------------------------------------####

list.of.packages <- c("shiny", "mlbench", "plotly", "shinythemes", "dplyr", 
                      "shinyWidgets", "ggplot2", "shinydashboard", "data.table",
                      "optimx", "shinycssloaders", "DT")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(new.packages) 

lapply(as.list(list.of.packages), require, character.only = TRUE)

#### Load required source files --------------------------------------------####
source("code/getbeta.R")
source("code/function_IEMAGSEIR_Dublin.R")
source('code/1_loadData.r')

# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("The Effect of Interventions"),
  
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
      #fluidRow(
      # Output: Data file ----
      #tableOutput("contents")
      #fluidRow(
      # column(6, 
      #       plotlyOutput("Plot")),
      #column(6, 
      #      plotlyOutput("PlotRt")),
      splitLayout(cellWidths = c("50%", "50%"),
                  
                  plotlyOutput("Plot")%>% withSpinner(color = "#0dc5c1", size = 2),
                  plotlyOutput("PlotRt")%>% withSpinner(color = "#0dc5c1", size = 2)
      ),
      br(),#br(),#br(),#br(),#br(),br(),
      DT::dataTableOutput("Estimates")#tableOutput("Estimates")#,
      #)
    )
    
  )
)

# Define server logic to read selected file ----
server <- function(input, output) {
  
  count_data <- reactive({
    req(input$CountData)
    County <- "Dublin"
    read.csv(input$CountData$datapath, sep = ",") %>% 
      filter(CountyName == County) %>%
      select(TimeStamp,ConfirmedCovidCases) %>% 
      rename(date = TimeStamp , cases = ConfirmedCovidCases) %>%
      mutate(date = as.Date(date)) 
    
  })
  
  pop_data <- reactive({
    req(input$PopData)
    
    read.csv(input$PopData$datapath, sep = ",")
    
  })
  Interventions_data <- reactive({
    req(input$InterInfo)
    
    read.csv(input$InterInfo$datapath, sep = ",") %>% 
      mutate(start = as.Date(start, "%d/%m/%y"),end = as.Date(end, "%d/%m/%y")) 
    
  })
  Ests <- eventReactive(input$go, {
    
    req(count_data())
    req(pop_data())
    req(Interventions_data())
    
    nlm_fun <- function(scalars, Actual_Cc, Pop_Struc, Inter_info) {
      
      scalars_ests <- exp(scalars)  
      
      names(scalars_ests) <- unique(Inter_info$policy)
      Base <- simulation_SEIR_model(R0t = 3.4,
                                    POP = Pop_Struc,
                                    contacts_ireland = contacts,
                                    interventions = Inter_info, 
                                    dt = 1,
                                    # = interventions$start[1],
                                    dateEnd = count_data()$date[dim(count_data())[1]],# Time step (days)
                                    scalars = scalars_ests)
      
      Est_Cc <- rowSums(Base$sol_out[grepl('Cc_',names(Base$sol_out))])
      
      sum((Actual_Cc - (Est_Cc[-1]))^2) 
      
    }
    
    scalars_init <- c(1.839589e+00, 0.510116460, 0.09050230, 0.008395608, 0.122922232,
                      0.207667229, 0.308520907, 0.101967661, 0.10010276)#c(6.858719e-02, 1.839589e+00, 9.594844e-02, 4.283682e-04, 2.721637e-06,
    # 3.610250e-01, 2.099779e-01, 1.736524e-01, 9.436696e-02)
    
    #scalars <- scalars_init
    ests <- nlm(nlm_fun, log(scalars_init), Actual_Cc = count_data()$cases,
                Pop_Struc = pop_data(), Inter_info = Interventions_data(),
                stepmax = 0.5, iterlim = 1)$est
    
    names(ests) <- unique(Interventions_data()$policy)
    
    exp(ests)
    
  })
  
  Sim <- reactive({
    req(Ests())
    req(pop_data())
    req(Interventions_data())
    groups <- dim(contacts[[1]])[1]
    sim <- simulation_SEIR_model(R0t = 3.4,
                                 POP = pop_data(),
                                 contacts_ireland = contacts,
                                 interventions = Interventions_data(), 
                                 dt = 1,
                                 scalars = Ests())
    # 
    # S <- sim$sol_out[grepl('S_',names(sim$sol_out))]
    # Ev <- sim$sol_out[grepl('Ev_',names(sim$sol_out))]
    # Ip <- sim$sol_out[grepl('Ip_',names(sim$sol_out))]
    # IA <- sim$sol_out[grepl('IA_',names(sim$sol_out))]
    # Ii <- sim$sol_out[grepl('Ii_',names(sim$sol_out))]
    # It <- sim$sol_out[grepl('It_',names(sim$sol_out))]
    # Iti <- sim$sol_out[grepl('Iti_',names(sim$sol_out))]
    # Iq <- sim$sol_out[grepl('Iq_',names(sim$sol_out))]
    # R <- sim$sol_out[grepl('R_',names(sim$sol_out))]
    Cc <- sim$sol_out[grepl('Cc_',names(sim$sol_out))]
    
    tstart_intervention <- sim$tstart_intervention  
    tend_intervention <- sim$tend_intervention
    
    scalars <- sim$scalars
    
    constraint <- vector()
    Rt <- vector()
    
    for ( j in 1:length(sim$sol_out$time)) {
      
      t <- sim$sol_out$time[j]
      
      if (t == 0) {
        INTERVENTION <- noquote(names(tstart_intervention)[1])
        constraint[j]  <- as.numeric(scalars[INTERVENTION])
      } else { for( i in 1:length(tstart_intervention)){
        if( (t > ((tstart_intervention[[i]]) - 1)) & (t <= (tend_intervention[[i]])))
        {
          INTERVENTION <- noquote(names(tstart_intervention)[i])
          constraint[j]  <- as.numeric(scalars[INTERVENTION])
        } 
      }
      } 
      
      
      Const <- list(home = diag(constraint[j], groups),
                    work = diag(constraint[j], groups), 
                    school = diag(constraint[j], groups),
                    others = diag(constraint[j], groups))
      
      ## Estimating R0t
      
      pars <- c(4.9, 5.9, 7.0, 0.25, 0.05, 0.05, 0.5, 0.75, 0.13, 3.6)
      names(pars) <- c("L","Cv","Dv","h","i","j","f","tv","q","TT")
      
      Beta0 <- getbeta(R0t = 3.4, pars = pars, p_age =  pop_data()$propage, CONTACTMATRIX = contacts)
      
      Rt[j] <- getR0t(beta =  Beta0, pars = pars, p_age = pop_data()$propage, constraints = Const,
                      CONTACTMATRIX = contacts)
    }
    
    as.data.frame(cbind("time" = sim$sol_out$time, 
                        "Cc" = rowSums(Cc), "constraint" = constraint, "Rt" = Rt))
    
  })
  
  output$Plot <- renderPlotly({  
    req(Sim())
    req(count_data())
    p <- plot_ly(data = Sim(),  x = ~time, y = ~(Cc), name = 'Estimates',
                 type = "scatter", mode = 'lines', # height = 480,
                 line = list(color = 'tomato')) %>%
      layout(yaxis = list(title = "Daily cumulative no. of cases"),
             xaxis = list(title ="Time(days)")) %>% 
      add_trace(data = count_data(), x = ~seq_along(date), y = ~cases, 
                name = 'Reported',
                mode = 'lines', 
                line = list(color = 'black') )
    
    p
  })
  
  output$PlotRt <- renderPlotly({  
    req(Sim())
    R0t <- 3.4
    p <- plot_ly(data = Sim(),  x = ~time, y = ~Rt,#(constraint*R0t), 
                 #name = 'The Effective reproduction Number',
                 type = "scatter", mode = 'lines',  #height = 480,
                 line = list(color = 'tomato')) %>%
      layout(yaxis = list(title = "The Effective Reproduction Number"),
             xaxis = list(title ="Time(days)")) 
    
    p
  })
  
  output$Estimates <- DT::renderDataTable(
    DT::datatable({
      req(Ests())
      req(Interventions_data())
      
      #ests <- c(mean(Ests()[1],Ests()[1]), ((Ests()[6]+ Ests()[7])/2),Ests()[7],Ests()[8],((Ests()[8] + Ests()[9])/2),Ests()[9])
      
      #names(ests) <- c(paste("Lockdown lvl",seq(0,5)))
      as.data.frame(cbind(Interventions_data()[[1]],Interventions_data()[[2]], 
                          Ests(), Interventions_data()[[3]]))}, 
      colnames = c("Start Date", " End Date", "Intervention Effect", "Intervention"),
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