server <- function(input, output, session) {
  
  # Process Population Data ----
  popSize <- reactive({
    if (is.null(input$pop_file)) { return(population) }
    popdat <- read.csv(input$pop_file$datapath)
    popsize <- popdat[[2]]
    data.frame(age_label = as.character(popdat[[1]]), popage = popsize, propage = popsize/sum(popsize))
  })
  
  # Update Checkbox Input when Population Size Changes ----
  observe({
    updateCheckboxGroupInput(session, inputId = "age_sel", label =  "Select Ages to Display:",
                             choices = popSize()[[1]], selected = popSize()[[1]])
  })
  
  selGroups <- reactive({
    popSize()[[1]] %in% input$age_sel
  })
  
  xstart <- reactive({
    
    num_group <- nrow(popSize())
    num_exp <- input$Estart/num_group
    num_inf <- input$Istart/num_group
    num_rec <- input$Rstart/num_group
    
    xstart <- c(popSize()$popage - num_exp - num_inf - num_rec,
                rep(num_exp, num_group),
                rep(num_inf, num_group),
                rep(0, num_group),
                rep(0, num_group),
                rep(0, num_group),
                rep(0, num_group),
                rep(0, num_group),
                rep(num_rec, num_group))
    
    names(xstart) <- c(paste0('S_',1:num_group),
                       paste0('Ev_',1:num_group),
                       paste0('Ip_',1:num_group),
                       paste0('IA_',1:num_group),
                       paste0('Ii_',1:num_group),
                       paste0('It_',1:num_group),
                       paste0('Iti_',1:num_group),
                       paste0('Iq_',1:num_group),
                       paste0('R_',1:num_group))
    
    xstart
    
  })
  
  # Process Contact Matrices Input ----
  input_h_CM <- reactive({
    if (is.null(input$House_CM)) { return(contacts$home) }
    as.matrix( read.csv(input$House_CM$datapath, header = FALSE) )
  })
  
  input_w_CM <- reactive({
    if (is.null(input$Work_CM)) { return(contacts$work) }
    as.matrix( read.csv(input$Work_CM$datapath, header = FALSE) )
  })
  
  input_s_CM <- reactive({
    if (is.null(input$School_CM)) { return(contacts$school) }
    as.matrix( read.csv(input$School_CM$datapath, header = FALSE) )
  })
  
  input_o_CM <- reactive({
    if (is.null(input$Other_CM)) { return(contacts$others) }
    as.matrix( read.csv(input$Other_CM$datapath, header = FALSE) )
  })
  
  input_CM <- reactive({
    all <- input_h_CM() + input_w_CM() + input_s_CM() + input_o_CM()
    list(home = input_h_CM(), work = input_w_CM(), school = input_s_CM(), 
         others = input_o_CM(), all = all)
  })
  
  # Process Lockdown Info ----
  linfo <- reactive({
    if (is.null(input$Lockdown_Info)) { return(input$Lockdown_Info) }
    read_csv(input$Lockdown_Info$datapath, col_types = "DDn")
  })
  
  # Process Dates ----
  tmax <- reactive({
    as.numeric( difftime(input$dates[2], input$dates[1], units = "days") )
  })
  
  # Get Beta Value ----
  gotbeta <- reactive({
    pars <- c(input$L, input$Cv, input$Dv, input$h, 
              input$f, input$tv, input$q, input$TT)
    names(pars) <- c("L","Cv","Dv","h","f","tv","q","TT")
    getbeta_new(input$R0, pars = pars, p_age = popSize()$propage, 
                CONTACTMATRIX = input_CM())
  })
  
  # Simulate SEIR Model ----
  sol <- reactive({
    SEIR_model_simulation(pars = c(input$L, input$Cv, input$Dv, input$h, 
                                   input$f, input$tv, input$q, input$TT),
                          contacts_ireland = input_CM(),
                          dateStart = input$dates[1],
                          startval = xstart(),
                          POP = popSize(),
                          beta = gotbeta(),
                          tmax = tmax(), 
                          lockdown_information = linfo())$solution 
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
    input$dates[1] + sol()$time
  })
  
  output$summary_plot <- renderPlotly({
    
    # Extract Compartments
    # browser()
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
    
    # Sum across desired age groups
    S_draw <- rowSums(S[selGroups()])
    E_draw <- rowSums(E[selGroups()])
    I_draw <- rowSums(I[selGroups()])
    R_draw <- rowSums(R[selGroups()])
    
    # Draw the Plot
    plot_ly(x = ~xval(), y = ~S_draw, name = 'Susceptible', type = 'scatter', mode = 'lines',
            line = list(color = 'rgb(69, 95, 245)')) %>%
      add_trace(y = ~E_draw, name = 'Exposed', mode = 'lines', line = list(color = 'rgb(214, 122, 17)')) %>% 
      add_trace(y = ~I_draw, name = 'Infected', mode = 'lines', line = list(color = 'rgb(186, 24, 19)')) %>%
      add_trace(y = ~R_draw, name = 'Removed', mode = 'lines', line = list(color = 'rgb(23, 191, 26)')) %>% 
      layout(shapes = shaded_regions()) %>%
      layout(xaxis = list(title = "Time"), 
             yaxis = list(title = "Compartment Size", 
                          range = c(0, rowSums(sol()[1, ]) * 1.01)))
    
  })
  
  output$S_age_plot <- renderPlotly({
    
    S <- sol()[grepl('S_',names(sol()))][selGroups()]
    
    ### Draw the Plot
    p_S <- plot_ly(x = ~xval(), y = NA, type = 'scatter', mode = 'lines')
    
    for ( i in seq_along(input$age_sel) ) {
      p_S <- p_S %>% add_trace(y = S[[i]], name = input$age_sel[i], mode = 'lines')
    }
    
    p_S %>% layout(shapes = shaded_regions()) %>%
      layout(xaxis = list(title = "Time"), 
             yaxis = list(title = list(text = "Susceptible",
                                       font = list(color = 'rgb(69, 95, 245)')),
                          range = c(0, max(S) * 1.01))
      )
    
  })
  
  output$E_age_plot <- renderPlotly({
    
    E <- sol()[grepl('Ev_',names(sol()))][selGroups()]
    
    ### Draw the Plot
    p_E <- plot_ly(x = ~xval(), y = NA, type = 'scatter', mode = 'lines')
    
    for ( i in seq_along(input$age_sel) ) {
      p_E <- p_E %>% add_trace(y = E[[i]], name = input$age_sel[i], mode = 'lines')
    }
    
    p_E %>% layout(shapes = shaded_regions()) %>%
      layout(xaxis = list(title = "Time"), 
             yaxis = list(title = list(text = "Exposed",
                                       font = list(color = 'rgb(214, 122, 17)')),
                          range = c(0, max(E) * 1.01))
      )
    
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
    p_I <- plot_ly(x = ~xval(), y = NA, type = 'scatter', mode = 'lines')
    
    for ( i in seq_along(input$age_sel) ) {
      p_I <- p_I %>% add_trace(y = I[[i]], name = input$age_sel[i], mode = 'lines')
    }
    
    p_I %>% layout(shapes = shaded_regions()) %>%
      layout(xaxis = list(title = "Time"), 
             yaxis = list(title = list(text = paste0("Infected: ", input$I_type),
                                       font = list(color = 'rgb(186, 24, 19)')),
                          range = c(0, max(I) * 1.01))
      )
    
  })
  
  output$R_age_plot <- renderPlotly({
    
    R <- sol()[grepl('R_',names(sol()))]
    
    ### Draw the Plot
    p_R <- plot_ly(x = ~xval(), y = NA, type = 'scatter', mode = 'lines')
    
    for ( i in seq_along(input$age_sel) ) {
      p_R <- p_R %>% add_trace(y = R[[i]], name = input$age_sel[i], mode = 'lines')
    }
    
    p_R %>% layout(shapes = shaded_regions()) %>%
      layout(xaxis = list(title = "Time"), 
             yaxis = list(title = list(text = "Removed",
                                       font = list(color = 'rgb(23, 191, 26)')),
                          range = c(0, max(R) * 1.01))
      )
    
  })
  
  ################################################################################ 
  ## Second tab
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
    
    nlm_fun <- function(scalars, Actual_Cc, population_df, contacts_list, interventions_df) {
      
      scalars_est <- exp(scalars)  
      
      names(scalars_est) <- unique(interventions_df$policy)
      
      Simulation <- try(simulation_SEIR_model(R0t = 3.4,
                                       POP = population_df,
                                       contacts_ireland = contacts_list,
                                       interventions = interventions_df, 
                                       dt = 1,
                                       dateStart = interventions_df$start[1],
                                       dateEnd = cumulative_cases$date[dim(cumulative_cases)[1]],# Time step (days)
                                       scalars = scalars_est), silent = TRUE)
      
      if (class(Simulation) != "try-error") {
        
        Est_Cc <- rowSums(Simulation$sol_out[grepl('Cc_',names(Simulation$sol_out))])
        
        ss <- try(sum((Actual_Cc - Est_Cc)^2), silent = TRUE)
        
        if ( !is.na(ss) & class(ss) != "try-error"){
          tss <- ss
        } else {
          tss <- Inf}
      } else {
        tss <- Inf
      }
      tss
    }
    
    scalars_init <- c(2.27619341, 0.94894146, 0.21845530, 0.01503771, 0.31684429, 0.40807037, 0.49130460,
                      0.34598300, 0.20890798)
    
    ests <- nlm(nlm_fun, log(scalars_init), Actual_Cc = count_data()$cases,
                population_df = pop_data(), interventions_df = Interventions_data(),
                contacts_list = contacts,
                # stepmax = 0.5, 
                iterlim = 1)$est
    
    names(ests) <- unique(Interventions_data()$policy)
    
    exp(ests)
    
  })
  
  Sim <- reactive({
    
    groups <- dim(contacts[[1]])[1]
    sim <- simulation_SEIR_model(R0t = 3.4,
                                 POP = pop_data(),
                                 contacts_ireland = contacts,
                                 interventions = Interventions_data(), 
                                 dt = 1,
                                 dateStart = Interventions_data()$start[1], 
                                 dateEnd = Interventions_data()$end[dim(Interventions_data())[1]],# 
                                 scalars = Ests())
    
    Cc <- sim$sol_out[grepl('Cc_',names(sim$sol_out))]
    
    tstart_intervention <- sim$tstart_intervention  
    tend_intervention <- sim$tend_intervention
    
    scalars <- sim$scalars
    
    constraint <- vector()
    Rt <- vector()
    
    for ( j in 1:length(sim$sol_out$time)) {
      
      t <- sim$sol_out$time[j]
      
      intervention_ind <- (t >= tstart_intervention) & (t < (tend_intervention + 1))  ## rounding t or not rounding t?
      
      intervention_label <- noquote(names(intervention_ind[intervention_ind == 1]))
      
      constraint[j] <- ifelse( !any(intervention_ind), scalars["No Intervention"], scalars[intervention_label])
      
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
  
  output$PlotCc <- renderPlotly({  
    
    p <- plot_ly(data = Sim(),  x = ~time, y = ~(Cc), name = 'Estimates',
                 type = "scatter", mode = 'lines', 
                 line = list(color = 'tomato')) %>%
      layout(yaxis = list(title = "Daily cumulative no. of cases"),
             xaxis = list(title ="Time(days)")) %>% 
      add_trace(data = count_data(), x = ~seq_along(date), y = ~cases, 
                name = 'Reported',
                mode = 'lines', 
                line = list(color = 'black') )
    
    p <- p %>% layout(legend = list(x = 0.1, y = 0.9))
    
    p
    
  })
  
  output$PlotRt <- renderPlotly({  
    R0t <- 3.4
    p <- plot_ly(data = Sim(),  x = ~time, y = ~Rt,
                 type = "scatter", mode = 'lines',  
                 line = list(color = 'tomato')) %>%
      layout(yaxis = list(title = "The Effective Reproduction Number"),
             xaxis = list(title ="Time(days)")) 
    
    p
  })
  
  output$Estimates <- DT::renderDataTable(
    DT::datatable({
      req(Interventions_data())
      
      as.data.frame(cbind(Interventions_data(), Ests()))},
      colnames = c("Start Date", " End Date", "Intervention", "Intervention Effect"),
      caption = 'The effect of different interventions on overall contacts.',
      rownames = FALSE,
      extensions = 'Buttons',
      
      options = list(
        paging = FALSE,
        searching = FALSE,
        fixedColumns = TRUE,
        autoWidth = TRUE,
        ordering = FALSE,
        dom = 'tB',
        buttons = c('copy', 'csv')
      ),
      
      class = "display"))
  
}