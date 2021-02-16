server <- function(input, output, session) {
  
  ####-- 2. Model Dashboard ------------------------------------------------####
  
  # Selected Groups ----
  selGroups <- reactive({
    dub_population[[1]] %in% input$age_sel
  })
  
  # Process Dates ----
  tmax <- reactive({
    as.numeric( difftime(input$end_date, as.Date("2020-02-29"), units = "days") )
  })
  
  # Intervention correction
  intinfo2 <- reactive({
    dat <- intervention_adjust(interventions_info, input$end_date)
    merge(dat, optim_res, by = "policy")[-1]
  })
  
  # Get Beta Value ----
  dub_gotbeta <- reactive({
    pars <- c(input$L, input$Cv, input$Dv, input$h, 
              input$f, input$tv, input$q, input$k,
              input$TT)
    names(pars) <- c("L","Cv","Dv","h","f","tv","q", "k", "TT")
    getbeta(input$R0, pars = pars, p_age = dub_population$propage,
            CONTACTMATRIX = contacts)
  })
  
  irl_gotbeta <- reactive({
    pars <- c(input$L, input$Cv, input$Dv, input$h, 
              input$f, input$tv, input$q, input$k, 
              input$TT)
    names(pars) <- c("L","Cv","Dv","h","f","tv","q", "k", "TT")
    getbeta(input$R0, pars = pars, p_age = irl_population$propage,
            CONTACTMATRIX = contacts)
  })
  
  # Simulate SEIR Model ----
  dub_sol <- reactive({
    SEIR_model_simulation(pars = c(input$L, input$Cv, input$Dv, input$h, 
                                   input$f, input$tv, input$q, input$k, 
                                   input$TT),
                          contacts_ireland = contacts,
                          dateStart = as.Date("2020-02-29"),
                          startval = dub_xstart,
                          POP = dub_population,
                          beta = dub_gotbeta(),
                          tmax = tmax(), 
                          lockdown_information = intinfo2())$solution 
  })
  
  irl_sol <- reactive({
    SEIR_model_simulation(pars = c(input$L, input$Cv, input$Dv, input$h, 
                                   input$f, input$tv, input$q, input$k,
                                   input$TT),
                          contacts_ireland = contacts,
                          dateStart = as.Date("2020-02-29"),
                          startval = irl_xstart,
                          POP = irl_population,
                          beta = irl_gotbeta(),
                          tmax = tmax(), 
                          lockdown_information = intinfo2())$solution 
  })
  
  # Define Shaded Regions for Plots ----
  shaded_regions <- reactive({
    
    if (!input$shade) { return(NULL) }
    
    raw_scales <- intinfo2()[[3]]
    scales <- 1 - (raw_scales - min(raw_scales)) / (max(raw_scales) - min(raw_scales))
    N <- nrow(intinfo2()) 
    shaded_list <- vector("list", length = N)
    
    for (i in 1:N) {
      
      shaded_list[[i]] <- list(type = "rect",
                               fillcolor = "gray", line = list(color = "gray"), 
                               opacity = 0.5 * scales[i], 
                               x0 = intinfo2()[[1]][i], 
                               x1 = intinfo2()[[2]][i], xref = "x",
                               ysizemode = "pixel",
                               yanchor = 0,
                               y0 = 0, y1 = 350, yref = "y",
                               layer = "below")
      
    }
    
    shaded_list
    
  })
  
  # Create Plots ----
  xval <- reactive({
    as.Date("2020-02-29") + dub_sol()$time
  })
  
  output$summary_plot_dub <- renderPlotly({
    
    summary_plot(dub_sol(), xval(), comp_vec, selGroups(), input$I_type) %>%
      layout(shapes = shaded_regions())
    
  })
  
  output$summary_plot_irl <- renderPlotly({
    
    summary_plot(irl_sol(), xval(), comp_vec, selGroups(), input$I_type) %>%
      layout(shapes = shaded_regions())
    
  })
  
  output$S_age_plot_dub <- renderPlotly({
    
    comp_plot(dub_sol(), xval(), "S_", selGroups(), input) %>% 
      layout(shapes = shaded_regions()) %>%
      layout(xaxis = list(title = ""), 
             title = list(text = "Age Breakdown"),
             yaxis = list(title = list(text = "Susceptible",
                                       font = list(color = 'rgb(69, 95, 245)')))
             )
  })
  
  output$S_age_plot_irl <- renderPlotly({
    
    comp_plot(irl_sol(), xval(), "S_", selGroups(), input) %>% 
      layout(shapes = shaded_regions()) %>%
      layout(xaxis = list(title = ""), 
             title = list(text = "Age Breakdown"),
             yaxis = list(title = list(text = "Susceptible",
                                       font = list(color = 'rgb(69, 95, 245)')))
      )
    
  })
  
  output$E_age_plot_dub <- renderPlotly({
    
    comp_plot(dub_sol(), xval(), "Ev_", selGroups(), input) %>% 
      layout(shapes = shaded_regions()) %>%
      layout(xaxis = list(title = ""), 
             title = list(text = "Age Breakdown"),
             yaxis = list(title = list(text = "Exposed",
                                       font = list(color = 'rgb(214, 122, 17)')))
      )
    
  })
  
  output$E_age_plot_irl <- renderPlotly({
    
    comp_plot(irl_sol(), xval(), "Ev_", selGroups(), input) %>% 
      layout(shapes = shaded_regions()) %>%
      layout(xaxis = list(title = ""), 
             title = list(text = "Age Breakdown"),
             yaxis = list(title = list(text = "Exposed",
                                       font = list(color = 'rgb(214, 122, 17)')))
      )
    
  })
  
  output$I_age_plot_dub <- renderPlotly({
    
    Sy <- c("Ii_", "It_", "Iti_", "Iq_")
    
    g <- function(x) {
      res <- numeric(0)
      for (i in length(x)) {
        res <- rbind(res, dub_sol()[grepl(x[i],names(dub_sol()))])
      }
      res
    }
    
    I_Sy <- g(Sy)
    I_Pr <- g('Ip_')
    I_As <- g('IA_')
    I_Al <- I_Sy + I_Pr + I_As
    I <- get( paste0("I_", substr(input$I_type, start = 1, stop = 2)) )
    
    ### Draw the Plot
    p_I <- plot_ly(x = ~xval(), y = NA, type = 'scatter', mode = 'lines')
    
    for ( i in seq_along(input$age_sel) ) {
      p_I <- p_I %>% add_trace(y = I[[i]], name = input$age_sel[i], mode = 'lines')
    }
    
    p_I %>% layout(shapes = shaded_regions()) %>%
      layout(xaxis = list(title = "Time"), 
             yaxis = list(title = list(text = paste0("Infected: ", input$I_type),
                                       font = list(color = 'rgb(186, 24, 19)')))
      )
    
  })
  
  output$I_age_plot_irl <- renderPlotly({
    
    Sy <- c("Ii_", "It_", "Iti_", "Iq_")
    
    g <- function(x) {
      res <- numeric(0)
      for (i in length(x)) {
        res <- rbind(res, irl_sol()[grepl(x[i],names(irl_sol()))])
      }
      res
    }
    
    I_Sy <- g(Sy)
    I_Pr <- g('Ip_')
    I_As <- g('IA_')
    I_Al <- I_Sy + I_Pr + I_As
    I <- get( paste0("I_", substr(input$I_type, start = 1, stop = 2)) )
    
    ### Draw the Plot
    p_I <- plot_ly(x = ~xval(), y = NA, type = 'scatter', mode = 'lines')
    
    for ( i in seq_along(input$age_sel) ) {
      p_I <- p_I %>% add_trace(y = I[[i]], name = input$age_sel[i], mode = 'lines')
    }
    
    p_I %>% layout(shapes = shaded_regions()) %>%
      layout(xaxis = list(title = "Time"), 
             yaxis = list(title = list(text = paste0("Infected: ", input$I_type),
                                       font = list(color = 'rgb(186, 24, 19)')))
      )
    
  })
  
  output$R_age_plot_dub <- renderPlotly({
    
    comp_plot(dub_sol(), xval(), "R_", selGroups(), input) %>% 
      layout(shapes = shaded_regions()) %>%
      layout(xaxis = list(title = "Time"), 
             yaxis = list(title = list(text = "Removed",
                                       font = list(color = 'rgb(23, 191, 26)')))
      )
    
  })
  
  output$R_age_plot_irl <- renderPlotly({
    
    comp_plot(irl_sol(), xval(), "R_", selGroups(), input) %>%
      layout(shapes = shaded_regions()) %>%
      layout(xaxis = list(title = "Time"), 
             yaxis = list(title = list(text = "Removed",
                                       font = list(color = 'rgb(23, 191, 26)')))
      )
    
  })
  
  ####-- 3. Forecast Settings ----------------------------------------------####
  
  # tmax
  tmax_forecast <- reactive({
    as.numeric( difftime(input$start_date + 55, as.Date('2020-02-29'), units = "days") )
  })
  
  # Intervention correction
  intinfo3 <- reactive({
    intervention_adjust(interventions_info, input$start_date - 1)
  })
  
  # Lockdown Info
  linfo_forecast <- reactive({
    date_start_seq <- seq.Date(input$start_date, input$start_date + 55, 14)
    chosen_levs <- c(input$res1, input$res2, input$res3, input$res4)
    lockdown_forecast <- rbind(intinfo3(),
                               data.frame(start = date_start_seq, 
                                          end = date_start_seq + 13,
                                          policy = chosen_levs))
    optim_dat <- merge(lockdown_forecast, optim_res, by = "policy")
    merge(optim_dat, boot_scales_t, by = "policy")
  })
  
  # Simulate SEIR Model
  dub_forecasters <- eventReactive(input$fit_forecast, {
    
    tmax <- tmax_forecast()
    linfo <- linfo_forecast()
    
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    clusterExport(
               cl, varlist=c("SEIR_model_simulation", "def_pars", 
                             "contacts", "dub_xstart", "dub_population",
                             "dub_def_beta", "lsoda", "SEIR_model_D", 
                             "tmax", "linfo"),
               envir=environment())
    
    All_Runs <- foreach(i = 5:ncol(linfo), .combine = rbind) %dopar% {
      SEIR_model_simulation(pars = def_pars,
                            contacts_ireland = contacts,
                            dateStart = as.Date('2020-02-29'),
                            startval = dub_xstart,
                            POP = dub_population,
                            beta = dub_def_beta,
                            tmax = tmax, 
                            lockdown_information = linfo[, c(2:3, i)])$solution
    }
    
    All_Runs <- split(All_Runs, All_Runs$time)
    UL <- foreach(x = All_Runs, .combine = rbind,
                  .final = as.data.frame) %dopar% {
      sapply(x, quantile, probs = 0.975, names = FALSE)
    }
    LL <- foreach(x = All_Runs, .combine = rbind,
                  .final = as.data.frame) %dopar% {
      sapply(x, quantile, probs = 0.025, names = FALSE)
    }
    rm(All_Runs)
    stopCluster(cl)
    
    MID <- SEIR_model_simulation(pars = def_pars,
                                 contacts_ireland = contacts,
                                 dateStart = as.Date('2020-02-29'),
                                 startval = dub_xstart,
                                 POP = dub_population,
                                 beta = dub_def_beta,
                                 tmax = tmax_forecast(), 
                                 lockdown_information = linfo_forecast()[, 2:4])$solution
    
    list(MID = MID, UL = as.data.frame(UL), LL = as.data.frame(LL))
    
  })

  dub_out <- reactive({
    Map(comp_sel, x = dub_forecasters(),
        MoreArgs = list(y = input$Disp_comp))
  })
  
  date_place <- eventReactive(input$fit_forecast, {
    input$start_date + seq(1-N_known, 55)
  })
  
  segs <- eventReactive(input$fit_forecast, {
    seq(input$start_date, input$start_date + 42, 14)
  })
  
  output$forecast_plot_dub <- renderPlotly({
    
    plot_ly(x = ~date_place(), y = ~dub_out()$UL, name = 'Upper Limit', type = 'scatter', 
            mode = 'lines', line = list(color = 'blue', dash = "dash")) %>%
      add_trace(y = ~dub_out()$LL, name = "Lower Limit", mode = 'lines', 
                fill = 'tonexty', fillcolor='rgba(70, 136, 242, 0.2)',
                line = list(color = 'blue')) %>%
      add_trace(y = ~dub_out()$MID, name = "Forecast", mode = 'lines',
                line = list(color = 'red')) %>%
      add_trace(x = ~date_place()[1:N_known], y = ~dub_out()$MID[1:N_known], name = "Forecast", 
                mode = 'lines',
                line = list(color = 'red', dash = "solid")) %>%
      add_trace(x = ~date_place()[1:N_known], y = ~dub_out()$LL[1:N_known], name = "Lower Limit", 
                mode = 'lines', fill = "none",
                line = list(color = 'blue', 
                            dash = "solid")) %>%
      add_trace(x = ~date_place()[1:N_known], y = ~dub_out()$UL[1:N_known], name = "Upper Limit", 
                mode = 'lines', fill = "none",
                line = list(color = 'blue', 
                            dash = "solid")) %>%
      add_segments(x = segs(), xend = segs(), y = min(dub_out()$LL), yend = max(dub_out()$UL),
                   line = list(color = 'black', dash = "solid", width = 0.5)) %>%
      layout(showlegend = FALSE, xaxis = list(title = "Date"), yaxis = list(title = input$Disp_comp),
             paper_bgcolor='rgb(255,255,255)', plot_bgcolor='rgb(229,229,229)')
    
  })
  
}

