server <- function(input, output, session) {
  
  ####-- 3. Model Dashboard ------------------------------------------------####
  xstart_dub <- reactive({
    
    num_group <- nrow(dub_population)
    num_exp <- input$Estart/num_group
    num_inf <- input$Istart/num_group
    num_rec <- input$Rstart/num_group
    
    xstart <- c(dub_population$popage - num_exp - num_inf - num_rec,
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
  
  xstart_irl <- reactive({
    
    num_group <- nrow(irl_population)
    num_exp <- input$Estart/num_group
    num_inf <- input$Istart/num_group
    num_rec <- input$Rstart/num_group
    
    xstart <- c(irl_population$popage - num_exp - num_inf - num_rec,
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
  
  # Selected Groups ----
  selGroups <- reactive({
    dub_population[[1]] %in% input$age_sel
  })
  
  # Process Dates ----
  tmax <- reactive({
    as.numeric( difftime(input$dates[2], input$dates[1], units = "days") )
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
                          dateStart = input$dates[1],
                          startval = xstart_dub(),
                          POP = dub_population,
                          beta = dub_gotbeta(),
                          tmax = tmax(), 
                          lockdown_information = linfo_mean)$solution 
  })
  
  irl_sol <- reactive({
    SEIR_model_simulation(pars = c(input$L, input$Cv, input$Dv, input$h, 
                                   input$f, input$tv, input$q, input$k,
                                   input$TT),
                          contacts_ireland = contacts,
                          dateStart = input$dates[1],
                          startval = xstart_irl(),
                          POP = irl_population,
                          beta = irl_gotbeta(),
                          tmax = tmax(), 
                          lockdown_information = linfo_mean)$solution 
  })
  
  # Define Shaded Regions for Plots ----
  shaded_regions <- reactive({
    
    if (!input$shade) { return(NULL) }
    
    raw_scales <- linfo_mean[[3]]
    scales <- 1 - (raw_scales - min(raw_scales)) / (max(raw_scales) - min(raw_scales))
    N <- nrow(linfo_mean)
    shaded_list <- vector("list", length = N)
    
    for (i in 1:N) {
      
      shaded_list[[i]] <- list(type = "rect",
                               fillcolor = "gray", line = list(color = "gray"), 
                               opacity = 0.5 * scales[i], 
                               x0 = linfo_mean[[1]][i], 
                               x1 = linfo_mean[[2]][i], xref = "x",
                               ysizemode = "pixel",
                               yanchor = 0,
                               y0 = 0, y1 = 350, yref = "y",
                               layer = "below")
      
    }
    
    shaded_list
    
  })
  
  # Create Plots ----
  xval <- reactive({
    input$dates[1] + dub_sol()$time
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
  
  ####-- 4. Forecast Settings ----------------------------------------------####
  
  # Placeholder Upper and Lower limits
  N_known <- 13
  N_full <- 68
  
  date_place <- reactive({
    input$start_date + seq(1-N_known, 55)
  })
  segs <- reactive({
    seq(input$start_date, input$start_date + 42, 14)
  })
  text_places <- reactive({
    c(input$start_date + 7, seq(input$start_date + 21, input$start_date + 49, 14))
  })
    
  UL <- rnorm(N_full, 450, 25)
  MID <- rnorm(N_full, 300, 25)
  LL <- rnorm(N_full, 150, 25)
  
  output$forecast_plot_dub <- renderPlotly({
    
    p <- plot_ly(x = ~date_place(), y = ~UL, name = 'Upper Limit', type = 'scatter', 
            mode = 'lines', line = list(color = 'blue', dash = "dash")) %>%
      add_trace(y = ~LL, name = "Lower Limit", mode = 'lines', 
                fill = 'tonexty', fillcolor='rgba(70, 136, 242, 0.2)',
                line = list(color = 'blue')) %>%
      add_trace(y = ~MID, name = "Forecast", mode = 'lines',
                #fill = 'tonexty', fillcolor='rgba(70, 136, 242, 0.2)', 
                line = list(color = 'red')) %>%
      add_trace(x = ~date_place()[1:N_known], y = ~MID[1:N_known], name = "Forecast2", 
                mode = 'lines',
                line = list(color = 'red', dash = "solid")) %>%
      add_trace(x = ~date_place()[1:N_known], y = ~LL[1:N_known], name = "LL2", 
                mode = 'lines', fill = "none",
                line = list(color = 'blue', 
                            dash = "solid")) %>%
      add_trace(x = ~date_place()[1:N_known], y = ~UL[1:N_known], name = "UL2", 
                mode = 'lines', fill = "none",
                line = list(color = 'blue', 
                            dash = "solid")) %>%
      add_segments(x = segs(), xend = segs(), y = min(LL), yend = max(UL),
                   line = list(color = 'black', dash = "solid", width = 0.5)) %>%
      layout(showlegend = FALSE,
             paper_bgcolor='rgb(255,255,255)', plot_bgcolor='rgb(229,229,229)')
    
    for (i in 1:4) {
      p <- p %>% add_text(x = text_places()[i], y = max(UL)*1.1, 
                          text = input[[paste0("res", i)]],
                          textfont = list(color = '#000000', size = 14))
    }
    
    p
    
  })
  
  output$forecast_plot_irl <- renderPlotly({
    
    p <- plot_ly(x = ~date_place, y = ~UL, name = 'Upper Limit', type = 'scatter', 
                 mode = 'lines', line = list(color = 'blue', dash = "dash")) %>%
      add_trace(y = ~LL, name = "Lower Limit", mode = 'lines', 
                fill = 'tonexty', fillcolor='rgba(70, 136, 242, 0.2)',
                line = list(color = 'blue')) %>%
      add_trace(y = ~MID, name = "Forecast", mode = 'lines',
                #fill = 'tonexty', fillcolor='rgba(70, 136, 242, 0.2)', 
                line = list(color = 'red')) %>%
      add_trace(x = ~date_place[1:N_known], y = ~MID[1:N_known], name = "Forecast2", 
                mode = 'lines',
                line = list(color = 'red', dash = "solid")) %>%
      add_trace(x = ~date_place[1:N_known], y = ~LL[1:N_known], name = "LL2", 
                mode = 'lines', fill = "none",
                line = list(color = 'blue', 
                            dash = "solid")) %>%
      add_trace(x = ~date_place[1:N_known], y = ~UL[1:N_known], name = "UL2", 
                mode = 'lines', fill = "none",
                line = list(color = 'blue', 
                            dash = "solid")) %>%
      add_segments(x = segs, xend = segs, y = min(LL), yend = max(UL),
                   line = list(color = 'black')) %>%
      layout(showlegend = FALSE,
             paper_bgcolor='rgb(255,255,255)', plot_bgcolor='rgb(229,229,229)')
    
    for (i in 1:4) {
      p <- p %>% add_text(x = text_places[i], y = max(UL)*1.1, 
                          text = input[[paste0("res", i)]],
                          textfont = list(color = '#000000', size = 14))
    }
    
    p
    
  })
  
}

