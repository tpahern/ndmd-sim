# Server code for simulation of non-differential misclassification of disease under (near)-perfect specificity and imperfect yet (expected) non-differential sensitivity

library(shiny)
library(dplyr)
library(ggplot2)
library(ggExtra)
#library(ggpubr)


# Define server logic required to draw a histogram
#shinyServer(
function(input, output) {
  
  observeEvent(input$run, {

    # Set random or fixed seed based on user selection
    seed <- ifelse(input$randseed==1, Sys.time(), 2718)
    set.seed(seed)
    
    i.row <- 0

      # retrieve parameter values from UI
      niter  <- as.numeric(input$niter)
      size   <- as.numeric(input$studyn)
      ip     <- as.numeric(input$ip_set)
      pe     <- as.numeric(input$pe_set)
      r      <- as.numeric(input$r_set)
      se     <- as.numeric(input$se_set)
      sp     <- as.numeric(input$sp_set)
      
      # Initialize results matrix based on selected 'niter' value
      ifelse(niter==1e4, sim.data <- matrix(NA, nrow=1e4, ncol=27),
             ifelse(niter==1e5, sim.data <- matrix(NA, nrow=1e5, ncol=27),
                           sim.data <- matrix(NA, nrow=1e6, ncol=27)))
      
      # give up to 2x niter to achieve niter valid records
      niter.max <- niter*2
            
      # calculate risk in the exposed
      prob   <- min(ip*r,1)
      
      # compute number of exposed (N1) and unexposed (N0) cases
      N1 <- round(size*pe, digits=0)
      N0 <- size-N1
      
      # record number of successes and counter
      nsuccess <- 0
      ctr <- 0
      
      # want up to 'niter' iterations, but stop after trying 'niter.max' times
      while((nsuccess < niter)&(ctr < niter.max)) {
    
      # uptick the counter
        ctr <- ctr+1
        
        # CORRECTLY CLASSIFIED DATA
        # generate table cells for correctly-classified scenario
        A <- rbinom(1,N1,prob)
        B <- rbinom(1,N0,ip)
        C <- N1-A
        D <- N0-B
        # compute risk ratio for correctly-classified data
        ipr_c <- (A/N1)/(B/N0)
        
        # STOCHASTIC MISCLASSIFICATION
        # false-negative cases
        FN1s <- A-rbinom(1,A,se)
        FN0s <- B-rbinom(1,B,se)
        # false-positive cases
        FP1s <- rbinom(1,N1-A,1-sp)
        FP0s <- rbinom(1,N0-B,1-sp)
        # compute cell frequencies
        a_s <- A - FN1s + FP1s
        b_s <- B - FN0s + FP0s
        c_s <- N1 - a_s
        d_s <- N0 - b_s
        # compute risk ratio for stochastically misclassified data
        iprm_sto <- (a_s / N1) / (b_s / N0)
        # calculate var[ln(iprm_sto)]
        v_sto <- 1 / a_s - 1 / N1 + 1 / b_s - 1 / N0
        # calculate actual values of Se and Sp
        se_exp_s <- (A - FN1s) / A
        se_unexp_s <- (B - FN0s) / B
        sp_exp_s <- (C - FP1s) / C
        sp_unexp_s <- (D - FP0s) / D
        # calculate degree of non-differentiality for Se and Sp
        delta_sp_sto <- sp_exp_s - sp_unexp_s
        delta_se_sto <- se_exp_s - se_unexp_s
        # calculate overall bias factor
        obf_sto <- iprm_sto / ipr_c
        
        
        # DETERMINISTIC MISCLASSIFICATION
        # false-negative cases
        FN1d <- A-round(se*A, digits=0)
        FN0d <- B-round(se*B,digits=0)
        # false-positive cases
        FP1d = round((N1 - A) * (1 - sp),digits=0)
        FP0d = round((N0 - B) * (1 - sp),digits=0)
        # compute cell frequencies
        a_d <- A - FN1d + FP1d
        b_d <- B - FN0d + FP0d
        c_d <- N1 - a_d
        d_d <- N0 - b_d
        # compute risk ratio for deterministically misclassified data
        iprm_det <- (a_d / N1) / (b_d / N0)
        # calculate var[ln(iprm_det)]
        v_det <- 1 / a_d - 1 / N1 + 1 / b_d - 1 / N0
        # calculate actual values of Se and Sp
        se_exp_d <- (A - FN1d) / A
        se_unexp_d <- (B - FN0d) / B
        sp_exp_d <- (C - FP1d) / C
        sp_unexp_d <- (D - FP0d) / D
        # calculate degree of non-differentiality for Se and Sp
        delta_sp_det <- sp_exp_d - sp_unexp_d
        delta_se_det <- se_exp_d - se_unexp_d
        # calculate overall bias factor
        obf_det <- iprm_det / ipr_c
        
        
        # check for NaN or Inf results (invalid iterations)
        flag <- min(is.finite(c(ipr_c, iprm_sto, iprm_det, 
                                se_exp_s, se_unexp_s, sp_exp_s, sp_unexp_s, 
                                se_exp_d, se_unexp_d, sp_exp_d, sp_unexp_d, 
                                v_sto, v_det, 
                                obf_sto, obf_det)))
        
        # flag==0 means there is Nan or Inf, flag==1 means we keep it
        if(flag == 1){
          nsuccess <- nsuccess+1
          i.row <- i.row+1
          
          # add valid results to sim.data
          result <- matrix(c(size, niter, se, sp, pe, ip, r, 
                             flag,
                             se_exp_s, se_unexp_s, sp_exp_s, sp_unexp_s, 
                             se_exp_d, se_unexp_d, sp_exp_d, sp_unexp_d, 
                             delta_sp_sto, delta_se_sto, delta_sp_det, delta_se_det, 
                             ipr_c, iprm_sto, iprm_det, 
                             v_sto, v_det, 
                             obf_sto, obf_det), nrow=1,ncol=27)
          
          sim.data[i.row,] <- result
      }
      
    }
    
    # remove unfilled rows from sim.data
    i.row <- i.row+1
    if(i.row <= niter) {
      sim.data <- sim.data[-(i.row:niter),]
    }
    
    # convert sim.data to a data frame
    # add column names to sim.data
    colnames(sim.data) <- c('size', 'niter', 'se_set', 'sp_set', 'pe_set', 'ip_set', 'r_set',
                            'flag',
                            'se1_s', 'se0_s', 'sp1_s', 'sp0_s',
                            'se1_d', 'se0_d', 'sp1_d', 'sp0_d', 
                            'delta_sp_s', 'delta_se_s', 'delta_sp_d', 'delta_se_d',
                            'iprc', 'iprm_s', 'iprm_d',
                            'v_s', 'v_d',
                            'obf_s', 'obf_d'
                            )
    
    sim.data.df <- data.frame(sim.data)
    
    records <- nrow(sim.data.df)
    badrecs <- niter-records
    output$recs <- renderText({records})
    output$badrecs <- renderText({badrecs})

        
    # generate derived variables
    sim.data.df <- sim.data.df |>
      mutate(
        
        # iprm_s=as.numeric(iprm_s),
        # v_s=as.numeric(v_s),
        
        # confidence limits for misclassified risk ratio
        ul_s = exp(log(iprm_s) + 1.96 * sqrt(v_s)),
        ll_s = exp(log(iprm_s) - 1.96 * sqrt(v_s)),
        ul_d = exp(log(iprm_d) + 1.96 * sqrt(v_d)),
        ll_d = exp(log(iprm_d) - 1.96 * sqrt(v_d)),

        # log bias factors (base 10)
        logobf_s = log10(obf_s),
        logobf_d = log10(obf_d),

        # indicator: observed true risk ratio within interval of misclassified risk ratio
        iprc_in_ci_s = ifelse(ll_s <= iprc & iprc <= ul_s,1,0),
        iprc_in_ci_d = ifelse(ll_d <= iprc & iprc <= ul_d,1,0),

        # indicator: set true risk ratio within interval of misclassified risk ratio
        iprt_in_ci_s = ifelse(ll_s <= r_set & r_set <= ul_s,1,0),
        iprt_in_ci_d = ifelse(ll_d <= r_set & r_set <= ul_d,1,0),

        # indicator: misclassified risk ratio greater than set true risk ratio
        iprm_gt_iprt_s = ifelse(iprm_s > r_set,1,0),
        iprm_gt_iprt_d = ifelse(iprm_d > r_set,1,0),

        # indicator: misclassified risk ratio greater than observed true risk ratio
        iprm_gt_iprc_s = ifelse(iprm_s > iprc,1,0),
        iprm_gt_iprc_d = ifelse(iprm_d > iprc,1,0),

        # indicator: misclassified risk ratio between null and set true risk ratio
        iprm_btw1_iprt_s = ifelse(1 <= iprm_s & iprm_s <= r_set,1,0),
        iprm_btw1_iprt_d = ifelse(1 <= iprm_d & iprm_d <= r_set,1,0),

        # indicator: misclassified risk ratio between null and observed true risk ratio
        iprm_btw1_iprc_s = ifelse(1 <= iprm_s & iprm_s <= iprc,1,0),
        iprm_btw1_iprc_d = ifelse(1 <= iprm_d & iprm_d <= iprc,1,0),

        # indicator: misclassified risk ratio less than null
        iprm_lt_1_s = ifelse(iprm_s < 1,1,0),
        iprm_lt_1_d = ifelse(iprm_d < 1,1,0),

        # define set specificity as difference from 1
        sp_setm1=round(1-sp_set,3),

        # re-scale deltas to unit=thousandths for regression
        delta_se_s_1k=delta_se_s*1000,
        delta_sp_s_1k=delta_sp_s*1000,

        # checksums for mutually exclusive indicator sets
        checksum_s=iprm_lt_1_s+iprm_btw1_iprt_s+iprm_gt_iprt_s,
        checksum_d=iprm_lt_1_d+iprm_btw1_iprt_d+iprm_gt_iprt_d
      )
    
    
    # Summarizing portions of the distribution
    pct_iprm_lt1_s=format(round(100*mean(sim.data.df$iprm_lt_1_s, na.rm=TRUE), digits=3), nsmall=2)
    pct_iprm_btw1_iprt_s=format(round(100*mean(sim.data.df$iprm_btw1_iprt_s, na.rm=TRUE), digits=3), nsmall=2)
    pct_iprm_gt_iprt_s=format(round(100*mean(sim.data.df$iprm_gt_iprt_s, na.rm=TRUE), digits=3), nsmall=2)
    
    ###temp$Se <- format(round(temp$se_set, digits=2), nsmall=2)
    
### Make simulation data downloadable
    output$simdat <- renderTable({
      sim.data.df()
    })
    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("ndmd-sim-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(sim.data.df, file)
      }
    )
    
      
### Text output

    # summary of simulated true risk ratios
    iprc_sum <- summary(sim.data.df$iprc, digits=3)
    output$iprc_sum <- renderPrint(iprc_sum)
    
    # summary of misclassified risk ratios
    iprm_s_sum <- summary(sim.data.df$iprm_s, digits=3)
    output$iprm_s_sum <- renderPrint(iprm_s_sum)
    
    # summary of areas of misclassified risk ratio distribution
    output$iprm_areas <- 
      renderPrint(cat(paste("Summary of misclassified risk ratio distribution:", "\n",
                            "\u2022", pct_iprm_gt_iprt_s,"%","are greater than the true risk ratio", "\n",
                            "\u2022", pct_iprm_btw1_iprt_s,"%", "are between the null and the true risk ratio", "\n",
                            "\u2022", pct_iprm_lt1_s,"%", "are less than the null", "\n")
                      )
                  )
    
    
### Visual output
    
    # scatter plot of Se_exp, Se_unexp with marginal densities
    sejd <- ggplot(sim.data.df, aes(x=se1_s, y=se0_s)) +
      geom_point(color="#1f78b4", alpha=1/20) +
      scale_x_continuous(limits=c(0.5,1.0)) +
      scale_y_continuous(limits=c(0.5,1.0)) +
      xlab("Sensitivity (exposed)") +
      ylab("Sensitivity (unexposed)") +
      theme_minimal(base_size = 16)
      
        output$se_jointdist <- 
        renderPlot({
          ggMarginal(sejd, type="density")
        })

      # scatter plot of Sp_exp, Sp_unexp with marginal densities
      spjd <- ggplot(sim.data.df, aes(x=sp1_s, y=sp0_s)) +
        geom_point(color="#fd8d3c", alpha=1/20) +
        scale_x_continuous(limits=c(0.9,1.0)) +
        scale_y_continuous(limits=c(0.9,1.0)) +
        xlab("Specificity (exposed)") +
        ylab("Specificity (unexposed)") +
        theme_minimal(base_size = 16)
      
          output$sp_jointdist <- 
            renderPlot({
              ggMarginal(spjd, type="density")
            })
        
    
    # density plot of simulated true risk ratios
    dens_iprc <- 
      ggplot(data=sim.data.df, aes(x=iprc)) +
      geom_density() +
      scale_x_continuous(breaks=c(0.5,1.0,2.0,3.0,4.0,5.0), limits=c(0.5,r+2)) +
      geom_vline(aes(xintercept=r),
                 color = "#1f78b4",
                 linetype = "dashed",
                 size = 1) +
      geom_text(aes(x = r + 0.08,
                    label = "Set true risk ratio",
                    y = 0.1),
                color = "#1f78b4",
                angle = 90,
                #vjust = "bottom",
                hjust = "left") +
      xlab("Simulated true risk ratio") +
      ylab("Density") +
      theme_minimal(base_size = 16)
    
    output$dens_iprc <- renderPlot({dens_iprc})
    
    # density plot of misclassified risk ratios
      # data frame with annotation info
      # annotation <- data.frame(
      #   x = c(r+0.53, r+0.57, r+0.57, r+0.57),
      #   y = c(1.7,1.5,1.3,1.1),
      #   label = c(
      #     paste0("Percent of misclassified risk ratios:"),
      #     paste0("\u2022", " greater than the true risk ratio = ", pct_iprm_gt_iprt_s,"%"),
      #     paste0("\u2022", " between the null and the true risk ratio = ", pct_iprm_btw1_iprt_s,"%"),
      #     paste0("\u2022", " less than the null = ", pct_iprm_lt1_s,"%")
      #   )
      #   )
      
      dens_iprms <- 
        ggplot(data=sim.data.df, aes(x=iprm_s)) +
        geom_density() +
        scale_x_continuous(breaks=c(0.5,1.0,2.0,3.0,4.0,5.0), limits=c(0.5,r+2)) +
        geom_vline(aes(xintercept=r),
                   color = "#1f78b4",
                   linetype = "dashed",
                   size = 1) +
        geom_text(aes(x = r + 0.08,
                      label = "Set true risk ratio",
                      y = 0.1),
                  color = "#1f78b4",
                  angle = 90,
                  #vjust = "bottom",
                  hjust = "left") +
        geom_vline(aes(xintercept=1), 
                   color = "grey", 
                   linetype = "dashed",
                   size = 1) +
        geom_text(aes(x = 1.08,
                      label = "Null",
                      y = 0.1),
                  color = "grey",
                  angle = 90,
                  hjust = "left") +
        xlab("Misclassified risk ratio") +
        ylab("Density") +
        theme_minimal(base_size = 16)
        # geom_text(data=annotation, aes(x=x, y=y, label=label),
        #            color="black",
        #            size=5, 
        #            fontface="plain",
        #            hjust = "left") +
        
      
    output$dens_iprms <- renderPlot({dens_iprms})

    # scatter plot of delta_se, delta_sp with points colored by bias factor
    deltascat <- ggplot(sim.data.df, aes(x=delta_se_s, y=delta_sp_s)) +
      geom_point(aes(color=logobf_s), size=3, alpha=0.5) +
      labs(caption = bquote("log(BF)=log"(RR[misclassified]/RR[true]))) +
      scale_x_continuous(name=bquote(Delta[Se]~" = "~Se[exposed] - Se[unexposed])) +
      scale_y_continuous(name=bquote(Delta[Sp]~" = "~Sp[exposed] - Sp[unexposed])) +
      theme(axis.text.x = element_text(size=30)) +
      scale_color_gradient2(low="#150254", 
                            mid="#f7f7f7",
                            high="#540402",
                            midpoint=0,
                            limits=c(min(sim.data.df$logobf_s), max(sim.data.df$logobf_s)),
                            space="Lab", 
                            na.value="black", 
                            guide = "colorbar",
                            aesthetics = "color",
                            name=bquote("log(BF)")) +
      theme_minimal(base_size = 16)
    
    output$deltascat <- renderPlot({deltascat})
    })

}
#)
