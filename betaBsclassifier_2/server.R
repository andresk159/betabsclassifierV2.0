suppressMessages(if(!require(shiny)){install.packages("shiny");library(shiny)}else{library(shiny)})
suppressMessages(if(!require(tidyverse)){install.packages("tidyverse");library(tidyverse)}else{library(tidyverse)})
suppressMessages(if(!require(rlang)){install.packages("rlang");library(rlang)}else{library(rlang)})
suppressMessages(if(!require(shinycssloaders)){install.packages("shinycssloaders");library(shinycssloaders)}else{library(shinycssloaders)})
suppressMessages(if(!require(pROC)){install.packages("pROC");library(pROC)}else{library(pROC)})
suppressMessages(if(!require(googleVis)){install.packages("googleVis");library(googleVis)}else{library(googleVis)})
suppressMessages(if(!require(highcharter)){install.packages("highcharter");library(highcharter)}else{library(highcharter)})
suppressMessages(if(!require(RColorBrewer)){install.packages("RColorBrewer");library(RColorBrewer)}else{library(RColorBrewer)})
suppressMessages(if(!require(collapsibleTree)){install.packages("collapsibleTree");library(collapsibleTree)}else{library(collapsibleTree)})
suppressMessages(if(!require(shinydashboardPlus)){install.packages("shinydashboardPlus");library(shinydashboardPlus)}else{library(shinydashboardPlus)})
suppressMessages(if(!require(shinydashboard)){install.packages("shinydashboard");library(shinydashboard)}else{library(shinydashboard)})
suppressMessages(if(!require(shinyWidgets)){install.packages("shinyWidgets");library(shinyWidgets)}else{library(shinyWidgets)})
suppressMessages(if(!require(googleway)){install.packages("googleway");library(googleway)}else{library(googleway)})
suppressMessages(if(!require(bsplus)){install.packages("bsplus");library(bsplus)}else{library(bsplus)})
suppressMessages(if(!require(stringr)){install.packages("stringr");library(stringr)}else{library(stringr)})
suppressMessages(if(!require(stringi)){install.packages("stringi");library(stringi)}else{library(stringi)})
suppressMessages(if(!require(shinyjs)){install.packages("shinyjs");library(shinyjs)}else{library(shinyjs)})
suppressMessages(if(!require(pROC)){install.packages("pROC");library(pROC)}else{library(pROC)})
suppressMessages(if(!require(readxl)){install.packages("readxl");library(readxl)}else{library(readxl)})
suppressMessages(if(!require(gtools)){install.packages("gtools");library(gtools)}else{library(gtools)})
suppressMessages(if(!require(rlang)){install.packages("rlang");library(rlang)}else{library(rlang)})
suppressMessages(if(!require(shinybusy)){install.packages("shinybusy");library(shinybusy)}else{library(shinybusy)})
suppressMessages(if(!require(caret)){install.packages("caret");library(caret)}else{library(caret)})
suppressMessages(if(!require(broom)){install.packages("broom");library(broom)}else{library(broom)})
suppressMessages(if(!require(epiR)){install.packages("epiR");library(epiR)}else{library(epiR)})





short_info <- function(input, title, place){
  
  input %>% shinyInput_label_embed(shiny_iconlink() %>% bs_embed_tooltip( title = title, placement = place))
}

vectors <- function(n.sint){
  vectors <- gtools::permutations(n=2, r=n.sint, v=c(0,1), repeats.allowed = T)
  vectors <- vectors[nrow(vectors):1, ]
  return(data.frame(vectors))
}


calc_params <- function(x, ...){
  if(length(x[!is.na(x)]) > 1){
    m <- mean(x, ...  )
    desv <- sd(x, ...)
  }else{
    m <- mean(x, ...)
    desv <- m*5
  }
  ret <- paste0(m,";", desv)
  
  return(ret)
}

source("www/helpers.R")



server <- function(input, output, session){
  rv_list <- reactiveValues()
  rv_list$symp_names <- character()
  rv_list$dx_names <- character()
  rv_list$hemo_names <- character()
  rv_list$n.sint    <- numeric()
  rv_list$comb     <- data.frame()
  rv_list$global_parms <- data.frame()# parametros globales de norm dist para hemo vars
  rv_list$thres_vals <- data.frame()# puntos de corte para hemo vars
  rv_list$arbol_only_symp <- data.frame()#arbol de solo sintomas
  rv_list$arbol_hemo_vars <- data.frame() #arbol con variables del hemograma
  rv_list$res_mtrs_only_symp <- data.frame()#metricas globales de solo sintomas
  rv_list$res_hemo_clas <- data.frame()#metricas globales de hemo vars
  rv_list$test_only_symp <- data.frame()#Validacion externa solo sint
  rv_list$test_hemo_vars <- data.frame()#validacion externa hemo vars
  
  
  
  data_in <- reactive({
    req(input$file_in)
    
 
   ext <- str_extract(input$file_in$name, ".[a-zA-Z]{3}$")
   
   if(ext == ".csv"){
     df <- read_csv(input$file_in$datapath)
   }else if(ext ==".xlsx"){
     df <- readxl::read_excel(input$file_in$datapath, sheet = 1)
   }else if(ext == ".xls"){
     df <- readxl::read_xls(input$file_in$datapath, sheet = 1)
   }else{
    
     df <- "error"
   }

   if( "data.frame" %in% class(df) ){
     
     rv_list$dx_names <- names(df)[grepl(names(df), pattern = "dxdengue")]
     dx_name <- as.character(rv_list$dx_names)
     
     
     rv_list$n.sint <- df %>% 
       dplyr::select(1:rv_list$dx_names) %>% 
       ncol()-1
     
     n.sint <- rv_list$n.sint
     
     rv_list$symp_names <- names(df)[1:rv_list$n.sint]
     
     if(ncol(df) >= (rv_list$n.sint+2)){
       rv_list$hemo_names <- names(df)[(rv_list$n.sint+2):ncol(df)]
     }else{
       rv_list$hemo_names <- NA_character_
     }
    
     
     rv_list$comb <- vectors(n.sint) %>%  
       tidyr::unite(., comb, sep = "", remove = F)
     
     df_final <- df %>% 
       mutate_at(vars(1:n.sint), .funs = function(i){ifelse(i >1, 1, i )}) %>% 
       tidyr::unite(., comb, 1:n.sint, sep = "", remove = F ) %>% 
       mutate(id = as.character( 1:nrow(.)))
    
     rv_list$global_parms <- df_final %>%
       group_by(dxdengue0) %>% 
       summarise_at(vars(rv_list$hemo_names), list("mean", "sd"), na.rm = T)
     
   }else{
     df_final <- "error"
   }
   return(df_final)
  
  })
  
  
  output$data_in_prev <- renderDataTable( {
   data_in()
  }, options = list(searching = FALSE, pageLength = 5, scrollX = TRUE))
  
  
  
  output$info_boxes <- renderUI({
   
     
    df <- data_in()
    if("data.frame" %in% class(df)) {
      dx_var <- rlang::sym(rv_list$dx_names)
      prev <- df %>%
        filter(!!dx_var == 1) %>% 
        nrow(.)/nrow(df)
      
      na_count <- df %>% 
        drop_na %>% 
        nrow(.)/nrow(df)
      
      
      fluidRow(
        valueBox("Prevalence", value = paste0(round(prev*100,1),"%"), icon = icon("fas fa-virus") ),
        valueBox("Missing values", value = paste0(round((1-na_count)*100,1),"%"), icon = icon("fas fa-exclamation-triangle"), color = "yellow" ),
        valueBox("symptoms", value = rv_list$n.sint, icon = icon("fas fa-list-ol"), color = "purple" )
      )
    }else{
      sendSweetAlert(
        session = session,
        title = "Error...",
        text = "Invalid input file !",
        type = "error"
      )
    }
    
    
  })
  
  output$hemo_fields <- renderUI({
    df <- data_in()
    if("data.frame" %in% class(df)){
      if(input$hemo_vrs == "Reclassify"){
        
        hemo_nms <- rv_list$hemo_names
        n<- length(hemo_nms)
        fields <- lapply(1:n, function(i){
         fluidRow(
           column(6, 
                  numericInput(inputId = paste0("hemo_vals",i),
                               label   = paste0(hemo_nms[i],":"),
                               value   = 0,
                               min = 0,
                               max = 1000000)),
           column(6, 
                  selectInput(inputId = paste0("symb_", i), 
                              label = "Condition:",
                              choices = c("Less than (<)", "Greather than (>)" ))
                  
                  )
         )
            
         
          
           
         
          
        })
        tagList(fields)
        
      }
      
    }else{
      sendSweetAlert(
        session = session,
        title = "Error...",
        text = "No valid file is loaded !",
        type = "error"
      )
    }
    
    
  })
  
  observeEvent(input$calc_probs, {
 
     
    df <- data_in()
    
    if("data.frame" %in% class(df)){
      show_modal_spinner(spin = "double-bounce",
                         color = "#112446",
                         text = "Caulculating probs...") # show the modal window
      
     
      symp_nms <- rv_list$symp_names 
     
      in_df <- df
      comb <-  rv_list$comb
      names(comb) <- c("comb", symp_nms)
      dx_name <- rlang::sym(rv_list$dx_names ) 
      
      arbol_only_symp <- in_df %>%
        dplyr::select(-id) %>%
        group_by(comb) %>%
        dplyr::summarise(., dx_pos = sum(!!dx_name == 1),
                         dx_neg = sum(!!dx_name == 0),
                         total  = n()) %>%
        ungroup() %>%
        full_join(comb, .) %>%
        dplyr::mutate_all(.funs = function(i){ifelse(is.na(i), 0, i)}) %>%
        dplyr::mutate(score_pos = purrr::pmap_dbl(.l = list(x = dx_pos,y = dx_neg , z = total), .f = function(x,y,z){
          mean(rbeta(10000, 1+x, 1+y))
        }),
        score_neg = purrr::pmap_dbl(.l = list(x = dx_pos,y = dx_neg , z = total), .f = function(x,y,z){
          mean(rbeta(10000, 1+y, 1+x))
        }),
        test_res = case_when(
          score_pos > score_neg ~ "1",
          score_neg > score_pos ~ "0",
          TRUE ~ "2"
        ) )

      rv_list$arbol_only_symp <- arbol_only_symp
      

      predictions <- left_join(in_df, arbol_only_symp %>% 
                                 dplyr::select(comb, score_pos, score_neg, test_res))
      
        
       # global_mtrs <-  table( predictions$test_res,  predictions$dxdengue0) %>% 
       #   confusionMatrix(.) %>% 
       #   broom::tidy(., by_class = T)
       # 
      global_mtrs <- epiR::epi.tests(table( predictions$test_res,  predictions$dxdengue0) , conf.level = 0.95)
     
      croc <- pROC::roc(predictions$dxdengue0, predictions$score_pos)
      ci <- ci.auc(croc)
      
       
      show_res <- data.frame(var = c("Apparent prevalence" ,
                                     "True prevalence" ,
                                     "Sensitivity",
                                     "Specificity",
                                     "Accuracy",
                                     "odds ratio",
                                     "number needed to diagnose",
                                     "Youden index",
                                     "Positive predictive value",
                                     "Negative predictive value",
                                     "Positive likelihood ratio",
                                     "Negative likelihood ratio ") , global_mtrs$rval[1:12] %>% bind_rows) %>% 
        bind_rows(., data.frame(var = "AUC", est = ci[2], lower= ci[1], upper = ci[3]))
      
      names(show_res) <- c("Metric", "Estimated", "IC.lower", "IC.upper")
      
      rv_list$res_mtrs_only_symp <- show_res
      
      k <- 5
      cv_res <- lapply(1:k, function(k){
        
        test <- in_df %>% 
          sample_n(size = round(0.2*nrow(.)))
        
        train <- in_df %>%
          filter(!id %in% test$id) 
        
        only_symp <- train %>%
          select(-id) %>% 
          group_by(comb) %>% 
          dplyr::summarise(., dx_pos = sum(dxdengue0 == 1),
                           dx_neg = sum(dxdengue0 == 0),
                           total  = n()) %>% 
          ungroup() %>% 
          full_join(comb, .) %>%
          dplyr::mutate_all(.funs = function(i){ifelse(is.na(i), 0, i)}) %>% 
          dplyr::mutate(score_pos = purrr::pmap_dbl(.l = list(x = dx_pos,y = dx_neg , z = total), .f = function(x,y,z){
            mean(rbeta(1000, 1+x, 1+y))
          }),
          score_neg = purrr::pmap_dbl(.l = list(x = dx_pos,y = dx_neg , z = total), .f = function(x,y,z){
            mean(rbeta(1000, 1+y, 1+x))
          }),
          test_res = case_when(
            score_pos > score_neg ~ "1",
            score_neg > score_pos ~ "0",
            TRUE ~ "2"
          ) )
        
        test_mts <- left_join(test, only_symp %>% 
                                select(comb, score_pos, score_neg)) %>% 
          dplyr::mutate(test_res = case_when(
            score_pos > score_neg ~ "1",
            score_neg > score_pos ~ "0",
            TRUE ~ "2"
          ) ) 
        
        
        global_mtrs <- epiR::epi.tests(table( test_mts$test_res,  test_mts$dxdengue0) , conf.level = 0.95)
        
        croc <- pROC::roc(test_mts$dxdengue0, test_mts$score_pos)
        ci <- ci.auc(croc)
        
        
        show_res <- data.frame(var = c("Apparent prevalence" ,
                                       "True prevalence" ,
                                       "Sensitivity",
                                       "Specificity",
                                       "Accuracy",
                                       "odds ratio",
                                       "number needed to diagnose",
                                       "Youden index",
                                       "Positive predictive value",
                                       "Negative predictive value",
                                       "Positive likelihood ratio",
                                       "Negative likelihood ratio ") , global_mtrs$rval[1:12] %>% bind_rows) %>% 
          bind_rows(., data.frame(var = "AUC", est = ci[2], lower= ci[1], upper = ci[3])) %>% 
          dplyr::select(var, est)
        
        return(show_res)
      })
      
     
      test_vars_only_symp <- purrr::reduce(cv_res, left_join, by = c("var"="var")) %>% 
        mutate(.,  mean = rowMeans(dplyr::select(.,matches("est"))),  
                   sd = apply(dplyr::select(.,matches("est")),1, sd),
                   lower = mean-1.96*sd,
                   upper =mean+1.96*sd) %>% 
        dplyr::select(-matches("est"), -sd)
      names(test_vars_only_symp) <- c("Metric", "Estimated", "IC.lower", "IC.upper")
      
      rv_list$test_only_symp <- test_vars_only_symp
      
      if(length(na.omit(rv_list$hemo_names)) > 1  ){
        
        hemo_nms <- as.character(rv_list$hemo_names)
        n<- length(hemo_nms)
       
        if(n != 0 | !any(is.na(hemo_nms))){
          if(input$hemo_vrs == "Reclassify"){
            
           
            thres <- lapply(1:n, function(i){
              vals <- input[[ paste0("hemo_vals",i)]]
              varname <- hemo_nms[i]
              return(vals)
            }) %>% 
              bind_cols
            
            names(thres) <- hemo_nms
            
            for(i in 1:n){
              vals <- as.numeric(thres[i])
              varname <- hemo_nms[i]
              to_ev <- input[[paste0("symb_",i)]]
              
              in_df <- in_df %>% 
                dplyr::mutate_at(vars(varname), function(k){
                  
                  ifelse(is.na(k), 
                         NA, 
                         ifelse(
                           (to_ev == "Less than (<)" & k <= vals)|(to_ev == "Greather than (>)" & k >= vals), 
                           "1",
                           "0"
                  ))
               
                  })
            }
           
            
            rv_list$thres_vals <- thres
            
            new_in_df <- in_df %>% 
              dplyr::select(-one_of("comb"), -one_of("id")) %>% 
              dplyr::select(1:rv_list$n.sint, hemo_nms, dx_name)
            
            n.sint <-  ncol(new_in_df)-1
            
            comb <- vectors(n.sint) %>%  
              tidyr::unite(., comb, sep = "", remove = F)
            
            names(comb) <- c("comb",names(new_in_df[1:n.sint]))
            
            new_in_df <- new_in_df %>% 
              mutate_at(vars(1:n.sint), .funs = function(i){ifelse(i >1, 1, i )}) %>% 
              tidyr::unite(., comb, 1:n.sint, sep = "", remove = F ) %>% 
              mutate(id = as.character( 1:nrow(.))) %>% 
              drop_na()
         
            arbol_hemo_vars <- new_in_df %>%
              dplyr::select(-id) %>%
              group_by(comb) %>%
              dplyr::summarise(., dx_pos = sum(!!dx_name == 1),
                               dx_neg = sum(!!dx_name == 0),
                               total  = n()) %>%
              ungroup() %>%
              full_join(comb, .) %>%
              dplyr::mutate_all(.funs = function(i){ifelse(is.na(i), 0, i)}) %>%
              dplyr::mutate(score_pos = purrr::pmap_dbl(.l = list(x = dx_pos,y = dx_neg , z = total), .f = function(x,y,z){
                mean(rbeta(1000, 1+x, 1+y))
              }),
              score_neg = purrr::pmap_dbl(.l = list(x = dx_pos,y = dx_neg , z = total), .f = function(x,y,z){
                mean(rbeta(1000, 1+y, 1+x))
              }),
              test_res = case_when(
                score_pos > score_neg ~ "1",
                score_neg > score_pos ~ "0",
                TRUE ~ "2"
              ) )
           
            rv_list$arbol_hemo_vars <- arbol_hemo_vars
            
            predictions <- left_join(new_in_df, arbol_hemo_vars %>% 
                                       dplyr::select(comb, score_pos, score_neg, test_res))
            
            
            # global_mtrs <-  table( predictions$test_res,  predictions$dxdengue0) %>% 
            #   confusionMatrix(.) %>% 
            #   broom::tidy(., by_class = T)
            # 
            global_mtrs <- epiR::epi.tests(table( predictions$test_res,  predictions$dxdengue0) , conf.level = 0.95)
            
            croc <- pROC::roc(predictions$dxdengue0, predictions$score_pos)
            ci <- ci.auc(croc)
           
            
            show_res_hemo <- data.frame(var = c("Apparent prevalence" ,
                                                "True prevalence" ,
                                                "Sensitivity",
                                                "Specificity",
                                                "Accuracy",
                                                "odds ratio",
                                                "number needed to diagnose",
                                                "Youden index",
                                                "Positive predictive value",
                                                "Negative predictive value",
                                                "Positive likelihood ratio",
                                                "Negative likelihood ratio ") , global_mtrs$rval[1:12] %>% bind_rows) %>% 
              bind_rows(., data.frame(var = "AUC", est = ci[2], lower= ci[1], upper = ci[3]))
            
            names(show_res_hemo) <- c("Metric", "Estimated", "IC.lower", "IC.upper")
            
            rv_list$res_hemo_clas <- show_res_hemo
           
            
            k <- 5
            cv_res_hemo <- lapply(1:k, function(k){
              
              
              test <- new_in_df %>% 
                sample_n(size = round(0.2*nrow(.)))
              
              train <- in_df %>%
                filter(!id %in% test$id) 
              
              only_symp <- train %>%
                select(-id) %>% 
                group_by(comb) %>% 
                dplyr::summarise(., dx_pos = sum(dxdengue0 == 1),
                                 dx_neg = sum(dxdengue0 == 0),
                                 total  = n()) %>% 
                ungroup() %>% 
                full_join(comb, .) %>%
                dplyr::mutate_all(.funs = function(i){ifelse(is.na(i), 0, i)}) %>% 
                dplyr::mutate(score_pos = purrr::pmap_dbl(.l = list(x = dx_pos,y = dx_neg , z = total), .f = function(x,y,z){
                  mean(rbeta(1000, 1+x, 1+y))
                }),
                score_neg = purrr::pmap_dbl(.l = list(x = dx_pos,y = dx_neg , z = total), .f = function(x,y,z){
                  mean(rbeta(1000, 1+y, 1+x))
                }),
                test_res = case_when(
                  score_pos > score_neg ~ "1",
                  score_neg > score_pos ~ "0",
                  TRUE ~ "2"
                ) )
              
              test_mts <- left_join(test, only_symp %>% 
                                      select(comb, score_pos, score_neg)) %>% 
                dplyr::mutate(test_res = case_when(
                  score_pos > score_neg ~ "1",
                  score_neg > score_pos ~ "0",
                  TRUE ~ "2"
                ) ) 
              
              
              global_mtrs <- epiR::epi.tests(table( test_mts$test_res,  test_mts$dxdengue0) , conf.level = 0.95)
              
              croc <- pROC::roc(test_mts$dxdengue0, test_mts$score_pos)
              ci <- ci.auc(croc)
              
              
              show_res <- data.frame(var = c("Apparent prevalence" ,
                                             "True prevalence" ,
                                             "Sensitivity",
                                             "Specificity",
                                             "Accuracy",
                                             "odds ratio",
                                             "number needed to diagnose",
                                             "Youden index",
                                             "Positive predictive value",
                                             "Negative predictive value",
                                             "Positive likelihood ratio",
                                             "Negative likelihood ratio ") , global_mtrs$rval[1:12] %>% bind_rows) %>% 
                bind_rows(., data.frame(var = "AUC", est = ci[2], lower= ci[1], upper = ci[3])) %>% 
                dplyr::select(var, est)
              
              return(show_res)
            })
          
            test_hemo_vars <- purrr::reduce(cv_res_hemo, left_join, by = c("var"="var")) %>% 
              mutate(.,  mean = rowMeans(dplyr::select(.,matches("est"))),  
                     sd = apply(dplyr::select(.,matches("est")),1, sd),
                     lower = mean-1.96*sd,
                     upper =mean+1.96*sd) %>% 
              dplyr::select(-matches("est"), -sd)
            
            names(test_hemo_vars) <- c("Metric", "Estimated", "IC.lower", "IC.upper")
            
            rv_list$test_hemo_vars <- test_vars_only_symp
            
            
          
            }else{
            
            
            only_hemo_vars <- purrr::map( .x = list("1", "0") , .f = function(.x ){
              nms <- ifelse(.x == "1", "pos", 'neg')
              
              
              in_df %>% 
                dplyr::select(comb, dxdengue0:ncol(.)) %>% 
                mutate(dxdengue0 = as.character(dxdengue0), id = as.character(id)) %>% 
                filter(dxdengue0 == .x) %>% 
                group_by(comb) %>% 
                summarise_if(is.numeric, calc_params , na.rm = T) %>% 
                rename_if(is.character, .funs = function(i){
                  ifelse(!grepl("comb", i), paste0(i,"_hemo_", nms), i)
                })
              
            }) %>% purrr::reduce(., full_join)
            
            
            arbol_hemo_vars <- arbol_only_symp %>% 
              full_join(only_hemo_vars) %>% 
              tibble
            
            arbol_hemo_vars_complete <-  arbol_hemo_vars
            
            for(i in 1:length(hemo_nms)){
              
              params_pos <- rv_list$global_parms %>% 
                filter(dxdengue0 == 1) %>% 
                dplyr::select(matches(hemo_nms[i]))
              
              params_neg <- rv_list$global_parms %>% 
                filter(dxdengue0 == 0) %>% 
                dplyr::select(matches(hemo_nms[i]))
              
              arbol_hemo_vars_complete <-  arbol_hemo_vars_complete %>% 
                mutate_at(.vars = vars(paste0(hemo_nms[i],"_hemo_pos")), .funs = function(k){
                  ret <- sapply(k, function(i){
                    if(is.na(i)|grepl("NaN",i)){
                      ret<- paste0(params_pos[1],";", params_pos[2] )
                    }else{
                      ret <- i
                    }
                    return(ret)
                  }, simplify = T)
                  return(ret)
                }) %>% 
                mutate_at(.vars = vars(paste0(hemo_nms[i],"_hemo_neg")), .funs = function(k){
                  ret <- sapply(k, function(i){
                    if(is.na(i)|grepl("NaN",i)){
                      ret<- paste0(params_neg[1],";", params_neg[2] )
                    }else{
                      ret <- i
                    }
                    return(ret)
                  }, simplify = T)
                  return(ret)
                })
            }
            
            rv_list$arbol_hemo_vars <- arbol_hemo_vars_complete
            
            #validacion algoritmo con variables del hemograma
            
            to_val <- in_df %>% 
              left_join(., arbol_hemo_vars, by = c("comb"="comb")) %>%
              mutate_if(.predicate = function(i){ length(table(i))==2 }, as.character) %>% 
              dplyr::select(comb, score_pos, score_neg,  
                            !!hemo_nms, 
                            !!paste0(hemo_nms,"_hemo_pos"),
                            !!paste0(hemo_nms,"_hemo_neg"),
                            dxdengue0)
            
            
            scores <- lapply(hemo_nms, function(i){
              x <- to_val %>% 
                dplyr::select(matches(i)) %>%
                apply(., 1, function(k){
                  
                  q <-  gsub(" ", "", k[1]) %>% as.double()
                  
                  if(!is.na(q)){
                    
                    if(!is.na(k[2]) & !grepl(pattern = "NaN", k[2])){
                      params_dist_pos <- str_split(k[2], ";", simplify = T)
                      
                      m_pos <- as.numeric(params_dist_pos[1])
                      sd_pos <-  as.numeric(params_dist_pos[2])
                      
                      #mx_dens <- dnorm(m_pos, m_pos, sd_pos) 
                      
                      res_pos <- dnorm(q, m_pos, sd_pos, log = F)/dnorm(m_pos, m_pos, sd_pos, log = F)
                      
                    }else if(is.na(k[2])){
                      params <- rv_list$global_parms %>% 
                        filter(dxdengue0 == 1) %>% 
                        dplyr::select(matches(i))
                      
                      res_pos <- dnorm(q, as.numeric(params[1]), as.numeric(params[2]))/dnorm(as.numeric(params[1]), as.numeric(params[1]), as.numeric(params[2]))
                      res_pos <- ifelse(res_pos > 1, 1, res_pos)
                      
                    }else{
                      res_pos<- 0.5
                    }
                    
                    
                    if(!is.na(k[3]) & !grepl(pattern = "NaN", k[3])){
                      params_dist_neg <- str_split(k[3], ";", simplify = T)
                      
                      m_neg <- as.numeric(params_dist_neg[1])
                      sd_neg <-  as.numeric(params_dist_neg[2])
                      
                      #mx_dens <- dnorm(m_neg, m_neg, sd_neg) 
                      
                      res_neg <- dnorm(q, m_neg, sd_neg)/dnorm(m_neg, m_neg, sd_neg)
                      
                      
                    }else if(is.na(k[3])){
                      params <- rv_list$global_parms %>% 
                        filter(dxdengue0 == 0) %>% 
                        dplyr::select(matches(i))
                      
                      res_neg <- dnorm(q, as.numeric(params[1]), as.numeric(params[2]))/dnorm(as.numeric(params[1]), as.numeric(params[1]), as.numeric(params[2]))
                      res_neg <- ifelse(res_neg > 1, 1, res_neg)
                      
                    }else{
                      res_neg<- 0.5
                    }
                    
                    
                  }else{
                    res_pos <- 0.5
                    res_neg<- 0.5
                  }
                  
                  
                  return(data.frame( res_pos, res_neg))
                }) %>% bind_rows()
              names(x) <- c(paste0(i, "_score_pos"), paste0(i, "_score_neg"))
              return(x)
            }) %>% 
              bind_cols()%>% 
              tibble
            
            # names(scores) <- paste0(hemo_vars, "_score")
            
            to_val <- scores %>% 
              dplyr::select(matches("_pos")) %>%
              mutate(score_hemo_pos = apply(., 1, prod)) %>% 
              bind_cols(to_val, .)
            
            to_val <- scores %>% 
              dplyr::select(matches("_neg")) %>% 
              mutate(score_hemo_neg = apply(., 1, prod, na.rm = T)) %>% 
              bind_cols(to_val, .)
            
            
            to_val <- to_val %>% 
              mutate(score_final_pos  = score_hemo_pos * score_pos,
                     score_final_neg  = score_hemo_neg * score_neg,
                     test_res         = ifelse(score_final_pos > score_final_neg, "1", "0") )  
            
            global_mtrs_hemo <- epiR::epi.tests(table( to_val$test_res,  to_val$dxdengue0) , conf.level = 0.95)
            
            croc_hemo <- pROC::roc(to_val$dxdengue0, to_val$score_final_pos)
            ci_hemo <- ci.auc(croc)
            
            
            show_res_hemo <- data.frame(var = c("Apparent prevalence" ,
                                                "True prevalence" ,
                                                "Sensitivity",
                                                "Specificity",
                                                "Accuracy",
                                                "odds ratio",
                                                "number needed to diagnose",
                                                "Youden index",
                                                "Positive predictive value",
                                                "Negative predictive value",
                                                "Positive likelihood ratio",
                                                "Negative likelihood ratio ") , global_mtrs_hemo$rval[1:12] %>% bind_rows) %>% 
              bind_rows(., data.frame(var = "AUC", est = ci_hemo[2], lower= ci_hemo[1], upper = ci_hemo[3]))
            
            names(show_res_hemo) <- c("Metric", "Estimated", "IC.lower", "IC.upper")
            
            rv_list$res_hemo_clas <- show_res_hemo
            
            
            k <-5
            
            cv_res_hemo <- lapply(1:k, function(i){
              
              test <- in_df %>% 
                sample_n(size = round(0.2*nrow(.)))
              
              train <- in_df %>%
                filter(!id %in% test$id) 
              
              only_symp <- train %>%
                dplyr::select(-id) %>% 
                group_by(comb) %>% 
                dplyr::summarise(., dx_pos = sum(dxdengue0 == 1),
                                 dx_neg = sum(dxdengue0 == 0),
                                 total  = n()) %>% 
                ungroup() %>% 
                full_join(comb, ., by = c("comb" = "comb")) %>%
                dplyr::mutate_all(.funs = function(i){ifelse(is.na(i), 0, i)}) %>% 
                dplyr::mutate(score_pos = purrr::pmap_dbl(.l = list(x = dx_pos,y = dx_neg , z = total), .f = function(x,y,z){
                  mean(rbeta(1000, 1+x, 1+y))
                }),
                score_neg = purrr::pmap_dbl(.l = list(x = dx_pos,y = dx_neg , z = total), .f = function(x,y,z){
                  mean(rbeta(1000, 1+y, 1+x))
                }),
                test_res = case_when(
                  score_pos > score_neg ~ "1",
                  score_neg > score_pos ~ "0",
                  TRUE ~ "2"
                ) )
              
              
              only_hemo_vars <- purrr::map( .x = list("1", "0") , .f = function(.x ){
                nms <- ifelse(.x == "1", "pos", 'neg')
                
                train %>% 
                  dplyr::select(comb, dxdengue0:ncol(.)) %>% 
                  mutate(dxdengue0 = as.character(dxdengue0), id = as.character(id)) %>% 
                  filter(dxdengue0 == .x) %>% 
                  group_by(comb) %>% 
                  summarise_if(is.numeric, calc_params , na.rm = T) %>% 
                  rename_if(is.character, .funs = function(i){
                    ifelse(!grepl("comb", i), paste0(i,"_hemo_", nms), i)
                  })
                
              }) %>% purrr::reduce(., full_join)
              
              
              arbol_hemo_vars <- only_symp %>% 
                full_join(only_hemo_vars, by = c("comb" = "comb")) %>% 
                tibble
              
              to_val <- test %>% 
                left_join(., arbol_hemo_vars, by = c("comb" = "comb")) %>%
                mutate_if(.predicate = function(i){ length(table(i))==2 }, as.character) %>% 
                dplyr::select(comb, score_pos, score_neg,  
                              !!hemo_nms, 
                              !!paste0(hemo_nms,"_hemo_pos"),
                              !!paste0(hemo_nms,"_hemo_neg"),
                              dxdengue0)
              
              
              scores <- lapply(hemo_nms, function(i){
                x <- to_val %>% 
                  dplyr::select(matches(i)) %>%
                  apply(., 1, function(k){
                    
                    q <-  gsub(" ", "", k[1]) %>% as.double()
                    
                    if(!is.na(q)){
                      
                      if(!is.na(k[2]) & !grepl(pattern = "NaN", k[2])){
                        params_dist_pos <- str_split(k[2], ";", simplify = T)
                        
                        m_pos <- as.numeric(params_dist_pos[1])
                        sd_pos <-  as.numeric(params_dist_pos[2])
                        
                        #mx_dens <- dnorm(m_pos, m_pos, sd_pos) 
                        
                        res_pos <- dnorm(q, m_pos, sd_pos, log = F)/dnorm(m_pos, m_pos, sd_pos, log = F)
                        
                      }else if(is.na(k[2])){
                        params <- rv_list$global_parms %>% 
                          filter(dxdengue0 == 1) %>% 
                          dplyr::select(matches(i))
                        
                        res_pos <- dnorm(q, as.numeric(params[1]), as.numeric(params[2]))/dnorm(as.numeric(params[1]), as.numeric(params[1]), as.numeric(params[2]))
                        res_pos <- ifelse(res_pos > 1, 1, res_pos)
                        
                      }else{
                        res_pos<- 0.5
                      }
                      
                      
                      if(!is.na(k[3]) & !grepl(pattern = "NaN", k[3])){
                        params_dist_neg <- str_split(k[3], ";", simplify = T)
                        
                        m_neg <- as.numeric(params_dist_neg[1])
                        sd_neg <-  as.numeric(params_dist_neg[2])
                        
                        #mx_dens <- dnorm(m_neg, m_neg, sd_neg) 
                        
                        res_neg <- dnorm(q, m_neg, sd_neg)/dnorm(m_neg, m_neg, sd_neg)
                        
                        
                      }else if(is.na(k[3])){
                        params <-  rv_list$global_parms %>% 
                          filter(dxdengue0 == 0) %>% 
                          dplyr::select(matches(i))
                        
                        res_neg <- dnorm(q, as.numeric(params[1]), as.numeric(params[2]))/dnorm(as.numeric(params[1]), as.numeric(params[1]), as.numeric(params[2]))
                        res_neg <- ifelse(res_neg > 1, 1, res_neg)
                        
                      }else{
                        res_neg<- 0.5
                      }
                      
                      
                    }else{
                      res_pos <- 0.5
                      res_neg<- 0.5
                    }
                    
                    
                    return(data.frame( res_pos, res_neg))
                  }) %>% bind_rows()
                names(x) <- c(paste0(i, "_score_pos"), paste0(i, "_score_neg"))
                return(x)
              }) %>% 
                bind_cols()%>% 
                tibble
              
              
              to_val <- scores %>% 
                dplyr::select(matches("_pos")) %>% 
                mutate(score_hemo_pos = apply(., 1, prod)) %>% 
                bind_cols(to_val, .)
              
              to_val <- scores %>% 
                dplyr::select(matches("_neg")) %>% 
                mutate(score_hemo_neg = apply(., 1, prod)) %>% 
                bind_cols(to_val, .)
              
              
              to_val <- to_val %>% 
                mutate(score_final_pos  = score_hemo_pos * score_pos,
                       score_final_neg  = score_hemo_neg * score_neg,
                       test_res         = ifelse(score_final_pos > score_final_neg, "1", "0"),
                       test_res         = as.factor(test_res),
                       dxdengue0        = as.factor(dxdengue0))  
              
              global_mtrs_hemo <- epiR::epi.tests(table( to_val$test_res,  to_val$dxdengue0) , conf.level = 0.95)
              
              croc_hemo <- pROC::roc(to_val$dxdengue0, to_val$score_final_pos)
              ci_hemo <- ci.auc(croc)
              
              
              show_res_hemo <- data.frame(var = c("Apparent prevalence" ,
                                                  "True prevalence" ,
                                                  "Sensitivity",
                                                  "Specificity",
                                                  "Accuracy",
                                                  "odds ratio",
                                                  "number needed to diagnose",
                                                  "Youden index",
                                                  "Positive predictive value",
                                                  "Negative predictive value",
                                                  "Positive likelihood ratio",
                                                  "Negative likelihood ratio ") , global_mtrs_hemo$rval[1:12] %>% bind_rows) %>% 
                bind_rows(., data.frame(var = "AUC", est = ci_hemo[2], lower= ci_hemo[1], upper = ci_hemo[3]))%>% 
                dplyr::select(var, est)
              
              
              return(show_res_hemo)
            })
            
            test_hemo_vars <- purrr::reduce(cv_res_hemo, left_join, by = c("var"="var")) %>% 
              mutate(.,  mean = rowMeans(dplyr::select(.,matches("est"))),  
                     sd = apply(dplyr::select(.,matches("est")),1, sd),
                     lower = mean-1.96*sd,
                     upper =mean+1.96*sd) %>% 
              dplyr::select(-matches("est"), -sd)
            
            names(test_hemo_vars) <- c("Metric", "Estimated", "IC.lower", "IC.upper")
            
            rv_list$test_hemo_vars <- test_hemo_vars
          }
          
          
        }
      } 
      
      remove_modal_spinner() 
      
    }else{
      sendSweetAlert(
        session = session,
        title = "Error...",
        text = "No valid file is loaded !",
        type = "error"
      )
    }
    
    
    
  })
  
  
  output$only_symp_res<- renderTable({
    req(rv_list$res_mtrs_only_symp )
    rv_list$res_mtrs_only_symp 
  })
  
  output$only_symp_res_val <- renderTable({
    req(rv_list$test_only_symp)
    rv_list$test_only_symp
  })
  
  output$hemo_vars_res <- renderUI({
    req(rv_list$hemo_names)
    hemo_nms <- rv_list$hemo_names
    if(!any(is.na(hemo_nms))){
      fluidRow(title = "Metrics with hemogram variables added",
        column(6,
               tags$h6("Global results"),
              renderTable({rv_list$res_hemo_clas}) 
               ),
        column(6,
               tags$h6("Validation Results"),
               renderTable({rv_list$test_hemo_vars })
               )
        
      )
    }
    
    
  })
  
  
  output$downl_res <- downloadHandler(
    
    filename = function(){
      paste0("arbol_file", Sys.Date(), ".csv")
    }, 
    content = function(con){
      hemo_nms <- rv_list$hemo_names
      if(length(na.omit(hemo_nms)) == 0 ){
       
        write_csv(rv_list$arbol_only_symp, con)
      }else{
       
        write_csv(rv_list$arbol_hemo_vars, con)
      }
    })
  
#####################################
##### ingresar y clasificar individuos
####################################  
  rv_list$n.fields <- data.frame()
  
  arbol_in <- reactive({
    
    ext <- str_extract(input$arbol_in$name, ".[a-zA-Z]{3}$")
    
    if(ext == ".csv"){
      df <- read_csv(input$arbol_in$datapath)
    }else if(ext ==".xlsx"){
      df <- readxl::read_excel(input$arbol_in$datapath, sheet = 1)
    }else if(ext == ".xls"){
      df <- readxl::read_xls(input$arbol_in$datapath, sheet = 1)
    }else{
      
      df <- "error"
    }
    
  
     symp_names <- df %>% 
      dplyr::select(1:dx_pos) %>% 
      dplyr::select(-comb, -dx_pos) %>% 
      names()
    hemo_names <- df %>%
      select(matches("hemo_pos")) %>% 
      names() %>% 
      str_extract(., "[a-zA-Z0-9]+")
     
    if(length(hemo_names) == 0){
      rv_list$n.fields <- data.frame(names = symp_names, type = "symp" )
    }else{
      rv_list$n.fields <- data.frame(names = c(symp_names, hemo_names), 
                                     type = c(rep("symp", length(symp_names)), rep("hemo", length(hemo_names))))
    }
    
    return(df)
  })
  
  
  observe({req(input$arbol_in)
          print(arbol_in())} )
  
  
  output$input_fields <- renderUI({
    req(input$arbol_in)
   
    df <- arbol_in()
   
    if("data.frame" %in% class(df)){
      arbol_names <- rv_list$n.fields
      n <- length(arbol_names$names)
      
      fields <- lapply(1:(n+3), function(i){
        nms <- c(arbol_names$names, "Positive score", "Negative score", "Test outcome")
        if(i <= n & arbol_names$type[i] == "symp"){
          div(style="display:inline-block",
              numericInput(inputId = paste0("input_fields_", i),label = nms[i], min = 0, max =1, step = 1, value = 0, width = "100px")
          )
        }else if(i <= n & arbol_names$type[i] == "hemo"){
          div(style="display:inline-block",
              numericInput(inputId = paste0("input_fields_", i),label = nms[i], min = 0, max =10000000 , value = 0, width = "100px")
          )
        }else{
          div(style="display:inline-block",
              shinyjs::disabled(
                numericInput(inputId = paste0("res_fields_", i), label = nms[i], value = 0, width = "100px")
              )
              
          )
        }
       
       
      })
      
    }
  
    
  })
  
  observeEvent(input$classify_ind, {
   
  arbol <- arbol_in()
  
  if("data.frame" %in% class(arbol)){
   
    arbol_names <- rv_list$n.fields
    n <- length(arbol_names$names)
    
    val <- c()
    for(i in 1:n){
      val[i] <- input[[paste0("input_fields_", i)]]
    }
    
    arbol_vals <- arbol_names %>% 
      mutate(value = val,
             status = ifelse(type == "hemo" & value == 0, "error", "ok")) 
   
    if(any(arbol_vals$status == "error")){
      show_alert(
        title = "Error !!",
        text = "Invalid input values for hemogram variables...",
        type = "error"
      )
    }else{
      
      symp_comb <- arbol_vals %>% 
        dplyr::filter(type == "symp") %>% 
         dplyr::pull(value) %>% 
        paste(., collapse = "")
      

      
      if(any(arbol_vals$type == "hemo")){
        hemo_vals <- arbol_vals %>% 
          dplyr::filter(type == "hemo") %>% 
          dplyr::pull(value)
      }
    
      
      
    }
    
    if(!any(arbol_vals$type == "hemo")){
      sc_pos <- arbol %>% 
        filter(comb == symp_comb) %>% 
        pull(score_pos)
      
      sc_neg <- arbol %>% 
        filter(comb == symp_comb) %>% 
        pull(score_neg)
      
      test_res <- ifelse(sc_pos > sc_neg, 1, 0)
      #output[[paste0("res_fields_", n+1)]] <- renderText("me lapelan")
      updateNumericInput(inputId = paste0("res_fields_", n+1), value =  sc_pos, label = "Score pos")
      updateNumericInput(inputId = paste0("res_fields_", n+2), value =  sc_neg, label = "Score neg")
      updateNumericInput(inputId = paste0("res_fields_", n+3), value =  test_res, label = "Test outcome")
      
     
      
    }else if(any(arbol_vals$type == "hemo")){
      
      score_pos_cont <- arbol %>% 
        filter(comb == symp_comb) 
      
      scores_pos <- c()
      scores_neg <- c()
      hemo_nms <- arbol_vals %>% 
        dplyr::filter(type == "hemo") %>% 
        dplyr::pull(names)
      
      qs <- arbol_vals %>% 
        dplyr::filter(type == "hemo") %>% 
        dplyr::pull(value)
      
      for(i in 1:length(hemo_nms)){
        params_pos <- score_pos_cont %>% 
          pull(sym(paste0(hemo_nms[i],"_hemo_pos"))) %>% 
          str_split(., pattern = ";", simplify = T)
        
        params_neg <- score_pos_cont %>% 
          pull(sym(paste0(hemo_nms[i],"_hemo_neg"))) %>% 
          str_split(., pattern = ";", simplify = T)
       
        scores_pos[i] <- dnorm(qs[i]  ,mean = as.numeric(params_pos[1]), sd = as.numeric(params_pos[2]))/dnorm( as.numeric(params_pos[1])  ,mean = as.numeric(params_pos[1]), sd = as.numeric(params_pos[2]))
        scores_neg[i] <- dnorm(qs[i]  ,mean = as.numeric(params_neg[1]), sd = as.numeric(params_neg[2]))/dnorm( as.numeric(params_neg[1])  ,mean = as.numeric(params_neg[1]), sd = as.numeric(params_neg[2]))
        
      }
      print(scores_pos)
      sc_pos <- arbol %>% 
        filter(comb == symp_comb) %>% 
        pull(score_pos)
      
      sc_neg <- arbol %>% 
        filter(comb == symp_comb) %>% 
        pull(score_neg)
      
      test_res <- ifelse(sc_pos > sc_neg, 1, 0)
     
      updateNumericInput(inputId = paste0("res_fields_", n+1), value =  sc_pos*prod(scores_pos), label = "Score pos")
      updateNumericInput(inputId = paste0("res_fields_", n+2), value =  sc_neg*prod(scores_neg), label = "Score neg")
      updateNumericInput(inputId = paste0("res_fields_", n+3), value =  test_res, label = "Test outcome")
      
      
    }
  
    
  }
    
  })
  
  
############################################
#### ingresar y clasificar base de datos ##
##########################################
  rv_list$n.fields2 <- data.frame()
  
  arbol_in2 <- reactive({
    ext <- str_extract(input$arbol_in2$name, ".[a-zA-Z]{3}$")
    
    if(ext == ".csv"){
      df <- read_csv(input$arbol_in2$datapath)
    }else if(ext ==".xlsx"){
      df <- readxl::read_excel(input$arbol_in2$datapath, sheet = 1)
    }else if(ext == ".xls"){
      df <- readxl::read_xls(input$arbol_in2$datapath, sheet = 1)
    }else{
      
      df <- "error"
    }
    
    
    symp_names <- df %>% 
      dplyr::select(1:dx_pos) %>% 
      dplyr::select(-comb, -dx_pos) %>% 
      names()
    hemo_names <- df %>%
      select(matches("hemo_pos")) %>% 
      names() %>% 
      str_extract(., "[a-zA-Z0-9]+")
    
    if(length(hemo_names) == 0){
      rv_list$n.fields2 <- data.frame(names = symp_names, type = "symp" )
    }else{
      rv_list$n.fields2 <- data.frame(names = c(symp_names, hemo_names), 
                                     type = c(rep("symp", length(symp_names)), rep("hemo", length(hemo_names))))
    }
   
    return(df)
  })
  
 
  
  bd_in <- reactive({
    
    ext <- str_extract(input$bd_in$name, ".[a-zA-Z]{3}$")
    
    if(ext == ".csv"){
      df <- read_csv(input$bd_in$datapath)
    }else if(ext ==".xlsx"){
      df <- readxl::read_excel(input$bd_in$datapath, sheet = 1)
    }else if(ext == ".xls"){
      df <- readxl::read_xls(input$bd_in$datapath, sheet = 1)
    }else{
      
      df <- "error"
    }
    
  
    
    return(df)
  })
  
 
observeEvent(input$classify_db, {
    #req(rv_list$n.fields)
  
   
   
    bd <-  bd_in()
    arbol <- arbol_in2()
    names_control <- rv_list$n.fields2
    
    
    
    symp_names <-  names_control %>% 
      filter(type == "symp") %>% 
      pull(names)
    
    
    bd <- bd %>% 
      tidyr::unite(., comb, 1:length(symp_names), sep = "", remove = F ) 
   
   
    
    if("data.frame" %in% class(arbol) & "data.frame" %in% class(bd)){
      show_modal_spinner(spin = "double-bounce",
                         color = "#112446",
                         text = "Classifying indvs...")
      
      if(!any(names_control$type == "hemo")){
        
        classified <- left_join(bd, arbol %>% 
                    dplyr::select(comb, score_pos, score_neg), by = c("comb" = "comb")) %>% 
          dplyr::mutate(test_res = ifelse(score_pos >= score_neg, "1", "0") )
          
        rv_list$classified_db <- classified
        
      }else if(any(names_control$type == "hemo")){
        
        to_class <- left_join(bd, arbol %>% 
                    dplyr::select(comb, score_pos, score_neg, matches("hemo_pos"), matches("hemo_neg")), by = c("comb" = "comb")) 
          
         var_nms <- names_control %>% 
           filter(type == "hemo") %>% 
           pull(names)
         
        scores <- lapply(var_nms, function(i){
          x <- to_class %>% 
            dplyr::select(matches(i)) %>%
            apply(., 1, function(k){
              
              q <-  gsub(" ", "", k[1]) %>% as.double()
              
              if(!is.na(q)){
                
                if(!is.na(k[2]) & !grepl(pattern = "NaN", k[2])){
                  params_dist_pos <- str_split(k[2], ";", simplify = T)
                  
                  m_pos <- as.numeric(params_dist_pos[1])
                  sd_pos <-  as.numeric(params_dist_pos[2])
                  
                  #mx_dens <- dnorm(m_pos, m_pos, sd_pos) 
                  
                  res_pos <- dnorm(q, m_pos, sd_pos, log = F)/dnorm(m_pos, m_pos, sd_pos, log = F)
                  
                }else{
                  res_pos<- 0.5
                }
                
                
                if(!is.na(k[3]) & !grepl(pattern = "NaN", k[3])){
                  params_dist_neg <- str_split(k[3], ";", simplify = T)
                  
                  m_neg <- as.numeric(params_dist_neg[1])
                  sd_neg <-  as.numeric(params_dist_neg[2])
                  
                  #mx_dens <- dnorm(m_neg, m_neg, sd_neg) 
                  
                  res_neg <- dnorm(q, m_neg, sd_neg)/dnorm(m_neg, m_neg, sd_neg)
                  
                  
                }else{
                  res_neg<- 0.5
                }
                
                
              }else{
                res_pos <- 0.5
                res_neg<- 0.5
              }
              
              
              return(data.frame( res_pos, res_neg))
            }) %>% bind_rows()
          names(x) <- c(paste0(i, "_score_pos"), paste0(i, "_score_neg"))
          return(x)
          
        })%>% 
          bind_cols()%>% 
          tibble
        
        
        to_class <- scores %>% 
          dplyr::select(matches("_pos")) %>% 
          mutate(score_hemo_pos = apply(., 1, prod)) %>% 
          bind_cols(to_class, .)
        
        to_class <- scores %>% 
          dplyr::select(matches("_neg")) %>% 
          mutate(score_hemo_neg = apply(., 1, prod)) %>% 
          bind_cols(to_class, .)
        
        
        to_class <- to_class %>% 
          mutate(score_final_pos  = score_hemo_pos * score_pos,
                 score_final_neg  = score_hemo_neg * score_neg,
                 test_res         = ifelse(score_final_pos > score_final_neg, "1", "0"),
                 test_res         = as.factor(test_res),
                 dxdengue0        = as.factor(dxdengue0))
        
        rv_list$classified_db <- to_class
        
      }
      print(rv_list$classified_db)
      remove_modal_spinner() 
       
    }else{
      show_alert(
        title = "Error !!",
        text = "Error reading database...",
        type = "error"
      )
    }
    
    
  })
 
output$class_fields <- renderUI({
  
  req(rv_list$classified_db)
  tagList(
    dataTableOutput("classified1"),
    shinyWidgets::downloadBttn(outputId = "downl_class",
                               style = "simple",
                               color = "primary")
  )
})


 
output$classified1 <- renderDataTable({
  req(rv_list$classified_db)
  rv_list$classified_db
}, options = list(searching = FALSE, pageLength = 5, scrollX = TRUE))  
  
output$downl_class <- downloadHandler(
  
  filename = function(){
    paste0("classified_bd_", Sys.Date(), ".csv")
  }, 
  content = function(con){
    write_csv(rv_list$classified_db, con)
  })
  
}#end server




