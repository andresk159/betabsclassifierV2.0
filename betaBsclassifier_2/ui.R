# Shiny web application for classify denge disease patients
# Author: Andres Camilo Mendez

#cargar paquetes

suppressMessages(if(!require(shiny)){install.packages("shiny");library(shiny)}else{library(shiny)})
suppressMessages(if(!require(tidyverse)){install.packages("tidyverse");library(tidyverse)}else{library(tidyverse)})
suppressMessages(if(!require(rlang)){install.packages("rlang");library(rlang)}else{library(rlang)})
suppressMessages(if(!require(sf)){install.packages("sf");library(sf)}else{library(sf)})
suppressMessages(if(!require(shinycssloaders)){install.packages("shinycssloaders");library(shinycssloaders)}else{library(shinycssloaders)})
suppressMessages(if(!require(googleVis)){install.packages("googleVis");library(googleVis)}else{library(googleVis)})
suppressMessages(if(!require(XML)){install.packages("XML");library(XML)}else{library(XML)})
suppressMessages(if(!require(highcharter)){install.packages("highcharter");library(highcharter)}else{library(highcharter)})
suppressMessages(if(!require(RColorBrewer)){install.packages("RColorBrewer");library(RColorBrewer)}else{library(RColorBrewer)})
suppressMessages(if(!require(collapsibleTree)){install.packages("collapsibleTree");library(collapsibleTree)}else{library(collapsibleTree)})
suppressMessages(if(!require(shinydashboardPlus)){install.packages("shinydashboardPlus");library(shinydashboardPlus)}else{library(shinydashboardPlus)})
suppressMessages(if(!require(shinydashboard)){install.packages("shinydashboard");library(shinydashboard)}else{library(shinydashboard)})
suppressMessages(if(!require(shinyWidgets)){install.packages("shinyWidgets");library(shinyWidgets)}else{library(shinyWidgets)})
suppressMessages(if(!require(googleway)){install.packages("googleway");library(googleway)}else{library(googleway)})
suppressMessages(if(!require(bsplus)){install.packages("bsplus");library(bsplus)}else{library(bsplus)})
suppressMessages(if(!require(shinyjs)){install.packages("shinyjs");library(shinyjs)}else{library(shinyjs)})


short_info <- function(input, title, place){
  
  input %>% shinyInput_label_embed(shiny_iconlink() %>% bs_embed_tooltip( title = title, placement = place))
}

source("www/helpers.R", local = TRUE)



ui <- dashboardPage(
  dashboardHeader(title = "Calculadora Dengue V.2.0", titleWidth = 400),
  dashboardSidebar(disable = T),
    dashboardBody(
      useSweetAlert(),
      useShinyjs(),
      fluidRow(
        column(4,
               box(title = "Tool box",
                   solidHeader = F,
                   width = 12,
                 fileInput("file_in", "Load file",
                           multiple = FALSE,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv",
                                      ".xlsx")),
                 prettyRadioButtons(
                   inputId = "hemo_vrs",
                   label = "Hemogram variables:", 
                   choices = c("Use continous values", "Reclassify"),
                   inline = TRUE, 
                   status = "danger",
                   fill = TRUE
                 ),
                 uiOutput("hemo_fields"),
                 actionBttn(
                   inputId = "calc_probs",
                   label = "Calculate",
                   style = "simple",
                   color = "primary"
                 )
               )
              
               ),
        column(8, 
               tabBox(
                 id = "tabPanel_1",
                 width = 12,
                 height = "100%",
                 tabPanel(title = "Calculate probs", value = "tab1",
                          fluidRow(
                            box(title = "Preview database",
                                collapsible = T,
                                width = 12,
                                dataTableOutput("data_in_prev"),
                                uiOutput("info_boxes")     
                            ),
                            box(title = "results",
                                collapsible = F,
                                width =12,
                                tags$h4("Only symptoms results metrics"),
                                fluidRow(
                                  column(6,
                                         tags$h6("Global results"),
                                         tableOutput("only_symp_res")
                                         ),
                                  column(6,
                                         tags$h6("Validation Results"),
                                         tableOutput("only_symp_res_val")
                                         )
                                 
                                ),
                                tags$h4("Hemogram variables results"),
                                uiOutput("hemo_vars_res"),
                                shinyWidgets::downloadBttn(outputId = "downl_res",
                                                           style = "simple",
                                                            color = "primary")
                                )
                          )
                         
                             
                 ),
                 tabPanel(title = "Classify indvs", value = "tab2",
                          fluidRow(
                            box(title = "Classify inndividuals",
                                width = 12,
                                collapsible = T,
                                clollapsed = T,
                                
                                fileInput("arbol_in", "Load file",
                                          multiple = FALSE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv",
                                                     ".xlsx")),
                                
                                
                                uiOutput("input_fields"),
                                actionBttn(inputId = "classify_ind", 
                                           label = "Classify", 
                                           style = "simple", 
                                           color = "primary")
                                
                                
                            )
                          )
                          
                 ),
                 tabPanel(title = "Classify data base", value = "tab3",
                          fluidRow(
                            box(title= "Classify data base",
                                width = 12,
                                collapsible = T,
                                clollapsed = T,
                                fileInput("arbol_in2", "Load Arbol file",
                                          multiple = FALSE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv",
                                                     ".xlsx")),
                                fileInput("bd_in", "Load database file",
                                          multiple = FALSE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv",
                                                     ".xlsx")),
                                actionBttn(inputId = "classify_db", 
                                           label = "Classify", 
                                           style = "simple", 
                                           color = "primary"),
                                uiOutput("class_fields")
                            )
                            )
                          )
               )
      )
     
    
       
     )
      
) 
    
  
  
  
)



