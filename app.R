options(warn=-1)
#lib_dir <- paste0(getwd(),"/libs")
#.libPaths( c( .libPaths(),lib_dir ) )
# p.R175H,p.R248Q,p.R273H,p.R248W, p.R273C, p.R282W, p.G245S
# R175H,R248Q,R273H,R248W,R273C,R282W,G245S
#p.R175H|p.R248Q|p.R273H|p.R248W|p.R273C|p.R282W|p.G245S

print("Loading libraries required...")
list.of.packages <- c("shiny","shinyWidgets","shinyjs","igraph","stringr", "dplyr")
invisible(lapply(list.of.packages, library, character.only = TRUE))

working_directory <- getwd()


source(paste0(working_directory,'//CCLE2.R'))
inputs_dir = paste0(working_directory,"//inputs")
TCGA_TP53_mutations <- as.data.frame(read.table(file = paste0(inputs_dir,"/TCGA_mutations.tsv"), sep = '\t', header = TRUE))
#TCGA_disease <- unique(TCGA_TP53_mutations$Cancer_Type)


TCGA_types <- as.data.frame(read.table(file = paste0(inputs_dir,"/subtype.txt"), sep = '\t', header = TRUE))
subtypes <- sort(unique(TCGA_types$Subtype))
#TCGA_disease <- subtypes


for (i in 1: length(subtypes)) {
  if (grepl( "/", subtypes[i], fixed = TRUE)) {
    print("Adjusting subtype name to be used for folder creation...")
    subtypes[i] <- str_replace_all(subtypes[i], "/", "_")
  } 
}
print(subtypes)
GRN_dir = paste0(inputs_dir, "//", "GRN.graphml")
GRN <- read_graph(GRN_dir, format = "graphml")

# ******************************************************
genes_HGNC  = V(GRN)$name
genes_HGNC <- c(genes_HGNC,"TP53")

genes_HGNC <- unique(genes_HGNC)
genes_HGNC <- str_sort(genes_HGNC)
ui <- fluidPage(
  tags$head(
    tags$style(HTML('#run{background-color:red}'))
  ),
  useShinyjs(),
  setBackgroundColor(
    color = c("#F7FBFF", "#2171B5"),
    gradient = "linear",
    direction = "right"
  ),
  titlePanel("CCLE & TCGA Data Analysis Tool 2021 v.1.0"),
  sidebarLayout(
    sidebarPanel( 
      actionButton("run_ccle", "RUN", icon("play"), 
                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
      #actionButton('run_ccle', 'Run the analysis'),
      shinyjs::hidden(p(id = "running_text", "Processing...")),
      checkboxInput("GLM_signature_user", label = "Use only GLM derived signature?", FALSE),
      checkboxInput("load_other_user", label = "Use TCGA data instead of CCLE (only for GLM and CARNIVAL)", FALSE),
      
      selectInput("violin_gene", "Select gene of interest (usually a TF):", choice = genes_HGNC, selected = "TP53"),
      
      selectInput('selectfile','Select disease from CCLE:',choice = tools::file_path_sans_ext(list.files('./inputs/CELL_LINES/')),selected = "TNBC"),
      radioButtons(inputId="choice", label="Cell line dataset selection:", 
                   choices=list("Run only for the selected disease" = 1, "Run for all cancer cell line samples in single run" = 2, 
                                "Run for all cancer types (cell line .csv files)" = 3),selected = 1),
      
      tags$hr(style="border-color: blue;"),
      checkboxInput("CARNIVAL_flag", label = "Run CARNIVAL optimization", TRUE),
      textInput("top", "Number of top measurements from Dorothea to use (0 = all):", 0),
      checkboxInput("FEM_flag", label = "Use whole expression matrix?", FALSE),
      checkboxInput("FC_flag", label = "Use fold-change (WT control) in expression input?", FALSE),
      textInput("milp_gap", "CPLEX MILP OPTIMALITY GAP:", 0),
      textInput("cpu", "Number of threads to use with CPLEX:", 8),
      textInput("score_user", "Set cut-off similarity score to save networks [0-1]:", 0.8),
      textInput("hotspots_user", "Input the mutation hotspots (protein change, watch naming as different in CCLE/TCGA):", "p.R175H,p.R248Q,p.R273H,p.R248W, p.R273C, p.R282W, p.G245S"),
      tags$hr(style="border-color: blue;"),
      checkboxInput("GLM_flag", label = "Perform Generalized Linear Regression", FALSE),
      radioButtons(inputId="GLM", label="GLM options:", 
                   choices=list("Run Binomial Regression" = 1,"Run Multinomial Regression" = 2),selected = 1),
      
      textInput("choice3", "Enter binomial keyword to classify:", "Missense"),
      checkboxInput("GLM_predict", label = "Use previous best model to predict (no training)", FALSE),
      checkboxInput("GLM_all_u", label = "Use all genes for GLM", TRUE),
      tags$hr(style="border-color: blue;")
      
      
    ),
    mainPanel(
      
      tags$hr(style="border-color: blue;"),
      
      
      radioButtons(inputId="choice2", label="Select mode for TCGA:", 
                   choices=list("Run for only for the selected subtype" = 1, 
                                "Run pancancer in single run" = 2,"Run per disease type (Carnival)" = 3),selected = 3),
      
      
      #checkboxInput("disease_based", label = "Run TCGA pancancer based on disease", FALSE),
      #checkboxInput("subtype_based", label = "Run TCGA pancancer based on subtype", FALSE),
      #checkboxInput("run_all_TCGA", label = "Run TCGA pancancer", FALSE),
      
      #selectInput('TCGA_disease','Select disease from TCGA to run CARNIVAL:',choice = TCGA_disease,selected = "BRCA"),
      selectInput('TCGA_subtype','OR select subtype from TCGA to run CARNIVAL:',choice = subtypes, selected = "BRCA_Basal"),
      
      
      
      
      
      
      tags$hr(style="border-color: blue;"),
      
      
      #"CDC20,CENPA,KIF2C,PLK1"
      
      textInput("signature_user", "Input your custom signature (HGNC symbols separated by commas) for the given TF:", ""),
      
      tags$hr(style="border-color: blue;"),
      checkboxInput("merge_flag", label = "Merge MAPK pathway with the Regulon", FALSE),
      checkboxInput("merge_flag_2", label = "Only use regulons for analysis", TRUE),
      
      tags$hr(style="border-color: blue;"),
      textInput("dataset_version_u", "CCLE_version:", "21Q2"),
      checkboxInput("new_version_u", label = "Read new CCLE version?", FALSE),
      
      img(src='logo.png', height = '149px', width = '763px', align = "center")
      
    )
  )
)

######################################################################################

server <- function(input, output,session) {
  observe({
    shinyjs::toggleState("GLM", input$GLM_flag == "TRUE")
    
    shinyjs::toggleState("choice2", input$load_other_user == "TRUE")
    
    
    shinyjs::toggleState("TCGA_subtype", input$choice2 == 1)
    #shinyjs::toggleState("TCGA_disease", input$choice2 == 1)
    
    
    shinyjs::toggleState("choice3", input$GLM == 1 & input$GLM_flag == TRUE)
    
    shinyjs::toggleState("merge_flag", input$merge_flag_2 == "FALSE")
    shinyjs::toggleState("merge_flag_2", input$merge_flag == "FALSE")
    
    
    shinyjs::toggleState("selectfile", input$choice == 1)
    shinyjs::toggleState("GLM_all_u", input$GLM_flag == "TRUE")
    
    shinyjs::toggleState("top", input$CARNIVAL_flag == "TRUE")
    shinyjs::toggleState("milp_gap", input$CARNIVAL_flag == "TRUE")
    shinyjs::toggleState("hotspots_user", input$CARNIVAL_flag == "TRUE")
    shinyjs::toggleState("FC_flag", input$CARNIVAL_flag == "TRUE")
    
    shinyjs::toggleState("cpu", input$CARNIVAL_flag == "TRUE")
    shinyjs::toggleState("score_user", input$CARNIVAL_flag == "TRUE")
    
    #shinyjs::toggleState("load_other_user", input$CARNIVAL_flag == "TRUE")
    
    
    #shinyjs::toggleState("TCGA_disease", input$TCGA_subtype == "FALSE")
    #shinyjs::toggleState("TCGA_subtype", input$TCGA_disease == "FALSE")
    
    
  })
  v <- reactiveValues()
  plotReady <- reactiveValues(ok = FALSE)
  
  observeEvent(input$run_ccle, {
    
    shinyjs::disable("run_ccle")
    shinyjs::show("running_text")
    sig <- input$signature_user
    
    
    run_all_TCGA <- FALSE
    
    
    top_score <- input$score_user
    TCGA_choice <- 1
    
    if (input$choice2 == 1) {
      TCGA_disease <- sort(c(input$TCGA_subtype))
      TCGA_choice <- 1
    }
    else if (input$choice2 == 2) {
      run_all_TCGA <- TRUE
      TCGA_disease <- c("pancancer")
      TCGA_choice <- 2
      
    }
    else if (input$choice2 == 3) {
      TCGA_disease <- subtypes
      TCGA_choice <- 3
      
    }
    key <- input$choice3
    
    hotspots <- input$hotspots_user
    if (input$load_other_user == TRUE) {
      hotspots_default <- c("R175H","R248Q","R273H","R248W", "R273C", "R282W", "G245S")
    }
    else{ 
      hotspots_default <- c("p.R175H","p.R248Q","p.R273H","p.R248W", "p.R273C", "p.R282W", "p.G245S")
    }
    
    dataset_version <- input$dataset_version_u
    GAP <- as.numeric(input$milp_gap)
    threads <-  as.numeric(input$cpu)
    load_secondary <- input$load_other_user
    top_user <- as.integer(input$top)
    plotReady$ok <- FALSE
    load_other <- FALSE
    c_flag = FALSE  # run CARNIVAL
    M_flag = FALSE # merge initial user pathway with the regulon provided for all analysis
    GLM_all = FALSE # flag to denote whether we use ALL genes for GLM
    new_version = FALSE # if true, then raw csv files are directly read from inputs directory instead of the already saved CCLE_version.Rdata
    M2_flag = FALSE 
    GLM_signature = FALSE
    FC_user =  FALSE
    FEM_user = FALSE
    GLM_signature <- input$GLM_signature_user
    GLM_predict_user <- FALSE
    GLM_predict_user <- input$GLM_predict
    
    if (input$merge_flag == TRUE){M_flag = TRUE}
    if (input$merge_flag_2 == TRUE){M2_flag = TRUE}
    if (input$GLM_all_u == TRUE){GLM_all = TRUE}
    if (input$CARNIVAL_flag == TRUE) {c_flag = TRUE}
    if (input$new_version_u == TRUE) {new_version=TRUE}
    load_other_GLM = FALSE
    if (input$load_other_user == TRUE) {load_other_GLM <- TRUE}
    if (input$FC_flag == TRUE) {FC_user <- TRUE}
    if (input$FEM_flag == TRUE) {FEM_user <- TRUE}
    if (input$GLM_predict == TRUE) {GLM_predict_user <- TRUE}
    
    sig <- stringr::str_replace_all(sig, fixed(" "), "")
    sig <- as.list(unlist(strsplit(sig, ',')))
    sig <- unlist(sig)
    
    hotspots <- stringr::str_replace_all(hotspots, fixed(" "), "")
    hotspots <- as.list(unlist(strsplit(hotspots, ',')))
    hotspots <- unlist(hotspots)
    
    
    if (load_secondary == TRUE) {load_other = TRUE}
    
    if (is.null((sig)))
    {
      sig_user <- NULL
    }
    else 
    {
      sig_user <- sig
    }
    
    
    if (is.null((hotspots)))
    {
      hotspot_user <- hotspots_default
    }
    else 
    {
      hotspot_user <- hotspots
    }
    
    if (input$choice == 1)
    {
      if (input$GLM_flag == TRUE) {
        if (input$GLM == 1) {
          v$res <- CCLE2(tools::file_path_sans_ext(input$selectfile), 
                         input$violin_gene,
                         TRUE,
                         "binomial",
                         c_flag,
                         M_flag,
                         M2_flag,
                         GLM_all,
                         new_version,
                         dataset_version,
                         GAP,
                         threads,
                         sig_user,
                         load_other,top_user,TCGA_disease,run_all_TCGA,hotspot_user,TCGA_choice,key,top_score,load_other_GLM,GLM_signature,FC_user,FEM_user,GLM_predict_user)
          
          setwd(working_directory)
          stopApp(returnValue = invisible())
        }
        else {
          v$res <- CCLE2(tools::file_path_sans_ext(input$selectfile),
                         input$violin_gene,
                         TRUE,
                         "multinomial",
                         c_flag,
                         M_flag,
                         M2_flag,
                         GLM_all,
                         new_version,
                         dataset_version,
                         GAP,
                         threads,
                         sig_user,
                         load_other,top_user,TCGA_disease,run_all_TCGA,hotspot_user,TCGA_choice,key,top_score,load_other_GLM,GLM_signature,FC_user,FEM_user,GLM_predict_user)
          
          setwd(working_directory)
          stopApp(returnValue = invisible())
          
        }
      }
      else {
        v$res <- CCLE2(tools::file_path_sans_ext(input$selectfile),
                       input$violin_gene,
                       FALSE,
                       "binomial",
                       c_flag,
                       M_flag,
                       M2_flag,
                       GLM_all,
                       new_version,
                       dataset_version,
                       GAP,
                       threads,
                       sig_user,
                       load_other,top_user,TCGA_disease,run_all_TCGA,hotspot_user,TCGA_choice,key,top_score,load_other_GLM,GLM_signature,FC_user,FEM_user,GLM_predict_user)
        
        setwd(working_directory)
        stopApp(returnValue = invisible())
      }
      
      
      
    }
    else if (input$choice == 2)
    {
      
      
      if (input$GLM_flag == TRUE) {
        if (input$GLM == 1) {
          v$res <- CCLE2( "all_cell_lines", input$violin_gene, TRUE, "binomial",c_flag,M_flag, M2_flag,GLM_all,
                          new_version,dataset_version,GAP,threads,sig_user,load_other,top_user,TCGA_disease,
                          run_all_TCGA,hotspot_user,TCGA_choice,key,top_score,load_other_GLM,GLM_signature,FC_user,FEM_user,GLM_predict_user)
          setwd(working_directory)
          stopApp(returnValue = invisible())
        }
        else {
          v$res <- CCLE2( "all_cell_lines", input$violin_gene, TRUE, "multinomial",c_flag,M_flag, M2_flag,GLM_all,
                          new_version,dataset_version,GAP,threads,sig_user,load_other,top_user,TCGA_disease,
                          run_all_TCGA,hotspot_user,TCGA_choice,key,top_score,load_other_GLM,GLM_signature,FC_user,FEM_user,GLM_predict_user)
          setwd(working_directory)
          stopApp(returnValue = invisible())
          
        }
      }
      else {
        v$res <- CCLE2( "all_cell_lines", input$violin_gene, FALSE, "multinomial",c_flag,M_flag, M2_flag,GLM_all,
                        new_version,dataset_version,GAP,threads,sig_user,load_other,top_user,TCGA_disease,
                        run_all_TCGA,hotspot_user,TCGA_choice,key,top_score,load_other_GLM,GLM_signature,FC_user,FEM_user,GLM_predict_user)
        setwd(working_directory)
        stopApp(returnValue = invisible())
      }
      
    }
    else if (input$choice == 3)
    {
      
      if (input$GLM_flag == TRUE) {
        if (input$GLM == 1) {
          v$res <- CCLE2( "all", input$violin_gene, TRUE, "binomial",c_flag,M_flag, M2_flag,GLM_all,
                          new_version,dataset_version,GAP,threads,sig_user,load_other,top_user,TCGA_disease,
                          run_all_TCGA,hotspot_user,TCGA_choice,key,top_score,load_other_GLM,GLM_signature,FC_user,FEM_user,GLM_predict_user)
          setwd(working_directory)
          stopApp(returnValue = invisible())
        }
        else {
          v$res <- CCLE2( "all", input$violin_gene, TRUE, "multinomial",c_flag,M_flag, M2_flag,GLM_all,
                          new_version,dataset_version,GAP,threads,sig_user,load_other,top_user,TCGA_disease,
                          run_all_TCGA,hotspot_user,TCGA_choice,key,top_score,load_other_GLM,GLM_signature,FC_user,FEM_user,GLM_predict_user)
          setwd(working_directory)
          stopApp(returnValue = invisible())
          
        }
      }
      else {
        v$res <- CCLE2( "all", input$violin_gene, FALSE, "binomial",c_flag,M_flag, M2_flag,GLM_all,
                        new_version,dataset_version,GAP,threads,sig_user,load_other,top_user,TCGA_disease,
                        run_all_TCGA,hotspot_user,TCGA_choice,key,top_score,load_other_GLM,GLM_signature,FC_user,FEM_user,GLM_predict_user)
        setwd(working_directory)
        stopApp(returnValue = invisible())
      }
      
    }
    
  })
  
  Sys.sleep(2)
  plotReady$ok <- TRUE
  
}

######################################################################################

shinyApp(ui = ui, server = server)
