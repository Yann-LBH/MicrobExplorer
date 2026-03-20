#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

# Define server logic required to draw a histogram
function(input, output, session) {
  # Logique pour changer de page au clic du bouton dans la card
  observeEvent(input$go_reads, {
    bslib::nav_select(
      id = "main_navbar",    # L'ID que vous avez mis dans page_navbar
      selected = "reads_page" # Le 'value' de votre nav_panel
    )
  })

  # --- LOGIQUE MÉTIER (Exemple) ---
  output$plot_reads <- renderPlot({
    # Simulation d'un plot
    hist(rnorm(100), col = "#18bc9c", main = paste("Analysis for", input$dataset))
  })
}

# options(shiny.maxRequestSize = 100 * 1024^2)
# 
# library(shiny)
# library(arrow) # Required for parquet files
# 
# # Define server logic required to draw a histogram
# function(input, output, session) {
#   
#   # Create a reactive value to store the loaded data
#   uploaded_data <- reactive({
#     req(input$upload_data) # Wait until a file is uploaded
#     
#     ext <- tools::file_ext(input$upload_data$name)
#     filepath <- input$upload_data$datapath
#     
#     # Logic to handle different extensions
#     switch(ext,
#            rds = readRDS(filepath),
#            parquet = arrow::read_parquet(filepath),
#            validate("Invalid file; Please upload a .rds or .parquet file")
#     )
#   })
#   
#   # Example: Display a preview table of the uploaded data
#   output$preview_table <- renderDT({
#     req(uploaded_data())
#     head(uploaded_data())
#   })