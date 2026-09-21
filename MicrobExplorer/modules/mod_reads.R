mod_reads_ui <- function(id) {
  ns <- NS(id)
  
  nav_panel(
    title = "Reads",
    value = "reads_page",
    layout_sidebar(
      sidebar = sidebar(
        title = "Reads Controls",
        # 1. Added ns() to inputs
        selectInput(ns("tax_level"), "Niveau taxonomique :", 
                    choices = c("Phylum", "Class", "Order", "Family", "Genus")),
        sliderInput(ns("top_n"), "Top N :", min = 5, max = 50, value = 10),
        hr(),
        sliderInput(ns("title_size"), "Taille titre :", min = 8, max = 24, value = 14),
        sliderInput(ns("axis_size"), "Taille axes :", min = 6, max = 18, value = 10),
        hr(),
        helpText("Talk with your personal agent"),
        bs_icon("robot", size = "2em"),
        # 2. Fixed duplicate ID 'dataset' and added ns()
        selectInput(ns("llm_model"), "LLM :", choices = c("GPT-4", "Claude 3.5"))
      ),
      
      # 3. Plot rendered inside the main card
      card(
        full_screen = TRUE, 
        card_header("Graphique Reads"), 
        plotOutput(ns("plot"), height = "500px")
      )
    )
  )
}

mod_reads_server <- function(id, data_master) {
  moduleServer(id, function(input, output, session) {
    output$plot <- renderPlot({ plot(1:10, main = "Reads Plot") })
    
    return(reactive({
      list(
        taxonomic_level = input$tax_level,
        top_n           = input$top_n,
        title_size      = input$title_size,
        axis_size       = input$axis_size
      )
    }))
  })
}