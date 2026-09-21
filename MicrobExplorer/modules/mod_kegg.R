mod_kegg_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "KEGG", icon = icon("chart-pie"),
    layout_sidebar(
      sidebar = sidebar(
        selectInput(ns("tax_level"), "Niveau fonctionnel :", choices = c("Pathways", "Modules")),
        sliderInput(ns("top_n"), "Top N :", min = 5, max = 50, value = 10),
        sliderInput(ns("title_size"), "Taille titre :", min = 8, max = 24, value = 14),
        sliderInput(ns("axis_size"), "Taille axes :", min = 6, max = 18, value = 10)
      ),
      card(card_header("Graphique KEGG"), plotOutput(ns("plot")))
    )
  )
}

mod_kegg_server <- function(id, data_master) {
  moduleServer(id, function(input, output, session) {
    output$plot <- renderPlot({ plot(1:10, main = "KEGG Plot") })
    
    return(reactive({
      list(
        functional_level = input$tax_level,
        top_n            = input$top_n,
        title_size       = input$title_size,
        axis_size        = input$axis_size
      )
    }))
  })
}