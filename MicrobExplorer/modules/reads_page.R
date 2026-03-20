reads_ui <- function() {
nav_panel(
  title = "Reads",
  value = "reads_page",
  layout_sidebar(
    sidebar = sidebar(
      title = "Reads Controls",
      selectInput("dataset", "Select Data", choices = c("Sample A", "Sample B")),
      hr(),
      helpText("Talk with your personal agent"),
      bs_icon("robot", size = "2em"),
      selectInput("dataset", "LLM", choices = c("Sample A", "Sample B")),
      hr()
    ),
    card(full_screen = TRUE, card_header("Results"), plotOutput("plot_reads"))
  )
)
}