mod_benchmarks_ui <- function(id) {
  ns <- NS(id)
  nav_panel(title = "Benchmarks", icon = icon("tachometer-alt"), card(p("Statistiques de calcul.")))
}

mod_benchmarks_server <- function(id) {
  moduleServer(id, function(input, output, session) {})
}