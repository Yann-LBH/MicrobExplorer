mod_qc_ui <- function(id) {
  ns <- NS(id)
  nav_panel(title = "QC", icon = icon("check-circle"), card(p("Métrique de qualité.")))
}

mod_qc_server <- function(id) {
  moduleServer(id, function(input, output, session) {})
}