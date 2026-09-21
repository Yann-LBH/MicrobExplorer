# Main Landing Page
mod_home_ui <- function(id) {
  ns <- NS(id)
nav_panel(
  title = "Home",
  value = "home_page",
  icon = bs_icon("house"),
  
  # Header Section
  layout_column_wrap(
    width = 1,
    card(
      card_header("Welcome to MicrobExplorer"),
      p("Select a module below to start visualizing your post-processed omics data.\n
          For more informations visit the github page of the project :"), 
      tags$a(href="https://github.com/Yann-LBH/MicrobExplorer", target="_blank", "Click here!"),
      full_screen = FALSE
    ),
    card(
      card_header("Data Import"),
      fileInput(
        inputId = "upload_data", 
        label = "Choose RDS or Parquet file",
        multiple = FALSE,
        accept = c(".rds", ".parquet")
      ),
      helpText("Upload your post-processed omics data.")
    )
  ),
)
}

mod_home_server <- function(id) {
  moduleServer(id, function(input, output, session) {})
}