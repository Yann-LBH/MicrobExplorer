# Other Navigation Tabs (Placeholders)
export_ui <- function() {
nav_panel(
  "Export Style", 
  icon = bs_icon("floppy"), 
  
  # Header Section
  layout_column_wrap(
    width = 1,
    card(
      card_header("Save your Style !"),
      p("You can save your graph style and add it to the pipeline to get perfect plot from the next run.\n
                  For more informations visit the github page of the project :"), 
      tags$a(href="https://github.com/Yann-LBH/MicrobExplorer", target="_blank", "Click here!"),
      full_screen = FALSE
    ),
    card(
      card_header("Style Export"),
      icon = bs_icon("download"),
      fileInput(
        inputId = "download_style", 
        label = "Choose RDS or Parquet file",
        multiple = FALSE,
        accept = c(".rds", ".parquet")
      ),
      helpText("Download your style.")
    )
  )
)
}