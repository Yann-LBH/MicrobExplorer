# Other Navigation Tabs (Placeholders)
export_ui <- function(id) {
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

mod_export_ui <- function(id) {
  ns <- NS(id)
  downloadButton(ns("btn_export"), "Exporter la config (pipeline_config.json)", class = "btn-success")
}

mod_export_server <- function(id, reads_config, contigs_config, kegg_config) {
  moduleServer(id, function(input, output, session) {
    
    output$btn_export <- downloadHandler(
      filename = function() {
        paste0("pipeline_config_", Sys.Date(), ".json")
      },
      content = function(file) {
        config_data <- list(
          pipeline_version = "1.0",
          export_date      = as.character(Sys.Date()),
          parameters       = list(
            reads   = reads_config(),
            contigs = contigs_config(),
            kegg    = kegg_config()
          )
        )
        write_json(config_data, file, pretty = TRUE, auto_unbox = TRUE)
      }
    )
  })
}