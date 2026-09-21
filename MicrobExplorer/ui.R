ui <- page_navbar(
  title = tags$a(
    href = "https://github.com/Yann-LBH/MicrobExplorer", 
    target = "_blank",
    style = "text-decoration: none;",
    tagList(
      img(
        src = "logo.png", 
        height = "100px", # Adjusted height to fit standard navbar
        style = "margin-right: 10px; vertical-align: middle;"
      ),
      tags$span(
        "MicrobExplorer", 
        style = "font-size: 22px; font-weight: bold; color: white; vertical-align: middle;"
      )
    )
  ),
  theme = bs_theme(
    version = 5, 
    bootswatch = "flatly", 
    heading_font = font_google("Inter")
  ),
  id = "main_navbar",
  
  # Navigation panels (tabs) using bslib structure
  mod_home_ui("home"),
  mod_benchmarks_ui("benchmarks"),
  mod_qc_ui("qc"),
  mod_reads_ui("reads"),
  mod_contigs_ui("contigs"),
  mod_kegg_ui("kegg"),
  
  # Footer for global settings/export
  footer = tagList(
    hr(),
    div(
      class = "p-3 text-center",
      mod_export_ui("export_config")
    )
  )
)