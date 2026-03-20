# Main Landing Page
home_ui <- function() {
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
  
  br(),
  
  # Choice Blocks (Action Cards)
  layout_column_wrap(
    width = 1/3, # 3 cards per row
    fixed_height = 200,
    
    # Card 1: Reads
    card(
      card_header(bs_icon("text-paragraph", size = "2em"), p("Reads")),
      card_body(
        p("Take a look at your clean and trimm data for reads")
      ),
      card_footer(
        actionButton("go_reads", "Explore", class = "btn-primary w-100")
      )
    ),
    
    # Card 2: Contigs
    card(
      card_header(bs_icon("dash-lg", size = "2em"), p("Contigs")),
      card_body(
        p("Take a look at your clean and trimm data for contigs")
      ),
      card_footer(
        actionButton("go_contigs", "Explore", class = "btn-primary w-100")
      )
    ),
    
    # Card 3: IA Model
    card(
      card_header(bs_icon("robot", size = "2em"), p("IA Model")),
      card_body(
        p("Try our random forest model")
      ),
      card_footer(
        actionButton("go_ia", "Use", class = "btn-primary w-100")
      )
    )
  )
)
}