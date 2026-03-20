#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

ui <- page_navbar(
  id = "main_navbar",
  title = tags$a(href="https://github.com/Yann-LBH/MicrobExplorer", target="_blank",
  tagList(
    img(
      src = "logo.png", 
      height = "125px", 
      style = "margin-right: 20px; vertical-align: middle;"
    ),
    tags$span("MicrobExplorer", 
              style = "font-size: 32px; font-weight: bold; color: white; vertical-align: middle; margin-right: 20px;")
  )
    ),
  
  #Invisible tab
  header = tags$head(
    tags$style(HTML("
      /* Masquer les pages de module dans la barre du haut */
      .navbar-nav .nav-link[data-value='reads_page'],
      .navbar-nav .nav-link[data-value='contigs_page'],
      .navbar-nav .nav-link[data-value='ia_page'] {
        display: none !important;
      }
      /* Optionnel : Rendre la barre latérale plus élégante */
      .sidebar { background-color: #f8f9fa !important; }
    "))
  ),
  
  theme = bs_theme(version = 5, 
                   bootswatch = "flatly", 
                   heading_font = "sans") %>% 
  bs_add_variables(
    "nav-link-font-size" = "1.2rem",   # Taille des onglets (Home, etc.)
    "navbar-padding-y" = "1rem" # Espacement interne de la barre
  ),

  home_ui(),    # Appel de la fonction définie dans modules/home_page.R
  reads_ui(),   # Appel de la fonction définie dans modules/reads_page.R
  export_ui()   # Appel de la fonction définie dans modules/export_page.R
)