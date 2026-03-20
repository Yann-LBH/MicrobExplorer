library(ggplot2)
library(ggimage) # Pour intégrer votre logo dans le ggplot
library(hexSticker)

# 1. Création du subplot avec ggplot2
# Cela permet un contrôle total sur l'interligne et les positions
p <- ggplot() +
  
  # AJOUT DU LOGO
  # On décale légèrement le logo vers la gauche (x = 0.5) 
  # et on le remonte un peu pour l'équilibre (y = 1.1)
  geom_image(aes(x = 4, y = 4.5, image = "MicrobExplorer/www/logo_shiny_ME.png"), size = 0.9) +
  
  # AJOUT DU TEXTE
  # hjust = 0 : Aligne le texte à GAUCHE sur la coordonnée x (sur la première lettre 'M')
  # x = 0.8 : Le point de départ horizontal du texte
  # y = 1.1 : Aligné verticalement avec le logo
  # lineheight = 0.4 : Interligne très serré, comme demandé
  geom_text(aes(x = 5.5, y = 6, label = "Microb\nExplorer"), 
            size = 14,             # Légèrement réduit pour l'équilibre
            family = "sans", 
            fontface = "bold", 
            lineheight = 0.4,
            hjust = 0) + 
  
  theme_void() + 
  
  # --- LA CORRECTION CLÉ N°1 ---
  # On donne BEAUCOUP plus d'espace (un cadre de 10x10) pour éviter que 
  # les éléments ne touchent les bords de la zone de dessin.
  xlim(0, 10) + ylim(0, 10) # Cadre de travail agrandi

# --- 2. Génération du sticker final ---
sticker(
  subplot = p,                # On utilise l'objet ggplot créé au-dessus
  package = "",               # Texte déjà inclus dans le ggplot
  s_x = 1, s_y = 1,           # Centre le subplot
  
  # --- LA CORRECTION CLÉ N°2 ---
  # On agrandit la zone du subplot pour qu'elle remplisse mieux l'hexagone.
  # Cela réduit la zone de "retrait" autour des éléments.
  s_width = 1.7, s_height = 1.7,
  
  h_fill = "white",
  h_color = "#18bc9c",
  filename = "MicrobExplorer/www/logo.png"
)