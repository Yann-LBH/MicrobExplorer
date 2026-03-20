# Créer une liste contenant tous les objets
global_list <- list(
  ventes = readRDS("data/ventes.rds"),
  clients = readRDS("data/clients.rds"),
  coords = readRDS("data/map_coords.rds")
)

# Sauvegarder la liste entière
saveRDS(global_list, "./././results/rds/app_data.rds")