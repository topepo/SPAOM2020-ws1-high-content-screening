library(tidymodels)
library(readr)
library(janitor)
library(magick)
library(learntidymodels)
library(embed)
library(ggforce)
library(plotly)

theme_set(theme_bw())

# ------------------------------------------------------------------------------

unneeded <- c('image_number', 'object_number', 'metadata_file_location',
              'metadata_frame', 'metadata_series', 
              'file_name_raw', 'path_name_raw', 'number_object_number',
              'euler_number', 'image_path')

process_image <- function(x) {
  x <- image_read(x)
  x <- image_modulate(x, brightness = 180, saturation = 180)
  x <- image_resize(x, "500x500")
  as.raster(x)
}

high_content_screen <-
  read_csv("dataset2D/output/objects.csv") %>%
  clean_names() %>%
  rename(id = metadata_colrow, treatment = metadata_treatment) %>% 
  relocate(id) %>% 
  mutate(image_path = paste0("dataset2D/output/", id, "--outline.png")) %>% 
  rename_with(~ gsub("^area_shape_", "", .x)) %>% 
  rename_with(~ gsub("^intensity_", "", .x)) %>% 
  rename_with(~ gsub("^radial_distribution_", "", .x)) %>% 
  mutate(images = map(image_path, process_image)) %>% 
  mutate(
    treatment = case_when(
      treatment == "matrigel5" ~ "matrigel",
      treatment == "pre2i" ~ "pre2",
      treatment == "pre2i" ~ "pre2 (cleaned)",
      TRUE ~ "control"
    ),
    treatment = factor(treatment)
    ) %>% 
  select(-all_of(unneeded))

names(high_content_screen$images) <- high_content_screen$id
images <- high_content_screen$images
high_content_screen$images <- NULL

# View(high_content_screen)

# ------------------------------------------------------------------------------

cell_rec <- 
  recipe(treatment ~ ., data = high_content_screen) %>% 
  update_role(id, new_role = "id") %>% 
  step_YeoJohnson(all_numeric_predictors()) %>% 
  step_normalize(all_numeric_predictors()) %>% 
  prep()

# ------------------------------------------------------------------------------

scatter_plots <- function(object, plotly = TRUE, title = NULL) {
  if (plotly) {
    p <- 
      bake(object, new_data = NULL) %>% 
      ggplot(aes(x = .panel_x, y = .panel_y, fill = treatment, colour = treatment)) +
      facet_matrix(vars(-treatment, -id), layer.diag = 2) + 
      geom_point(alpha = 0.5, position = 'auto', aes(text = id)) + 
      ggtitle(title)
    p <- ggplotly(p) 
  } else {
    p <- 
      bake(object, new_data = NULL) %>% 
      ggplot(aes(x = .panel_x, y = .panel_y, fill = treatment, colour = treatment)) +
      facet_matrix(vars(-treatment, -id), layer.diag = 2) + 
      geom_point(alpha = 0.5, position = 'auto') + 
      ggtitle(title)
  }
  p
}

# ------------------------------------------------------------------------------

pca_rec <- 
  cell_rec %>% 
  step_pca(all_numeric_predictors(), num_comp = 3) %>% 
  prep()

pca_rec %>% scatter_plots("PCA")

plot_loadings(pca_rec, component_number <= 3)

pca_rec %>% 
  bake(new_data = NULL) %>% 
  filter(PC1 > 5)
plot(images[["06A"]])

pca_rec %>% 
  bake(new_data = NULL) %>% 
  filter(PC2 < -5)
plot(images[["05E"]])

# ------------------------------------------------------------------------------


pls_rec <- 
  cell_rec %>% 
  step_pls(all_numeric_predictors(), outcome = vars(treatment), num_comp = 3) %>% 
  prep()

pls_rec %>% scatter_plots("PLS")

# ------------------------------------------------------------------------------

ica_rec <- 
  cell_rec %>% 
  step_ica(all_numeric_predictors(), num_comp = 3) %>% 
  prep()

ica_rec %>% scatter_plots("ICA")

# ------------------------------------------------------------------------------

mds_rec <- 
  cell_rec %>% 
  step_isomap(all_numeric_predictors(), num_terms = 3) %>% 
  prep()

mds_rec %>% scatter_plots("MDS")

# ------------------------------------------------------------------------------

umap_rec <- 
  cell_rec %>% 
  step_umap(all_numeric_predictors(), num_comp = 3) %>% 
  prep()

umap_rec %>% scatter_plots("UMAP") 


umap_sup_rec <- 
  cell_rec %>% 
  step_umap(all_numeric_predictors(), outcome = vars(treatment), num_comp = 3) %>% 
  prep()

umap_sup_rec %>% scatter_plots("Supervised UMAP")

umap_sup_rec %>% 
  bake(new_data = NULL) %>% 
  filter(umap_3 > -5 & treatment == "control")

plot(images[["01D"]])
plot(images[["02F"]])
plot(images[["01H"]])

plot(images[["05E"]])
plot(images[["06A"]])
