# Environment preparation ####
rm(list = ls()) # Reset R`s brain
# Set working directory
setwd("~/git/inaturalist_ua_paper")

# Load libraries
library(tidyverse)
library(sf)
library(knitr)
library(patchwork)
library(treemap)

# Download Ukraine regions as an `sf` polygon
# Uncomment following lines and run it only for the first time
# library(remotes)
# remotes::install_github("dickoa/rgeoboundaries")
# library(rgeoboundaries)
# ukraine_poly <- rgeoboundaries::gb_adm1("UKR")
# save(ukraine_poly, file = "Ukraine_poly_adm1.Rdata")

load(file = "Ukraine_poly_adm1.Rdata")   # regions (oblast)

# Set custom theme for ggplot2 plots
mytheme <- theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# For colorizing density scatterplots
library(RColorBrewer)
# Define palette
paletteForUse <- c('#d10000', '#ff6622', '#ffda21', '#33dd00', '#1133cc', '#220066', '#330044')
colors <- colorRampPalette(paletteForUse)(256)

# Load data ####
# Specify previously saved iNaturalist data, as *.csv file.
# All fields must be checked before downloading.
inat_data <- "observations-391636.csv" # downloaded 2024-01-03

# Load data and convert it to `sf`
fungi <- read.csv(inat_data) %>% 
  # select required columns
  select(id,
         observed_on,
         user_id,
         user_name,
         quality_grade,
         url,
         num_identification_agreements,
         num_identification_disagreements,
         latitude,
         longitude,
         coordinates_obscured,
         place_admin1_name,
         place_admin2_name,
         scientific_name,
         taxon_id,
         taxon_phylum_name,
         taxon_subphylum_name,
         taxon_class_name,
         taxon_order_name,
         taxon_family_name,
         taxon_genus_name,
         taxon_species_name,
         taxon_subspecies_name,
         taxon_variety_name,
         taxon_form_name) %>% 
  mutate_at("id", as.character) %>% 
  mutate_at("observed_on", as.Date) %>% 
  filter(observed_on <= "2023-12-31") %>% 
  mutate_at(c(
    "user_name",
    "quality_grade",
    "coordinates_obscured",
    "place_admin1_name",
    "place_admin2_name",
    "scientific_name",
    "taxon_phylum_name",
    "taxon_subphylum_name",
    "taxon_class_name",
    "taxon_order_name",
    "taxon_family_name",
    "taxon_genus_name",
    "taxon_species_name",
    "taxon_subspecies_name",
    "taxon_variety_name",
    "taxon_form_name"
  ), as.factor) %>% 
  mutate(place_admin1_name = recode(place_admin1_name,
                                    "Crimea" = "Autonomous Republic of Crimea", # "old one" = "new one"
                                    "Cherkasy" = "Cherkasy Oblast",
                                    "Chernihiv" = "Chernihiv Oblast",
                                    "Chernivtsi" = "Chernivtsi Oblast",
                                    "Dnipropetrovs'k" = "Dnipropetrovsk Oblast",
                                    "Donets'k" = "Donetsk Oblast",
                                    "Kharkiv" = "Kharkiv Oblast",
                                    "Kherson" = "Kherson Oblast",
                                    "Khmel'nyts'kyy" = "Khmelnytskyi Oblast",
                                    "Kirovohrad" = "Kirovohrad Oblast",
                                    "Kiev" = "Kyiv Oblast",
                                    "Kiev City" = "Kyiv",
                                    "Luhans'k" = "Luhansk Oblast",
                                    "L'viv" = "Lviv Oblast",
                                    "Ivano-Frankivs'k" =  "Ivano-Frankivsk Oblast",
                                    "Mykolayiv" = "Mykolaiv Oblast",
                                    "Odessa" = "Odessa Oblast",
                                    "Poltava" = "Poltava Oblast",
                                    "Rivne" = "Rivne Oblast",
                                    "Sevastopol'" = "Sevastopol",
                                    "Sumy" = "Sumy Oblast",
                                    "Ternopil'" = "Ternopil Oblast",
                                    "Vinnytsya" = "Vinnytsia Oblast",
                                    "Volyn" = "Volyn Oblast",
                                    "Transcarpathia" = "Zakarpattia Oblast",
                                    "Zaporizhzhya" = "Zaporizhia Oblast",
                                    "Zhytomyr" = "Zhytomyr Oblast"
  )
  ) %>%
  # convert to `simple features` spatial object
  st_as_sf(dim = "XY", remove = FALSE, na.fail = F, 
           coords = c("longitude", "latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs") %>% 
  # clip by country polygon
  st_filter(ukraine_poly) %>% 
  # drop records without coordinates
  filter(!st_is_empty(.))

# How many records?
nrow(fungi)

# Input statistics ####
num_obs <- nrow(fungi)                                   # total number of observations
num_observers <- nlevels(fungi$user_name)                # how many observers?
# obscured_coords <- nrow(fungi[fungi$coordinates_obscured == "true",]) # obscured coordinates
# scientificNames <- length(unique(fungi$scientific_name)) # scientific names
species <- length(unique(fungi$taxon_species_name))      # Number of species
genera <- length(unique(fungi$taxon_genus_name))         # Number of genera
families <- length(unique(fungi$taxon_family_name))      # Number of families
classes <- length(unique(fungi$taxon_class_name))        # Number of classes
orders <- length(unique(fungi$taxon_order_name))         # Number of orders
sybphylums <- length(unique(fungi$taxon_subphylum_name)) # Number of subphylum
phylums <- length(unique(fungi$taxon_phylum_name))       # Number of phylum

basic_stat <- data.frame(Metric = c(
  "Observations",
  "Observers",
  # "Obscured coordinates",
  # "Scientific names",
  "Species",
  "Genera",
  "Families",
  "Orders",
  "Classes",
  "Subphyla",
  "Phyla"
),
Number = c(
  num_obs,
  num_observers,
  # obscured_coords,
  # scientificNames,
  species,
  genera,
  families,
  orders,
  classes,
  sybphylums,
  phylums
))

kable(basic_stat, caption = "Basic statistics, iNaturalist data")

# Spatial distribution ####
# Density scatterplot
fungi$dens <- col2rgb(densCols(fungi[['longitude']], fungi[['latitude']]))[1,] + 1L
# Map densities to colors
fungi$colors = colors[fungi$dens]

ggplot() +
  geom_sf(data = ukraine_poly, fill = "transparent", colour = "gray") +
  geom_point(data = fungi, aes(x = longitude, y = latitude),
             size = 0.2,
             colour = fungi$colors,
             alpha = 0.5) +
  # scale_color_manual(values = rev(brewer.pal(6, "Paired"))) +
  labs(
    fill = "",
    title = "a) Observation density"
  ) +
  theme_void() -> all_inat_obs_map

# all_inat_obs_map

# https://www.statology.org/order-bars-ggplot2-bar-chart/
fungi %>% 
  ggplot(
    aes(x = reorder(place_admin1_name, place_admin1_name, function(x) length(x)))
  ) +
  geom_bar() +
  coord_flip() +
  scale_y_reverse(position = "right", breaks = c(0, 10000)) +
  scale_x_discrete(position = "top") +
  labs(x = "", y = "",
       title = "b) Number of observations per region") +
  mytheme +
  theme(panel.border = element_blank()) -> num_obs_barplot

# Combine both plots
all_inat_obs_map + num_obs_barplot +
  # the panel area of the first column is three time bigger that of the second column
  plot_layout(widths = c(4, 1)) -> distribution_plot

png("./figures/distribution_plot.png", width = 20, height = 12, units = "cm", res = 300)
distribution_plot
dev.off()

distribution_plot

# How many observations observers make?
fungi %>% 
  st_drop_geometry() %>% 
  count(user_id) %>%
  ggplot(aes(x = n)) +
  geom_histogram(fill="red", colour="red", alpha = 0.3) +
  scale_y_log10() +
  labs(
    x = "Number of observation",
    y = "Number of users, log-transformed",
    title = "a)"
  ) +
  mytheme -> observers_p

# Identification agreements
fungi %>% 
  mutate(num_identification_disagreements = as.factor(num_identification_disagreements)) %>% 
  ggplot(aes(x = num_identification_agreements,
             fill = num_identification_disagreements,
             colour = num_identification_disagreements)) +
  geom_bar(alpha = 0.8) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7)) +
  labs(
    x = "Number of identification agreements",
    y = "Number of observations",
    fill = "Num. of identification disagreements",
    title = "b)"
  ) +
  guides(colour = FALSE) +
  mytheme +
  theme(legend.title = element_text(size = 10),
        legend.position = "top") -> id_agreements_p

png("./figures/observers_and_id_agreements.png", width = 22, height = 12, units = "cm", res = 300)
observers_p + id_agreements_p
dev.off()

observers_p + id_agreements_p

# Identification precision ####
# How many observations identified to that level?
num_species <- fungi %>% filter(!taxon_species_name == "") %>% 
  nrow() # ID to species level

num_genus <- fungi %>% filter(taxon_species_name == "",
                              !taxon_genus_name == "") %>% 
  nrow() # ID only to genus level

num_family <- fungi %>% filter(taxon_genus_name == "",
                               !taxon_family_name == "") %>% 
  nrow() # ID only to family level

num_order <- fungi %>% filter(taxon_family_name == "",
                              !taxon_order_name == "") %>% 
  nrow() # ID only to order level

num_class <- fungi %>% filter(taxon_order_name == "",
                              !taxon_class_name == "") %>% 
  nrow() # ID only to class level

num_subphylum <- fungi %>% filter(taxon_class_name == "",
                                  !taxon_subphylum_name == "") %>% 
  nrow() # ID only to subphylum level

num_phylum <- fungi %>% filter(taxon_subphylum_name == "",
                               !taxon_phylum_name == "") %>% 
  nrow() # ID only to phylum level

num_kingdom <- fungi %>% filter(taxon_phylum_name == "") %>% 
  nrow() # ID only as 'Fungi'


# Merge numbers together to make labels for the barplot
labels <- tibble(level = c(
  "species",
  "genus",
  "family",
  "order",
  "class",
  "subphylum",
  "phylum",
  "kingdom"
),
values = c(
  num_species,
  num_genus,
  num_family,
  num_order,
  num_class,
  num_subphylum,
  num_phylum,
  num_kingdom
)) %>% 
  mutate(level = factor(level, levels = c(
    "species",
    "genus",
    "family",
    "order",
    "class",
    "subphylum",
    "phylum",
    "kingdom"
  ), ordered = TRUE)
  )


# Barplot
to_species <- fungi %>% filter(!taxon_species_name == "") %>% 
  mutate(id_level = "species") # ID to species level

to_genus <- fungi %>% filter(taxon_species_name == "",
                             !taxon_genus_name == "") %>% 
  mutate(id_level = "genus") # ID only to genus level

to_family <- fungi %>% filter(taxon_genus_name == "",
                              !taxon_family_name == "") %>% 
  mutate(id_level = "family") # ID only to family level

to_order <- fungi %>% filter(taxon_family_name == "",
                             !taxon_order_name == "") %>% 
  mutate(id_level = "order") # ID only to order level

to_class <- fungi %>% filter(taxon_order_name == "",
                             !taxon_class_name == "") %>% 
  mutate(id_level = "class") # ID only to class level

to_subphylum <- fungi %>% filter(taxon_class_name == "",
                                 !taxon_subphylum_name == "") %>% 
  mutate(id_level = "subphylum") # ID only to subphylum level

to_phylum <- fungi %>% filter(taxon_subphylum_name == "",
                              !taxon_phylum_name == "") %>% 
  mutate(id_level = "phylum") # ID only to phylum level

to_kingdom <- fungi %>% filter(taxon_phylum_name == "") %>% 
  mutate(id_level = "kingdom") # ID only as 'Fungi'

bind_rows(to_species,
          to_genus,
          to_family,
          to_order,
          to_class,
          to_subphylum,
          to_phylum,
          to_kingdom) %>% 
  select(quality_grade, id_level) %>% 
  ggplot(aes(x = id_level, 
             fill = quality_grade,
             colour = quality_grade
  )) +
  geom_bar(alpha = 0.8) +
  scale_x_discrete(limits = c("species",
                              "genus",
                              "family",
                              "order",
                              "class",
                              "subphylum",
                              "phylum",
                              "kingdom")) +
  annotate("text",
           x = labels$level,
           y = labels$values + 1300,
           label = labels$values) +
  scale_fill_manual(values = c("#EC754A",
                               "#EACF65", "#3C8D53")) +
  scale_colour_manual(values = c("#EC754A",
                                 "#EACF65", "#3C8D53")) +
  labs(x = "Most precise level of identification",
       y = "Number of observations"
  ) +
  guides(fill = FALSE,
         colour = FALSE) +
  mytheme -> p1

# How many observations are per grade category?
research_grade <- nrow(fungi[fungi$quality_grade == "research",])
casual <- nrow(fungi[fungi$quality_grade == "casual",])
needs_id <- nrow(fungi[fungi$quality_grade == "needs_id",])

# Pie-chart
p2 <- tibble(grade = c("Research Grade",
                       "Needs ID",
                       "Casual"),
             value = c(
               research_grade,
               needs_id,
               casual
             )) %>% 
  ggplot(aes(x = "", y = value, fill = grade)) +
  geom_col(color = "black", alpha = 0.8) +
  geom_label(aes(label = value),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c("#EC754A",
                               "#EACF65", "#3C8D53")) +
  guides(fill = guide_legend(title = "")) +
  theme_void() +
  theme(legend.position = "top")

# p2
# Combine both sub-plots
p1 + annotation_custom(ggplotGrob(p2), 
                       xmin = "order", xmax = "phylum",
                       ymin = max(length(unique(fungi$taxon_species_name))/2),
                       ymax = Inf) -> id_precision_p


# Save combined plot
png("./figures/identification_depth.png", width = 16, height = 12, units = "cm", res = 300)
id_precision_p
dev.off()

id_precision_p

# Taxonomic structure ####
fungi %>% 
  filter(taxon_species_name != "") %>% 
  select(taxon_species_name,
         taxon_genus_name,
         taxon_family_name,
         taxon_order_name,
         taxon_class_name,
         taxon_phylum_name) %>% 
  st_drop_geometry() %>% 
  unique() -> inat_checklist

table(inat_checklist$taxon_phylum_name, inat_checklist$taxon_class_name) %>% 
  as.data.frame() %>% 
  filter(Var1 != "") %>% 
  filter(Freq != 0) %>% 
  # as_tibble() %>% 
  mutate(Var2 = recode(Var2,
                       " " = "Incertae sedis")
  ) %>% 
  rename(Phylum = Var1,
         Class = Var2) -> taxa


png("./figures/treemap_inat.png", width = 16, height = 12, units = "cm", res = 300)
treemap(taxa,
        index = c("Phylum", "Class"),
        vSize = "Freq",
        vColor = "Freq",
        title = "",
        title.legend="number of species and infraspecific taxa",
        format.legend = list(scientific = FALSE, big.mark = " ")
)
dev.off()

treemap(taxa,
        index = c("Phylum", "Class"),
        vSize = "Freq",
        vColor = "Freq",
        title = "",
        title.legend="number of species and infraspecific taxa",
        format.legend = list(scientific = FALSE, big.mark = " ")
)

# Teporal structure ####
# Accumulation by date
fungi %>%
  filter(!is.na(observed_on)) %>%
  count(observed_on) %>%
  mutate(cum_sum = cumsum(n)) -> stacked.df

stacked.df %>%
  ggplot(aes(x = observed_on, y = cum_sum)) +
  geom_area(alpha = .8) +
  scale_x_date(limits = c(as.Date("2000-01-01"), as.Date("2023-12-31"))) +
  annotate("rect", xmin = as.Date("2022-02-24"), xmax = as.Date("2023-12-31"),
           ymin = -1, ymax = Inf, alpha = .2, fill = "red") +
  annotate("rect", xmin = as.Date("2000-06-01"), xmax = as.Date("2001-06-01"),
           ymin = 48000, ymax = 52000, alpha = .2, fill = "red") +
  annotate("text", x = as.Date("2004-09-01"), y = 50000,
           label = "Full-scale invasion"
  ) +
  annotate("rect", xmin = as.Date("2019-03-01"), xmax = as.Date("2022-02-24"),
           ymin = -1, ymax = Inf, alpha = .2, fill = "blue") +
  annotate("rect", xmin = as.Date("2000-06-01"), xmax = as.Date("2001-06-01"),
           ymin = 42000, ymax = 46000, alpha = .2, fill = "blue") +
  annotate("text", x = as.Date("2003-05-01"), y = 44000,
           label = "COVID-19"
  ) +
  labs(x = "Time",
       y = "Cumulative number of observations"
  ) +
  mytheme -> accumulation_plot


png("./figures/accumulation_plot.png", width = 16, height = 12, units = "cm", res = 300)
accumulation_plot
dev.off()

accumulation_plot

# Choroplets by years
# Oblast level
fungi %>% 
  mutate(year = year(observed_on)) %>% 
  filter(!is.na(year)) %>% 
  filter(year >= "2020") %>% 
  count(place_admin1_name, year) %>% 
  rename(shapeName = place_admin1_name) %>% 
  st_drop_geometry() %>% 
  left_join(ukraine_poly, by = "shapeName") %>% 
  st_as_sf() %>% 
  ggplot() +
  geom_sf(aes(fill = n)) +
  scale_fill_viridis_c() +
  labs(
    fill = "Number of\nobservations"
  ) +
  facet_wrap(vars(year)) +
  theme_void() -> regions_choroplet_map

png("./figures/regions_choroplet_map.png", width = 16, height = 12, units = "cm", res = 300)
regions_choroplet_map
dev.off()

regions_choroplet_map

fungi %>% 
  mutate(month = month(observed_on),
         year = year(observed_on)) %>% 
  filter(!is.na(year)) %>% 
  filter(year >= "2020") %>% 
  ggplot(aes(x = month)) +
  geom_bar() +
  coord_polar() +
  scale_x_continuous(
    breaks = c(1:12),
    labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Now", "Dec")
  ) +
  labs(x = "", y = "") +
  facet_wrap(vars(year)) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) -> season_dynamics_plot

png("./figures/season_dynamics_plot.png", width = 16, height = 12, units = "cm", res = 300)
season_dynamics_plot
dev.off()

season_dynamics_plot

# Rare and alien species ####
fungi %>% 
  filter(taxon_species_name %in% c("Hericium coralloides",
                                   "Polyporus umbellatus",
                                   "Morchella steppicola",
                                   "Clathrus archeri"
  )
  ) %>% 
  mutate(taxon_species_name = factor(taxon_species_name, ordered = TRUE,
                                     levels = c("Hericium coralloides",
                                                "Polyporus umbellatus",
                                                "Morchella steppicola",
                                                "Clathrus archeri"
                                     )
  )
  ) %>% 
  mutate(taxon_species_name = recode(taxon_species_name,
                                     "Hericium coralloides" = "a) Hericium coralloides", # "old one" = "new one"
                                     "Polyporus umbellatus" = "b) Polyporus umbellatus",
                                     "Morchella steppicola" = "c) Morchella steppicola",
                                     "Clathrus archeri" = "d) Clathrus archeri"
  )
  ) %>% 
  ggplot() +
  geom_sf(data = ukraine_poly, fill = "transparent", colour = "gray") +
  geom_sf(alpha = 0.5) +
  facet_wrap(vars(taxon_species_name)) +
  theme_void() -> redlisted_plot

png("./figures/redlisted_plot.png", width = 16, height = 12, units = "cm", res = 300)
redlisted_plot
dev.off()

redlisted_plot

# End of the script ####
