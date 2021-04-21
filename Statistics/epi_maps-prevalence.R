#!/usr/bin/Rscript

# Mercè Garí
# 210317

library(ggmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rgdal)
library(spdep)
library(stringr)
library(sf)
library(sp)
library(colorspace)
library(cowplot)
library(geojson)
library(grid)
library(spatial)
library(ggspatial)
library(xlsx)

#######################################################
### Follows https://juanitorduz.github.io/germany <- plots/
### Data from http://insideairbnb.com/get-the-data.html
### file http://data.insideairbnb.com/germany/bv/munich/2020-06-20/visualisations/neighbourhoods.geojson

# Input and output locations used in this script
here_r = function (...) here:here("Statistics", ...)
here_output = function (...) here::here("Output", ...)
here_maps_data = function (...) here::here("maps_data", ...)

# Polygon of districts
m.sh <- readOGR(here_maps_data("neighbourhoods.geojson"))
class(m.sh)
head(m.sh)

districts <- m.sh
d.sf <- st_as_sf(districts)
d.sf <- d.sf %>%
  rename(District = neighbourhood)

# Polygon of constituencies
S.original <- st_read(here_maps_data("KW_2020_Stimmbezirke_V2.shp"))
constituencies <- S.original

# Depending on the data, we would need 4 or 5 digits constituency ID (cid)
# 4 digits refers to: last two digits are the constituencies and first one
#  or first two are the districts
# 5 digits refers to: 2+3, 2 for the district, and 3 with the leading 0 for
#  the constituencies
cid <- S.original %>%
  as_tibble() %>%
  select(KW_SB_2020, KW_SB_20_2) %>%
  unique() %>%
  rename(cid.4 = KW_SB_2020,  # 4 digits
         cid.5 = KW_SB_20_2)  # 5 digits

# Polygon of constituencies clean
S <- S.original %>%
  mutate(cid.4 = KW_SB_2020,
         cid.5 = KW_SB_20_2,
         District_ID = as.numeric(str_sub(cid.5, 1, 2)))

# Match district with district-id
rd <- tibble(District = attributes(m.sh)$data$neighbourhood) %>%
  mutate(District_ID = (1:(n()))) 
rd
# Polygons of districts
rd1 <- tibble(District = attributes(m.sh)$data$neighbourhood) %>%
  mutate(region = as.character(0:(n()-1)))
rd1
# Joining polygon constituencies with districts
S <- S %>%
  left_join(rd)  # joining by District_ID

# Create and write table of 4/5 CIDs with district IDs and district names
dc <- S %>%
  as_tibble() %>%
  select(cid.4, cid.5, District_ID, District) %>%
  unique()
#write.csv(dc, here_maps_out("District-Constituency-IDs.csv"))

# Create table with longitude, latitude, group, region
map.m <- map_data(m.sh)

# Group by district, with longitude and latitude as the means
map.centroids <- map.m %>%
  left_join(rd1) %>%
  group_by(District) %>%
  summarize(long = mean(long), lat = mean(lat))

# Load population density
pop.d <- read.xlsx(here_maps_data("Muc_popn_density_Wiki.xlsx", 1))
head(pop.d)

pop.d <- pop.d %>%
  rename(District_ID = dist) %>%
  select(District_ID, density) %>%
  rename(`Pop. Density` = density)

# Load prevalence for districts
prev.d <- read.xlsx(here_maps_data("District_Seroprev.xlsx", 1))
head(prev.d)

prev.d <- prev.d %>%
  filter(Outcome_dist %in% "Binomial") %>%
  filter(smoothing %in% "Spatial AR - Global") %>%
  filter(prev_type %in% "Not adjusted") %>%
  filter(outcome_type %in% "Weighted") %>%
  rename(District = Borough) %>%
  mutate(neighbourhood = District) %>%
  mutate(`Prevalence (%)` = SIR.50*100,
         `Prevalence (%)\nCI low` = SIR.025*100,
         `Prevalence (%)\nCI high` = SIR.975*100)#

# The number of individuals should be the ones in sheet 2 (for crude analysis)
num.ind <- read.xlsx(here_maps_data("District_Seroprev.xlsx", 2))
head(num.ind)  
num.ind <- num.ind %>%
  rename(`N. Participants` = Individuals.Tested_per_dist) 

prev.d <- left_join(prev.d, num.ind)

prev.d.S <- left_join(S, prev.d) 
prev.d.S <- left_join(prev.d.S, pop.d)

prev.d.C <- left_join(d.sf, prev.d)
prev.d.C <- left_join(prev.d.C, dc) %>%
  left_join(pop.d)

library(viridis)
map.prev.mean <- ggplot(data=prev.d.C, size=0.1, 
                        aes(fill=`Prevalence (%)`)) + 
  geom_sf(color="white") +
  coord_sf() + theme_bw() +
  scale_fill_continuous_sequential(palette = "Viridis",
                                   breaks=c(3, 3.5, 4)) +
  theme(panel.background = element_blank(),
        axis.text=element_blank(), panel.grid.major = element_blank()) +
  annotation_scale()
map.prev.mean  

map.prev.low <- ggplot(data=prev.d.C, size=0.1, 
                       aes(fill=`Prevalence (%)\nCI low`)) + 
  geom_sf(color="white") +
  coord_sf() + theme_bw() +
  scale_fill_continuous_sequential(palette = "Viridis",
                                   breaks=c(1.75, 2.25, 2.75)) +
  theme(panel.background = element_blank(),
        axis.text=element_blank(), panel.grid.major = element_blank()) +
  annotation_scale()

map.prev.low

map.prev.high <- ggplot(data=prev.d.C, size=0.1, 
                        aes(fill=`Prevalence (%)\nCI high`)) + 
  geom_sf(color="white") +
  coord_sf() + theme_bw() +
  scale_fill_continuous_sequential(palette = "Viridis",
                                   breaks=c(4.5, 5.75, 7)) +
  theme(panel.background = element_blank(),
        axis.text=element_blank(), panel.grid.major = element_blank()) +
  annotation_scale()
map.prev.high

#districts.nind <- merge(districts, prev.d)
map.centroids <- left_join(map.centroids, prev.d)

map.dens.ind <- ggplot() +
  theme_bw() +
  geom_sf(data = prev.d.C, 
          aes(fill=`Pop. Density`),
          colour="white") +
  geom_point(data=map.centroids, aes(x=long, y=lat,
                                     size=`N. Participants`),
             color="orange") +
  scale_size_continuous(breaks=c(50, 150, 300),
                        labels=c("<100", "100-200", "250-500"),
                        range=c(1,5)) +
  xlab("") + ylab("") +
  coord_sf() + theme_bw() +
  scale_fill_continuous_sequential(palette = "Viridis",
                                   name=bquote(atop("Pop. Density","("~persons/km^2~")"))) +
                                   #name="Pop. density\n(persons/km2)") +
                                     #expression(paste("Pop. Density", person/km^2))) +
  theme(panel.background = element_blank(),
        axis.text=element_blank(), panel.grid.major = element_blank()) +
  annotation_scale()
map.dens.ind

plot_grid(
  plot_grid(map.dens.ind, map.prev.mean, nrow=1, labels=c("A", "B")),
  plot_grid(map.prev.low, map.prev.high,
          labels=c("C", "D"), nrow=1),  nrow=2)
ggsave(here_output("Epi2-PrevMap-Bin-Global-NA-Weighted.pdf"), height=10, width=12)
ggsave(here_output("Epi2-PrevMap-Bin-Global-NA-Weighted.png"), height=10, width=12)

  