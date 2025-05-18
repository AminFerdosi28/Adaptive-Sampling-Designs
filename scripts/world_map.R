library(readxl)
library(dplyr)
library(ggplot2)
library(maps)
library(here)

data <- read_xlsx(here('data/data_plot.xlsx'))

######## Grid tile map
c1<- summary(data$New_Cases)[['3rd Qu.']]
covid_index<- rep(0, nrow(data))
covid_index[which(data$New_Cases >= c1)] = 1
covid_index<- as.factor(covid_index)
data<- cbind(data, covid_index)
c2<- 8
temp_index<- rep(0, nrow(data))
temp_index[which(data$Temperature > c2)] = 1
temp_index<- as.factor(temp_index)
data<- cbind(data, temp_index)

# create worldmap
worldmap <- ggplot(data)

# covid map
worldC19<- worldmap + 
  geom_rect(aes(xmin = mapx, ymin = mapy, 
                xmax = mapx + 0.96, ymax = mapy + 0.96,
                fill = covid_index)) +
  scale_fill_manual(values = c('0' = 'limegreen', '1' = 'indianred3'), label = c('Low', 'High'), name = expression(bold('Status'))) + 
  ggtitle(expression(bold("New cases of COVID-19 in Dec 2021"))) +
  geom_text(aes(x = mapx, y = mapy, 
                label = Country_Code),
            size = 3, 
            nudge_x = 0.5, nudge_y = -0.5,
            vjust = 0.5, hjust = 0.5) +
  scale_y_reverse() + 
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text = element_blank(), axis.title = element_blank(), 
        plot.title = element_text(hjust = 0.5))
worldC19


# temp map
worldTemp<- worldmap + 
  geom_rect(aes(xmin = mapx, ymin = mapy, 
                xmax = mapx + 0.96, ymax = mapy + 0.96,
                fill = temp_index)) +
  scale_fill_manual(values = c('0' = 'steelblue', '1' = 'indianred3'), label = c('Cold', 'Tropic'), name = expression(bold('Climate'))) + 
  ggtitle(expression(bold("Cold and Tropical Countries"))) +
  geom_text(aes(x = mapx, y = mapy, 
                label = Country_Code),
            size = 3, 
            nudge_x = 0.5, nudge_y = -0.5,
            vjust = 0.5, hjust = 0.5) +
  scale_y_reverse() + 
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text = element_blank(), axis.title = element_blank(), 
        plot.title = element_text(hjust = 0.5))
worldTemp



######## world map
world <- map_data("world")
data_world <- read_xlsx(here("data/data_worldmap.xlsx"))
c1<- summary(data_world$New_Cases)[['3rd Qu.']]
covid_index<- rep(0, nrow(data_world))
covid_index[which(data_world$New_Cases >= c1)] = 1
covid_index<- as.factor(covid_index)
data_world<- cbind(data_world, covid_index)
c2<- 8
temp_index<- rep(0, nrow(data_world))
temp_index[which(data_world$Temperature > c2)] = 1
temp_index<- as.factor(temp_index)
data_world<- cbind(data_world, temp_index)

worldSubset <- inner_join(world, data_world, by = "region")

dff <- worldSubset %>%
  group_by(Country_Code) %>%
  summarize(long = mean(long, na.rm = T), lat = mean(lat, na.rm = T), group = group)


plain <- theme(
  axis.text = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.title = element_blank(),
  panel.background = element_rect(fill = "white"),
  plot.title = element_text(hjust = 0.5)
)

worldC19_map <- ggplot(data = worldSubset, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  geom_polygon(aes(fill = log(New_Cases)), show.legend=F) +
  scale_fill_distiller(palette = 'Reds', direction = 1) + 
  ggtitle(expression(bold("New cases of COVID-19 in Dec 2021"))) +
  borders('world', colour = 'white') + plain
worldC19_map


worldTemp_map <- ggplot(data = worldSubset, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) +
  geom_polygon(aes(fill = Temperature), show.legend = F) +
  scale_fill_distiller(palette = 'Blues', direction = -1) + 
  ggtitle(expression(bold("Cold and Tropical Countries"))) +
  borders('world', colour = 'white') + plain
worldTemp_map



