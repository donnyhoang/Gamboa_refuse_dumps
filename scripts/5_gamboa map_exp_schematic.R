library(ggmap)
library(dplyr)
library(ggplot2)
library(cowplot)
library(magick)
library(tiff)
library(jpeg)


#get map first

register_stadiamaps("13f7981e-5a38-4f91-bb4e-e23c4244e12f")

#read in site coordinates
data <-read.csv("site_coord.csv")

#R stripped the '-' sign on import, add it back so my coordinates are on the correct side of the world
data$Longitude <- data$Longitude * -1


#only first 3 sites are relevant
#data <- data[1:3,]

#set map limits
myLocation <- c(left=-79.76, bottom=9.105,right=-79.69,top=9.165)
#make map
myMap <- get_stadiamap(myLocation, maptype= "stamen_terrain", zoom = 14)
#ggmap(myMap)

#add coordinates

map <- ggmap(myMap) +
  geom_point(data=data, aes(x=Longitude,y=Latitude, color = Type), size=5) +
  scale_color_manual(values = c("black", "red")) +
  labs(y="Latitude", x="Longitude") +
  theme_classic(base_size = 20) +
  #theme(legend.position = "bottom") +
  theme(legend.position = "none") +
  guides(color=guide_legend(title = "Sampling Sites")) +
  ggtitle("Sampling Sites in Gamboa, Panama")

map

#ggsave("gamboa_map.tiff", device = "tiff", dpi = 700)
#ggsave("gamboa_map.png", device = "png", dpi = 700)
#ggsave("gamboa_map.pdf", device = "pdf", dpi = 700)


#add dump iamge
dump <- readJPEG("dump_example.jpg")
schematic <- readJPEG("schematic2.jpg")
#schematic <- readTIFF("schematic2.tif")


p1 <- ggdraw() +
  draw_image(dump, scale = 0.9) +
  theme(plot.margin = unit(c(0,0,0,0), "cm") )

p2 <- map

p3 <- ggdraw() +
  draw_image(schematic) +
  theme(plot.margin = unit(c(0,0,0,0), "cm") )


figtemp <- plot_grid(p1, p2, rel_heights = c(8,10), rel_widths = c(4, 6), labels = "AUTO", label_size = 25)
figtemp


fig <- plot_grid(figtemp, p3, ncol = 1, rel_heights = c(10,8), labels = c("","C"), label_size = 25, align = "l")
fig

ggsave("map_schematic.tiff", device = "tiff", dpi = 700)
ggsave("map_schematic.png", device = "png", dpi = 700)
ggsave("map_schematic.pdf", device = "pdf", dpi = 700)
