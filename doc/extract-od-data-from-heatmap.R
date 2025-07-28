library(jpeg)
library(data.table)
library(ggplot2)

arr <- readJPEG("out-004.jpg")
dt <- rgb(arr[,,1], arr[,,2], arr[,,3]) |> matrix(nrow = nrow(arr), ncol = ncol(arr)) |> as.data.table()

dt[, row := rev(1:nrow(arr))]
dt <- melt(dt, id.vars = "row", variable.name = "column")
dt[, column := as.numeric(gsub(pattern = "^V([[:digit:]]+)$", replacement = "\\1", x = dt$column))]

# grid 22 * 21 / 2 = 231 OD pairs
bus_x <- seq(from  = 242.2, to = 3195, length.out = 520)
length(bus_x)

# stop_y <- seq(from = 106, to = 730, by = 2.295)
stop_y <- seq(from = 194.5, to = 828.5, length.out = 231)
length(stop_y)

dots <- CJ(stop_y, bus_x)
dots[, bus_id := .GRP, by = bus_x]
dots[, stop_id := .GRP, by = stop_y]
dots[, stop_id := max(stop_id) + 1L - stop_id]
setkey(dots, bus_x)

# OD
ggplot() +
  geom_raster(data = dt[row < 1000 & column < 3250], 
              mapping = aes(x = column, y = row, fill = value)) +
  scale_fill_identity() +
  coord_equal() +
  geom_text(data = dots, mapping = aes(x = bus_x, y = stop_y, label = bus_id), color = "black", size = 1, angle = 0)

ggsave("test.png", width = 16*3, height = 9*3, dpi = 300)

