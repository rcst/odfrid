library(jpeg)
library(data.table)
library(ggplot2)

arr <- readJPEG("doc/out-003.jpg")
dt <- rgb(arr[,,1], arr[,,2], arr[,,3]) |> matrix(nrow = nrow(arr), ncol = ncol(arr)) |> as.data.table()

dt[, row := rev(1:nrow(arr))]
dt <- melt(dt, id.vars = "row", variable.name = "column")
dt[, column := as.numeric(gsub(pattern = "^V([[:digit:]]+)$", replacement = "\\1", x = dt$column))]

# grid
bus_x <- seq(from  = 207, to = 1675, by = 2.848)
length(bus_x)

stop_y <- seq(from = 200, to = 730, by = 25)
length(stop_y)

dots <- CJ(stop_y, bus_x)
dots[, bus_id := .GRP, by = bus_x]
dots[, stop_id := .GRP, by = stop_y]
dots[, stop_id := max(stop_id) + 1L - stop_id]
setkey(dots, bus_x)

# boardings
# ggplot() +
#   geom_raster(data = dt[row > 120 & row < 740 & column > 200 & column < 1680], 
#               mapping = aes(x = column, y = row, fill = value)) +
#   scale_fill_identity() +
#   coord_equal() +
#   geom_text(data = dots, mapping = aes(x = bus_x, y = stop_y, label = bus_id), color = "black", size = 1, angle = 90)
# 
# ggsave("test.png", width = 16*2, height = 9*2, dpi = 300)

boarding_data <- dt[dots, on = c(row = "stop_y", column = "bus_x"), roll = "nearest"]

rgb2grey <- Vectorize(FUN = function(hex) {
  x <- col2rgb(hex)
  return(0.1 * x[1] + 0.2 * x[2] + 0.7 * x[3])
}, vectorize.args = "hex")

boarding_data[, grey := rgb2grey(value)]

# darkest blue - 78
rgb2grey("#102e60")
rgb2grey("#ffffff")
boarding_data[grey < 70][, .SD[1], by = bus_id]
boarding_data[, day := "Fri"]
boarding_data[bus_id < 409, day := "Thu"]
boarding_data[bus_id < 308, day := "Wed"]
boarding_data[bus_id < 204, day := "Tue"]
boarding_data[bus_id < 101, day := "Mon"]
boarding_data[grey < 70, grey := NA]
boarding_data[is.na(grey), day := NA]
boarding_data[, pax := 20  / (78 - 255) * grey + (((20 * 78) / (255 - 78)) + 20)]

bus_x <- seq(from  = 1817, to = 3287, by = 2.848)
length(bus_x)

stop_y <- seq(from = 200, to = 730, by = 25)
length(stop_y)

dots <- CJ(stop_y, bus_x)
dots[, bus_id := .GRP, by = bus_x]
dots[, stop_id := .GRP, by = stop_y]
dots[, stop_id := max(stop_id) - stop_id + 1L]
setkey(dots, bus_x)

# alightings
# ggplot() +
#   geom_raster(data = dt[row > 120 & row < 740 & column > 1810 & column < 3287], 
#               mapping = aes(x = column, y = row, fill = value)) +
#   scale_fill_identity() +
#   coord_equal() +
#   geom_text(data = dots, mapping = aes(x = bus_x, y = stop_y, label = bus_id), color = "black", size = 1, angle = 90)
# 
# ggsave("test-alighting.png", width = 16*2, height = 9*2, dpi = 300)

alighting_data <- dt[dots, on = c(row = "stop_y", column = "bus_x"), roll = "nearest"]

rgb2grey <- Vectorize(FUN = function(hex) {
  x <- col2rgb(hex)
  return(0.1 * x[1] + 0.2 * x[2] + 0.7 * x[3])
}, vectorize.args = "hex")

alighting_data[, grey := rgb2grey(value)]

# darkest blue - 78
rgb2grey("#102e60")
rgb2grey("#ffffff")
alighting_data[grey < 70][, .SD[1], by = bus_id]
alighting_data[, day := "Fri"]
alighting_data[bus_id < 409, day := "Thu"]
alighting_data[bus_id < 308, day := "Wed"]
alighting_data[bus_id < 204, day := "Tue"]
alighting_data[bus_id < 101, day := "Mon"]
alighting_data[grey < 70, grey := NA]
alighting_data[is.na(grey), day := NA]
alighting_data[, pax := 20  / (78 - 255) * grey + (((20 * 78) / (255 - 78)) + 20)]

# ggplot(data = data, mapping = aes(x = bus_id, y = stop_id, fill = pax)) +
#   geom_raster()

boarding_data[, pax_type := "boarding"]
alighting_data[, pax_type := "alighting"]

data <- na.omit(rbind(boarding_data, alighting_data))
data[, pax_round := round(pax, 0L)]
data[, pax_floor := floor(pax)]

fwrite(x = data, file = "../data/csv/chen-sample-afc.csv")

# check load
ddt <- dcast(data = data, formula = bus_id + stop_id ~ pax_type, value.var = "pax_round")
setnames(ddt, c("boarding", "alighting"), c("pax_on", "pax_off"))
ddt[, pax_diff := pax_on - pax_off]
ddt[, load := cumsum(pax_diff), by = bus_id]

ggplot(data = ddt[bus_id < 17], mapping = aes(x = stop_id, y = load)) +
  geom_step() +
  facet_wrap(~bus_id)

ddt[, .N, by = bus_id][N != 22]

# ggplot(data = data, mapping = aes(x = bus_id, y = stop_id, fill = pax)) +
#   geom_tile() +
#   scale_fill_continuous(low = "white", high = "blue") +
#   facet_wrap(~pax_type)
