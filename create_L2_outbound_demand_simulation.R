# install.packages(c("osmdata", "sf", "dplyr", "ggplot2"))

library(osmdata)
library(sf)
library(glue)
library(dplyr)
library(plyr)
library(ggplot2)
library(data.table)

set.seed(42)

get_route_stops <- function(route_id) {
        db <- DBI::dbConnect(drv = odbc::odbc(), "mot")
        x <- glue_sql("
                SELECT seq_id,
                ST_X(ST_Transform(geom, 4326)) AS lon,
                ST_Y(ST_Transform(geom, 4326)) AS lat
                FROM route_stops
                JOIN stops USING (stop_id)
                WHERE route_id = {route_id}
                AND inbound = false
                ORDER BY seq_id
                ", .con = db, route_id = route_id, inbound = !!inbound) |>
                DBI::dbGetQuery(conn = db, statement = _)
        DBI::dbDisconnect(db)
        return(x)
}

stops <- get_route_stops(22) |> setDT()
stops_sf <- st_as_sf(stops, coords = c("lon", "lat"), crs = 4326)

# 2. Fetch OSM POIs
bbox <- st_bbox(stops_sf)
q <- opq(bbox = bbox) |>
        add_osm_features(list(amenity = NULL, shop = NULL, office = NULL, tourism = NULL))
pois <- osmdata_sf(q)$osm_points %>% st_as_sf()

# 3. Count POIs within 300m buffer
buffers <- st_buffer(st_transform(stops_sf, 3857), 300)
pois_t <- st_transform(pois, 3857)
poi_counts <- colSums(st_within(pois_t, buffers, sparse = FALSE))

stops$poi_counts <- poi_counts

# Fetch population in catchement area around stops
pop_raster <- system.file("data/tif",
        "GHS_POP_E2030_GLOBE_R2023A_4326_3ss_V1_0_R5_C20.tif",
        package = "odfrid",
        mustWork = TRUE
) |>
        terra::rast()

# r <- terra::project(pop_raster, "epsg:4326")
sext <- terra::ext(stops_sf)
rc <- terra::crop(pop_raster, sext)
pop_values <- terra::extract(pop_raster, buffers, fun = sum, na.rm = TRUE)
stops$pop_catch <- pop_values[, 2]

stops[, production := pop_catch / sum(pop_catch)]
stops[, attraction := poi_counts / sum(poi_counts)]

# stops$production <- round(runif(nrow(stops), 10, 50)) # synthetic population

# 4. Distance matrix for gravity model
n <- nrow(stops)
dist_matrix <- st_distance(st_transform(stops_sf, 3857)) |> units::drop_units()
dist_decay <- exp(-0.001 * dist_matrix)

# ---- Per-Bus-Trip Simulation ----

# Time periods with start times
# time_periods <- data.frame(
#   name = c("AM_peak", "Midday", "PM_peak", "Evening"),
#   start_time = as.POSIXct(c("2025-08-28 07:00:00",
#                             "2025-08-28 10:00:00",
#                             "2025-08-28 16:00:00",
#                             "2025-08-28 20:00:00"))
# )
#
# scaling <- list(
#   AM_peak = list(prod = 1.5, attr = 2.0),
#   Midday  = list(prod = 0.8, attr = 1.0),
#   PM_peak = list(prod = 2.0, attr = 1.5),
#   Evening = list(prod = 0.5, attr = 0.8)
# )

# passengers_per_period <- c(300, 150, 350, 100)
# buses_per_period <- c(6, 3, 7, 2)  # evenly spaced departures

# ---- Extract Schedule Times --- 
schedule_arrival_times <- function() {
	db <- DBI::dbConnect(odbc::odbc(), "mot")
	x <- glue_sql("SELECT arrival_time 
		      FROM schedule 
		      WHERE route_id  = 2 
		      AND season = 'before-summer' 
		      AND day_type = 'normal' 
		      AND direction = 'outbound' 
		      AND year = 2025 
		      ORDER BY arrival_time;", .con = db) |>
		      DBI::dbGetQuery(db, statement = _) |>
		      setDT()
	DBI::dbDisconnect(db)
	# x[, arrival_time := as.POSIXct(x = arrival_time, "%H:%M:%S")]
	      x
}

xt <- schedule_arrival_times()
xt[, arrival_time := as.ITime(arrival_time)]
xt <- CJ(dd = seq(from = as.POSIXct("2025-09-01 00:00:00"),
		  to = as.POSIXct("2025-09-05 00:00:00"), 
		  by = "1 day"), 
	 at = xt$arrival_time)
xt[, arrival_t := as.integer(at)]
xt[, arrival_t := arrival_t - min(arrival_t)]
xt[, arrival_time := dd + at]

scalings <- data.table(name = c("am", "ip", "pm", "ap"),
		       period_begin = c(as.ITime("05:00:00"), 
					as.ITime("10:00:00"),
					as.ITime("14:00:00"),
					as.ITime("18:00:00")),
		       period_end = c(as.ITime("10:00:00"), 
				      as.ITime("14:00:00"),
				      as.ITime("18:00:00"),
				      as.ITime("23:59:59")),
		       prod = c(1.5, 0.8, 2.0, 0.5), 
		       attr = c(2.0, 1.0, 1.5, 0.8))

xt <- scalings[xt, on = c("period_begin<=at", "period_end>at")]
xt[, trip_id := .I]

# ---- Compute OD Matrix ----
S <- nrow(stops)
D <- as.integer(S * (S-1) / 2)
y <- integer(0) # OD vector
b <- integer(0L) # boarders
a <- integer(0L) # alighters
# ---- Simulate APC ----
for(trip in 1:nrow(xt)) {
	od_matrix <- outer(stops$production * xt[trip, prod], stops$attraction * xt[trip, attr]) * dist_decay
	od_matrix[lower.tri(od_matrix, diag = TRUE)] <- NA
	od_matrix <- od_matrix / rowSums(od_matrix, na.rm = TRUE)
	nop <- 500L # no. passengers of whole trip
	# OD of passenger numbers
	pod <- matrix(data = rep(NA_integer_, nrow(stops)*nrow(stops)), 
		      nrow = nrow(stops), 
		      ncol = nrow(stops))
	for (i in 1:(nrow(stops) - 1L)) {
		u <- c(rep(NA, i), 
		       rmultinom(n = 1, 
				 size = round(stops$production[i] * xt[trip, prod]  * nop), 
				 prob = na.omit(od_matrix[i, ])))
		pod[i,] <- u 
	}
	# convert to compressed format
	y <- cbind(y, t(pod) |> as.vector() |> na.omit())
	b <- cbind(b, rowSums(pod, na.rm = TRUE))
	a <- cbind(a, colSums(pod, na.rm = TRUE))
}

x <- rbind(b, a)

# ---- Visualize Load ----
ld <- matrix(data = b[1, ], nrow = 1L, ncol = ncol(b)) 
for(i in 2:nrow(b))
	ld <- rbind(ld, ld[i-1L,] + b[i,] - a[i,])

dt_ld <- data.table(load = as.vector(ld), 
		    stop_id = rep(1:nrow(b), ncol(b)), 
		    trip_id = rep(1:ncol(b), each = nrow(b)))
dt_ld <- xt[dt_ld, on = .(trip_id)]

ggplot(data = dt_ld, mapping = aes(x = stop_id, y = load, color = trip_id, group = trip_id)) +
	geom_step() +
	facet_wrap(~name)
