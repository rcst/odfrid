S <- 25
S * (S-1) / 2
inp <- generate_fake_pc(S = S, N = 100, boarders = 20)
inp$dest
out <- model_sample(ax = inp$x,
                    dep_time = inp$t, 
                    sample = 100, 
                    warmup = 100, 
                    D = 4, 
                    print_n = 1) 

save(inp, out, file = "model_check_S25_data.rds")
load(file = "model_check_S25_data.rds")

# calculate median of matrices
library(data.table)
library(ggplot2)
library(plyr)
dt <- data.table(apply(X = out$y, MARGIN = c(1, 2), FUN = median))
dt[, od_index := .I]
mdt <- melt(dt, id.vars = "od_index", value.name = "flow")
mdt[, type := "estimate"]

dt2 <- data.table(inp$y)
dt2[, od_index := .I]
mdt2 <- melt(dt2, id.vars = "od_index", value.name = "flow")
mdt2[, type := "original"]

dt <- rbind(mdt, mdt2)
dt[, n := as.integer(gsub("^V(.*)$", "\\1", dt$variable))]
ddt <- dcast(data = dt, od_index + variable ~ type, value.var = "flow")

ggplot(data = ddt, mapping = aes(x = original, y = estimate)) +
  geom_point(alpha = I(1/2)) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal() +
  xlab("original [passenger flow]") +
  ylab("estimate [passenger flow]") +
  ggtitle("OD estimates simulated data set of a 6 stops route and 100 trips")

ggsave(filename = "model-estimate-comparision-scatter-S6.png")


# dt <- data.table(apply(X = out$Y, MARGIN = c(1, 2), FUN = median))
# dt[, od_index := .I]
# mdt <- melt(dt, id.vars = "od_index", value.name = "flow")
# mdt[, n := as.integer(gsub("^V(.*)$", "\\1", mdt$variable))]
# mdt[, variable := NULL]

ggplot(data = dt, mapping = aes(x = n, y = od_index, fill = flow)) +
  geom_raster() +
  # coord_equal() +
  scale_fill_continuous(type = "viridis", name = "OD flow") +
  facet_grid(type~.) +
  xlab("bus trip ID") +
  ylab("OD index") +
  ggtitle("OD estimates vs. original simulated data set - 6 stops route and 100 trips")

ggsave(filename = "model-estimate-comparison-S6.png", width = 20, height = 10)

dt <- data.table(likelihood = out$lq)
dt[, index := .I]
setnames(dt, colnames(dt)[1], "likelihood")
ggplot(data = dt, mapping = aes(y = likelihood, x = index)) +
  geom_line() +
  ylab("log likelihood")

dt_psi <- adply(.data = out$psi, 
                .margins = c(2, 3), 
                .id = c("psi", "smpl"),
                .fun = function(x) x, 
                .progress = "text") |> setDT() |>
melt(id.vars = c("psi", "smpl"), variable.name = "tt")
dt_psi[, t := as.integer(tt)]
dt_psi[, tt := NULL]
dt_psi[, s := as.integer(smpl)]
dt_psi[, smpl := NULL]
dt_psi[, group := .GRP, by = .(psi, s)]

ggplot(data = dt_psi[s %in% sample(dt_psi[, max(s)], 100)], mapping = aes(y = value, x = t, color = psi, group = group)) +
  geom_line(alpha = I(1/2))


all(out$lambda[300,,] == 1)
inp$dest

dt_lbd <- adply(.data = out$lambda[,,1],
                .margins = 2,
                .fun = function(x) x[1:(S-1)],
                .progress = "text",
                .id = c("n")) |> setDT() |>
melt(id.vars = c("n"), variable.name = "index")

dt_lbd[, nbus := as.integer(n)]
dt_lbd[, idx := factor(as.integer(index))]
dt_lbd[, c("n", "index") := NULL] 
dt_lbd[, cvalue := cumsum(value), keyby = .(nbus)]
dt_lbd[, mvalue := shift(x = cvalue, fill = 0, type = "lag"), keyby = .(nbus)]

require(scales)
n <- length(levels(factor(dt_lbd$idx))) # number of colors
cols <- hue_pal(h = c(0, 360) + 15, 
                c = 100, l = 65, 
                h.start = 0, direction = 1)(n)[order(sample(1:n, n))] # color palette in random orde

ggplot(dt_lbd, 
       mapping = aes(x = nbus, 
                     ymin = mvalue, 
                     ymax = cvalue, 
                     fill = as.factor(idx), 
                     group = idx)) + 
  geom_ribbon(color = "black", linewidth = 0.2) +
  # scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = cols) +
  theme(legend.position = "none")
  # scale_x_continuous(limits = c(0, 10))

