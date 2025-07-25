S <- 25
S * (S-1) / 2
inp <- generate_fake_pc(S = S, N = 100, boarders = 20)
inp$dest
out <- model_sample(ax = inp$x,
                    dep_time = inp$t, 
                    sample = 2000, 
                    warmup = 2000, 
                    D = 4, 
                    print_n = 1) 

plot(out$lq, type = 'l')


all(inp$x == routing_matrix(S) %*% out$y[,,4])
all(out$y[,,10] == 0)

names(out)
save(out, file = "model_check_S25_data.rds")

# calculate median of matrices
library(data.table)
library(ggplot2)
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

inp$dest

S <- 25
M <- S * (S - 1) / 2
K <- matrix(nrow = S, ncol = S)
K[lower.tri(K)] <- 1:M
K <- t(K)
K[,inp$dest[2]][1:(inp$dest[2]-1)]
out$lambda[K[inp$dest[1],][(inp$dest[1]+1):S],,]

i <- function(j) K[,j][1:(j-1)]
j <- function(i) K[i,][(i+1):S]

test()
