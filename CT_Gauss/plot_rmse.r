library(ggplot2)
library(ggthemes)
library(scales)

# rmse_gauss <- as.matrix(read.csv("out/rmse_t.csv"))

steps <- 1000

dist <- c(rep("Gaussian Distributions", 5*steps))

comp <- c(rep("Position RMSE (m)",   5*steps - 1), 5)

Filter <- c(rep("HOUSE", steps),
                rep("UKF",   steps),
                rep("CUT-4", steps),
                rep("CUT-6", steps),
                rep("CUT-8", steps)) 

rmse_gauss   <- as.matrix(read.csv("R_testing/rmse_t.csv"))

t <- rep(rmse_gauss[,1], 5)

RMSE <- c(rmse_gauss[,2:6])

filter <- c("HOUSE", "UKF", "CUT-4", "CUT-6","CUT-8")



data <- data.frame(t, RMSE, Filter, comp)

data$Filter <- factor(data$Filter, levels=filter)

colors <- c("#cc3311", "#0077bb", "#ee7733", "#009988", "#ee3377")
dashes <- c("solid", "dotted", "dotdash", "longdash", "twodash")

ggplot(data, aes(x=t, y=RMSE)) + geom_line(aes(col=Filter, linetype=Filter)) +
    scale_color_manual(values=colors) +
    # facet_grid(rows=vars(comp), cols=vars(dist),
    #     scales="free_y", switch="y", shrink=FALSE) +
    theme_bw(base_size = 9) +
    theme(panel.grid.major = element_line(),
          legend.title=element_blank()) +
    xlab("Time-of-Flight (%)") + ylab(NULL) +
    scale_y_log10(breaks = 10^(-10:10),
                  labels = trans_format("log10", math_format(10^.x)))

ggsave("ct_gauss_rmse.eps", width=6, height=4.5, units="in")

