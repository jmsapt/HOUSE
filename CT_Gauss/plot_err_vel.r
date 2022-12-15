library(ggplot2)
library(scales)

pos_err_gauss <- as.matrix(read.csv("R_testing/vel_err.csv"))


Ng <- nrow(pos_err_gauss)

err <- c(pos_err_gauss)

dst <- c(rep("Gaussian Distributions", 5*Ng))

fil <- c(rep("HOUSE", Ng),
            rep("UKF",   Ng),
            rep("CUT-4", Ng),
            rep("CUT-6", Ng),
            rep("CUT-8", Ng))

cmp <- c(rep("Velocity Error (m/s)",   5*Ng))

data <- data.frame(err, dst, fil, cmp)

filter <- c("HOUSE", "UKF", "CUT-4", "CUT-6", "CUT-8")
data$fil <- factor(data$fil, levels=filter)

colors <- c("#cc3311", "#0077bb", "#ee7733", "#009988", "#ee3377", "#37c920")

ggplot(data, aes(x=fil, y=err, fill=fil)) +
    scale_fill_manual(values=colors) +
    geom_violin(trim=FALSE) + 
    geom_boxplot(width=0.1, outlier.shape=NA, fill="white") +
    theme_bw(base_size = 9) + 
    theme(panel.grid.major = element_line(),
          legend.position = "none") +
    facet_grid(rows=vars(cmp), cols=vars(dst),
        scales="free_y", switch="y", shrink=FALSE) +
    xlab(NULL) + ylab(NULL) +
    scale_y_log10(breaks = 10^(-10:10),
                  labels = trans_format("log10", math_format(10^.x)))

ggsave("velocity_err.eps", width=6, height=6, units="in")

