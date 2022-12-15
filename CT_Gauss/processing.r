library(matrixStats)

process_rmse <- function() {

    print("Post-processing RMSE", quote=FALSE)
    print(dist, quote=FALSE)

    trials = 100
    steps = 1000

    house_r_se <- matrix(, nrow=steps, ncol=trials)
    ukf_r_se   <- matrix(, nrow=steps, ncol=trials)
    cut4_r_se  <- matrix(, nrow=steps, ncol=trials)
    cut6_r_se  <- matrix(, nrow=steps, ncol=trials)
    cut8_r_se  <- matrix(, nrow=steps, ncol=trials)


    house_v_se <- matrix(, nrow=steps, ncol=trials)
    ukf_v_se   <- matrix(, nrow=steps, ncol=trials)
    cut4_v_se  <- matrix(, nrow=steps, ncol=trials)
    cut6_v_se  <- matrix(, nrow=steps, ncol=trials)
    cut8_v_se  <- matrix(, nrow=steps, ncol=trials)

    state_interp <- function(filename) {
        t <- as.matrix(read.csv(filename))[,1]
        x <- as.matrix(read.csv(filename))[,2:7]
        xi <- matrix(, nrow=steps, 6)
        for (i in 1:5) {
            xi[,i] <- approx(t, x[,i], n=steps, method="linear")$y
        }
        return (xi)
    }
    
    
    tru_file <- paste("out/ct_true_1.csv")
    for (j in 1:trials) {

        print(paste("   Trial ", j), quote=FALSE)


        house_file <- paste("out/ct_est_house_1_", j, ".csv", sep="")
        ukf_file   <- paste("out/ct_est_ukf_1_", j, ".csv", sep="")
        cut4_file  <- paste("out/ct_est_cut4_1_", j, ".csv", sep="")
        cut6_file  <- paste("out/ct_est_cut6_1_", j, ".csv", sep="")
        cut8_file  <- paste("out/ct_est_cut8_1_", j, ".csv", sep="")


        tru <- state_interp(tru_file)

        house_est <- state_interp(house_file)
        ukf_est   <- state_interp(ukf_file)
        cut4_est  <- state_interp(cut4_file)
        cut6_est  <- state_interp(cut6_file)
        cu86_est  <- state_interp(cut8_file)


        house_r_se[,j] <- ((house_est[,1] - tru[,1])^2 + (house_est[,3] - tru[,3])^2)
        ukf_r_se  [,j] <- ((ukf_est  [,1] - tru[,1])^2 + (ukf_est  [,3] - tru[,3])^2)
        cut4_r_se [,j] <- ((cut4_est [,1] - tru[,1])^2 + (cut4_est [,3] - tru[,3])^2)
        cut6_r_se [,j] <- ((cut6_est [,1] - tru[,1])^2 + (cut6_est [,3] - tru[,3])^2)
        cut8_r_se [,j] <- ((cut6_est [,1] - tru[,1])^2 + (cut6_est [,3] - tru[,3])^2)


        house_v_se[,j] <- ((house_est[,2] - tru[,2])^2 + (house_est[,2] - tru[,2])^2)
        ukf_v_se  [,j] <- ((ukf_est  [,2] - tru[,2])^2 + (ukf_est  [,2] - tru[,2])^2)
        cut4_v_se [,j] <- ((cut4_est [,2] - tru[,2])^2 + (cut4_est [,2] - tru[,2])^2)
        cut6_v_se [,j] <- ((cut6_est [,2] - tru[,2])^2 + (cut6_est [,2] - tru[,2])^2)
        cut8_v_se [,j] <- ((cut6_est [,2] - tru[,2])^2 + (cut6_est [,2] - tru[,2])^2)


    }

    # --- RMSE

    rmse_pos <- c(sqrt(mean(house_r_se)),
                  sqrt(mean(ukf_r_se)),
                  sqrt(mean(cut4_r_se)),
                  sqrt(mean(cut6_r_se)),
                  sqrt(mean(cut8_r_se)))

    rmse_vel <- c(sqrt(mean(house_v_se)),
                  sqrt(mean(ukf_v_se)),
                  sqrt(mean(cut4_v_se)),
                  sqrt(mean(cut6_v_se)),
                  sqrt(mean(cut8_v_se)))

    filter <- c("HOUSE", "UKF", "CUT-4", "CUT-6", "CUT-8")

    rmse <- data.frame(filterno = 1:5, filter, rmse_pos, rmse_vel)

    file <- paste("R_testing/rmse",  ".csv", sep="")
    write.csv(rmse, file, row.names=FALSE, quote=FALSE)

    # --- RMSE vs. time

    time_percent <- seq(0, 100, length.out=steps)

    rmse_t <- data.frame(time_percent, 
          house_pos = sqrt(rowMeans(house_r_se)),
          house_vel = sqrt(rowMeans(house_v_se)),
          ukf_pos   = sqrt(rowMeans(ukf_r_se)),
          ukf_vel   = sqrt(rowMeans(ukf_v_se)),
          cut4_pos  = sqrt(rowMeans(cut4_r_se)),
          cut4_vel  = sqrt(rowMeans(cut4_v_se)),
          cut6_pos  = sqrt(rowMeans(cut6_r_se)),
          cut6_vel  = sqrt(rowMeans(cut6_v_se)))

    file <- paste("R_testing/rmse_t", ".csv", sep="")
    write.csv(rmse_t, file, row.names=FALSE, quote=FALSE)

    # --- Convergence

    convmse <- function(se) {
        cmse <- vector(, length = trials)
        for (j in 1:trials) {
            cmse[j] <- mean(se[,1:j])        
        }
        return(cmse)
    }

    convsse <- function(se) {
        csse <- vector(, length = trials)
        for (j in 1:trials) {
            csse[j] <- sd(se[,1:j])
        }
        return(csse)
    }

    conv <- data.frame(ntrials = 1:trials,
                house_r_mse = convmse(house_r_se), house_r_sse = convsse(house_r_se),
                house_v_mse = convmse(house_v_se), house_v_sse = convsse(house_v_se),
                ukf_r_mse = convmse(ukf_r_se), ukf_r_sse = convsse(ukf_r_se),
                ukf_v_mse = convmse(ukf_v_se), ukf_v_sse = convsse(ukf_v_se),
                cut4_r_mse = convmse(cut4_r_se), cut4_r_sse = convsse(cut4_r_se),
                cut4_v_mse = convmse(cut4_v_se), cut4_v_sse = convsse(cut4_v_se),
                cut6_r_mse = convmse(cut6_r_se), cut6_r_sse = convsse(cut6_r_se),
                cut6_v_mse = convmse(cut6_v_se), cut6_v_sse = convsse(cut6_v_se))

    file <- paste("R_testing/conv", ".csv", sep="")
    write.csv(conv, file, row.names = FALSE, quote=FALSE)

}

process_box <- function() {

    print("Post-processing for boxplot", quote=FALSE)

    print("Post-processing", quote=FALSE)

    trialno <- 1:100

    trials = length(trialno)

    f = 10
    tmin = 1
    kstart = f * tmin + 1

    house_r_err <- c()
    ukf_r_err   <- c()
    cut4_r_err  <- c()
    cut6_r_err  <- c()
    cut8_r_err  <- c()


    house_v_err <- c()
    ukf_v_err   <- c()
    cut4_v_err  <- c()
    cut6_v_err  <- c()
    cut8_v_err  <- c()

    house_w_err <- c()
    ukf_w_err   <- c()
    cut4_w_err  <- c()
    cut6_w_err  <- c()
    cut8_w_err  <- c()

    # replace with 1:trails
    # for testing purposes
    for (j in 1:trials) {

        k <- trialno[j]

        print(paste("   Trial ", k), quote=FALSE)

        tru_file <- paste("out/ct_true_1", ".csv", sep="")

        house_file <- paste("out/ct_est_house_1_", k, ".csv", sep="")
        ukf_file   <- paste("out/ct_est_ukf_1_",   k, ".csv", sep="")
        cut4_file  <- paste("out/ct_est_cut4_1_",  k, ".csv", sep="")
        cut6_file  <- paste("out/ct_est_cut6_1_",  k, ".csv", sep="")
        cut8_file  <- paste("out/ct_est_cut8_1_",  k, ".csv", sep="")


        tru <- as.matrix(read.csv(tru_file))[,2:7]

        kstop <- length(tru[,1])

        tru <- tru[kstart:kstop,]

        house_est <- as.matrix(read.csv(house_file))[kstart:kstop,2:7]
        ukf_est   <- as.matrix(read.csv(ukf_file))  [kstart:kstop,2:7]
        cut4_est  <- as.matrix(read.csv(cut4_file)) [kstart:kstop,2:7]
        cut6_est  <- as.matrix(read.csv(cut6_file)) [kstart:kstop,2:7]
        cut8_est  <- as.matrix(read.csv(cut8_file)) [kstart:kstop,2:7]


        house_r_err <- c(sqrt((house_est[,1] - tru[,1])^2 + (house_est[,3] - tru[,3])^2), house_r_err)
        ukf_r_err   <- c(sqrt((ukf_est[,1] - tru[,1])^2 + (ukf_est[,3] - tru[,3])^2), ukf_r_err)
        cut4_r_err  <- c(sqrt((cut4_est[,1] - tru[,1])^2 + (cut4_est[,3] - tru[,3])^2), cut4_r_err)
        cut6_r_err  <- c(sqrt((cut6_est[,1] - tru[,1])^2 + (cut6_est[,3] - tru[,3])^2), cut6_r_err)
        cut8_r_err  <- c(sqrt((cut8_est[,1] - tru[,1])^2 + (cut8_est[,3] - tru[,3])^2), cut8_r_err)

        house_v_err <- c(sqrt((house_est[,2] - tru[,2])^2 + (house_est[,4] - tru[,4])^2), house_v_err)
        ukf_v_err   <- c(sqrt((ukf_est[,2] - tru[,2])^2 + (ukf_est[,4] - tru[,4])^2), ukf_v_err)
        cut4_v_err  <- c(sqrt((cut4_est[,2] - tru[,2])^2 + (cut4_est[,4] - tru[,4])^2), cut4_v_err)
        cut6_v_err  <- c(sqrt((cut6_est[,2] - tru[,2])^2 + (cut6_est[,4] - tru[,4])^2), cut6_v_err)
        cut8_v_err  <- c(sqrt((cut8_est[,2] - tru[,2])^2 + (cut8_est[,4] - tru[,4])^2), cut8_v_err)

        house_w_err <- c(abs(house_est[,5] - tru[,5]), house_w_err)
        ukf_w_err   <- c(abs(ukf_est[,5] - tru[,5]), ukf_w_err)
        cut4_w_err  <- c(abs(cut4_est[,5] - tru[,5]), cut4_w_err)
        cut6_w_err  <- c(abs(cut6_est[,5] - tru[,5]), cut6_w_err)
        cut8_w_err  <- c(abs(cut8_est[,5] - tru[,5]), cut8_w_err)

    }

    names <- c("HOUSE", "UKF", "CUT-4", "CUT-6", "CUT-8")

    r_err <- data.frame(HOUSE=house_r_err, UKF=ukf_r_err, CUT4=cut4_r_err, CUT6=cut6_r_err, CUT8=cut8_r_err)
    v_err <- data.frame(HOUSE=house_v_err, UKF=ukf_v_err, CUT4=cut4_v_err, CUT6=cut6_v_err,CUT8=cut8_v_err)
    w_err <- data.frame(HOUSE=house_w_err, UKF=ukf_w_err, CUT4=cut4_w_err, CUT6=cut6_w_err,CUT8=cut8_w_err)

    pos_file <- paste("R_testing/pos_err", ".csv", sep="")
    vel_file <- paste("R_testing/vel_err", ".csv", sep="") 
    omega_file <- paste("R_testing/omega_err", ".csv", sep="") 


    write.csv(r_err, pos_file, row.names = FALSE, quote=FALSE)
    write.csv(v_err, vel_file, row.names = FALSE, quote=FALSE)
    write.csv(w_err, omega_file, row.names = FALSE, quote=FALSE)


}

# process_rmse("gauss")
process_rmse()

# process_box()
# process_box("pearson")
