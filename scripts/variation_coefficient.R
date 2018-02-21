
rep_data = read.csv("/Volumes/Seagate Expansion Drive/simulations/reGRiE2/Kr_60sec_repression_on/rep_variances_related_to_mean.txt", header = F, sep = " ")
norep_data = read.csv("/Volumes/Seagate Expansion Drive/simulations/reGRiE2/Kr_60sec_repression_off/variances_related_to_mean.txt", header = F, sep = " ")

rep_data$values = sqrt(rep_data$V1) / rep_data$V2
norep_data$values = sqrt(norep_data$V1) / norep_data$V2

joint_df = data.frame(time = seq(5,60,5), values = rep_data$values, type = "С репрессией")
joint_df = rbind(joint_df, data.frame(time = seq(5,60,5), values = norep_data$values, type = "Без репрессии"))

library(ggplot2)
ggplot(joint_df, aes(x = time, y = values, colour = type)) +
    geom_line() + geom_point() +
    scale_x_continuous(breaks = seq(5,60,5)) +
    scale_y_continuous(limits = c(0.05,0.13), breaks = seq(0.05,0.13,0.01)) +
    xlab("Время, с") +
    ylab("Значения коэффициента вариации") +
    guides(colour = guide_legend(title = "Эксперименты:"))
