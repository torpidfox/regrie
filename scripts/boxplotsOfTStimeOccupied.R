# read premodified csv
tsData = read.csv("/Users/dmitrav/Politech/Laboratory/StoсhasticModelling/reGRiE/versions/v.1.2/simulations/3 серия/analysis/modified_target_site_5607590001078610445.csv", header = T)
library(reshape2)

ts.melted = melt(tsData[,-(2:4)], id.vars = "specie")
ts.melted[,2] = ""

require(ggplot2)

# non zero values
ts.melted1 = ts.melted[ts.melted[,"specie"] %in% c("bcd","cad","hb","hkb","Kr"),]
ggplot(data = ts.melted1, aes(x=variable, y=value / 10 ** 4)) + # simulations were for 10^4
    geom_boxplot(aes(fill=specie)) +
    xlab("ТФ") +
    ylab("Время занятости СПТФ, с") + 
    scale_y_continuous(breaks = seq(0,16,1)) +
    guides(fill=guide_legend(title="Типы ТФ"))

# another df with values close to 0
ts.melted2 = ts.melted[ts.melted[,"specie"] %in% c("gt","kni","tll"),]
ggplot(data = ts.melted2, aes(x=variable, y=value / 10 ** 4)) + # simulations were for 10^4
    geom_boxplot(aes(fill=specie)) +
    xlab("ТФ") +
    ylab("Время занятости СПТФ, с") +
    scale_y_continuous(breaks = seq(0,0.5,0.05)) +
    guides(fill=guide_legend(title="Типы ТФ"))