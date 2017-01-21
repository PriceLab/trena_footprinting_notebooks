plot.roc.curve <- function(stat.dataframe, printAUC=FALSE) {
    
    if (printAUC == TRUE){
        specif.sense.df <- stat.dataframe[,c("specificity","sensitivity")]
        ordered.specif.sense.df <- specif.sense.df[with(specif.sense.df, order(specificity)), ]
        ss.x <- ordered.specif.sense.df$specificity
        ss.y <- ordered.specif.sense.df$sensitivity
        AUC <- trapz(ss.x,ss.y)
    }
    
    p <- ggplot(stat.dataframe) +
    geom_line(aes(x=specificity, y=sensitivity, color=Model.Name)) +
    geom_abline(intercept = 1, slope = 1, color='lightgrey') +
    scale_x_reverse(lim=c(1,0)) + scale_y_continuous(lim=c(0,1)) +
    coord_fixed(ratio=1) +
    labs(x="Specificity", y="Sensitivity") +
    theme_minimal()
    if (printAUC == TRUE) {
        p <- p + annotate("text", x=0.25, y=0.2, label=paste('AUC = ', format(AUC,digits=5), sep=''))
    }
    p
}

plot.precrecall.curve <- function(stat.dataframe) {
    ggplot(stat.dataframe) + 
    geom_line(aes(x=sensitivity, y=ppv, color=Model.Name)) + 
    scale_x_continuous(lim=c(0,1)) + scale_y_continuous(lim=c(0,1)) +
    coord_fixed(ratio=1) +
    labs(x="Precision", y="Recall") +
    theme_minimal()
}

plot.mattcc.curve <- function(stat.dataframe) {
    ggplot(stat.dataframe) +
    geom_line(aes(x=threshold, y=MattCC, color=Model.Name)) +
    scale_x_continuous(lim=c(0,1)) + scale_y_continuous(lim=c(0,1)) +
    coord_fixed(ratio=1) +
    labs(x="Threshold", y="Matthews correlation coefficient") +
    theme_minimal()
}

plot.all.curve <- function(stat.dataframe) {
    stat.dataframe.tidy <- gather(stat.dataframe, 'Statistic', 'Value', -threshold, -Model.Name)    
    ggplot(stat.dataframe.tidy) +
    geom_line(aes(x=threshold, y=Value, color=Statistic)) +
    scale_x_continuous(lim=c(0,1)) + scale_y_continuous(lim=c(0,1)) +
    coord_fixed(ratio=1) +
    labs(x="Threshold") +
    theme_minimal()
}


