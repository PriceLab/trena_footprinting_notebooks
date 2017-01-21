mattcc <- function(true_classes, predicted_classes) {
        
    pos <- true_classes == 1
    neg <- true_classes == 0
    
    TP <- as.numeric(sum( predicted_classes[pos] == true_classes[pos] ))
    TN <- as.numeric(sum( predicted_classes[neg] == true_classes[neg] ))
    FP <- as.numeric(sum( predicted_classes[pos] != true_classes[pos] ))
    FN <- as.numeric(sum( predicted_classes[neg] != true_classes[neg] )) 
    
    N1 <- TP*TN
    N2 <- FP*FN
    D1 <- TP+FP
    D2 <- TP+FN
    D3 <- TN+FP
    D4 <- TN+FN
    
    N <- N1-N2
    D <- D1*D2*D3*D4
    
    if (D==0) {
        M <- 0
    } else {
        M <- N/sqrt(D)
    }
    
    return(M)
}

mattcc.curve <- function(true.values, pred.probs, num.points=101) {
    threshholds <- seq(from=0, to=1, length.out=num.points)
    thresh.inds <- 1:length(threshholds)

    mattcc.data <- sapply(thresh.inds,
                      function(t) mattcc(true.values,as.numeric(pred.probs>threshholds[t]))
                   )
    return(mattcc.data)
}

# make dataframe for preditions vs ChIPseq
make.pred.df.from.model <- function(gbdt.model,X.testdata,y.testdata) {
    
    preds.test <- predict(gbdt.model, X.testdata, missing=NA)
    pred.df <- data.frame("ChIPseq.bound"=y.testdata,"Prediction"=preds.test)
    
    pred.df$Model.Name = gbdt.model$Model.Name
    return(pred.df)
    
}

make.pred.df.from.glm <- function(glm.model,df.testdata) {
    
    preds.test <- predict.glm(object=glm.model, newdata=df.testdata, type="response")
    
    pred.df <- data.frame("ChIPseq.bound"=df.testdata$CHIPseq.bound,"Prediction"=preds.test)
    
    pred.df$Model.Name = glm.model$Model.Name
    return(pred.df)
    
}

# compute roc curve and prec/recall-like stats
make.stats.df.from.preds <- function(pred.df, num.roc.points=101) {
    
    preds <- pred.df$Prediction
    cs.vals <- pred.df$ChIPseq.bound

    roc.tf <- roc(cs.vals, preds,
                  ret=c("threshold", "sens", "spec", "ppv", "npv", "acc"))
    stats.tf <- coords((roc.tf),
       seq(from=0, to=1, length.out=num.roc.points),
       ret=c("threshold", "sens", "spec", "ppv", "npv", "acc"))
    mattcc.tf <- mattcc.curve(cs.vals, preds, num.points=num.roc.points)
    df.stats <- as.data.frame(t(stats.tf))
    df.stats$MattCC <- mattcc.tf
    df.stats$Model.Name <- pred.df$Model.Name[1]
    
    return(df.stats)
}

