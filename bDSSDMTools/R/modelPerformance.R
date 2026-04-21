#' Calculate Model Performance Metrics
#'
#' Evaluates model predictions against observed presence/absence values using
#' various threshold-dependent and independent metrics. It utilizes functions
#' from the 'sdm' and 'ecospat' packages and calculates AIC and R-squared.
#' Allows selection of metrics based on a specific threshold criterion.
#'
#'
#' @export
#' @importFrom sdm evaluates
#' @importFrom stats glm binomial AIC lm na.omit var cor complete.cases kmeans dist as.matrix sd quantile aggregate family formula coefficients
#' @importFrom ecospat ecospat.boyce
#'

modelPerformance <- function(observed,predicted,index) {

  library(ggplot2)
  library(grid)

  precision_gain = function(TP,FN,FP,TN) {
    n_pos = TP+FN
    n_neg = FP+TN
    prec_gain = 1-(n_pos*FP)/(n_neg*TP)
    prec_gain[TN+FN==0] = 0
    return(prec_gain)
  }

  recall_gain = function(TP,FN,FP,TN) {
    n_pos = TP+FN
    n_neg = FP+TN
    rg = 1-(n_pos*FN)/(n_neg*TP)
    rg[TN+FN==0] = 1
    return(rg)
  }

  # create a table of segments
  .create.segments = function(labels,pos_scores,neg_scores) {
    # reorder labels and pos_scores by decreasing pos_scores, using increasing neg_scores in breaking ties
    new_order = order(pos_scores,-neg_scores,decreasing=TRUE)
    labels = labels[new_order]
    pos_scores = pos_scores[new_order]
    neg_scores = neg_scores[new_order]
    # create a table of segments
    segments = data.frame(pos_score=NA,neg_score=NA,pos_count=0,neg_count=rep(0,length(labels)))
    j = 0
    for (i in seq_along(labels)) {
      if ((i==1)||(pos_scores[i-1]!=pos_scores[i])||(neg_scores[i-1]!=neg_scores[i])) {
        j = j + 1
        segments$pos_score[j] = pos_scores[i]
        segments$neg_score[j] = neg_scores[i]
      }
      segments[j,4-labels[i]] = segments[j,4-labels[i]] + 1
    }
    segments = segments[1:j,]
    return(segments)
  }

  .create.crossing.points = function(points,n_pos,n_neg) {
    n = n_pos+n_neg
    points$is_crossing = 0
    # introduce a crossing point at the crossing through the y-axis
    j = min(which(points$recall_gain>=0))
    if (points$recall_gain[j]>0) { # otherwise there is a point on the boundary and no need for a crossing point
      delta = points[j,,drop=FALSE]-points[j-1,,drop=FALSE]
      if (delta$TP>0) {
        alpha = (n_pos*n_pos/n-points$TP[j-1])/delta$TP
      } else {
        alpha = 0.5
      }
      new_point = points[j-1,,drop=FALSE] + alpha*delta
      new_point$precision_gain = precision_gain(new_point$TP,new_point$FN,new_point$FP,new_point$TN)
      new_point$recall_gain = 0
      new_point$is_crossing = 1
      points = rbind(points,new_point)
      points = points[order(points$index),,drop=FALSE]
    }   
    # now introduce crossing points at the crossings through the non-negative part of the x-axis
    crossings = data.frame()
    x = points$recall_gain
    y = points$precision_gain
    for (i in which((c(y,0)*c(0,y)<0)&(c(1,x)>=0))) {
      cross_x = x[i-1]+(-y[i-1])/(y[i]-y[i-1])*(x[i]-x[i-1])
      delta = points[i,,drop=FALSE]-points[i-1,,drop=FALSE]
      if (delta$TP>0) {
        alpha = (n_pos*n_pos/(n-n_neg*cross_x)-points$TP[i-1])/delta$TP
      } else {
        alpha = (n_neg/n_pos*points$TP[i-1]-points$FP[i-1])/delta$FP
      }
      new_point = points[i-1,,drop=FALSE] + alpha*delta
      new_point$precision_gain = 0
      new_point$recall_gain = recall_gain(new_point$TP,new_point$FN,new_point$FP,new_point$TN)
      new_point$is_crossing = 1
      crossings = rbind(crossings,new_point)
    }
    # add crossing points to the 'points' data frame
    points = rbind(points,crossings)
    points = points[order(points$index,points$recall_gain),2:ncol(points),drop=FALSE]
    rownames(points) = NULL
    points$in_unit_square = 1
    points$in_unit_square[points$recall_gain<0] = 0
    points$in_unit_square[points$precision_gain<0] = 0
    return(points)
  }

  create_prg_curve = function(labels,pos_scores,neg_scores=-pos_scores) { 
    create_crossing_points = TRUE
    n = length(labels)
    n_pos = sum(labels)
    n_neg = n - n_pos
    # convert negative labels into 0s
    labels = 1*(labels==1) 
    segments = .create.segments(labels,pos_scores,neg_scores)
    # calculate recall gains and precision gains for all thresholds
    points = data.frame(index=1:(1+nrow(segments)))
    points$pos_score=c(Inf,segments$pos_score)
    points$neg_score=c(-Inf,segments$neg_score)
    points$TP = c(0,cumsum(segments$pos_count))
    points$FP = c(0,cumsum(segments$neg_count))
    points$FN = n_pos-points$TP
    points$TN = n_neg-points$FP
    points$precision_gain = precision_gain(points$TP,points$FN,points$FP,points$TN)
    points$recall_gain = recall_gain(points$TP,points$FN,points$FP,points$TN)
    if (create_crossing_points) {
      points = .create.crossing.points(points,n_pos,n_neg)
    } else {
      points = points[,2:ncol(points)]
    }
    return(points)
  }

  calc_auprg = function(prg_curve) {
    area = 0
    for (i in 2:nrow(prg_curve)) {
      if (!is.na(prg_curve$recall_gain[i-1]) && (prg_curve$recall_gain[i-1]>=0)) {
        width = prg_curve$recall_gain[i]-prg_curve$recall_gain[i-1]
        height = (prg_curve$precision_gain[i]+prg_curve$precision_gain[i-1])/2
        area = area + width*height
      }
    }
    return(area)
  }

  prg_convex_hull = function(prg_curve) {
    y = prg_curve$precision_gain
    x = prg_curve$recall_gain
    m = length(x)
    y[is.na(x)] = NA
    y_peak = max(which(y==max(y,na.rm=TRUE)),na.rm=TRUE)
    ch = !is.na(y) & ((1:m)>=y_peak)
    ch[(c(Inf,x[1:(m-1)])==x)] = 0
    chi = which(ch==1)
    while (length(chi)>=3) {
      changed = FALSE
      for (i in 3:length(chi)) {
        s1 = (y[chi[i-1]]-y[chi[i-2]]) / (x[chi[i-1]]-x[chi[i-2]])
        s2 = (y[chi[i]]-y[chi[i-1]]) / (x[chi[i]]-x[chi[i-1]])
        if (s1<=1.00001*s2) {
          chi = chi[-(i-1)]
    changed = TRUE
    break
        }
      }
      if (!changed) {
        break
      }
    }
    convex_hull = prg_curve[chi,c("pos_score","neg_score","precision_gain","recall_gain")]
    convex_hull = rbind(c(Inf,-Inf,y[y_peak],-Inf),convex_hull)
    y = convex_hull$precision_gain
    x = convex_hull$recall_gain
    slopes = (c(0,y)-c(y,0))/(c(0,x)-c(x,0))
    convex_hull$f_calibrated_score = 1/(1-slopes[1:nrow(convex_hull)])
    return(convex_hull)
  }

  plot_prg = function(prg_curve,show_convex_hull=TRUE,show_f_calibrated_scores=TRUE) {
    d = prg_curve
    d = d[(!is.na(d$precision_gain))&(!is.na(d$recall_gain)),]
    d2 = d
    d2$precision_gain[d2$in_unit_square==0] = NA
    d3 = d[(d$is_crossing==0)&(d$in_unit_square==1),]
    p = ggplot2::ggplot(d)
    p = p + ggplot2::geom_segment(x=-0.015,xend=1,y=0.00,yend=0.00,color="grey",size=0.1)
    p = p + ggplot2::geom_segment(x=-0.015,xend=1,y=0.25,yend=0.25,color="grey",size=0.1)
    p = p + ggplot2::geom_segment(x=-0.015,xend=1,y=0.50,yend=0.50,color="grey",size=0.1)
    p = p + ggplot2::geom_segment(x=-0.015,xend=1,y=0.75,yend=0.75,color="grey",size=0.1)
    p = p + ggplot2::geom_segment(x=-0.015,xend=1,y=1.00,yend=1.00,color="grey",size=0.1)
    p = p + ggplot2::geom_segment(y=-0.015,yend=1,x=0.00,xend=0.00,color="grey",size=0.1)
    p = p + ggplot2::geom_segment(y=-0.015,yend=1,x=0.25,xend=0.25,color="grey",size=0.1)
    p = p + ggplot2::geom_segment(y=-0.015,yend=1,x=0.50,xend=0.50,color="grey",size=0.1)
    p = p + ggplot2::geom_segment(y=-0.015,yend=1,x=0.75,xend=0.75,color="grey",size=0.1)
    p = p + ggplot2::geom_segment(y=-0.015,yend=1,x=1.00,xend=1.00,color="grey",size=0.1)
    p = p + ggplot2::geom_rect(xmin=0,xmax=1,ymin=0,ymax=1,fill="transparent",color="black",size=0.3)
    p = p + ggplot2::geom_line(ggplot2::aes(x=recall_gain,y=precision_gain),color="lightblue",size=1.5)
    p = p + ggplot2::geom_line(data=d2,ggplot2::aes(x=recall_gain,y=precision_gain),color="blue",size=1.5,na.rm=TRUE)
    p = p + ggplot2::geom_point(data=d3,ggplot2::aes(x=recall_gain,y=precision_gain),color="blue",size=3)
    p = p + ggplot2::xlab("Recall Gain")
    p = p + ggplot2::ylab("Precision Gain")
    p = p + ggplot2::coord_cartesian(xlim=c(-0.015,1.015),ylim=c(-0.015,1.015))
    p = p + ggplot2::theme_bw()
    p = p + ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.ticks.margin = grid::unit(-0.1,"lines"),
      axis.ticks = ggplot2::element_blank()
      )
    convex_hull = prg_convex_hull(prg_curve)
    if (show_convex_hull) {
      p = p + ggplot2::geom_line(data=convex_hull,ggplot2::aes(x=recall_gain,y=precision_gain),color="red",linetype=2)
    }
    if (show_f_calibrated_scores) {
      y = convex_hull$precision_gain
      x = convex_hull$recall_gain
      convex_hull$ya = 0.5*(y+c(0,y[1:(length(y)-1)]))
      convex_hull$xa = 0.5*(x+c(0,x[1:(length(x)-1)]))
      if (nrow(convex_hull)>=3) {
        p = p + ggplot2::geom_text(data=convex_hull[3:nrow(convex_hull),],ggplot2::aes(x=xa,y=ya,label=round(f_calibrated_score,2)),color="red",hjust=0,vjust=0)
      }
    }
    return(p)
  }

  .test = function() {
    dir.create("prg_r_test_fig_auto",showWarnings=FALSE)
    dir.create("prg_r_test_fig_unit",showWarnings=FALSE)
    dir.create("prg_r_test_tab",showWarnings=FALSE)
    dir.create("prg_r_test_auprg",showWarnings=FALSE)
    test_labels = c("0","1","00","01","10","11","001","010","011","100","101","110","0001","0010","0011","0100","0101","0110","0111","1000","1001","1010","1011","1100","1101","1110","01100","0011000","010110100101","000010101100101","1111110010101100101") 
    for (i in seq_along(test_labels)) {
      labels = as.numeric(strsplit(test_labels[i],"")[[1]])
      pos_scores = rev(1:length(labels))
      prg_curve = create_prg_curve(labels,pos_scores)
      p = plot_prg(prg_curve,xlim=c(-0.6,1.1),ylim=c(-0.6,1.1))
      ggplot2::ggsave(paste("prg_r_test_fig_auto/prg_r_test_",i,".pdf",sep=""),p)
      p = plot_prg(prg_curve)
      ggplot2::ggsave(paste("prg_r_test_fig_unit/prg_r_test_",i,".pdf",sep=""),p)
      prg_curve = round(prg_curve,5)
      auprg = round(calc_auprg(prg_curve),5)
      write.csv(prg_curve,file=paste("prg_r_test_tab/prg_r_test_",i,".csv",sep=""),row.names=FALSE)
      write.csv(data.frame(auprg=auprg),file=paste("prg_r_test_auprg/prg_r_test_",i,".csv",sep=""),row.names=FALSE)
    }
  }


  if( is.null(predicted)) {

    predicted.accuracy <- data.frame(threshold=NA, auc=NA, omission.rate=NA, sensitivity=NA, specificity=NA, prop.correct=NA, Kappa=NA, aicc=NA, deviance=NA, boyce=NA, tss=NA, prevalence=NA)

  }

  if( ! is.null(predicted)) {

    keepData <- which(!is.na(observed) & !is.na(predicted))
    observed <- observed[keepData]
    predicted <- predicted[keepData]

    if( length(unique(predicted)) == 1 & length(observed) > 1 ) { predicted[which(observed == 1)] <- predicted[which(observed == 1)] + 0.01 }

    prg_curve <- create_prg_curve(labels = observed, pos_scores = predicted)
    auprg_score <- calc_auprg(prg_curve)

    predicted.accuracy <- sdm::evaluates( x=observed , p=predicted )
    predicted.accuracy <- data.frame(criteria=predicted.accuracy@threshold_based$criteria,
                                     threshold=predicted.accuracy@threshold_based$threshold,
                                     auc=predicted.accuracy@statistics$AUC,
                                     auc.prg=auprg_score,
                                     omission.rate=predicted.accuracy@threshold_based$ommission,
                                     sensitivity=predicted.accuracy@threshold_based$sensitivity,
                                     specificity=predicted.accuracy@threshold_based$specificity,
                                     prop.correct=predicted.accuracy@threshold_based$ommission + predicted.accuracy@threshold_based$commission - 1,
                                     Kappa=predicted.accuracy@threshold_based$Kappa,
                                     aicc=NA,
                                     deviance=NA,
                                     boyce=NA,
                                     tss=predicted.accuracy@threshold_based$sensitivity + predicted.accuracy@threshold_based$specificity - 1,
                                     prevalence.obs=sum(observed == 1) / sum(observed == 0),
                                     prevalence.pred=predicted.accuracy@threshold_based$prevalence)
    
    predicted.accuracy$aicc <- AIC(glm(observed~predicted, family = binomial))
    predicted.accuracy$deviance <- summary(lm(observed~predicted))$r.squared

    boyce.a <- predicted
    boyce.a <- boyce.a[!is.na(boyce.a)]
    boyce.p <- predicted[which(observed==1)]
    boyce.p <- boyce.p[!is.na(boyce.p)]

    if( length(boyce.p) == 1) { boyce.p <- c(boyce.p-0.01,boyce.p-0.005,boyce.p,boyce.p+0.005,boyce.p+0.01) }

    tryCatch( boyce.value <- ecospat::ecospat.boyce(boyce.a , boyce.p ,PEplot=FALSE)$Spearman.cor, error=function(e) { Error <<- TRUE })
    tryCatch( boyce.value <- ecospat::ecospat.boyce(boyce.a , boyce.p ,PEplot=FALSE)$cor, error=function(e) { Error <<- TRUE })
    if( is.null(boyce.value) ) { boyce.value <- NA }

    val <- mean(boyce.p,na.rm=T)

    while( is.na(boyce.value) ) {

      val <- val - 0.001

      tryCatch( boyce.value <- ecospat::ecospat.boyce(boyce.a , c(boyce.p, val, val) ,PEplot=FALSE)$Spearman.cor, error=function(e) { Error <<- TRUE })
      tryCatch( boyce.value <- ecospat::ecospat.boyce(boyce.a , c(boyce.p, val, val) ,PEplot=FALSE)$cor, error=function(e) { Error <<- TRUE })
      if( is.null(boyce.value) ) { boyce.value <- NA }

      if( val < 0 ) { boyce.value <- 0; break }

    }

    predicted.accuracy$boyce <- boyce.value

    if( index == "tss" | index == "auc.prg" | index == "auc" | index == "boyce" ) { predicted.accuracy <- predicted.accuracy[which(predicted.accuracy$criteria == "max(se+sp)"),] }
    if( index == "minimumTrain" ) { predicted.accuracy <- predicted.accuracy[which(predicted.accuracy$criteria == "P0"),] }
    if( index == "minimumTrain95" ) { predicted.accuracy <- predicted.accuracy[which(predicted.accuracy$criteria == "P5"),] }
    if( index == "minimumTrain90" ) { predicted.accuracy <- predicted.accuracy[which(predicted.accuracy$criteria == "P10"),] }

    predicted.accuracy[predicted.accuracy$tss < 0,"tss"] <- 0
    # predicted.accuracy <- predicted.accuracy[,-1]

  }

  rownames(predicted.accuracy) <- NULL
  return(predicted.accuracy)

}
