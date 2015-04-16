################################################################################
## Functions to run a histogram based naive bayes
## ToDO: change inputs to X,Y,weights, ... or use formula interface

################################################################################

## naiveBayes
##==============================================================================
#'  Naive bayes classifier using histograms and shrinkage
#'
#'  Fits a naive bayes model to continous and categorical/factor predictors.
#'    Continous predictors are first binned, then estimates shrunk towards zero.
##  Inputs:
#'  @param vars the names or column numbers of specific predictors.
#'  @param X a data.frame of predictors, can include continuous and
#'    categorical/factors along with \code{X$type} (linked or unlinked) and
#'    \code{X$weight}
#'  @param df the degrees of freedom for each component density. if vector, each
#'    predictor can use a different df
#'  @param nbins the number of bins for continuous predictors
#'  @param partition for binning; indicates if breaks generated from quantiles
#'    or equal spacing
##  Outputs:
#'  @return BF a bayes factor object; list of component bayes factors
##  Notes:
#' @description After binning, this adds pseudo counts to each bin count to give
#'  df approximate degrees of freedom. If partition=quantile, this does not
#'  assume a continuous uniform prior over support, but rather a discrete uniform
#'  over all (unlabeled) observations points.
#'  @examples
#'  # See vignette: "Statistical Methods for Crime Series Linkage" for usage.  
#'  @export
##==============================================================================
naiveBayes <- function(vars,X,df=20,nbins=30,partition=c('quantile','width')){
  partition = match.arg(partition)
  nvars = length(vars)
  df = rep(df,length=nvars)
  BF = vector("list",nvars)
  for(j in 1:nvars){
    var = vars[j]
    var.type = class(X[,var])
    if(var.type %in% c('numeric','integer')){
      bks = make.breaks(X[,var],partition,nbins=nbins)
    } else  bks = NULL
    BF[[j]] = getBF(X,var,breaks=bks,df=df[j])
  }
  names(BF) = vars
  return(BF)
}


## predictNB
##==============================================================================
#'  Generate prediction (sum of log bayes factors) from a \code{naiveBayes} object
##  Inputs:
#'  @param NB a naive bayes object from \code{\link{naiveBayes}}
#'  @param X data.frame of new predictors, column names must match NB names
#'  @param components (logical) return the log bayes factors from each component
#'    or return the sum of log bayes factors
#'  @param vars the names or column numbers of specific predictors. If NULL, then
#'    all predictors will be used
##  Outputs:
#'  @return BF if \code{components = FALSE}, the sum of log bayes factors, if
#'    \code{components = TRUE} the component bayes factors (useful for plotting)
##  Notes:
#'  @description This does not include the log prior odds, so will be off by a
#'    constant
#'  @examples
#'  # See vignette: "Statistical Methods for Crime Series Linkage" for usage.     
#'  @export
##  Make into S3 class
##==============================================================================
predictNB <- function(NB,X,components=FALSE,vars=NULL){
  if(is.null(vars))  vars = names(NB)
  nvars = length(vars)
  BF = matrix(NA,nrow(X),nvars)
  for(j in 1:nvars){
    var = vars[j]
    BF[,j] = predictBF(NB[[var]],X[,var],log=TRUE)
  }
  colnames(BF) = vars
  if(components) return(BF)
  BF = rowSums(BF)
  return(BF)
}


## getBF
##==============================================================================
#'  Estimates the bayes factor for continous and categorical predictors.
#'
#'  Continous predictors are first binned, then estimates shrunk towards zero.
##  Inputs:
#'  @param X data.frame of predictors, can include continuous and
#'    categorical/factors along with X$type (linked or unlinked) and X$weight
#'  @param var the names or column number of one specific predictor
#'  @param breaks - set of break point for continuous predictors or NULL for
#'    categorical or discrete
#'  @param df the effective degrees of freedom for the cetegorical density
#'    estimates
##  Outputs:
#'  @return data.frame containing the levels/categories with estimated Bayes factor
##  Notes:
#'  @description This adds pseudo counts to each bin count to give df effective
#'    degrees of freedom. Must have all possible factor levels and must be of
#'    factor class.
#'  @note Give linked and unlinked a different prior according to sample size
#'  @examples
#'  # See vignette: "Statistical Methods for Crime Series Linkage" for usage.     
#'  @export
##==============================================================================
getBF <- function(X,var,breaks=NULL,df=5){
  replaceNA <- function(x,r=0) as.numeric(ifelse(is.na(x),r,x))
  linked = (X$type=='linked')
  var = as.character(var)
  var.type = class(X[,var])
  x = X[,var]
  weights = X$weight
  if(var.type %in% c('numeric','integer')){
    n.bks = length(breaks)
    x = cut(x,breaks=breaks,include.lowest=TRUE)
    x = addNA(x,ifany=TRUE)
    t.linked = tapply(weights[linked],x[linked],sum)
    t.unlinked = tapply(weights[!linked],x[!linked],sum)
    fromto = data.frame(from=c(breaks[-n.bks]),to=c(breaks[-1]))
    if(nrow(fromto)<length(t.linked)) fromto = rbind(fromto,c(NA,NA))
    E = data.frame(fromto,
                   value=levels(x),
                   N.linked=replaceNA(t.linked),
                   N.unlinked=replaceNA(t.unlinked))
        var.type = 'numeric'
  }
  if(var.type %in% c('factor','character')){
    x = as.factor(x)
    x = addNA(x,ifany=TRUE)
    t.linked = tapply(weights[linked],x[linked],sum)
    t.unlinked = tapply(weights[!linked],x[!linked],sum)
    E = data.frame(value=names(t.linked),
                   N.linked=replaceNA(t.linked),
                   N.unlinked=replaceNA(t.unlinked))
    attr(E,'levels') = levels(x)
    var.type = 'categorical'
  }
  #- get 'a' based on degrees of freedom
  df2a <- function(df,k,N)  (N/k) * ((k-df)/(df-1)) # k is # levels
  a2df <- function(a,k,N) k*(N+a)/(N+k*a)
  nlevs = nrow(E)
  df = min(df,nlevs-1e-8)
  a.linked = df2a(df,k=nlevs,N=sum(E$N.linked))
  a.unlinked = df2a(df,k=nlevs,N=sum(E$N.unlinked))

  getP <- function(N,a) (N+a)/sum(N+a)
  E$p.linked = getP(E$N.linked,a.linked)
  E$p.unlinked = getP(E$N.unlinked,a.unlinked)
  #E = transform(E,p.linked = getP(N.linked,a.linked),
  #                p.unlinked = getP(N.unlinked,a.unlinked))
  E$BF = with(E, p.linked/p.unlinked)
  E[is.na(E$BF),'BF'] = 1           # Set 0/0=1
  attr(E,'var') = var               # add variable name
  attr(E,'breaks') = breaks         # add breaks (or NULL)
  attr(E,'a') = c(linked=a.linked,unlinked=a.unlinked) # add shrinkage parameters
  attr(E,'df') = df
  attr(E,'df2') = c(linked=a2df(a.linked,k=nlevs,N=sum(E$N.linked)),
                                unlinked=a2df(a.unlinked,k=nlevs,N=sum(E$N.unlinked)))
  attr(E,'type') = var.type
  return(E)
}


## predictBF
##==============================================================================
#'  Generate prediction of a component bayes factor
##  Inputs:
#'  @param BF bayes factor data.frame from \code{\link{getBF}}
#'  @param x vector of new predictor values
#'  @param log (logical) if \code{TRUE}, return the \bold{log} bayes factor estimate
##  Outputs:
#'  @return estimated (log) bayes factor from a single predictor
#'  @description This does not include the log prior odds, so will be off by a constant
#'  @examples
#'  # See vignette: "Statistical Methods for Crime Series Linkage" for usage.     
#'  @export
##==============================================================================
predictBF <- function(BF,x,log=TRUE){
  breaks = attr(BF,'breaks')
  if(!is.null(breaks)){
    x = cut(x,breaks=breaks,include.lowest=TRUE)
  }
  ind = match(as.character(x),as.character(BF$value))
  bf = BF$BF[ind]
  bf[is.na(bf)] = 1    # if data in x is outside of breaks, then set bf to 1
  if(log) bf = log(bf)
  attr(bf,'log') = log
  return(bf)
}


## make.breaks
##==============================================================================
#'  Make break points for binning continuous predictors
##  Inputs:
#'  @param x observed sample
#'  @param type one of \code{width} (fixed width) or \code{quantile} binning
#'  @param nbins number of bins
#'  @param binwidth bin width; corresponds to quantiles if type='quantile'
##  Outputs:
#'  @return set of unique break points for binning
#'  @keywords internal
##==============================================================================
make.breaks <- function(x,type='quantile',nbins=NULL,binwidth=NULL){
  type <- match.arg(type, c("quantile", "width"))
  if ((!is.null(nbins) && !is.null(binwidth)) || (is.null(nbins) &&
                                                    is.null(binwidth))) {
    stop("Specify exactly one of nbins or width")
  }
  if (type == "width") {
    rng <- range(x, na.rm = TRUE, finite = TRUE)
    if (!is.null(binwidth)) bks = unique(c(rng[1],seq(rng[1],rng[2],by=binwidth),rng[2]))
    else                    bks = seq(rng[1], rng[2], length = nbins + 1)
  }
  if (type == "quantile") {
    if (!is.null(binwidth)) probs <- seq(0, 1, by = binwidth)
    else                    probs <- seq(0, 1, length = nbins + 1)
    bks = quantile(x, probs, na.rm = TRUE)
  }
  return(sort(unique(bks)))
}


#==============================================================================#
# Plotting functions
#==============================================================================#

## plotBKG
##==============================================================================
#'  Generate a background plot
#'
#'  This facilitates a common plot background
#'  @param xlim range of x-axis
#'  @param ylim range of y-axis
#'  @param background (logical) should a background color be used
#'  @param x.minor values for minor axis
#'  @param x.major values for major axis
#'  @param grid.lines (logical) should grid lines be added to plot
#'  @param glwd1 linewidth of gridlines
#'  @param glwd2 linewidth of gridlines
#'  @param bkg.col color of plot background
#'  @param boxed (logical) should plot be boxed
#'  @param \ldots other arguments passed to plot
#'  @return makes a plot
#'  @keywords internal
##==============================================================================
plotBKG <- function(xlim,ylim,background=TRUE,x.minor,x.major,grid.lines=TRUE,
                        glwd1=2,glwd2=1,bkg.col='grey90',boxed=TRUE,...){
  plot(xlim,ylim,typ='n',
       ylab='',xlab='',
       col.axis='grey50',cex.axis=.8,tcl=-.25,las=1,...)
  if(background){    # Add background color
    rng = par('usr')
    if(par('ylog')) rng[c(3,4)] = 10^(rng[c(3,4)])
    if(par('xlog')) rng[c(1,2)] = 10^(rng[c(1,2)])
    rect(rng[1],rng[3],rng[2],rng[4],col=bkg.col,border="transparent")
  }
  if(grid.lines){
    abline(h=axTicks(2),col = 'white',lty=1,lwd=glwd1)
    if(missing(x.major))  abline(v=axTicks(1),col = 'white',lty=1,lwd=glwd1)
    else                  abline(v=x.major   ,col = 'white',lty=1,lwd=glwd1)
    get.minor <- function(side=c(1,2)){  # return equal spaced minor axis
      xt = axTicks(side)
      nt = length(xt)
      axlog = ifelse(side==1,par('xlog'),par('ylog'))
      if(axlog) xt = log10(xt)
      minor = diff(xt[1:2])/2 + xt
      if(axlog) minor = 10^minor
    return(minor)
    }
    if(missing(x.minor)) abline(v=get.minor(1),col = 'white',lty=1,lwd=glwd2)
    else abline(v=x.minor,col='white',lty=1,lwd=glwd2)
    abline(h=get.minor(2),col = 'white',lty=1,lwd=glwd2)
    if(boxed) box(col='grey50')
  }
}


## plotBF
##==============================================================================
#'  plots 1D bayes factor
#'  @param BF Bayes Factor
#'  @param log.scale (logical)
#'  @param show.legend (logical)
#'  @param x.rng range of x-axis
#'  @param ylim range of y-axis
#'  @param cols Colors for plotting. First element is for linkage, second unlinked
#'  @param \ldots arguemnts passed into \code{\link{plotBKG}}
#'  @return plot of Bayes factor
#'  @examples
#'  # See vignette: "Statistical Methods for Crime Series Linkage" for usage.     
#'  @export
##==============================================================================
plotBF <- function(BF,log.scale=TRUE,show.legend=TRUE,x.rng,ylim,
  cols = c(color('darkred',alpha=.75),color('darkblue',alpha=.75)),...){
  if(missing(ylim)){
    ylim = range(BF$BF,na.rm=TRUE,finite=TRUE)
    if(log.scale) ylim = c(-1,1)*min(12,max(abs(log(ylim))))
  }
  if(attr(BF,'type') %in% "numeric"){
    BF = BF[!is.na(BF$to),]     # Remove NAs
    n = nrow(BF)
    xx = c(BF$from[1],BF$to)
    yy = c(BF$BF,BF$BF[n]) # Add extra line for step plotting
    if(log.scale) yy = log(yy)
    if(missing(x.rng))  x.rng = range(xx,na.rm=TRUE,finite=TRUE)
    plotBKG(x.rng,ylim,...)
    title(ylab=ifelse(log.scale,'log(BF)','BF'))
    baseline = ifelse(log.scale,0,1)
    segments(x.rng[1],baseline,x.rng[2],baseline)
    y.thres = ifelse(log.scale,0,1)
    rect(xx[-(n+1)], y.thres, xx[-1], yy[-(n+1)],
      col=ifelse(yy>y.thres,cols[1],cols[2]),border=NA)
  }
  if(attr(BF,'type') %in% "categorical"){
    if(nrow(BF)>10) BF = rankBF(BF,n=10)
    else BF = transform(BF,logBF=log(BF))
    mp = barplot(BF$logBF,names.arg=BF$value,plot=FALSE)
    xadd = diff(mp)[1]/2 # .2
    plotBKG(c(min(mp)-xadd,max(mp)+xadd),ylim,xaxt='n',yaxt='n',x.major=NULL,x.minor=mp)
    barplot(BF$logBF,names.arg=BF$value,
            cex.names=.9,
            col=ifelse(BF$logBF>0,cols[1],cols[2]),
            col.axis='grey50',cex.axis=.8,tcl=-.25,
            ylim=ylim, las=1, add=TRUE)
    title(ylab=ifelse(log.scale,'log(BF)','BF'))
  }
  if(show.legend){
    leg.names = c(expression(paste('Favors ',H[L])),expression(paste('Favors ',H[U])))
    legend('topright',leg.names,#col=1:2,lwd=2,
       fill = cols,
       bty = 'n',border = NA)
  }
}

## rankBF
##==============================================================================
##  Orders Category levels according to absolute value log bayes factor
##==============================================================================
rankBF <- function(BF,n=10,thres=NULL){
  n = min(n,nrow(BF))
  logBF = log(BF$BF)  # log-Bayes factor
  ord = order(abs(logBF),decreasing=TRUE)
  if(is.null(thres)) ind = ord[1:n]
  else ind = ord[abs(logBF[ord])>=thres]
  data.frame(BF[ind,],logBF=logBF[ind])
}

## color
##==============================================================================
#'  Creates transparent colors
#'  @param col Color that \R recognizes (names or number)
#'  @param alpha transparency value
#'  @seealso \code{\link{col2rgb}}
#'  @keywords internal
##==============================================================================
color <- function(col,alpha=NULL){
  col = col2rgb(col,alpha=TRUE)/255
  if(!is.null(alpha) &  alpha>=0 & alpha <=1){
    col[4] = alpha
  }
  rgb.col = rgb(col[1],col[2],col[3],alpha=col[4])
return(rgb.col)
}





