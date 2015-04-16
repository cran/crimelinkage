
################################################################################
##  compareCrimes.R
##  Various functions for creating evidence variables
################################################################################

#-- Required R Packages
#require(geosphere) # for compareSpatial() if lat/long coordinates used

##  compareCrimes
##==============================================================================
#'  Creates evidence variables by calculating \sQuote{distance} between crime pairs
#'
#'  Calculates spatial and temporal distance, difference in categorical, and
#'   absolute value of numerical crime variables
##  Inputs:
#'  @param Pairs (n x 2) matrix of crimeIDs
#'  @param crimedata data.frame of crime incident data. \code{crimedata} must have
#'    a column named: \code{crimeID}. Other column names must correspond to what 
#'    is given in \code{varnames} list.
#'  @param varnames list of column names corresponding to:
#'  \itemize{
#'    \item spatial - X,Y coordinates (in long,lat or X,Y)
#'    \item temporal - DT.FROM, DT.TO
#'    \item categorical - (optional) any categorical variables
#'    \item numerical - (optional) any numeric variables
#'  }
#'  @param binary (logical) match/no match or all combinations for categorical
#'    class
#'  @param longlat (logical) are spatial coordinates in (long,lat)?
#'  @param show.pb (logical) show the progress bar
#'  @param seed seed for random number generation 
##  Outputs:
#'  @return data.frame of various proximity measures between the two crimes
#'  @examples
#'  data(crimes)
#'  
#'  varnames = list(
#'    spatial = c("X", "Y"),
#'    temporal = c("DT.FROM","DT.TO"),
#'    categorical = c("MO1",  "MO2", "MO3"))
#'  pairs = t(combn(as.character(crimes$crimeID[1:4]),m=2))
#'
#'  compareCrimes(pairs,crimes,varnames,binary=TRUE)    
#'  
#'  @references
#'  Porter, M. D. (2014). A Statistical Approach to Crime Linkage.
#'    \emph{arXiv preprint arXiv:1410.2285.}.
#'  \url{http://arxiv.org/abs/1410.2285}
#'  @export
##==============================================================================
compareCrimes <- function(Pairs,crimedata,varnames,binary=TRUE,longlat=FALSE,
                          show.pb=FALSE,seed=NULL){
    if(!is.null(seed)) set.seed(seed)
#  i1 = Pairs[,1]
#  i2 = Pairs[,2] 
  i1 = match(Pairs[,1],crimedata$crimeID)
  i2 = match(Pairs[,2],crimedata$crimeID)

  spatial=temporal=categorical=numerical=NA
  if(!is.null(varnames$spatial)){
    spatial = compareSpatial(crimedata[i1,varnames$spatial],
                             crimedata[i2,varnames$spatial],
                              longlat=longlat) }

  if(!is.null(varnames$temporal)){
    temporal = compareTemporal(crimedata[i1,varnames$temporal],
                               crimedata[i2,varnames$temporal],show.pb=show.pb) }

  catNames = varnames$categorical
  if(!is.null(catNames)){
    categorical = data.frame(matrix(nrow=length(i1),ncol=length(catNames),
                                    dimnames=list(NULL,catNames)))
    for(cat in catNames){
      categorical[,cat] = compareCategorical(crimedata[i1,cat],crimedata[i2,cat],
                                             binary=binary)
    } }

  numNames = varnames$numerical
  if(!is.null(numNames)){
    numerical = data.frame(matrix(nrow=length(i1),ncol=length(numNames),
                                  dimnames=list(NULL,numNames)))
    for(num in numNames){
      numerical[,num] = compareNumeric(crimedata[i1,num],crimedata[i2,num])
    } }
  Dist = data.frame(spatial,temporal,categorical,numerical)
  Dist[sapply(Dist,function(x) all(is.na(x)))] <- list(NULL)  
  E = cbind(Pairs,Dist)  # combine crime pair info with evidence variables
  rownames(E) <- NULL  # remove rownames
return(E)
}


##  compareTemporal
##==============================================================================
#'  Make temporal evidence variable from (possibly uncertain) temporal info
#'
#'  Calculates the temporal distance between crimes
##  Inputs:
#'  @param DT1 (n x 2) data.frame of (DT.FROM,DT.TO) for the crimes
#'  @param DT2 (n x 2) data.frame of (DT.FROM,DT.TO) for the crimes
#'  @param niters number of iterations for Monte Carlo expectation
#'  @param show.pb (logical) show the progress bar
##  Outputs:
#'  @return data.frame of expected absolute differences:
#'    \itemize{
#'      \item temporal - overall difference (in days)  [0,max]
#'      \item tod - time of day difference (in hours)  [0,12]
#'      \item dow - fractional day of week difference (in days) [0,3.5]
#'    }
#'  @details Uses Monte Carlo expected value - aoristic like analysis
#'  @keywords internal
##==============================================================================
compareTemporal <- function(DT1,DT2,niters=2000,show.pb=TRUE){
  #-- make a DateTime object for further processing
  makeDTObj <- function(DT,origen){
    A = data.frame(day=as.numeric(DT-origen,units='days'),
                   tod=with(as.POSIXlt(DT), hour + min/60 + sec/3600),
                   dow=match(weekdays(DT,abbreviate=TRUE),
                             c('Sun','Mon','Tue','Wed','Thu','Fri','Sat')))
    A$dow.f = A$dow+A$tod/24 - 1 # set between [0,7]
    return(A)
  }

  #-- time of day difference
  # tod - [0,24]
  diff.tod <- function(tod1,tod2)  diff.circ(tod1,tod2,mod=24)

  #-- day of week (fractional) difference
  # dow - [0,7]
  diff.dow <- function(dow1,dow2)  diff.circ(dow1,dow2,mod=7)

  #-- difference between locations on a circle.
  # mod is circle circumference
  diff.circ <- function(a,b,mod){
    # mod=7: day of week (fractional) difference, dow - [0,7]
    # mod=24: time of day difference, tod - [0,24]
    a = a %% mod       # ensure between [0,mod]
    b = b %% mod       # ensure between [0,mod]
    circ.diff = a - b
    pmin(circ.diff %% mod ,-circ.diff %% mod)      # minimum difference
  }

  #-- Format DateTime Objects
  origen = min(DT1$DT.FROM,DT2$DT.FROM)
  A0 = makeDTObj(DT1$DT.FROM,origen)
#  A0$interval.hr = DT1$TIMERANGE
  A0$interval.hr = as.numeric(abs(difftime(DT1$DT.TO,DT1$DT.FROM,units='hours')))
  B0 = makeDTObj(DT2$DT.FROM,origen)
#  B0$interval.hr = DT2$TIMERANGE
  B0$interval.hr = as.numeric(abs(difftime(DT2$DT.TO,DT2$DT.FROM,units='hours')))

  #-- Monte Carlo Expected Absolute Differences
  n = nrow(A0)
  t.diff = A0$day - B0$day
  X = matrix(0,n,3)
  iters = 0L
  if(show.pb) pb = txtProgressBar(style=3,max=niters)
  while(iters < niters){
    iters = iters + 1L
    h.A = runif(n,min=0,max=A0$interval.hr) # random hours from DT.FROM
    h.B = runif(n,min=0,max=B0$interval.hr)
    temporal = abs(t.diff + (h.A-h.B)/24)   # Absolute Difference in Days
#    temporal = t.diff + (h.A-h.B)/24        # Difference in Days
    tod.A = (A0$tod + h.A)
    tod.B = (B0$tod + h.B)
    tod = diff.tod(tod.A %% 24,tod.B %% 24) # abs Difference in time of day
    dow.A = (A0$dow.f +  tod.A / 24) %% 7
    dow.B = (B0$dow.f +  tod.B / 24) %% 7
    dow = diff.dow(dow.A,dow.B)             # abs Difference in (fractional) day of week
    X = X + cbind(temporal,tod,dow)
    if(show.pb) setTxtProgressBar(pb,iters)
  }
#  X = abs(X)
  if(show.pb) close(pb)
  return(data.frame(X/niters))  # Average of absolute differences
}


##  compareSpatial
##==============================================================================
#'  Make spatial evidence variables
#'
#'  Calculates spatial distance between crimes (in km)
##  Inputs:
#'  @param C1 (n x 2) matrix of coordinates for the crimes
#'  @param C2 (n x 2) matrix of coordinates for the crimes
#'  @param longlat (logical) if true, the the coordinates are in (Long,Lat), else
##      assume a suitable project where euclidean distance can be applied
##  Outputs:
#'  @return numeric vector of distances between the crimes (in km)
#'  @keywords internal
##==============================================================================
#library(geosphere)
compareSpatial <- function(C1,C2,longlat=FALSE){
  if(longlat)  d = geosphere::distHaversine(C1,C2)
  else         d = sqrt(rowSums((C1 - C2)^2))
  return(d/1000)            # distance in km
}

##  compareNumeric
##==============================================================================
#' Make evidence variables from numeric crime data
#'
#' Calculates absolute difference between crimes variables
##  Inputs:
#'  @param C1 length n numerical values of crime attributes
#'  @param C2 length n numerical values of crime attributes
##  Outputs:
#'  @return numeric vector of absolute differences
#'  @keywords internal
##==============================================================================
compareNumeric <- function(C1,C2) abs(C1 - C2)


##  compareCategorical
##==============================================================================
#'  Make evidence variables from categorical crime data
#'
#'  Compares categorical crime data to check if they match.
##  Inputs:
#'  @param C1 length n categorical values of crime attributes
#'  @param C2 length n categorical values of crime attributes
#'  @param levs the levels of all possible values
#'  @param  binary (logical) match/no match or all combinations
##  Outputs:
#'  @return if binary=TRUE: 1 for match, 0 for non-matches;
#'  if binary=FALSE: factor vector of merged values (in form of f1:f2)
#'  @keywords internal
##==============================================================================
compareCategorical <- function(C1,C2,levs,binary=FALSE){
  if(binary){
    match.na = is.na(C1) & is.na(C2)   # counts matching NA's as a match
    match.value = as.character(C1)==as.character(C2)
    B = ifelse(match.value | match.na,1,0)
    B[is.na(B)] = 0    # set rest of NAs to non-match
    return(factor(B,levels=c(0,1)))
  }
  C1 = as.factor(C1); C2 = as.factor(C2)
  if(missing(levs)) levs = sort(unique(levels(C1),levels(C2)))
  A = data.frame(C1=factor(C1,levels=levs), C2=factor(C2,levels=levs))
  flip = which(is.na(A[,1]) | unclass(A[,1]) > unclass(A[,2]))
  A[flip,] = A[flip,2:1]    # Sort rows
  B = paste(A[,1],A[,2],sep=':')        # Merge values
  B = factor(B,levels=catLevels(levs))  # Make into factor
  return(B)
}

##  catLevels
##==============================================================================
#'  Make levels for merging category predictors
##  Inputs:
#'  @param levs levels of a catagorical variable (factor)
##  Outputs:
#'  @return levels for a new categorical variable of form f1:f2
#'  @keywords internal
##==============================================================================
catLevels <- function(levs){
  levs = unique(c(levs,NA))  # Add NA if not already included
  nlevs = length(levs)
  a = NULL
  for(i in 1:nlevs){
    b = cbind(levs[i],levs[i:nlevs])
    a = rbind(a,b)
  }
  levs2 = paste(a[,1],a[,2],sep=':')
  return(levs2)
}

