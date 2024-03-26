## ao = \u00e5
## ae = \u00e4
## oe= \u00f6 
## mu = \u03bc
## sigma = \u03c3
## alpha = \u03b1
## beta = \u03b2
## in = \u2208

## THESE AVAILABLE FROM pkg:evd
## gumbcdf.m
## gumbinv.m
## gumbpdf.m
## gumbrnd.m
##AMEND LAB FOR FITTING TO INCLUDE ALSO evd FUNCTIONS
## gumbfit.m
## gumbplot.m


#' QQ-plot for data anpassad till standardfordelningar.
#' 
#' @param x En vektor med data.
#' @param dist Vilken fordelning som ska anpassas. Se \link[MASS:fitdistr]{fitdistr}.
#' @param ... Ytterligare parametrar till \link[stats:qqplot]{qqplot}. De
#'     viktigaste ar: 
#' \describe{
#'   \item{main}{Text till titel}
#'   \item{xlab,ylab}{Text till x- och y-axlar}
#' }
#' @return Ingenting
#' @import graphics stats evd
#' @importFrom MASS fitdistr
#' @export
#' @examples
#' par(mfrow=c(2,2))
#' x <- rnorm(500,10,2)
#' qqdist(x)
#'
#' x <- rgamma(500, 2, 2)
#' qqdist(x)
#' qqdist(x, 'gamma')
#'
#' x <- rpois(500, 0.5)
#' qqdist(x, 'Poisson')
qqdist <- function(x, dist=c('normal', 'cauchy', 'exponential', 'gamma',
                             'geometric', 'gumbel', 'gev', 'lognormal', 'logistic',
                             'negative binomial', 'Poisson', 'weibull'),
                   ...){
    ##check for distributions, default to normal
    dist <- match.arg(dist)
    ##check additional arguments
    args <- list(...)
    ##title
    if(is.null(args$main)){ args$main <- paste(dist,'Q-Q plot')}
    ##xlab
    if(is.null(args$xlab)){ args$xlab <- 'Theoretical Quantiles' }
    ##xlab
    if(is.null(args$ylab)){ args$ylab <- 'Sample Quantiles' }

    if( dist=='gumbel'){
        ##Fit distribution (drop warnings)
        suppressWarnings( fit <- fgev(x, shape=0) )
    }else if( dist=='gev'){
        ##Fit distribution (drop warnings)
        suppressWarnings( fit <- fgev(x) )
    }else{
        ##Fit distribution (drop warnings)
        suppressWarnings( fit <- fitdistr(x, densfun=dist) )
    }

    ##compute quantiles
    x.seq <- ppoints(length(x))
    ##compute quantile functions
    if(dist=='normal'){
        fnc <- function(p){qnorm(p, mean=fit$estimate[1], sd=fit$estimate[2])}
    }else if(dist=='cauchy'){
        fnc <- function(p){qcauchy(p, location=fit$estimate[1],
                                   scale=fit$estimate[2])}
    }else if(dist=='exponential'){
        fnc <- function(p){qexp(p, rate=fit$estimate[1])}
    }else if(dist=='gamma'){
        fnc <- function(p){qgamma(p, shape=fit$estimate[1], rate=fit$estimate[2])}
    }else if(dist=='geometric'){
        fnc <- function(p){qgeom(p, prob=fit$estimate[1])}
    }else if(dist=='gumbel'){
        fnc <- function(p){qgumbel(p, loc=fit$estimate[1],
                                   scale=fit$estimate[2])}
    }else if(dist=='gev'){
        fnc <- function(p){qgev(p, loc=fit$estimate[1], scale=fit$estimate[2],
                                shape=fit$estimate[3])}
    }else if(dist=='lognormal'){
        fnc <- function(p){qlnorm(p, meanlog=fit$estimate[1],
                                  sdlog=fit$estimate[2])}
    }else if(dist=='logistic'){
        fnc <- function(p){qlogis(p, location=fit$estimate[1],
                                  scale=fit$estimate[2])}
    }else if(dist=='negative binomial'){
        fnc <- function(p){qnbinom(p, size=fit$estimate[1], mu=fit$estimate[2])}
    }else if(dist=='Poisson'){
        fnc <- function(p){qpois(p, lambda=fit$estimate[1])}
    }else if(dist=='weibull'){
        fnc <- function(p){qweibull(p, shape=fit$estimate[1], scale=fit$estimate[2])}
    }
    y <- fnc(x.seq)

    ##Plot the distribution (x=theoretical data; y=sample)
    do.call( qqplot, c(list(x=y, y=x), args) )
    ##add a reference line
    qqline(x, distribution = fnc)

    ##no return
    invisible(NULL)
}

#' Histrogram med anpassad tathetsfunktion
#' 
#' @param x En vektor med data.
#' @param dist Vilken fordelning som ska anpassas. Se \link[MASS:fitdistr]{fitdistr}.
#' @param ... Ytterligare parametrar till \link[graphics:hist]{hist}. De viktigaste ar:
#' \describe{
#'   \item{breaks}{Antalet staplar i histogramet}
#'   \item{main}{Text till titel}
#'   \item{xlab,ylab}{Text till x- och y-axlar}
#'   \item{freq}{Rita histogramet i frekvens eller sannolikhet}
#' }
#' @return De skattade parametrarna (se \link[MASS:fitdistr]{fitdistr}).
#' @import graphics stats
#' @importFrom MASS fitdistr
#' @export
#' @examples
#' par(mfrow=c(2,2))
#' x <- rnorm(100,10,2)
#' print( histfit(x) )
#' 
#' histfit(x, freq=TRUE, breaks=25)
#' histfit(x[x>0], 'gamma')
#' 
#' x <- rgamma(100, 2, 2)
#' print( histfit(x, 'gamma') )
histfit <- function(x, dist=c('normal', 'cauchy', 'exponential', 'gamma',
                              'geometric', 'gumbel', 'gev', 'lognormal', 'logistic',
                              'negative binomial', 'Poisson', 'weibull'),
                    ...){
    ##check for distributions, default to normal
    dist <- match.arg(dist)

    ##extract additional args and set defaults
    args <- list(...)
    ##ensure that we plot and do it with 'tight' axes
    args$plot <- TRUE
    args$yaxs <- 'i'
    ##xlim
    if(is.null(args$xlim)){
        args$xlim <- range(x)
        args$xlim <- args$xlim + 0.15*diff(args$xlim)*c(-1,1)
    }
    ##title
    if(is.null(args$main)){ args$main <- paste('Fit of',dist,'to data')}
    ##xlab
    if(is.null(args$xlab)){ args$xlab <- 'Data' }
    ##colour
    if(is.null(args$col)){ args$col <- 'lightblue' }
    ##check if plot is in frequency or probability scale
    is.freq <- TRUE
    if(!is.null(args$freq)){ is.freq <- args$freq }
    if(!is.null(args$probability)){ is.freq <- !args$probability }

    if( dist=='gumbel'){
        ##Fit distribution (drop warnings)
        suppressWarnings( fit <- fgev(x, shape=0) )
    }else if( dist=='gev'){
        ##Fit distribution (drop warnings)
        suppressWarnings( fit <- fgev(x) )
    }else{
        ##Fit distribution (drop warnings)
        suppressWarnings( fit <- fitdistr(x, densfun=dist) )
    }
    
    ##Strictly possitive data -> Set xlim[1] = max(xlim[1],0)
    if(dist %in% c('exponential', 'gamma', 'geometric', 'lognormal', 'logistic',
        'negative binomial', 'Poisson', 'weibull')){
        args$xlim[1] <- max(args$xlim[1],0)
    }
    ##compute density
    x.seq <- seq(args$xlim[1], args$xlim[2], length.out=1e3)
    ##compute relevant distribution
    if(dist=='normal'){
        y <- dnorm(x.seq, mean=fit$estimate[1], sd=fit$estimate[2])
    }else if(dist=='cauchy'){
        y <- dcauchy(x.seq, location=fit$estimate[1], scale=fit$estimate[2])
    }else if(dist=='exponential'){
        y <- dexp(x.seq, rate=fit$estimate[1])
    }else if(dist=='gamma'){
        y <- dgamma(x.seq, shape=fit$estimate[1], rate=fit$estimate[2])
    }else if(dist=='geometric'){
        y <- dgeom(x.seq, prob=fit$estimate[1])
    }else if(dist=='gumbel'){
        y <- dgumbel(x.seq, loc=fit$estimate[1], scale=fit$estimate[2])
    }else if(dist=='gev'){
        y <- dgev(x.seq, loc=fit$estimate[1], scale=fit$estimate[2],
                  shape=fit$estimate[3])
    }else if(dist=='lognormal'){
        y <- dlnorm(x.seq, meanlog=fit$estimate[1], sdlog=fit$estimate[2])
    }else if(dist=='logistic'){
        y <- dlogis(x.seq, location=fit$estimate[1], scale=fit$estimate[2])
    }else if(dist=='negative binomial'){
        y <- dnbinom(x.seq, size=fit$estimate[1], mu=fit$estimate[2])
    }else if(dist=='Poisson'){
        y <- dpois(x.seq, lambda=fit$estimate[1])
    }else if(dist=='weibull'){
        y <- dweibull(x.seq, shape=fit$estimate[1], scale=fit$estimate[2])
    }
    ##Plot histogram
    ##as a 2-step procedure first we compute histogram to get y-limits.
    args$plot <- FALSE ##don't plot
    args$warn.unused <- FALSE ##don't warn
    H <- do.call(hist, c(list(x=x),args))
    ##making sure that we check for scaling (freq=TRUE/FALSE)
    if( is.freq ){
        y <- y * fit$n * (H$mids[2]-H$mids[1]) 
        args$ylim <- c(0, 1.05*max(H$counts,y))
    }else{
        args$ylim <- c(0, 1.05*max(H$density,y))
    }
    ##plot the histogram with good ylim-scaling
    args$plot <- TRUE ##plot!
    do.call(hist, c(list(x=x),args))
    ##add the normal box around the plot
    box()
    ##and add line of the fitted distribution
    lines(x.seq, y, col='red', lwd=2)
    
    ##Distributions available in matlab but missing in MASS
    ## 'beta'
    ## 'birnbaumsaunders'                 Birnbaum-Saunders
    ## 'extreme value' or 'ev'            Extreme value
    ## 'generalized extreme value' 'gev'  Generalized extreme value
    ## 'generalized pareto' or 'gp'       Generalized Pareto (threshold 0)
    ## 'inverse gaussian'                 Inverse Gaussian
    ## 'loglogistic'                      Log logistic
    ## 'nakagami'                         Nakagami
    ## 'rayleigh'                         Rayleigh
    ## 'rician'                           Rician

    ##return fitting results invisible
    invisible(fit)
}

#' Illustrerar berakning av sannolikheter for normalfordelningar.
#' 
#' @param bnd Granser for sannolikheten P(bnd[1]<X<bnd[2]) beraknas. For
#'     ensidiga sannolikheter anvand Inf eller -Inf.
#' @param mu,sigma Parametrar i normalfordelningen
#' @return Ingenting
#' @import graphics stats
#' @export
#' @examples
#' par(mfrow=c(1,3))
#' plotnorm(c(-Inf,0))  ## P(X<0)
#' plotnorm(c(1,Inf))   ## P(X>1)
#' plotnorm(c(0,1))     ## P(0<X<1)
plotnorm <- function(bnd, mu=0, sigma=1){
    if( bnd[1]>bnd[2] ){
        stop('Values in bnd must be ordered: bnd[1]<=bnd[2]')
    }
    ##Compute probability
    p <- pnorm(bnd[2],mu,sigma)-pnorm(bnd[1],mu,sigma)
    ##setup edges of vector for plotting
    x <- c(mu-4*sigma, mu+4*sigma)
    if( is.finite(bnd[1]) ){
        x[1] <- min(x[1], bnd[1]-sigma, bnd[2]+sigma)
    }else{
        x[1] <- min(x[1], bnd[2]-sigma)
        bnd[1] <- x[1]
    }
    if( is.finite(bnd[2]) ){
        x[2] <- max(x[2], bnd[1]-sigma, bnd[2]+sigma)
    }else{
        x[2] <- max(x[2], bnd[1]+sigma)
        bnd[2] <- x[2]
    }
    ##construct vectors of x
    x <- seq(x[1], x[2], length.out=1000)
    x.bnd <- seq(bnd[1], bnd[2], length.out=1000)
    ##plot density
    plot(x, dnorm(x,mu,sigma), type='l', main=sprintf('Probability %.3g', p),
         xlab='Critical value', ylab='Density')
    ##add polygon
    polygon(c(x.bnd,max(x.bnd),min(x.bnd)), c(dnorm(x.bnd,mu,sigma),0,0),
            col='blue', border=NA) 
  ##no return
  invisible(NULL)
}


#' Illustrerar kritiskt omrade for hypotestest
#' 
#' Illustrerar kritiskt omrade (och styrka) for hypotestest av mu under
#' antagande om n observationer fran en normalfordelning, N(mu_0, sigma).
#' @param  mu0,sigma Parametrar i normalfordelningen under H0.
#' @param  n Antal observationer.
#' @param  alfa Signifikansniva for testet
#' @param  riktning Vilken mothypotes att använda. H0: mu=mu0 mot
#' \describe{
#'   \item{'<'}{H1: mu < mu0}
#'   \item{'>'}{H1: mu > mu0}
#'   \item{'!='}{H1: mu != mu0}
#' }
#' @param  mu.sant Varda pa mu att berakna styrkan for.
#' @return Ingenting
#' @import graphics stats
#' @export
#' @examples
#' hypotes(1, 10, 5, .05, '<')
#' hypotes(1, 10, 5, .05, '<', 4.52)
hypotes <- function(sigma, n, mu0, alfa, riktning=c('<','>','!='), mu.sant){
  riktning <- match.arg(riktning)

  ##compute standard deviation of estimate
  s <-  sigma/sqrt(n)
  ##and compute maximum of density function
  f0 <- dnorm(x=0, mean=0, sd=s)
  ##setup figures
  if(!missing(mu.sant)){
    par(mfrow=c(2,1), mar=c(3.1,3.1,3.1,0.5), mgp=c(1.5,0.5,0),
        xaxs='i', yaxs='i')
  }else{
    par(mfrow=c(1,1), mar=c(3.1,3.1,3.1,0.5), mgp=c(1.5,0.5,0),
        xaxs='i', yaxs='i')
  }

  ##plot the normal density
  x <- plotNormalDens(mu0, s)
  ##and plot the area (and compute k-values)
  k <- plotCritArea(mu0, s, alfa, riktning)
  ##construct x-labels
  if(riktning=='!='){
    xlab <- sprintf('k1 = %4.2f    k2 = %4.2f', k$k.un, k$k.ov)
  }else{
    xlab <- sprintf('k = %4.2f', k[[ which(sapply(k,is.finite)) ]])
  }

  ##figurtext till figure 1
  title(xlab=xlab, ylab='', main =
    sprintf('Kritiskt omr\u00e5de, H0: \u03bc=%4.1f, H1: \u03bc%s%4.1f',
            mu0, riktning, mu0))
  text(mu0+2*s, 1.3*f0, sprintf('H0: \u03bc = %4.1f', mu0), pos=4);
  text(mu0+2*s, 1.1*f0, sprintf('H1: \u03bc %s %4.1f', riktning, mu0), pos=4);
  text(mu0+2*s, 0.9*f0, sprintf('n = %4.0f',n), pos=4);
  text(mu0+2*s, 0.7*f0, sprintf('\u03c3 = %4.1f', sigma), pos=4);
  text(mu0+2*s, 0.5*f0, sprintf('\u03b1 = %4.3f', alfa), pos=4);

  ##figure 2
  if(!missing(mu.sant)){
    ##plot the normal density
    x <- plotNormalDens(mu0, s)
    ##plot the area(s)
    k <- plotCritArea(mu0, s, alfa, riktning)
    P.typ2 <- plotBetaArea(mu0, mu.sant, s, k)
    ##titels and labels
    title(xlab=sprintf('\u03b1=%4.3f (r\u00f6d), \u03b2=%4.3f (bl\u00e5); S(%4.1f)=1-\u03b2=%4.3f',
                      alfa, P.typ2, mu.sant, 1-P.typ2),
         ylab='',
         main=sprintf('Sannolikheter f\u00f6r fel av typ 1 och typ 2; Styrka d\u00e5 \u03bc=%4.1f',
                      mu.sant))
  }
  ##no return
  invisible(NULL)
}##hypotes <- function

#' Illustrerar skillnaden mellan känt och skattat sigma
#'
#' Illustrerar skillnaden mellan kant och skattat sigma da intervall gors for mu
#' i N(mu sigma). Givet ett varde pa sigma och en signifikansniva, alfa
#' simulerar funktionen normalfordelade slumptal och skattar sigma baserat pa
#' olika antal observationer (typiskt n=2 till n=100). En tabell med skattade
#' sigma, lambda- och t-kvantiler samt resulterande halv-bredd for tvasidiga
#' intervall skrivs ut.
#' @param sigma Standardavvikelse hos observationerna
#' @param alfa Signifikans niva for intervallen
#' @return Den resulterande tabellen.
#' @import stats
#' @export
#' @examples
#'   kvantilintervall(2,0.05)
kvantilintervall <- function(sigma, alfa){
  ##sigma should be positive
  stopifnot(sigma>0)

  n <- c(2:5, seq(10,30, by=5), 40, 50, 75, 100)
  res <- data.frame(n=n, sigma=sigma, s=NA, lambda=qnorm(1-alfa/2),
                    t=qt(1-alfa/2,n-1), bredd=2*qnorm(1-alfa/2)*sigma/sqrt(n),
                    bredd2=NA)
  for(i in 1:length(n)){
    varden <- rnorm(n[i], mean=0, sd=sigma)
    res$s[i] <- sd(varden)
  }
  res$bredd2 <- 2 * res$t * res$s / sqrt(res$n)
  names(res)[6] <- 'halv bredd sigma kand'
  names(res)[7] <- 'halv bredd sigma skattad'
  ##print results
  tmp <- format(unname(res),digits=3)
  names.1 <- substr(names(res),1,10)
  names.2 <- substr(names(res),12,max(nchar(names(res))))
  N <- pmax(nchar(tmp[1,])+1,nchar(names.1)+1,nchar(names.2)+1)
  fmt <- list(paste(sprintf("%%%.0fs",N),collapse=""))
  str <- paste(do.call(sprintf,c(fmt,as.list(names.1))),
               do.call(sprintf,c(fmt,as.list(names.2))),
               sep="\n")
  for(i in 1:dim(tmp)[1]){
  str <- paste(str, do.call(sprintf,c(fmt,as.list(tmp[i,]))),
               sep="\n")
  }
  cat( paste(str, "\n", sep="") )
  ##and return invisible
  invisible(res)
}##kvantilintervall <- function

#' Illustrerar direktmetoden vid test av mu i N(mu,sigma)
#'
#' Illustrerar direktmetoden och jamfor med kritiskt omrade for hypotestest av
#' mu under antagande om n observationer fran en normalfordelning, N(mu_0,
#' sigma).
#' @param  mu0,sigma Parametrar i normalfordelningen under H0.
#' @param  n Antal observationer.
#' @param  alfa Signifikansniva for testet
#' @param  riktning Vilken mothypotes att använda. H0: mu=mu0 mot
#' \describe{
#'   \item{'<'}{H1: mu < mu0}
#'   \item{'>'}{H1: mu > mu0}
#'   \item{'!='}{H1: mu != mu0}
#' }
#' @param  medel Observerat medelvarde av stickprovet.
#' @return Ingenting
#' @import graphics stats 
#' @export
#' @examples
#' Pvarde(2, 15, 5, .05, '<', 4)
#' 
#' x <- rnorm(mean=5, sd=2, n=15)
#' Pvarde(2, length(x), 5, .05, '<', mean(x))
Pvarde <- function(sigma, n, mu0, alfa, riktning=c('<','>','!='), medel){
  riktning <- match.arg(riktning)
  
  ##compute standard deviation of estimate
  s <-  sigma/sqrt(n)
  ##and compute maximum of density function
  f0 <- dnorm(x=0, mean=0, sd=s)
  ##setup figures
  par(mfrow=c(2,1), mar=c(3.1,3.1,3.1,0.5), mgp=c(1.5,0.5,0),
      xaxs='i', yaxs='i')

  ##plot the normal density
  x <- plotNormalDens(mu0, s)
  ##and plot the area (and compute k-values)
  k <- plotCritArea(mu0, s, alfa, riktning)
  lines(c(medel,medel), c(0,dnorm(medel,mu0,s)), col='black')
  ##construct x-labels
  if(riktning=='!='){
    xlab <- sprintf('k1 = %4.2f    k2 = %4.2f', k$k.un, k$k.ov)
  }else{
    xlab <- sprintf('k = %4.2f', k[[ which(sapply(k,is.finite)) ]])
  }

  ##figurtext till figure 1
  title(main ='Test med fixt \u03b1-v\u00e4rde',
        xlab=xlab, ylab='')
  ##Beskrivning av testet
  text(mu0-3.8*s, 1.3*f0, sprintf('H0: \u03bc = %4.1f', mu0), pos=4);
  text(mu0-3.8*s, 1.1*f0, sprintf('H1: \u03bc %s %4.1f', riktning, mu0), pos=4);
  text(mu0-3.8*s, 0.9*f0, sprintf('n = %4.0f',n), pos=4);
  text(mu0-3.8*s, 0.7*f0, sprintf('\u03c3 = %4.1f', sigma), pos=4);
  text(mu0-3.8*s, 0.5*f0, sprintf('\u03b1 = %4.3f', alfa), pos=4);
  ##text for medel
  text(mu0+2*s, 1.3*f0, sprintf('Medel = %4.3f',medel), pos=4)
  text(mu0+2.5*s, 0.9*f0,'\u21d4', pos=4)
  ##testa hur medel och k-varde forhaller sig
  HO.forkast <- medel<k$k.un || k$k.ov<medel
  if(riktning=='>'){
      if(HO.forkast){
          slutsats.medel <- 'Medel > k'
      }else{
          slutsats.medel <- 'Medel < k'
      }
  }else if(riktning=='<'){
      if(HO.forkast){
          slutsats.medel <- 'Medel < k'
      }else{
          slutsats.medel <- 'Medel > k'
      }
  }else{ #if(riktning=='!=')
      if(HO.forkast && medel<k$k.un){
          slutsats.medel <- 'Medel < k1'
      }else if(HO.forkast && medel>k$k.ov){
          slutsats.medel <- 'Medel > k2'
      }else{
          slutsats.medel <- 'k1 < Medel < k2'
      }
  }
  ##text for slutsats
  text(mu0+2*s, 1.1*f0, slutsats.medel, pos=4)
  if(HO.forkast){
      text(mu0+2*s, 0.7*f0, 'H0 f\u00f6rkastas', pos=4)
  }else{
      text(mu0+2*s, 0.7*f0, 'H0 f\u00f6rkastas ej', pos=4)
  }

  ##Second figure
  ##plot the normal density
  x <- plotNormalDens(mu0, s)
  ##area for sannolikhets berakning.
  if(riktning=='>'){
      p <- 1 - pnorm(medel, mu0, s)
  }else if(riktning=='<'){
      p <- pnorm(medel, mu0, s)
  }else{ ##if(riktning=='!=')
      p <- pnorm(medel, mu0, s)
      p <- 2*min(p, 1-p)
  }
  plotCritArea(mu0, s, alfa=p, riktning, col='darkred')

  title( main='Test med direktmetoden')
  ##text for medel och sannolikhet
  text(mu0+2*s, 1.3*f0, sprintf('Medel = %4.3f',medel), pos=4)
  text(mu0+2.5*s, 0.7*f0,'\u21d4', pos=4)
  if(riktning=='!='){
      text(mu0+2*s, 1.1*f0, sprintf('P = 2 * %4.3f', p/2), pos=4)
  }else{
      text(mu0+2*s, 1.1*f0, sprintf('P = %4.3f', p), pos=4)
  }
  ##text for slutsats
  if(HO.forkast){
      text(mu0+2*s, 0.9*f0, 'P < \u03b1', pos=4)
      text(mu0+2*s, 0.5*f0, 'H0 f\u00f6rkastas', pos=4)
  }else{
      text(mu0+2*s, 0.9*f0, 'P > \u03b1', pos=4)
      text(mu0+2*s, 0.5*f0, 'H0 f\u00f6rkastas ej', pos=4)
  }
}##Pvarde <- function(sigma, n, mu0, alfa, riktning, medel){

#'Illustrerar mu och sigma2-skattning samt konfidensintervall
#'
#' Ritar histogram for mu- och sigma^2-skattning samt illustrerar
#' konfidensintervall for mu. Skattningarna baseras pa n1 respektive n2
#' observationer fran en normalfordelning, N(mu, sigma).
#' 
#' @param mu,sigma Parametrar i normalfordelningen
#' @param n1,n2 Antal observationer att simulera i de tva stickproven. Gor det
#'     mojligt att jamfora hur skattningarna beter sig med olika antal
#'     observationer.
#' @param alternativ Textstrang som talar om vad som ska illustreras:
#' \describe{
#'   \item{'muskatt'}{Histogram for skattningar av mu}
#'   \item{'sigmaskatt'}{Histogram for skattningar av sigma^2}
#'   \item{'konfint'}{Illustrera konfidensintervall for mu}
#'   \item{'alla'}{Illustrera bade skattningar och intervall}
#' }
#' @return Ingenting
#' @import stats graphics
#' @importFrom grDevices dev.new
#' @export
#' @examples
#' skattningar(35, 2, 5, 25, 'muskatt')
skattningar <- function(mu, sigma, n1, n2, alternativ = c('alla','muskatt','sigmaskatt','konfint')){
    alternativ <- match.arg(alternativ)

    ##Antal simuleringar som gors
    n <- 1000;
    ##simulera tva sample
    x <- matrix(rnorm(n1*n, mu, sigma), n1, n)
    y <- matrix(rnorm(n2*n, mu, sigma), n2, n)

    ##mu skattningar
    if(alternativ %in% c('muskatt','alla')){
        ##intervallens bredd, for att satta axlar.
        width <- 3.5*sigma/sqrt(min(n1,n2))
        ##figure
        par(mfrow=c(2,1), mar=c(3.1,3.1,3.1,0.5), mgp=c(1.5,0.5,0),
            xaxs='i', yaxs='i')
        ##histograms
        hist(colMeans(x), 25, col='lightblue',
             xlim=c(mu-width, mu+width), ylim=c(0,200),
             main=sprintf('\u03bc-skattningens f\u00f6rdelning, n1=%2.0f', n1),
             ylab='', xlab='')
        abline(v=mu, col='red')

        hist(colMeans(y), 25, col='lightblue',
             xlim=c(mu-width, mu+width), ylim=c(0,200),
             main=sprintf('\u03bc-skattningens f\u00f6rdelning, n2=%2.0f', n2),
             ylab='', xlab='')
        abline(v=mu, col='red')
    }

    ##sigma2 skattningar
    if(alternativ %in% c('sigmaskatt','alla')){
        if(alternativ=='alla'){ dev.new() }
        ##intervallens bredd, for att satta axlar.
        width <- sigma^2 * max(qchisq(0.9995, n1-1)/(n1-1),
                               qchisq(0.9995, n2-1)/(n2-1))
        ##figure
        par(mfrow=c(2,1), mar=c(3.1,3.1,3.1,0.5), mgp=c(1.5,0.5,0),
            xaxs='i', yaxs='i')
        ##histograms
        hist(apply(x,2,var), 25, col='lightblue',
             xlim=c(0,width), ylim=c(0,200),
             main=sprintf('\u03c3^2-skattningens f\u00f6rdelning, n1=%2.0f', n1),
             ylab='', xlab='')
        abline(v=sigma^2, col='red')

        hist(apply(y,2,var), 25, col='lightblue',
             xlim=c(0,width), ylim=c(0,200),
             main=sprintf('\u03c3^2-skattningens f\u00f6rdelning, n2=%2.0f', n2),
             ylab='', xlab='')
        abline(v=sigma^2, col='red')
    }

    ##Simulering av konfidensintervall
    if(alternativ %in% c('konfint','alla')){
        if(alternativ=='alla'){ dev.new() }
        ##  berakna de 1000 konfidens intervallen
        CI.x <- cbind(colMeans(x)-qnorm(0.975)*sigma/sqrt(n1),
                      colMeans(x)+qnorm(0.975)*sigma/sqrt(n1))
        CI.y <- cbind(colMeans(y)-qnorm(0.975)*sigma/sqrt(n2),
                      colMeans(y)+qnorm(0.975)*sigma/sqrt(n2))
        I.x <- CI.x[,1]<mu & mu<CI.x[,2]
        I.y <- CI.y[,1]<mu & mu<CI.y[,2]
        ##compute intervall widths for plotting
        width <- 1.2*max( abs(c(CI.x[1:100,],CI.y[1:100,])-mu) )
        ##use only the first 100 for plotting
        ## IndGood = find(I_x(1:100));
        ## IndBad = find(~I_x(1:100));
        
        ##figure
        par(mfrow=c(1,2), mar=c(4.1,3.1,3.1,0.5), mgp=c(2.5,0.5,0),
            xaxs='i', yaxs='i')
        ##plots
        plot(0, 0, xlim=c(mu-width, mu+width), ylim=c(0,101), type='n',
             ylab='', main=sprintf('100 intervall f\u00f6r \u03bc, n1=%2.0f\nkonfidensgrad = 0.95',n1),
             xlab=sprintf('Andel av 1000 intervall som\nmissar sant \u03bc: %1.3f',1-mean(I.x)))
        for(i in 1:100){
            lines(CI.x[i,], c(i,i), col=c('red', 'blue')[I.x[i]+1])
        }
        abline(v=mu, col='black')

        plot(0, 0, xlim=c(mu-width, mu+width), ylim=c(0,101), type='n',
             ylab='', main=sprintf('100 intervall f\u00f6r \u03bc, n2=%2.0f\nkonfidensgrad = 0.95',n2),
             xlab=sprintf('Andel av 1000 intervall som\nmissar sant \u03bc: %1.3f',1-mean(I.y)))
        for(i in 1:100){
            lines(CI.y[i,], c(i,i), col=c('red', 'blue')[I.y[i]+1])
        }
        abline(v=mu, col='black')
    }
}##skattningar <- function(mu, sigma, n1, n2, alternativ)
    
#' Illustrerar styrekfunktion for hypotestest
#'
#' Illustrerar kritiskt omrade och styrekfunktion for hypotestest av mu under
#' antagande om n observationer fran en normalfordelning, N(mu_0, sigma).
#' @param  mu0,sigma Parametrar i normalfordelningen under H0.
#' @param  n Antal observationer.
#' @param  alfa Signifikansniva for testet
#' @param  riktning Vilken mothypotes att använda. H0: mu=mu0 mot
#' \describe{
#'   \item{'<'}{H1: mu < mu0}
#'   \item{'>'}{H1: mu > mu0}
#'   \item{'!='}{H1: mu != mu0}
#' }
#' @param  mu.sant Varda pa mu att berakna styrkan for.
#' @return Ingenting
#' @import stats graphics
#' @importFrom grDevices dev.new
#' @export
#' @examples
#' styrka(1, 10, 5, .05, '<', 4.52)
styrka <- function(sigma, n, mu0, alfa, riktning, mu.sant){
  riktning <- match.arg(riktning, c('<', '>', '!='))

  ##compute standard deviation of estimate
  s <-  sigma/sqrt(n)
  ##and compute maximum of density function
  f0 <- dnorm(x=0, mean=0, sd=s)
  ##setup figures
  if(!missing(mu.sant)){
    par(mfrow=c(2,1), mar=c(3.1,3.1,3.1,0.5), mgp=c(1.5,0.5,0),
        xaxs='i', yaxs='i')
  }else{
    par(mfrow=c(1,1), mar=c(3.1,3.1,3.1,0.5), mgp=c(1.5,0.5,0),
        xaxs='i', yaxs='i')
  }

  ##Beräkna och plota styrke funktionen.
  if(riktning=='<'){
      k <- qnorm(alfa, mu0, s)
      x <- seq(k-3*s, mu0+s, length.out=1000)
      y <- pnorm(k, x, s)
  }else if(riktning == '>'){
      k <- qnorm(1-alfa, mu0, s)
      x <- seq(mu0-s, k+3*s, length.out=1000)
      y <- 1-pnorm(k, x, s)
  }else{ #if(riktning=='!=')
      k1 <- qnorm(1-alfa/2, mu0, s)
      k2 <- qnorm(alfa/2, mu0, s)
      x <- seq(k2-3.5*s, k1+3.5*s, length.out=1000)
      y <- 1-pnorm(k1,x,s) + pnorm(k2,x,s)
  }

  ##If mu.sant given, start by plotting rejection regions
  if(!missing(mu.sant)){
      ##plot the normal density
      plotNormalDens(mu0, s, xlim=range(x))
      ##and plot the area (and compute k-values)
      k <- plotCritArea(mu0, s, alfa, riktning)
      ##and plot the beta-area
      P.beta <- plotBetaArea(mu0, mu.sant, s, k)
      title( main=sprintf('Sannolikheter f\u00f6r fel av typ 1 och typ 2; Styrka d\u00e5 \u03bc=%4.1f',
                          mu.sant))
  }
  ##plot styrka
  plot(x, y, type='l', col='blue', lwd=2,
       xlim=range(x), ylim=c(-0.1,1.1),
       xlab='\u03bc', ylab='', main='')
  ##add grid to plot
  grid(col='grey')
  ##and add reference line
  if(missing(mu.sant)){
      title( main='S(\u03bc) = P(f\u00f6rkasta H0)' )
  }else{
      lines(c(mu.sant, mu.sant, x[1]), c(0,1-P.beta,1-P.beta),
            col='red', lty=1, lwd=2)
      title(main = sprintf('S(\u03bc) = P(f\u00f6rkasta H0); S(%2.1f) = %2.2f',
                           mu.sant, 1-P.beta))
  }

  ## textning
  y <- par('usr')[4]
  if(riktning=='>'){
      x <- mu0 + 0.5*s
  }else if(riktning=='<'){
      x <- mu0 - s
  }else{ #if(riktning=='!=')
      x <- mu0
  }
  text(x, 0.9*y, sprintf('H0: \u03bc = %4.1f', mu0), pos=4)
  text(x, 0.8*y, sprintf('H1: \u03bc %s %4.1f', riktning, mu0), pos=4)
  text(x, 0.7*y, sprintf('n = %4.0f',n), pos=4)
  text(x, 0.6*y, sprintf('\u03c3 = %4.1f', sigma), pos=4)
  text(x, 0.5*y, sprintf('\u03b1 = %4.3f', alfa), pos=4)
}##styrkefkn <- function(sigma, n, mu0, alfa, riktning, mu.sant)

################################
## Internatl helper functions ##
################################
##fun that plots density curve
plotNormalDens <- function(mu0, s, xlim=c(mu0-4*s,mu0+4*s)){
  ##and compute maximum of density function
  f0 <- dnorm(x=0, mean=0, sd=s);
  ##density function
  x <- seq(xlim[1], xlim[2], length.out=1000)
  y <- dnorm(x, mu0, s)
  ##plot density function and help lines
  plot(x, y, type='l', col='black', lwd=2,
       xlim=range(x), ylim=c(-0.1,1.5)*f0,
       xlab='', ylab='', main='')
  lines(c(mu0,mu0), c(0,f0), col='black', lty=2)
  abline(a=0, b=0, lty=1, col='black')
  ##return the x vector used
  return(x)
}##function <- plotNormalDens

##fun that adds critical regions to density plot
plotCritArea <- function(mu0, s, alfa, riktning, col='red', border='black'){
  ##compute k-value and P.value
  if(riktning=='>'){
    k.ov <- qnorm(1-alfa, mu0, s)
    k.un <- -Inf
  }else if(riktning=='<'){
    k.ov <- Inf
    k.un <- qnorm(alfa, mu0, s)
  }else{ #if(riktning=='!=')
    k.ov <- qnorm(1-alfa/2, mu0, s)
    k.un <- qnorm(alfa/2, mu0, s)
  }
  ##plot the area(s)
  x <- par('usr') ##first find x (and y) limits
  ##ovre
  if( is.finite(k.ov) ){
    x.ov <- seq(k.ov, x[2], length.out=100)
    polygon(c(k.ov,x.ov), c(0,dnorm(x.ov,mu0,s)), col=col, border=border)
  }
  ##undre
  if( is.finite(k.un) ){
    x.un <- seq(x[1], k.un, length.out=100)
    polygon(c(x.un,k.un), c(dnorm(x.un,mu0,s),0), col=col, border=border)
  }
  return( list(k.un=k.un, k.ov=k.ov) )
}##function <- plotCritArea

##fun that adds rejection region
plotBetaArea <- function(mu0, mu.sant, s, k){
    f0 <- dnorm(0, 0, s)
    x <- par('usr') ##first find x (and y) limits
    x <- seq(x[1], x[2], length.out=1000)
    y <- dnorm(x, mu.sant, s)
    P <- pnorm(k$k.ov,mu.sant,s) - pnorm(k$k.un,mu.sant,s)

    ##plot the normal density
    lines(x, y, col='blue', lwd=2)
    ##plot the area(s)
    x.tmp <- seq(max(k$k.un, min(x)), min(k$k.ov,max(x)),
                 length.out=100)
    polygon(c(x.tmp[1],x.tmp,max(x.tmp)), c(0,dnorm(x.tmp,mu.sant,s),0),
            col='blue', border='blue')
    lines(c(mu.sant,mu.sant), c(0,f0), col='black', lty=2)
    ##return Prob ov rejection
    return(P)
}##plotBetaArea <- function(mu0, mu.sant, s, k)
