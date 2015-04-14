###################################################################
# sivipm R package
# Copyright INRA 2015
# INRA, UR1404, Research Unit MaIAGE
# F78352 Jouy-en-Josas, France.
#
# URL: http://cran.r-project.org/web/packages/sivipm
#
# This file is part of sivipm R package.
# sivipm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# See the GNU General Public License at:
# http://www.gnu.org/licenses/
#
###################################################################



##########################################################################
# Calcul des indices de sensibilite isivip et tsivip après regression PLS
##########################################################################
#  TSIVIP total or individual sensitivity indices or
# rule implemented in SIMCA-P to determine the significative components
# @title compute several TSIVIP results 
# @sivipm
# @param Y outputs  data.frame 
# @param XIndic, an object of class 'polyX'
# @param nc number of components
# @param options outputs selection. String vector.
# Valid values: 'tsivip', 'isivip', 'signif', 'all'
# @param graph if TRUE, display graph (for tsivip when alea =false)
# @param alea alea is added (only for tsivip)
# @param output additional results are returned
# @return
# \itemize{
# \item{tsivip}{ - when \code{options} include 'tsivip',  total sensitivity indices for each input variable and their percentage.}
# \item{isivip}{ - when \code{options} include 'isivip',
# first order individual sensitivity indices for each monome.
# When alea=TRUE, the non significant monomes are removed}
# \item{SIMCA}{ - when \code{options} include 'simca',
# the significative components (when no missing) calculated
# by SIMCA-P rule,
# \item{Lazraq}{ - when \code{options} include 'lazraq',
# the significative components
# (missing accepted if one response, only) calculated
# by the test de Lazraq).}
# \item{output}{ - required additional results.}
# }
# @export

sivipm <- function(Y, XIndic, nc = 2, options = c("fo.isivip", "tsivip", "simca", 
    "lazraq"), graph = FALSE, alea = FALSE, output = NULL) {
    # Initialisation des retours
    retour <- list()
    retisivip <- NULL
    rettsivip <- NULL
    retsimca <- NULL
    retlaz <- NULL
    regpls <- NULL
    
    # Acces aux donnees
    data.exp <- XIndic@dataX.exp
    if (!is.null(XIndic@Pindic)) {
        Indic <- XIndic@Pindic@indic
        nvar <- ncol(Indic)
    } else {
        Indic <- NULL
        nvar <- ncol(data.exp)
    }
    
    # Acces aux options
    opt.isivip <- any(grepl("isivip", options, ignore.case = TRUE)) || any(grepl("isivip", 
        output, ignore.case = TRUE))
    opt.tsivip <- any(grepl("tsivip", options, ignore.case = TRUE))
    opt.simca <- any(grepl("simca", options, ignore.case = TRUE))
    opt.laz <- any(grepl("lazraq", options, ignore.case = TRUE))
    
    
    if (!opt.tsivip && alea) {
        warning("sivipm : option alea ignored when tsivip not required")
    }
    if (!opt.tsivip && graph) {
        warning("sivipm : option graph ignored when tsivip not required")
    }
    
    # avec/sans missing values
    if (any(is.na(data.exp)) || any(is.na(Y))) 
        na.miss <- TRUE else na.miss <- FALSE
    
    # Faut-il calculer des sorties supplémentaires?
    if (!is.null(output)) 
        routput <- TRUE else routput <- FALSE
    
    if (opt.simca && na.miss) {
        warning("sivipm : option simca ignored when there are missing values.")
        opt.simca <- FALSE
    }
    if (opt.laz && na.miss && (ncol(as.matrix(Y)) > 1)) {
        warning("sivipm : option Lazraq ignored when there are missing values and several output variables.")
        opt.laz <- FALSE
    }
    
    if (opt.simca || opt.laz) {
        # Call regpls2 with all outputs because Q2 required
        regpls <- regpls2(Y, data.exp, nc, output = TRUE)
        if (opt.simca) 
            retsimca <- fsimcarule(regpls$Q2)
        if (opt.laz) 
            retlaz <- rlaz(regpls$x.scores, Y, nc)
    }
    
    
    # pas de signif et (isivip ou (tsivip sans alea))
    if (is.null(regpls) && (opt.isivip || (opt.tsivip && alea == FALSE))) {
        regpls <- regpls2(Y, data.exp, nc, output = routput)
    }
    
    
    
    if (opt.isivip) {
        retisivip <- isivip(regpls$VIP, data.exp)
    }
    
    if (opt.tsivip) 
        {
            if (all(is.na(Indic))) {
                warning("sivipm : option tsivip ignored because there is no polynome description")
            } else {
                if (alea) {
                  # cas avec alea
                  rett <- tsivipalea(Y, XIndic, nc, output = routput)
                  if (!is.null(retisivip)) 
                    {
                      # rajouter le isivip de l'alea parmi ceux-ci
                      nvar <- ncol(Indic)
                      mono <- names(retisivip)
                      nmono <- length(mono)
                      bb <- c(retisivip[1:nvar], rett$isivipalea, retisivip[(nvar + 
                        1):nmono])
                      names(bb) <- c(mono[1:nvar], "X_alea", mono[(nvar + 1):nmono])
                      retisivip <- bb
                    }  # fin !is.null(retisivip)
    # Oter des sorties ce qui concerne l'alea
    rett$isivipalea <- NULL
    nrett <- length(rett$tsivip)
    rett$tsivip <- rett$tsivip[-nrett]
    rett$percentage <- rett$percentage[-nrett]
    nmono <- length(rett$monosignif)
    rett$monosignif <- rett$monosignif[-nmono]

                  
                  regpls <- rett$regpls
                  rett$regpls <- NULL
                  rettsivip <- rett
                  if (graph == TRUE) 
                    tsivipgraph(rettsivip$tsivip, retisivip, Indic, nc, alea)
                } else {
                  # cas sans alea
                  rettsivip <- tsivipnoalea(regpls$VIP, Y, Indic, nc)
                  if (graph == TRUE) {
                    tsivipgraph(rettsivip$tsivip, retisivip, Indic, nc, alea, rettsivip$Evol)
                  }
                  rettsivip$Evol <- NULL
                }
            }  # end else
        }  # end options=='tsivip'
    
    # Retour
    
    if (any(grepl("isivip", options))) {
        # First order isivip
        retour$fo.isivip <- retisivip[1:nvar]
    }
    
    retour <- c(retour, rettsivip)
    retour <- c(retour, retsimca, retlaz)
    
    if (!is.null(regpls)) 
        {
            retour$output <- list()
            if (any(grepl("isivip", output, ignore.case = TRUE))) {
                retour$output$isivip <- retisivip
            }
            
            if (any(grepl("VIP", output))) {
                retour$output$VIP <- regpls$VIP
                retour$output$VIPind <- regpls$VIPind
                regpls$VIP <- NULL
                regpls$VIPind <- NULL
            }
            if (any(grepl("RSS", output, ignore.case = TRUE))) {
                retour$output$RSS <- regpls$RSS
                regpls$RSS <- NULL
            }
            if (any(grepl("PRESS", output, ignore.case = TRUE))) {
                retour$output$PRESS <- regpls$PRESS
                regpls$PRESS <- NULL
            }
            if (any(grepl("Q2", output, ignore.case = TRUE))) {
                retour$output$Q2 <- regpls$Q2
                retour$output$Q2cum <- regpls$Q2cum
                regpls$Q2 <- NULL
                regpls$Q2cum <- NULL
            }
            if (any(grepl("betaNat", output, ignore.case = TRUE))) {
                retour$output$betaNat <- regpls$betaNat
                retour$output$betaNat0 <- regpls$betaNat0
                regpls$betaNat <- NULL
                regpls$betaNat0 <- NULL
            }
 

            if (any(grepl("PLS", output, ignore.case = TRUE))) {
# ne  mettre dans PLS que les sorties spécifiques au cas ou on
# ne serait pas passé par les tests ci-dessus 
# (PLS demandé seulement, par exemple)
                regpls$VIP <- NULL
                regpls$VIPind <- NULL
                regpls$RSS <- NULL
		regpls$PRESS <- NULL
		regpls$Q2 <- NULL
                regpls$Q2cum <- NULL
		 regpls$betaNat <- NULL
                regpls$betaNat0 <- NULL
                retour$output$PLS <- regpls
            }
            if (length(retour$output) == 0) 
                retour$output <- NULL
        }  # fin !regpls
    
    if (any(grepl("isivip", output, ignore.case = TRUE))) 
        retour$output$isivip <- retisivip
    
    return(retour)
}  # end sivipm

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# @fsimcarule
# @title compute rule implemented in SIMCA-P to determine the number of significative components from Q2 (no missing values)
# Internal function
#@return  significative components

fsimcarule <- function(Q2) {
    
    a <- ncol(Q2)
    # signif <- as.integer(Q2[, a] >= 0.0975) s <- max(2, max(which(signif == 1)) +
    # 1) ret <- list(signifcomponents = s)
    signif <- (Q2[, a] >= 0.0975)
    names(signif) <- paste("t", 1:length(signif), sep = "")
    ret <- list(simca.signifcomponents = signif)
    return(ret)
}


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Calcul des indices globaux tsivip
# @title compute TSIVIP indices
# @tsivipf
# @param vecin VIP vector
# @param Indic indic matrix
# @return tsivip sensitivity indices for each input variable
# Internal function

tsivipf <- function(vecin, Indic) {
    a <- isivip(vecin, t(Indic))
    tsivip <- t(Indic) %*% a
    # AB 21/10/2013
    tri <- sort(tsivip, index.return = TRUE, decreasing = TRUE)
    if (length(tri$x) == 0) {
        rank <- 1:length(tsivip)
    } else {
        rank <- tri$ix
        tsivip <- tri$x
    }
    
    # b<-barplot(tsivip,xlab=as.character(rank))
    ret <- list(tsivip, rank)
    return(ret)
}  # end tsivipf

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# @tsivipgraph
# internal
# TSIVIP graphs

tsivipgraph <- function(tsivip, isivip, Indic, nc, alea, Evol = NULL) {
    taillept <- 0.8  # taille des annotations
    d <- seq(1, nc)
    grcomp <- (!is.null(Evol) && (nc > 1))
    if (grcomp) {
        # on fera le graphe des composantes
        nombrecomposantes <- d
        for (i in 1:(ncol(Indic) - 1)) {
            nombrecomposantes <- cbind(nombrecomposantes, d)
        }
        par(mfrow = c(2, 1))
    } else {
        # on ne fera pas le graphe des composantes
        par(mfrow = c(1, 1))
    }
    sortsivip <- sort(tsivip, decreasing = TRUE, index.return = TRUE)
    xsivip <- sortsivip$ix
    label <- names(tsivip)
    # space=0, pour que les barres soient jointes las=1, pour que les noms des
    # monomes soient horizontaux
    barplot(sortsivip$x, xlab = "TSIVIP (%)", names.arg = label[xsivip], horiz = TRUE, 
        space = 0, cex.axis = taillept, cex.names = taillept, cex.lab = taillept, 
        las = 1)
    if (!is.null(isivip)) 
        {
            # superposer les isivip des monomes unitaires
            visivip <- isivip[xsivip]
            barplot(visivip, add = TRUE, names.arg = rep("", length(visivip)), col = 1, 
                angle = 45, density = 20, legend.text = "First order ISIVIP", space = 0, 
                horiz = TRUE, args.legend = list(bty = "n", cex = taillept), 
                yaxt = "n", xaxt = "n")
        }  # fin (!is.null(isivip))
    if (grcomp) 
        {
            # Evol est nul dans le cas alea pour mettre les ticks-marks sur l'axe des x en
            # entier, on les met 'a la main'
            
            matplot(nombrecomposantes, Evol, type = "l", ylab = "TSIVIP", xlab = "number of components", 
                xaxt = "n", cex.axis = taillept, cex.lab = taillept)
            axis(1, at = 1:nc, labels = 1:nc, cex.axis = taillept)
            # pie(Evol[nrow(Evol),])
        }  # fin (!is.null(Evol)
    return(invisible())
}  # fin tsivipgraph


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# @tsivipnoalea
# @title compute TSIVIP total sensitivity indices from VIP
# when no alea
# Internal function
# @param VIP component VIP returned by regpls2
# @param Y outputs  data.frame 
# @param nc number of components
# @param Indic matrix computed by indic function
# @param graph  display graph
# @return tsivip total indices for each input variable

tsivipnoalea <- function(VIP, Y, Indic, nc = 2, graph = FALSE) {
    
    if (!is.null(colnames(Indic))) 
        varnames <- colnames(Indic) else varnames <- paste("X", 1:ncol(Indic), sep = "")
    
    
    TSIVIP <- matrix(nrow = nc, ncol = ncol(Indic))
    TSIVIPr <- matrix(nrow = nc, ncol = ncol(Indic))
    for (j in 1:nc) {
        VIPj <- VIP[, j]
        aa <- tsivipf(VIPj, Indic)
        TSIVIP[j, ] <- aa[[1]]
        TSIVIPr[j, ] <- aa[[2]]
    }
    r <- matrix(nrow = ncol(Indic), ncol = nc)
    for (i in 1:ncol(Indic)) {
        for (j in 1:nc) {
            r[i, j] <- which(TSIVIPr[j, ] == i)
        }
    }
    Evol <- matrix(ncol = ncol(Indic), nrow = nc)
    for (i in 1:ncol(Indic)) {
        for (k in 1:nc) {
            Evol[k, i] <- TSIVIP[k, r[i, k]]
        }
    }
    dimnames(Evol) <- list(1:nc, varnames)
    # Prendre la derniere ligne de Evol, i.e ce qui correspond a la derniere
    # composante
    tsivip <- Evol[nrow(Evol), ]
    names(tsivip) <- varnames
    
    percent <- sort((tsivip/sum(tsivip)) * 100, decreasing = TRUE)
    ret <- list(tsivip = tsivip, percentage = percent, Evol = Evol)
    return(ret)
    
}  # end tsivipmnoalea
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# @tsivipalea
# @title compute TSIVIP total sensitivity indices from VIP
# with alea
# Internal function
# @param Y outputs  data.frame 
# @param XIndic, an object of class 'polyX'
# @param nc number of components
# @param Indic matrix computed by indic function
# @return tsivip total indices for each input variable
#  The non significant monomes are removed
# from the tsivip calculation.


tsivipalea <- function(Y, XIndic, nc, output = FALSE) {
    dataX.exp <- XIndic@dataX.exp
    Indic <- XIndic@Pindic@indic
    
    if (!is.null(colnames(Indic))) 
        varnames <- colnames(Indic) else varnames <- paste("X", 1:ncol(Indic), sep = "")
    
    X_alea <- runif(nrow(dataX.exp))
    # Add a random variable into dataX.exp
    Xalea <- cbind(dataX.exp, X_alea)
    # Add a column and a row into Indic
    Indicalea <- cbind(Indic, rep(0, nrow(Indic)))
    varnames <- c(varnames, "X_alea")
    colnames(Indicalea) <- varnames
    Indicalea <- rbind(Indicalea, rep(0, ncol(Indicalea)))
    Indicalea[nrow(Indicalea), ncol(Indicalea)] <- 1
    regpls <- regpls2(Y, Xalea, nc, output = output)
    VIP <- regpls$VIP
    VIP <- VIP[, ncol(VIP)]
    isivipalea <- isivip(VIP, Xalea)
    leisivipalea <- isivipalea[length(isivipalea)]
    signif <- (isivipalea > leisivipalea)
    isivipalea <- isivipalea[signif]
    Indicalea <- Indicalea[signif, ]
    tsivipa <- t(Indicalea) %*% isivipalea
    tsivipa <- as.vector(tsivipa)
    names(tsivipa) <- varnames
    percent <- sort((tsivipa/sum(tsivipa)) * 100, decreasing = TRUE)
    correlation <- cor(X_alea, Y)
    retour <- list(tsivip = tsivipa, isivipalea = leisivipalea, percentage = percent, 
        monosignif = signif, correlalea = correlation)
    if (output) {
        retour$regpls <- regpls
    }
    
    return(retour)
    
    
}  # end tsivipalea
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# rlaz
# compute the Lazraq et Cléroux inferential test It is used to determine the
# significant components when there are missing values and a single response
# variable.  Input TT: x.scores. Matrix nc X nobs Y : response variable. Matrix
# nobs X 1 nc: number of components Internal function
rlaz <- function(TT, Y, nc) {
    YY <- scale(Y)  # centered-reduced response variable
    alpha <- 0.05  # alpha est fixe
    nobs <- nrow(YY)
    stu <- qt(1 - alpha, nobs)
    res <- apply(TT, 1, function(T, Y, nc, nobs, stu) {
        mat <- cbind(Y, as.matrix(T, ncol = 1))
        cr <- cor(mat, use = "na.or.complete")[1, 2]
        tlaz <- (sqrt(nobs) * cr)/sqrt(1 - cr * cr)
        if (tlaz > stu) {
            return(TRUE)
        } else {
            return(FALSE)
        }
    }, Y, nc, nobs, stu)
    return(list(lazraq.signifcomponents = res))
}  # end rlaz 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calcul des indices individuels isivip
# @title compute isivip individual sensitivity indices for each monome
# @isivip
# @param vecin VIP vector
# @param dataX.exp inputs data.frame
# @return isivip individual sensitivity indices for each monome

isivip <- function(vecin, dataX.exp) {
    sivip <- rep(0, ncol(dataX.exp))
    if (!is.null(colnames(dataX.exp)))
      names(sivip) <- colnames(dataX.exp)
    else
      names(sivip) <- paste("X", 1:ncol(dataX.exp))
    
    for (i in 1:ncol(dataX.exp)) {
        sivip[i] <- vecin[i] * vecin[i]/ncol(dataX.exp)
    }
    return(sivip)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# @title compute TSIVIP total sensitivity indices 
# @gosivip
# @param Y outputs  data.frame
# @param XIndic, an object of class 'polyX'
# @param nc number of components
# @param graph=FALSE display graph
# @return tsivip total indices for each input variable
# When alea=TRUE, the non significant monomes are removed
# from the tsivip calculation.
# It is a fast version of function of sivipm, because reduced to the
# case where only tsivip are required. Used in bootstrap.

gosivip <- function(Y, XIndic, nc=2, graph = FALSE, alea = FALSE) {

    dataX.exp <- XIndic@dataX.exp
    Indic <- XIndic@Pindic@indic
    if (!alea) 
        {
          regpls <- regpls2(Y, dataX.exp, nc,  output=FALSE)
          ret <- tsivipnoalea(regpls$VIP, Y,   Indic, nc, graph= graph)

        }  # end !alea
    
    if (alea) {
        ret <- tsivipalea(Y, XIndic, nc,  output=FALSE)
      } # end alea
    
    return(ret)  
    
} # end gosivip
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# @title realise des intervalles de confiance pour les tsivip par methode de bootstrap
# @sivipboot
# @param Y ouputs data.frame 
# @param XIndic, an object of class \code{\linkS4class{polyX}}
# @param nc number of components
# @param B number of bootstrap replicates
# @param alpha level of bootstrap confidence intervals
# @return tsivip percentile bootstrap confidence intervals

sivipboot <- function(Y, XIndic, B, nc=2, graph=FALSE, alea=FALSE,
                      alpha = 0.05) {
    dataX.exp <- XIndic@dataX.exp
    Indic <- XIndic@Pindic@indic
    if (!alea)
      boot <- matrix(nrow = B, ncol = ncol(Indic))
    else
      boot <- matrix(nrow = B, ncol = ncol(Indic)+1)
    # dans le cas aleatoire, on rajoute une variable
    
    # MAXTRY: maximal number of X generation at each loop
    MAXTRY = 5
    for (b in 1:B) {
        T <- sample(nrow(dataX.exp), nrow(dataX.exp), replace = TRUE)

        # AB: check scale is possible
        Xscale <- scale(as.matrix(dataX.exp[T, ]))
        ntry <- 1
        while (any(is.nan(Xscale)) && (ntry < MAXTRY) ) {
          T <- sample(nrow(dataX.exp), nrow(dataX.exp), replace = TRUE)
          Xscale <- scale(as.matrix(dataX.exp[T, ]))
          ntry <- ntry + 1
        } # end while

        
        if (any(is.nan(Xscale)))
            stop("sivipboot: unsuccessful generation of X")
        # End check scale

        XIndic@dataX.exp <- dataX.exp[T, ]
        ressivipm <- gosivip(Y[T, ], XIndic, nc, alea=alea,
                            graph = FALSE)$tsivip
          
        boot[b, ] <- ressivipm
         
    }
    varnames <- names(ressivipm)
    IC.inf <- rep(0, ncol(boot))
    IC.sup <- rep(0, ncol(boot))
    
    for (i in 1:ncol(boot)) {
        IC.inf[i] <- sort(boot[, i])[ceiling(B * (alpha/2))]
        IC.sup[i] <- sort(boot[, i])[ceiling(B * (1 - alpha/2))]
    }
    IC <- cbind(IC.inf, IC.sup)
    if (graph) {
      taillept <- 0.8 # taille des annotations
      colnames(boot) <- colnames(XIndic@Pindic@indic)
      # search for the median
         bb <- boxplot(boot, plot=FALSE)
      # tri selon la mediane
      bootsort <- boot[, order(bb$stats[3,])]
      boxplot(bootsort,cex.axis=taillept, cex.names=taillept, cex.lab=taillept, las=3)
       }
    rownames(IC) <- varnames
    return(IC)
} # end sivipboot
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


