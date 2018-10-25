\name{rem-package}
\alias{rem-package}
\alias{rem}
\alias{relational-event-model}
\alias{relational event model}
\docType{package}
\title{Fit Relational Event Models (REM)}
\description{The \pkg{rem} package uses a combination of event history and network analysis to test network dependencies in event sequences. If events in an event sequence depend on each other, network structures and patterns can be calculated and estimated using relational event models. The \code{rem}-package includes functions to calculate endogenous network statistics in (signed) one-, two- and multi-mode network event sequences. The statistics include inertia (\link{inertiaStat}), reciprocity (\link{reciprocityStat}), in- or outdegree statistics (\link{degreeStat}), closing triads (\link{triadStat}), closing four-cycles (\link{fourCycleStat}) or endogenous similarity statistics (\link{similarityStat}). The rate of event occurrence can then be tested using standard models of event history analysis, such as a stratified Cox model (or a conditional logistic regression). \link{createRemDataset} can be used to create counting process data sets with dynamic risk sets.
}
\details{
\tabular{ll}{
Package: \tab rem\cr
Type: \tab Package\cr
Version: \tab 1.3.1\cr
Date: \tab 2018-10-24\cr
}
}
\author{
Laurence Brandenberger \email{lbrandenberger@ethz.ch}
}
\references{
Lerner, Jurgen, Bussmann, Margit, Snijders, Tom. A., & Brandes, Ulrik. 2013. Modeling frequency and type of interaction in event networks. Corvinus Journal of Sociology and Social Policy, (1), 3-32.

Brandenberger, Laurence. 2018. Trading Favors - Examining the Temporal Dynamics of Reciprocity in Congressional Collaborations Using Relational Event Models. Social Networks, 54: 238-253.

Malang, Thomas, Laurence Brandeberger and Philip Leifeld. 2018. Networks and Social Influence in European Legislative Politics. British Journal of Political Science. DOI: 10.1017/S0007123417000217.
}

