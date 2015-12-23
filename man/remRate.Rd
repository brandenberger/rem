\name{remRate}
\alias{remRate}
\title{Rate function for relational event models}
%
\description{Estimate the rate function for relational event models.}
\usage{
remRate(formula, 
    dist = 'weibull', 
    var.names = NULL,  
    event.sequence.type = 'continuous', ...)
}
\arguments{
\item{formula}{ A formula expression in the conventional R syntax used for linear modelling. The response (i.e. the dependent variable) must be the event sequence. On the right-hand side, covariates are listed with a \code{+}. It is possible to use one of the endogenous statistics functions in the rem-package as independent (or control) variable. \link{flexsurvreg} is called to estimate the prameters of the model for now.}
\item{dist}{ The distribution specified for the parametric model. See \link{flexsurvreg} for more details.}
\item{var.names}{ An optional string vector containing variable names for the proposed independent variables in the \code{formula}. To be used if endogenous network effects are calculated within the \code{remRate}-function.}
\item{event.sequence.type}{ \code{'continuous'} or \code{'ordinal'}. Only 'continuous' is implemented as of now.}
\item{...}{ Optional arguments for \link{flexsurvreg}. }
}
\details{

The function \code{remRate()} allows users to estimate the rate function of the relational event model. After the variables are prepared, a \link{flexsurvreg} is called to estimate the parameters of the model. 
}
% \value{
% 
% }
% \references{
% 
% }
% \note{
% 
% }
\author{
Laurence Brandenberger \email{laurence.brandenberger@eawag.ch}
}
\seealso{
\link{rem-package}
}
\examples{
# create some data two-mode network event sequence data with
# a 'sender', 'target' and a 'time'-variable
sender <- c('A', 'B', 'A', 'C', 'A', 'D', 'F', 'G', 'A', 'B',
            'B', 'C', 'D', 'E', 'F', 'B', 'C', 'D', 'E', 'C', 
            'A', 'F', 'E', 'B', 'C', 'E', 'D', 'G', 'A', 'G', 
            'F', 'B', 'C')
target <- c('T1', 'T2', 'T3', 'T2', 'T1', 'T4', 'T6', 'T2', 
            'T4', 'T5', 'T5', 'T5', 'T1', 'T6', 'T7', 'T2', 
            'T3', 'T1', 'T1', 'T4', 'T5', 'T6', 'T8', 'T2',
            'T7', 'T1', 'T6', 'T7', 'T3', 'T4', 'T7', 'T8', 'T2')
time <- c('03.01.15', '04.01.15', '10.02.15', '28.02.15', '01.03.15', 
          '07.03.15', '07.03.15', '12.03.15', '04.04.15', '28.04.15',
          '06.05.15', '11.05.15', '13.05.15', '17.05.15', '22.05.15', 
          '09.08.15', '09.08.15', '14.08.15', '16.08.15', '29.08.15',
          '05.09.15', '25.09.15', '02.10.15', '03.10.15', '11.10.15', 
          '18.10.15', '20.10.15', '28.10.15', '04.11.15', '09.11.15', 
          '10.12.15', '11.12.15', '12.12.15')
type <- sample(c('con', 'pro'), 33, replace = TRUE)
important <- sample(c('important', 'not important'), 33,
                    replace = TRUE)
sender.country <- sample(c('country1', 'country2'), 33, 
                         replace = TRUE)

# combine them into a data.frame
dt <- data.frame(sender, target, time, type, important, sender.country)

# create event sequence and order the data
dt <- eventSequence(datevar = dt$time, dateformat = '\%d.\%m.\%y', 
                    data = dt, type = 'continuous', 
                    byTime = 'daily', returnData = TRUE,
                    sortData = TRUE)

# create some endogenous variables: 
dt$inertia <- inertiaStat(dt, dt$event.seq.cont, dt$sender, dt$target, 30)
dt$activitySender <- degreeStat(dt, dt$event.seq.cont, dt$sender, 30)
dt$targetPopularity <- degreeStat(dt, dt$event.seq.cont,dt$target, 30)

# check out their correlation coefficient
cor(subset(dt, select = c(inertia, activitySender, 
                          targetPopularity)), method="pearson") 

# fit event history model using inertia, actor activity and
# target popularity 
fit1 <- remRate(
  dt$event.seq.cont ~ 
    inertiaStat(dt, dt$event.seq.cont, dt$sender, dt$target, 30) +
    degreeStat(dt, dt$event.seq.cont, dt$sender, 30) +
    degreeStat(dt, dt$event.seq.cont,dt$target, 30),
  var.names = c('inertia', 'actorAct', 'targetPop'))
fit1
## Interpretation of the weibull-model coefficients:
# For positive coefficients, if the independent variable increases 
# in value, the survival time increases, meaning the event will occurr 
# at a lower rate. The coefficients are log odds. Calculate the odds 
# ratios by runnding exp(coefficient). This will give you the factor 
# by which the risk of the event occurring increases 
# (negative coefficient, OR < 1) or decreases 
# (positive coefficient, OR > 1).

# calculate closing four-cycles
dt$fourCyc <- fourCycleStat(dt, dt$event.seq.cont, 
                            dt$sender, dt$target, 30)
dt$fourCyc.support <- fourCycleStat(dt, dt$event.seq.cont, 
                                    dt$sender, dt$target, 30, 
                                    eventtypevar = dt$type, 
                                    eventtypevalue = 'positive')

# do actors from country1 engage more often in positive
# closing four-cycles?
library('ggplot2')
ggplot(dt, aes (event.seq.cont, fourCyc.support, 
                color = sender.country)) +
  geom_point()+ geom_smooth() 

# test above idea of closing four-cycles with an 
# interaction term with 'sender.country'
fit2 <- remRate(
  dt$event.seq.cont ~ 
    inertiaStat(dt, dt$event.seq.cont, dt$sender, dt$target, 30) +
    degreeStat(dt, dt$event.seq.cont, dt$sender, 30) +
    degreeStat(dt, dt$event.seq.cont,dt$target, 30) +
    dt$fourCyc.support*dt$sender.country,
  var.names = c('inertia', 'actorActivity', 
                'targetPopularity', 
                'fourCyc.support', 'country',
                'fourCyc.supportXcountry'))
fit2

# sender similarity: how likely is it that two
# senders are alike?
fit3 <- remRate(
  dt$event.seq.cont ~ 
    inertiaStat(dt, dt$event.seq.cont, dt$sender, dt$target, 30) +
    degreeStat(dt, dt$event.seq.cont, dt$sender, 30) +
    degreeStat(dt, dt$event.seq.cont,dt$target, 30) +
    similarityStat(dt, dt$event.seq.cont, dt$sender, dt$target,
                   halflife.last.event = 10),
  var.names = c('inertia', 'actorActivity', 
                'targetPopularity', 'simSender'))
fit3

# target similarity: how likely is it that two
# targets are used together (e.g. the current actor 
# copies others/learns from others)
fit4 <- remRate(
  dt$event.seq.cont ~ 
    inertiaStat(dt, dt$event.seq.cont, dt$sender, dt$target, 30) +
    degreeStat(dt, dt$event.seq.cont, dt$sender, 30) +
    degreeStat(dt, dt$event.seq.cont,dt$target, 30) +
    similarityStat(dt, dt$event.seq.cont, dt$sender, dt$target,
                   senderOrTarget = 'target',
                   halflife.last.event = 10),
  var.names = c('inertia', 'actorActivity', 
                'targetPopularity', 'simTarget'))
fit4

}
%\keyword{key}




