\name{fourCycleStat}
\alias{fourCycleStat}
\alias{fourCycle}
\title{Calculate four cycle statistics}
%
\description{Calculate the endogenous network statistic \code{fourCycle} that 
measures the tendency for events to close four cycles in two-mode event sequences.}
\usage{
fourCycleStat(data, time, sender, target, halflife, 
    weight = NULL, 
    eventtypevar = NULL, 
    eventtypevalue = 'standard', 
    eventfiltervar = NULL, 
    eventfilterAB = NULL, eventfilterAJ = NULL, 
    eventfilterIB = NULL, eventfilterIJ = NULL,
    eventvar = NULL,
    variablename = 'fourCycle',
    returnData = FALSE,
    dataPastEvents = NULL,
    showprogressbar = FALSE, 
    inParallel = FALSE, cluster = NULL
)
}
\arguments{
\item{data}{ A data frame containing all the variables.}

\item{time}{ Numeric variable that represents the event sequence. The variable
has to be sorted in ascending order.}

\item{sender}{ A string (or factor or numeric) variable that represents the sender of the event.}

\item{target}{ A string (or factor or numeric) variable that represents the target of the event.}

\item{halflife}{ A numeric value that is used in the decay function.
The vector of past events is weighted by an exponential decay function using the specified halflife. The halflife parameter determins after how long a period the event weight should be halved. E.g. if \code{halflife = 5}, the weight of an event that occured 5 units in the past is halved. Smaller halflife values give more importance to more recent events, while larger halflife values should be  used if time does not affect the sequence of events that much.}

\item{weight}{ An optional numeric variable that represents the weight of each event. If \code{weight = NULL} each event is given an event weight of \code{1}.
}

\item{eventtypevar}{ An optional variable that represents the type of the event. Use \code{eventtypevalue} to specify how the \code{eventtypevar} should be used to filter past events.
}

\item{eventtypevalue}{ An optional value (or set of values) used to specify how paste events should be filtered depending on their type. \code{'standard'}, \code{'positive'} or \code{'negative'} may be used. Default set to \code{'standard'}. \code{'standard'} referrs to closing four cylces where the type of the events is irrelevant. \code{'positive'} closing four cycles can be classified as reciprocity via the second mode. It indicates whether senders have a tendency to reciprocate or show support by engaging in targets that close a four cycle between two senders. \code{'negative'} closing four cycles represent opposition between two senders, where the current event is more likely if the two senders have opposed each other in the past. Support or opposition is represented by the \code{eventtypevar} value for each event.}

\item{eventfiltervar}{ An optinoal variable that allows filtering of past events using an event attribute. It can be a sender attribute, a target attribute, time or dyad attribute.
Use \code{eventfilterAB}, \code{eventfilterAJ}, \code{eventfilterIB} or \code{eventfilterIJ} to specify how the \code{eventfiltervar} should be used.}

\item{eventfilterAB}{ An optional value used to specify how 
paste events should be filtered depending on their attribute. Each distinct edge that form a four cycle can be filtered. \code{eventfilterAB} refers to the current event. \code{eventfilterAJ} refers to the event involving the current sender and target \code{j} that has been used by the current as well as the second actor in the past. \code{eventfilterIB} refers to the event involving the second sender and the current target. \code{eventfilterIJ} filters events that involve the second sender and the second target. See the four cycle formula in the \code{details} section for more information.}

\item{eventfilterAJ}{ see \code{eventfilterAB}.}

\item{eventfilterIB}{see \code{eventfilterAB}.}

\item{eventfilterIJ}{see \code{eventfilterAB}.}

\item{eventvar}{ An optional dummy variable with 0 values for null-events and 1 values for true events. If the \code{data} is in the form of counting process data, use the \code{eventvar}-option to specify which variable contains the 0/1-dummy for event occurrence. If this variable is not specified, all events in the past will be considered for the calulation of the four cycle statistic, regardless if they occurred or not (= are null-events). Misspecification could result in grievous errors in the calculation of the network statistic.}

\item{variablename}{ An optional value (or values) with the name the four cycle statistic variable should be given. To be used if \code{returnData = TRUE}.}

\item{returnData}{ \code{TRUE/FALSE}. Set to \code{FALSE} by default. The new  variable(s) are bound directly to the \code{data.frame} provided and the data frame is returned in full.}

\item{dataPastEvents}{ An optional \code{data.frame} with the following variables:
	column 1 = time variable, 
	column 2 = sender variable, 
	column 3 = target on other variable (or all "1"), 
	column 4 = weight variable (or all "1"), 
	column 5 = event type variable (or all "1"), 
	column 6 = event filter variable (or all "1"). Make sure that the data frame does not contain null events. Filter it out for true events only.}

\item{showprogressbar}{\code{TRUE/FALSE}. To be implemented.}

\item{inParallel}{ \code{TRUE/FALSE}. An optional boolean to specify if the loop should be run in parallel.}

\item{cluster}{ An optional numeric or character value that defines the cluster. By specifying a single number, the cluster option uses the provided number of nodes to parallellize. By specifying a cluster using the \code{makeCluster}-command in the \code{doParallel}-package, the loop can be run on multiple nodes/cores. E.g., \code{cluster = makeCluster(12, type="FORK")}.}

}
\details{
The \code{fourCycleStat()}-function calculates an endogenous statistic that measures whether events have a tendency to form four cycles.

The effect is calculated as follows:

\deqn{G_t = G_t(E) = (A, B, w_t), }{G_t = G_t(E) = (A, B, w_t),}

\eqn{G_t} represents the network of past events and includes all events \eqn{E}. These events consist
each of  a sender \eqn{a \in A}{a in A} and a target \eqn{b \in B}{b in B} and a weight function \eqn{w_t}:

\deqn{  w_t(i, j) = \sum_{e:a = i, b = j} | w_e | \cdot e^{-(t-t_e)\cdot\frac{ln(2)}{T_{1/2}}} \cdot \frac{ln(2)}{T_{1/2}}, }{ w_t(i, j) = \sum_{e:a = i, b = j} | w_e | * exp^{-(t-t_e)* (ln(2)/T_{1/2})} * (ln(2)/T_{1/2}),}

where \eqn{w_e} is the event weight (usually a constant set to 1 for each event), \eqn{t} is the current event time, \eqn{t_e} is the past event time and \eqn{T_{1/2}} is a halflife parameter.

For the four-cylce effect, the past events \eqn{G_t} are filtered to include only events
where the current event closes an open four-cycle in the past.

\deqn{fourCycle(G_t , a , b) = \sqrt[3]{\sum_{i \in A \& j \in B} w_t(a, j) \cdot w_t(i, b) \cdot w_t(i, j)}}{fourCycle(G_t , a , b) = (\sum_{i in A and j in B} w_t(a, j) * w_t(i, b) * w_t(i, j))^(1/3)}

An exponential decay function is used to model the effect of time on the endogenous statistics. The further apart the past event is from the present event, the less weight is given to this event. The halflife parameter in the \code{fourCycleStat()}-function determins at which rate the weights of past events should be reduced. Therefore, if the one (or more) of the three events in the four cycle have ocurred further in the past, less weight is given to this four cycle because it becomes less likely that the two senders reacted to each other in the way the four cycle assumes.

The \code{eventtypevar}- and \code{eventfiltervar}-options help filter the past events more specifically. How they are filtered depends on the \code{eventtypevalue}- and \code{eventfilter__}-option.
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

# combine them into a data.frame
dt <- data.frame(sender, target, time, type, important)

# create event sequence and order the data
dt <- eventSequence(datevar = dt$time, dateformat = '\%d.\%m.\%y', 
                    data = dt, type = 'continuous', 
                    byTime = "daily", returnData = TRUE,
                    sortData = TRUE)

# create counting process data set (with null-events) - conditional logit setting
dts <- createRemDataset(dt, dt$sender, dt$target, dt$event.seq.cont, 
                          eventAttribute = dt$type, 
                          atEventTimesOnly = TRUE, untilEventOccurrs = TRUE, 
						  returnInputData = TRUE)
## divide up the results: counting process data = 1, original data = 2
dtrem <- dts[[1]]
dt <- dts[[2]]
## merge all necessary event attribute variables back in
dtrem$type <- dt$type[match(dtrem$eventID, dt$eventID)]
dtrem$important <- dt$important[match(dtrem$eventID, dt$eventID)]
# manually sort the data set
dtrem <- dtrem[order(dtrem$eventTime), ]

# calculate closing four-cycle statistic
dtrem$fourCycle <- fourCycleStat(data = dtrem, 
                              time = dtrem$eventTime, 
                              sender = dtrem$sender, 
                              target = dtrem$target, 
                              eventvar = dtrem$eventDummy,
                              halflife = 20)

# plot closing four-cycles over time:
library("ggplot2")
ggplot(dtrem, aes (eventTime, fourCycle,
                group = factor(eventDummy), color = factor(eventDummy)) ) +
  geom_point()+ geom_smooth() 

# calculate positive closing four-cycles: general support
dtrem$fourCycle.pos <- fourCycleStat(data = dtrem, 
                                 time = dtrem$eventTime, 
                                 sender = dtrem$sender, 
                                 target = dtrem$target, 
                                 eventvar = dtrem$eventDummy,
                                 eventtypevar = dtrem$type, 
                                 eventtypevalue = 'positive',
                                 halflife = 20)

# calculate negative closing four-cycles: general opposition
dtrem$fourCycle.neg <- fourCycleStat(data = dtrem, 
                                     time = dtrem$eventTime, 
                                     sender = dtrem$sender, 
                                     target = dtrem$target, 
                                     eventvar = dtrem$eventDummy,
                                     eventtypevar = dtrem$type, 
                                     eventtypevalue = 'negative',
                                     halflife = 20)
}
%\keyword{key}




