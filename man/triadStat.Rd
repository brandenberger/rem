\name{triadStat}
\alias{triadStat}
\alias{triad}
\title{Calculate triad statistics}
%
\description{Calculate the endogenous network statistic \code{triads} that measures the tendency for events to close open triads.}
\usage{
triadStat(data, time, sender, target, halflife, 
    weight = NULL, 
    eventtypevar = NULL, 
    eventtypevalues = NULL, 
    eventfiltervar = NULL,
    eventfilterAI = NULL, 
    eventfilterBI = NULL,
    eventfilterAB = NULL, 
    eventvar = NULL,
    variablename = 'triad', 
    returnData = FALSE, 
    showprogressbar = FALSE, 
    inParallel = FALSE, cluster = NULL
)
}
\arguments{
\item{data}{ A data frame containing all the variables.}

\item{time}{ Numeric variable that represents the event sequence. The variable has to be sorted in ascending order.}

\item{sender}{ A string (or factor or numeric) variable that represents the sender of the event.}

\item{target}{ A string (or factor or numeric) variable that represents the target of the event.}

\item{halflife}{ A numeric value that is used in the decay function. The vector of past events is weighted by an exponential decay function using the specified halflife. The halflife parameter determins after how long a period the event weight should be halved. E.g. if \code{halflife = 5}, the weight of an event that occured 5 units in the past is halved. Smaller halflife values give more importance to more recent events, while larger halflife values should be used if time does not affect the sequence of events that much.}
 
\item{weight}{ An optional numeric variable that represents the weight of each event. If \code{weight = NULL} each event is given an event weight of \code{1}.}

\item{eventtypevar}{ An optional dummy variable that represents the type of the event. Use \code{eventtypevalues} to specify how the \code{eventtypevar} should be used to filter past events. Specifying the \code{eventtypevar} is needed to calculate effects of social balance theory, such as 'friend-of-friend' or 'enemy-of-enemy' statistics.}

\item{eventtypevalues}{ Two string values that represent the type of the past events. The first string value represents the eventtype that exists for all past events that include the current sender (either as sender or target) and a third actor. The second value represents the eventtype for all past events that include the target (either as sender or target) as well as the third actor.
An example: Let the \code{eventtypevar} indicate whether an event is of cooperative or hostile nature.
To test whether the hypothesis 'the friend of my friend is my friend' holds, both \code{eventtypevalues} must be the same and point to the cooperative type (e.g. \code{eventtypevalues = c("cooperation", "cooperation")}) depending on
how the \code{eventtypevar} is coded.
To test whether the hypothesis 'the friend of my enemy is my enemy' holds, 
the first value in \code{eventtypevalues} represents the hostile event between current sender and a third actor and the second value represents the cooperative event between the third actor and the target.
To test the hypothesis 'the enemy of my enemy is my friend', the first value represents the hostile events between current sender and a third actor and the second value represents the hostile event between the current target and the third actor.
For the fourth hypothesis, to test social balance theory 'the enemy of my friend is my enemy', the first value represents a cooperative event between the current sender and a third actor and the second value represents a hostile event between the current target and the third actor.}

\item{eventfiltervar}{ An optional string (or factor or numeric) variable that can be used to filter past and current events. Use \code{eventfilterAI}, \code{eventfilterBI} or \code{eventfilterAB} to specify which past events should be filtered and by what value.}

\item{eventfilterAI}{ An optional value used to specify how paste events should be filtered depending on their attribute. Each distinct edge that form a triad can be filtered. \code{eventfilterAI} refers to the past event involving the current sender (a) and a third actor (i). \code{eventfilterBI}referrs to past events involving target (b) and the third actor (i). \code{eventfilterAB} refers to the current event involving sender (a) and target (b).}

\item{eventfilterBI}{ see \code{eventfilterAI}.}

\item{eventfilterAB}{ see \code{eventfilterAI}.}

\item{eventvar}{ An optional dummy variable with 0 values for null-events and 1 values for true events. If the \code{data} is in the form of counting process data, use the \code{eventvar}-option to specify which variable contains the 0/1-dummy for event occurrence. If this variable is not specified, all events in the past will be considered for the calulation of the triad statistic, regardless if they occurred or not (= are null-events).}

\item{variablename}{ An optional value (or values) with the name the triad
statistic variable should be given. To be used if \code{returnData = TRUE}.}

\item{returnData}{ \code{TRUE/FALSE}. Set to \code{FALSE} by default. The new 
variable is bound directly to the \code{data.frame} provided and the
data frame is returned in full.}

\item{showprogressbar}{ \code{TRUE/FALSE}. Can only be set to TRUE if the function is not run in parallel.}

\item{inParallel}{ \code{TRUE/FALSE}. An optional boolean to specify if the loop should be run in parallel.}

\item{cluster}{ An optional numeric or character value that defines the cluster. By specifying a single number, the cluster option uses the provided number of nodes to parallellize. By specifying a cluster using the \code{makeCluster}-command in the \code{doParallel}-package, the loop can be run on multiple nodes/cores. E.g., \code{cluster = makeCluster(12, type="FORK")}.}
}
\details{

The \code{triadStat()}-function calculates an endogenous statistic that measures whether events have a tendency to form closing triads.

The effect is calculated as follows:

\deqn{G_t = G_t(E) = (A, B, w_t), }{G_t = G_t(E) = (A, B, w_t),}

\eqn{G_t} represents the network of past events and includes all events \eqn{E}. These events consist
each of  a sender \eqn{a \in A}{a in A} and a target \eqn{b \in B}{b in B} and a weight function \eqn{w_t}:

\deqn{  w_t(i, j) = \sum_{e:a = i, b = j} | w_e | \cdot e^{-(t-t_e)\cdot\frac{ln(2)}{T_{1/2}}} \cdot \frac{ln(2)}{T_{1/2}}, }{ w_t(i, j) = \sum_{e:a = i, b = j} | w_e | * exp^{-(t-t_e)* (ln(2)/T_{1/2})} * (ln(2)/T_{1/2}),}

where \eqn{w_e} is the event weight (usually a constant set to 1 for each event), \eqn{t} is the current event time, \eqn{t_e} is the past event time and \eqn{T_{1/2}} is a halflife parameter.

For the triad effect, the past events \eqn{G_t} are filtered to include only events
where the current event closes an open triad in the past.

\deqn{triad(G_t , a , b) = \sqrt{\sum_{i \in A} w_t(a, i) \cdot w_t(i, b)}}{triad(G_t , a , b) = (\sum_{i in A} w_t(a, i) * w_t(i, b))^(1/2)}

An exponential decay function is used to model the effect of time on the endogenous statistics. The further apart the past event is from the present event, the less weight is given to this event. The halflife parameter in the \code{triadStat()}-function determines at which rate the weights of past events should be reduced. Therefore, if the one (or more) of the two events in the triad have occurred further in the past, less weight is given to this triad because it becomes less likely that the sender and target actors reacted to each other in the way the triad assumes.

The \code{eventtypevar}- and \code{eventattributevar}-options help filter the past events more specifically. How they are filtered depends on the \code{eventtypevalue}- and \code{eventattributevalue}-option.


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
# create some data with 'sender', 'target' and a 'time'-variable
# (Note: Data used here are random events from the Correlates of War Project)
sender <- c('TUN', 'UNK', 'NIR', 'TUR', 'TUR', 'USA', 'URU', 
            'IRQ', 'MOR', 'BEL', 'EEC', 'USA', 'IRN', 'IRN', 
            'USA', 'AFG', 'ETH', 'USA', 'SAU', 'IRN', 'IRN',
            'ROM', 'USA', 'USA', 'PAN', 'USA', 'USA', 'YEM', 
            'SYR', 'AFG', 'NAT', 'UNK', 'IRN')
target <- c('BNG', 'RUS', 'JAM', 'SAU', 'MOM', 'CHN', 'IRQ', 
            'AFG', 'AFG', 'EEC', 'BEL', 'ITA', 'RUS', 'UNK',
            'IRN', 'RUS', 'AFG', 'ISR', 'ARB', 'USA', 'USA',
            'USA', 'AFG', 'IRN', 'IRN', 'IRN', 'AFG', 'PAL',
            'ARB', 'USA', 'EEC', 'IRN', 'CHN')
time <- c('800107', '800107', '800107', '800109', '800109', 
          '800109', '800111', '800111', '800111', '800113',
          '800113', '800113', '800114', '800114', '800114', 
          '800116', '800116', '800116', '800119', '800119',
          '800119', '800122', '800122', '800122', '800124', 
          '800125', '800125', '800127', '800127', '800127', 
          '800204', '800204', '800204')
type <- sample(c('cooperation', 'conflict'), 33,
               replace = TRUE)
important <- sample(c('important', 'not important'), 33,
                    replace = TRUE)

# combine them into a data.frame
dt <- data.frame(sender, target, time, type, important)

# create event sequence and order the data
dt <- eventSequence(datevar = dt$time, dateformat = "\%y\%m\%d", 
                    data = dt, type = "continuous", 
                    byTime = "daily", returnData = TRUE,
                    sortData = TRUE)

# create counting process data set (with null-events) - conditional logit setting
dts <- createRemDataset(dt, dt$sender, dt$target, dt$event.seq.cont, 
                          eventAttribute = dt$type, 
                          atEventTimesOnly = TRUE, untilEventOccurrs = TRUE, 
						  returnInputData = TRUE)
dtrem <- dts[[1]]
dt <- dts[[2]]
# manually sort the data set
dtrem <- dtrem[order(dtrem$eventTime), ]
# merge type-variable back in
dtrem$type <- dt$type[match(dtrem$eventID, dt$eventID)]

# calculate triad statistic
dtrem$triad <- triadStat(data = dtrem, time = dtrem$eventTime, 
                      sender = dtrem$sender, target = dtrem$target, 
                      eventvar = dtrem$eventDummy,
                      halflife = 2)

# calculate friend-of-friend statistic
dtrem$triad.fof <- triadStat(data = dtrem, time = dtrem$eventTime, 
                          sender = dtrem$sender, target = dtrem$target, 
                          halflife = 2, eventtypevar = dtrem$type, 
                          eventtypevalues = c("cooperation",
                                              "cooperation"), 
                          eventvar = dtrem$eventDummy)

# calculate friend-of-enemy statistic
dtrem$triad.foe <- triadStat(data = dtrem, time = dtrem$eventTime, 
                             sender = dtrem$sender, target = dtrem$target, 
                             halflife = 2, eventtypevar = dtrem$type, 
                             eventtypevalues = c("conflict",
                                                 "cooperation"), 
                             eventvar = dtrem$eventDummy)

# calculate enemy-of-friend statistic
dtrem$triad.eof <- triadStat(data = dtrem, time = dtrem$eventTime, 
                             sender = dtrem$sender, target = dtrem$target, 
                             halflife = 2, eventtypevar = dtrem$type, 
                             eventtypevalues = c("cooperation",
                                                 "conflict"), 
                             eventvar = dtrem$eventDummy)

# calculate enemy-of-enemy statistic
dtrem$triad.eoe <- triadStat(data = dtrem, time = dtrem$eventTime, 
                             sender = dtrem$sender, target = dtrem$target, 
                             halflife = 2, eventtypevar = dtrem$type, 
                             eventtypevalues = c("conflict",
                                                 "conflict"), 
                             eventvar = dtrem$eventDummy)
}
%\keyword{key}




