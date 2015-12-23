\name{inertiaStat}
\alias{inertiaStat}
\alias{inertia}
\title{Calculate inertia statistics}
%
\description{Calculate the endogenous network statistic \code{inertia} for relational event models. \code{inertia} measures the tendency for events to consist of the same sender and target (i.e. repeated events).}
\usage{
inertiaStat(data, time, sender, target, halflife,
    weight = NULL, 
    eventtypevar = NULL,
    eventtypevalue = "valuematch", 
    eventattributevar = NULL,
    eventattributevalue = "valuematch", 
    variablename = "inertia", 
    returnData = FALSE,
    showprogressbar = FALSE)
}
\arguments{
\item{data}{ A data frame containing all the variables.}

\item{time}{ Numeric variable that represents the event sequence. The variable has to be sorted in ascending order.}

\item{sender}{ A string (or factor or numeric) variable that represents the sender of the event.}

\item{target}{ A string (or factor or numeric) variable that represents the target of the event.}

\item{halflife}{ A numeric value that is used in the decay function. The vector of past events is weighted by an exponential decay function using the specified halflife. The halflife parameter determins after how long a period the event weigth should be halved. E.g. if \code{halflife = 5}, the weigth of an event that occured 5 units in the past is halved. Smaller halflife values give more importance to more recent events, while larger halflife values should be  used if time does not affect the sequence of events that much.} 

\item{weight}{ An optional numeric variable that represents the weigth of each event. If \code{weight = NULL} each event is given an event weigth of \code{1}.}

\item{eventtypevar}{ An optional variable that represents the type of the event. Use \code{eventtypevalue} to specify how the \code{eventtypevar} should be used to filter past events.}

\item{eventtypevalue}{ An optional value (or set of values) used to specify how  paste events should be filtered depending on their type. 
\code{eventtypevalue = "valuematch"} indicates that only past events that have the same type as the current event should be used to calculate the inertia statistic.
\code{eventtypevalue = "valuemix"} indicates that past and present events of specific types should be used for the inertia statistic. All the possible combinations of the eventtypevar-values will be used. E.g. if \code{eventtypevar} contains two unique values "a" and "b", 4 inertia statistics will be  calculated. The first variable calculates the inertia effect where the present  event is of type "a" and all the past events are of type "b". The next variable  calculates inertia for present events of type "b" and past events of type "a". Additionally, a variable is calculated, where present events as well as past events are of type "a" and a fourth variable calculates inertia for events with type "b" (i.e. valuematch on value "b").
\code{eventtypevalue = c(.., ..)} is similar to the \code{"nodmix"}-option, all different combinations of the values specified in \code{eventtypevalue} are used to create inertia statistics.}

\item{eventattributevar}{ An optional variable that represents an attribute of the event. It can be a sender attribute, a target attribute, time or dyad attribute. Use \code{eventattributevalue} to specify how the \code{eventattributevar} should be used.}

\item{eventattributevalue}{ An optional value (or set of values) used to specify how paste events should be filtered depending on their attribute. Similar to \code{eventtypevalue}.
\code{eventattributevalue = "valuematch"} indicates that only past events that have the same attribute should be used to calculate the inertia statistic.
\code{eventtypevalue = "valuemix"} indicates that past and present events of  specific attributes should be used for the inertia statistic. All the possible combinations of the eventattributevar-values will be used. 
\code{eventtypevalue = c(.., ..)} similar to the \code{"valuemix"}-option, all different combinations of the values specified in \code{eventtypevalue} are used to create inertia statistics.}

\item{variablename}{ An optional value (or values) with the name the inertia statistic variable should be given. To be used if \code{returnData = TRUE} or multiple inertia statistics are calculated.}

\item{returnData}{ \code{TRUE/FALSE}. Set to \code{FALSE} by default. The new variable(s) are bound directly to the \code{data.frame} provided and the data frame is returned in full.}

\item{showprogressbar}{\code{TRUE/FALSE}. To be implemented.}
}
\details{

The \code{inertiaStat()}-function calculates an endogenous statistic that  measures whether events have a tendency to be repeated with the same sender and target over the entire event sequence.

The effect is calculated as follows.

\deqn{G_t = G_t(E) = (A, B, w_t), }{G_t = G_t(E) = (A, B, w_t),}

\eqn{G_t} represents the network of past events and includes all events \eqn{E}. These events consist
each of  a sender \eqn{a \in A}{a in A} and a target \eqn{b \in B}{b in B} and a weight function \eqn{w_t}:

\deqn{  w_t(i, j) = \sum_{e:a = i, b = j} | w_e | \cdot e^{-(t-t_e)\cdot\frac{ln(2)}{T_{1/2}}} \cdot \frac{ln(2)}{T_{1/2}}, }{ w_t(i, j) = \sum_{e:a = i, b = j} | w_e | * exp^{-(t-t_e)* (ln(2)/T_{1/2})} * (ln(2)/T_{1/2}),}

where \eqn{w_e} is the event weight (usually a constant set to 1 for each event), \eqn{t} is the current event time, \eqn{t_e} is the past event time and \eqn{T_{1/2}} is a halflife parameter.

For the inertia effect, the past events \eqn{G_t} are filtered to include only events
where the senders and targets are identical to the current sender and target.

\deqn{inertia(G_t , a , b) = w_t(a, b)}

An exponential decay function is used to model the effect of time on the endogenous statistics. Each past event that contains the same sender and target and fulfills additional filtering options specivied via event type or event attributes is weighted with an exponential decay. The further apart the past event is from the present event, the less weight is given to this event. The halflife parameter in the \code{inertiaStat()}-function determins at which rate the weights of past events should be reduced. 

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
sender <- c('TUN', 'NIR', 'NIR', 'TUR', 'TUR', 'USA', 'URU', 
            'IRQ', 'MOR', 'BEL', 'EEC', 'USA', 'IRN', 'IRN', 
            'USA', 'AFG', 'ETH', 'USA', 'SAU', 'IRN', 'IRN',
            'ROM', 'USA', 'USA', 'PAN', 'USA', 'USA', 'YEM', 
            'SYR', 'AFG', 'NAT', 'NAT', 'USA')
target <- c('BNG', 'ZAM', 'JAM', 'SAU', 'MOM', 'CHN', 'IRQ', 
            'AFG', 'AFG', 'EEC', 'BEL', 'ITA', 'RUS', 'UNK',
            'IRN', 'RUS', 'AFG', 'ISR', 'ARB', 'USA', 'USA',
            'USA', 'AFG', 'IRN', 'IRN', 'IRN', 'AFG', 'PAL',
            'ARB', 'USA', 'EEC', 'BEL', 'PAK')
time <- c('800107', '800107', '800107', '800109', '800109', 
          '800109', '800111', '800111', '800111', '800113',
          '800113', '800113', '800114', '800114', '800114', 
          '800116', '800116', '800116', '800119', '800119',
          '800119', '800122', '800122', '800122', '800124', 
          '800125', '800125', '800127', '800127', '800127', 
          '800204', '800204', '800204')
type <- sample(c('cooperation', 'conflict'), 33,
               replace = TRUE)

# combine them into a data.frame
dt <- data.frame(sender, target, time, type)

# create event sequence and order the data
dt <- eventSequence(datevar = dt$time, dateformat = '\%y\%m\%d', 
                    data = dt, type = 'continuous', 
                    byTime = "daily", returnData = TRUE,
                    sortData = TRUE)

# calculate inertia statistics
dt$inertia <- inertiaStat(data = dt, time = dt$event.seq.cont, 
                          sender = dt$sender, target = dt$target, 
                          halflife = 2, returnData = FALSE)

# plot inertia over time
library('ggplot2')
ggplot(dt, aes ( event.seq.cont, inertia) ) +
  geom_point()+ geom_smooth() 

# inertia with typematch (e.g. for 'cooperation' events only 
# count past 'cooperation' events)             
dt$inertia.tm <- inertiaStat(data = dt, 
                             time = dt$event.seq.cont, 
                             sender = dt$sender, 
                             target = dt$target, 
                             halflife = 2, 
                             eventtypevar = dt$type, 
                             eventtypevalue = 'valuematch',
                             returnData = FALSE)

# inertia with valuemix: for each combination of types
# in the eventtypevar, create a variable
dt <- inertiaStat(data = dt, time = dt$event.seq.cont, 
                  sender = dt$sender, target = dt$target, 
                  halflife = 2, eventtypevar = dt$type, 
                  eventtypevalue = 'valuemix',
                  returnData = TRUE)
}
%\keyword{key}




