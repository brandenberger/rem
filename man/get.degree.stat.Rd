\name{get.degree.stat}
\alias{get.degree.stat}
\alias{degree}
\alias{outdegree}
\alias{indegree}
\alias{path-dependency}
\title{Calculate (in/out)-degree statistics}
%
\description{Calculate the endogenous network statistic \empth{indegree/outdegree} 
that measures the senders' tendency to be involved in events (sender activity,
sender out- or indegree) or the tendency of events surrounding a specific target
(target popularity, target in- or outdegree)}
\usage{
get.degree.stat(data, time, degreevar, halflife, weight = NULL,
  eventtypevar = NULL, eventtypevalue = "valuematch", eventattributevar = NULL, 
  eventattributevalue = "valuematch", degree.on.other.var = NULL,
  variablename = "degree", returnData = TRUE, showprogressbar = FALSE, ...)
}
\arguments{
\item{data}{ An optional data frame containing all the variables.}
\item{time}{ Numeric variable that represents the event sequence. The variable
has to be sorted in ascending order.}
\item{degreevar}{ A string (or factor or numeric) variable that represent the
sender or target of the event. The degree statistic will calculate how often 
in the past a given sender or target has been active by counting the number of
events in the past where the \code{degreevar} is repeated. See \code{details} for
more information on which variable to chose as \code{degreevar} for
one- and two-mode networks.}
\item{halflife}{ A numeric value that is used in the decay function.
The vector of past events is weighted by an exponential decay function using the
specified halflife. The halflife parameter determins after how long a period the
event weigth should be halved. E.g. if \code{halflife = 5}, the weigth of an
event that occured 5 units in the past is halved. Smaller halflife values give
more importance to more recent events, while larger halflife values should be 
used if time does not affect the sequence of events that much.}
\item{weight}{ An optional numeric variable that represents the weigth of each
event. If \code{weight = NULL} each event is given an event weigth of \code{1}.
}
\item{eventtypevar}{ An optional variable that represents the type of the event.
Use \code{eventtypevalue} to specify how the \code{eventtypevar} should be used
to filter past events.
}
\item{eventtypevalue}{ An optional value (or set of values) used to specify how 
paste events should be filtered depending on their type. 
\code{eventtypevalue = "`valuematch"'} indicates that only past events that have the
same type should be used to calculate the degree statistic.
\code{eventtypevalue = "`valuemix"'} indicates that past and present events of 
specific types should be used for the degree statistic. All the possible
combinations of the eventtypevar-values will be used. E.g. if \code{eventtypevar}
contains three unique values "`a"' and "`b"', 4 degree statistics will be 
calculated. The first variable calculates the degree effect where the present 
event is of type "`a"' and all the past events are of type "`b"'. The next variable 
calculates the degree statistic for present events of type "`b"' and past events of type "`a"'.
Additionally, a variable is calculated, where present events as well as past 
events are of type "`a"' and a fourth variable calculates the degree statistic for events with
type "`b"' (i.e. valuematch on value "`b"').
\code{eventtypevalue = c("..", "`.."')} similar to the "`\code{nodmix}"'-option, all 
different combinations of the values specified in \code{eventtypevalue} are used
to create the degree statistics.
}
\item{eventattributevar}{ An optinoal variable that represents an attribute of the event.
It can be a sender attribute, a target attribute, time or dyad attribute.
Use \code{eventattributevalue} to specify how the \code{eventattributevar} should be used.
}
\item{eventattributevalue}{ An optional value (or set of values) used to specify how 
paste events should be filtered depending on their attribute. Similar to \code{eventtypevalue}.
\code{eventattributevalue = "`valuematch"'} indicates that only past events that have the
same attribute should be used to calculate the degree statistic.
\code{eventtypevalue = "`valuemix"'} indicates that past and present events of 
specific attributes should be used for the degree statistic. All the possible
combinations of the eventattributevar-values will be used. 
\code{eventtypevalue = c("..", "`.."')} similar to the "`\code{valuemix}"'-option, all 
different combinations of the values specified in \code{eventtypevalue} are used
to create degree statistics.
}
\item{degree.on.other.var}{ A string (or factor or numeric) variable that 
represents the sender or target of the event. It can be used to calculate
target-outdegree or sender-indegree statistics in one-mode networks. For the
sender indegree statistic, fill the sender variable into the \code{degreevar} and
the target variable into the \code{degree.on.other.var}. For the target-outdegree
statistic, fill the target variable into the \code{degreevar} and the sender
variable into the \code{degree.on.other.var}.}
\item{variablename}{ An optional value (or values) with the name the degree
statistic variable should be given. Default 'degree' is used. To be used if
\code{returnData = TRUE} or multiple degree statistics are calculated.}
\item{returnData}{ \code{TRUE/FALSE}. Set to \code{TRUE} by default. The new 
variable(s) are bound directly to the \code{data.frame} provided and the
data frame is returned in full.}
\item{showprogressbar}{\code{TRUE/FALSE}. To be implemented.}
}
\details{

The \code{get.degree.stat()}-function calculates an endogenous statistic that
measures whether events have a tendency to include either the same sender or the
same target over the entire event sequence. The statistic follows the function

$ $, 

where $..$,

with decay function $ $.

An exponential decay function is used to model the effect of time on the
endogenous statistics. Each past event that contains the same sender or the same
target (depending on the variable specified in \code{degreevar}) and fulfills
additional filtering options (specivied via event type or event attributes)
is weighted with an exponential decay. The further apart the past event is
from the present event, the less weight is given to this event. The halflife
parameter in the \code{get.degree.stat()}-function determins at which rate
the weights of past events should be reduced. 

The \code{eventtypevar}- and \code{eventattributevar}-options help filter the
past events more specifically. How they are filtered depends on the
\code{eventtypevalue}- and \code{eventattributevalue}-option.

Depending on whether the degree statistic is measured on the sender variable or
the target variable, either activity or popularity effects are calculated. 

For one-mode networks: Four distinct statistics can be calculated: sender-indegree, 
sender-outdegree, target-indegree or target-outdegree. The sender-indegree measures
how often the current sender was targeted by other senders in the past (i.e. how
popular were current senders). The
sender-outedegree measures how often the current sender was involved in an event,
where they were also marked as sender (i.e. how active the current sender has 
been in the past). The target-indegree statistic measures 
how often the current targets were targeted in the past (i.e. how popular
were current targets). And the target-outdegree
measures how often the current targets were senders in the past (i.e. how active
were current targets in the past). 

For two-mode networks: Two distinct statistics can be calculated: sender-outdegree
and target-indegree. Sender-outdegree measures how often the current sender has
been involved in an event in the past (i.e. how active the sender has been up
until now). The target-indegree statistic measures how often the current target has
been involved in an event in the past (i.e. how popular a given target has been
before the current event).

}
\value{

}
\references{

}
\note{

}
\author{
Laurence Brandenberger \email{laurence.brandenberger@eawag.ch}
}
\seealso{
\code{\link{rem-package}} 
}
\examples{

}
\keyword{key}




