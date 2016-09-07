\name{degreeStat}
\alias{degreeStat}
\alias{degree}
\alias{outdegree}
\alias{indegree}
\title{Calculate (in/out)-degree statistics}
%
\description{Calculate the endogenous network statistic \code{indegree/outdegree} for relational event models. \code{indegree/outdegree} measures the senders' tendency to be involved in events (sender activity, sender out- or indegree) or the tendency of events to surround a specific target (target popularity, target in- or outdegree)}
\usage{
degreeStat(data, time, degreevar, halflife, 
    weight = NULL, 
    eventtypevar = NULL, 
    eventtypevalue = "valuematch", 
    eventfiltervar = NULL,
    eventfiltervalue = NULL, 
    eventvar = NULL,
    degreeOnOtherVar = NULL,
    variablename = "degree", 
    returnData = FALSE, 
    showprogressbar = FALSE, 
    inParallel = FALSE, cluster = NULL)
}
\arguments{
\item{data}{ A data frame containing all the variables.}
\item{time}{ Numeric variable that represents the event sequence. The variable has to be sorted in ascending order.}

\item{degreevar}{ A string (or factor or numeric) variable that represents the sender or target of the event. The degree statistic will calculate how often in the past, a given sender or target has been active by counting the number of events in the past where the \code{degreevar} is repeated. See \code{details} for more information on which variable to chose as \code{degreevar} for one- and two-mode networks.}

\item{halflife}{ A numeric value that is used in the decay function. The vector of past events is weighted by an exponential decay function using the specified halflife. The halflife parameter determines after how long a period the event weight should be halved. E.g. if \code{halflife = 5}, the weight of an event that occurred 5 units in the past is halved. Smaller halflife values give more importance to more recent events, while larger halflife values should be used if time does not affect the sequence of events that much.}

\item{weight}{ An optional numeric variable that represents the weight of each event. If \code{weight = NULL} each event is given an event weight of \code{1}.}

\item{eventtypevar}{ An optional variable that represents the type of the event. Use \code{eventtypevalue} to specify how the \code{eventtypevar} should be used to filter past events.}

\item{eventtypevalue}{ An optional value (or set of values) used to specify how paste events should be filtered depending on their type. 
\code{eventtypevalue = "valuematch"} indicates that only past events that have the same type should be used to calculate the degree statistic.
\code{eventtypevalue = "valuemix"} indicates that past and present events of specific types should be used for the degree statistic. All the possible combinations of the eventtypevar-values will be used. E.g. if \code{eventtypevar} contains two unique values "a" and "b", 4 degree statistics will be  calculated. The first variable calculates the degree effect where the present  event is of type "a" and all the past events are of type "b". The next variable  calculates the degree statistic for present events of type "b" and past events of type "a". Additionally, a variable is calculated, where present events as well as past  events are of type "a" and a fourth variable calculates the degree statistic for events with type "b" (i.e. valuematch on value "b").
\code{eventtypevalue = c("..", "..")} is similar to the \code{"nodemix"}-option, all different combinations of the values specified in \code{eventtypevalue} are used to create the degree statistics.}

\item{eventfiltervar}{ An optional numeric/character/or factor variable for each event. If \code{eventfiltervar} is specified, \code{eventfiltervalue} has to be provided as well.}

\item{eventfiltervalue}{ An optional character string that represents the value for which past events should be filtered. To filter the current events, use \code{eventtypevar}.}

\item{eventvar}{ An (optional) dummy variable with 0 values for null-events and 1 values for true events. If the \code{data} is in the form of counting process data, use the \code{eventvar}-option to specify which variable contains the 0/1-dummy for event occurrence. If this variable is not specified, all events in the past will be considered for the calulation of the degree statistic, regardless if they occurred or not (= are null-events).}

\item{degreeOnOtherVar}{ A string (or factor or numeric) variable that represents the sender or target of the event. It can be used to calculate target-outdegree or sender-indegree statistics in one-mode networks. For the sender indegree statistic, fill the sender variable into the \code{degreevar} and the target variable into the \code{degree.on.other.var}. For the target-outdegree statistic, fill the target variable into the \code{degreevar} and the sender variable into the \code{degree.on.other.var}.}

\item{variablename}{ An optional value (or values) with the name the degree statistic variable should be given. Default "degree" is used. To be used if \code{returnData = TRUE} or multiple degree statistics are calculated.}

\item{returnData}{ \code{TRUE/FALSE}. Set to \code{FALSE} by default. The new variable(s) are bound directly to the \code{data.frame} provided and the data frame is returned in full.}

\item{showprogressbar}{ \code{TRUE/FALSE}. Can only be set to TRUE if the function is not run in parallel.}

\item{inParallel}{ \code{TRUE/FALSE}. An optional boolean to specify if the loop should be run in parallel.}

\item{cluster}{ An optional numeric or character value that defines the cluster. By specifying a single number, the cluster option uses the provided number of nodes to parallellize. By specifying a cluster using the \code{makeCluster}-command in the \code{doParallel}-package, the loop can be run on multiple nodes/cores. E.g., \code{cluster = makeCluster(12, type="FORK")}.}

}

\details{

The \code{degreeStat()}-function calculates an endogenous statistic that measures whether events have a tendency to include either the same sender or the same target over the entire event sequence. 

The effect is calculated as follows.

\deqn{G_t = G_t(E) = (A, B, w_t), }{G_t = G_t(E) = (A, B, w_t),}

\eqn{G_t} represents the network of past events and includes all events \eqn{E}. These events consist
each of  a sender \eqn{a \in A}{a in A} and a target \eqn{b \in B}{b in B} (in one-mode networks \eqn{A = B}{A = B}) and a weight function \eqn{w_t}:

\deqn{  w_t(i, j) = \sum_{e:a = i, b = j} | w_e | \cdot e^{-(t-t_e)\cdot\frac{ln(2)}{T_{1/2}}} \cdot \frac{ln(2)}{T_{1/2}}, }{ w_t(i, j) = \sum_{e:a = i, b = j} | w_e | * exp^{-(t-t_e)* (ln(2)/T_{1/2})} * (ln(2)/T_{1/2}),}

where \eqn{w_e} is the event weight (usually a constant set to 1 for each event), \eqn{t} is the current event time, \eqn{t_e} is the past event time and \eqn{T_{1/2}} is a halflife parameter.

For the degree effect, the past events \eqn{G_t} are filtered to include only events
where the senders or targets are identical to the current sender or target.

\deqn{sender-outdegree(G_t , a , b) = \sum_{j \in B}  w_t(a, j)}{sender-outdegree(G_t , a , b) = \sum_{j in B} w_t(a, j)}

\deqn{target-indegree(G_t , a , b) = \sum_{i \in A}  w_t(i, b)}{target-indegree(G_t , a , b) = \sum_{i in A} w_t(i, b)}

\deqn{sender-indegree(G_t , a , b) = \sum_{i \in A}  w_t(i, a)}{sender-indegree(G_t , a , b) = \sum_{i in A} w_t(i, a)}

\deqn{target-outdegree(G_t , a , b) = \sum_{j \in B}  w_t(b, j)}{target-outdegree(G_t , a , b) = \sum_{j in B} w_t(b, j)}

Depending on whether the degree statistic is measured on the sender variable or the target variable, either activity or popularity effects are calculated. 

For one-mode networks: Four distinct statistics can be calculated: sender-indegree, sender-outdegree, target-indegree or target-outdegree. The sender-indegree measures how often the current sender was targeted by other senders in the past (i.e. how popular were current senders). The sender-outedegree measures how often the current sender was involved in an event, where they were also marked as sender (i.e. how active the current sender has been in the past). The target-indegree statistic measures how often the current targets were targeted in the past (i.e. how popular were current targets). And the target-outdegree measures how often the current targets were senders in the past (i.e. how active were current targets in the past). 

For two-mode networks: Two distinct statistics can be calculated: sender-outdegree and target-indegree. Sender-outdegree measures how often the current sender has been involved in an event in the past (i.e. how active the sender has been up until now). The target-indegree statistic measures how often the current target has been involved in an event in the past (i.e. how popular a given target has been before the current event).

An exponential decay function is used to model the effect of time on the endogenous statistics. Each past event that contains the same sender or the same target (depending on the variable specified in \code{degreevar}) and fulfills additional filtering options (specified via event type or event attributes) is weighted with an exponential decay. The further apart the past event is from the present event, the less weight is given to this event. The halflife parameter in the \code{degreeStat()}-function determines at which rate the weights of past events should be reduced. 

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
}
%\keyword{key}




