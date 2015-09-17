\name{get.fourCycle.stat}
\alias{get.fourCycle.stat}
\alias{fourCycle}
\alias{cycle}
\alias{path-dependency}
\title{Calculate four cycle statistics}
%
\description{Calculate the endogenous network statistic \empth{fourCycle} that 
measures the tendency for events close four cycles. Two-mode networks.}
\usage{
get.fourCycle.stat(data, time, sender, target, halflife, weight = NULL, 
  eventtypevar = NULL, eventtypevalue = "standard", eventattributevar = NULL, 
  eventattributeAB = NULL, eventattributeAJ = NULL, eventattributeIB = NULL, 
  eventattributeIJ = NULL, variablename = "fourCycle", returnData = TRUE,
  showprogressbar = FALSE, ...)
}
\arguments{
\item{data}{ An optional data frame containing all the variables.}
\item{time}{ Numeric variable that represents the event sequence. The variable
has to be sorted in ascending order.}
\item{sender}{ A string (or factor or numeric) variable that represents the sender
of the event.}
\item{target}{ A string (or factor or numeric) variable that represents the target
of the event.}
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
paste events should be filtered depending on their type. "`\code{standard}"',
"`\code{positive}"' or "`\code{negative}"' may be used. Default set to 
"`\code{standard}"'. "`\code{standard}"' referrs to closing four cylces where
the type of the events is irrelevant. "`\code{positive}"' closing four cycles 
can be classified as reciprocity via the second mode. It indicates whether 
senders have a tendency to reciprocate or show support by engaging in targets that
close a four cycle between two senders. "`\code{negative}"' closing four cycles
represent opposition between two senders, where the current event is more likely if
the two senders have opposed each other in the past.}
\item{eventattributevar}{ An optinoal variable that represents an attribute of the event.
It can be a sender attribute, a target attribute, time or dyad attribute.
Use \code{eventattributeAB}, \code{eventattributeAJ}, \code{eventattributeIB} or
\code{eventattributeIJ} to specify how the \code{eventattributevar} should be used.
}
\item{eventattributeAB}{ An optional value used to specify how 
paste events should be filtered depending on their attribute. Each distinct edge
that form a four cycle can be filtered. \code{eventattributeAB} refers to the
current event. \code{eventattributeAJ} refers to the event involving the current
sender and target $j$ that has been used by the current as well as the second 
actor in the past. \code{eventattributeIB} refers to the event involving the
second sender and the current target. \code{eventattributeIJ} filters events
that involve the second sender and the second target. See the four cycle formula
in the \link{details} for more information.}
\item{eventattributeAJ}{ see \link{eventattributeAB}.}
\item{eventattributeIB}{see \link{eventattributeAB}.}
\item{eventattributeIJ}{see \link{eventattributeAB}.}
\item{variablename}{ An optional value (or values) with the name the four cycle
statistic variable should be given. To be used if \code{returnData = TRUE}.}
\item{returnData}{ \code{TRUE/FALSE}. Set to \code{TRUE} by default. The new 
variable(s) are bound directly to the \code{data.frame} provided and the
data frame is returned in full.}
\item{showprogressbar}{\code{TRUE/FALSE}. To be implemented.}
}
\details{

The \code{get.fourCycle.stat()}-function calculates an endogenous statistic that 
measures whether events have a tendency to form four cycles.

function.

where $ $, 

with decay function $ $.

An exponential decay function is used to model the effect of time on the
endogenous statistics.
The further apart the past event is
from the present event, the less weight is given to this event. The halflife
parameter in the \code{get.fourCycle.stat()}-function determins at which rate
the weights of past events should be reduced. Therefore, if the one (or more) of
the three events in the four cycle have ocurred further in the past, less weight
is given to this four cycle because it becomes less likely that the two senders
reacted to each other in the way the four cycle assumes.

The \code{eventtypevar}- and \code{eventattributevar}-options help filter the
past events more specifically. How they are filtered depends on the
\code{eventtypevalue}- and \code{eventattributevalue}-option.

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




