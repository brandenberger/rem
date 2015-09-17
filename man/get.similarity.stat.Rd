\name{get.similarity.stat}
\alias{get.similarity.stat}
\alias{similarity}
\alias{path-dependency}
\title{Calculate similarity statistics}
%
\description{Calculate the endogenous network statistic \code{similarity} that 
measures the tendency for senders to adapt their behavior to that of their peers.}
\usage{
get.similarity.stat(data, time, sender, target, senderOrTarget = "sender",
  whichSimilarity = NULL, halflife.last.event = NULL, 
  halflife.time.between.events = NULL, eventtypevar = NULL, 
  eventattributevar = NULL, eventattributevalue = NULL, 
  variablename = "similarity", returnData = TRUE, showprogressbar = FALSE, ...)
}
\arguments{
\item{data}{ An optional data frame containing all the variables.}
\item{time}{ Numeric variable that represents the event sequence. The variable
has to be sorted in ascending order.}
\item{sender}{ A string (or factor or numeric) variable that represents the sender
of the event.}
\item{target}{ A string (or factor or numeric) variable that represents the target
of the event.}
\item{senderOrTarget}{ "`\code{sender}"' or "`\code{target}"'. Indicates on
which variable (sender or target) the similarity should be calculated.}
\item{whichSimilarity}{ "`\code{total}"' or "`\code{average}"'. Indicates how
the variable should be aggregated. "`\code{total}"' counts the number of
similar events there are in the past event history.  "`\code{average}"'
divides the count of similar events by the number of senders or the number of 
targets, depending on which mode of similarity is chosen.
}
\item{halflife.last.event}{ A numeric value that is used in the decay function. 
The vector of past events is weighted by an exponential decay function using the
specified halflife. The halflife parameter determins after how long a period the
event weigth should be halved. For sender similarity: The halflife determins the
weight of the count of targets that two actors have in common. The further back 
the second sender was active, the less weight is given the similarity between this
sender and the current sender. For target similarity: The halflife determins
the weight of the count of targets that have used both been used by other senders
in the past. The longer ago the current sender engaged in an event with the 
other target, the less weight is given the count.}
\item{halflife.time.between.events}{ A numeric value that is used in the decay 
function. Instead of counting each past event for the similarity statistic, each
event is reduced depending on the time that passed between the current event
and the past event. For sender similarity: Each target that two actors have in
common is weighted by the time that passed between the two events. For target 
similarity: Each sender that two targets have in common is weighted by the time 
that passed between the two events. }
\item{eventtypevar}{ An optional dummy variable that represents the type of the 
event. If specified, only past events are considered for the count that reflect the 
same type as the current event.
}
\item{eventattributevar}{ An optional variable that filters past events by the
\code{eventattributevalue} specified. }
\item{eventattributevalue}{ A string that represents an event attribute by which
all past events have to be filtered by.}
\item{variablename}{ An optional value (or values) with the name the similarity
statistic variable should be given. To be used if \code{returnData = TRUE}.}
\item{returnData}{ \code{TRUE/FALSE}. Set to \code{TRUE} by default. The new 
variable(s) are bound directly to the \code{data.frame} provided and the
data frame is returned in full.}
\item{showprogressbar}{\code{TRUE/FALSE}. To be implemented.}
}
\details{

The \code{get.inertia.stat()}-function calculates an endogenous statistic that 
measures whether sender (or targets) have a tendency to cluster together. Tow
distinct types of similarity measures can be calculated: sender similarity or
target similarity. 

Sender similarity: How many targets does the current sender have in common with
senders who used the current target in the past? How likely is it that two 
senders are alike?

$ $, 
 
where $ $,

Procedure of the sender similarity calcuation: 
list - Filter out all the targets that $a$ has used in the past
    - Filter out all the senders that also used the current target $b$
    - Filter out events that involve senders that have used $b$ in the past
    - Check how many targets the current sender and the past sender (that also
    used $b$) have in common

Target similarity: How many senders have used the same two concepts that
the current sender has used (in the past and is currently using)? For each
target that the current sender has used in the past, how many senders have also
used these past targets as well as the current target? How likely is it that two
targets are used together?

list - - Filter out all the targets that $a$ has used in the past
    - Filter out all the senders that also used the current target $b$
    - Filter out events that involve a target that $a$ has used in the past
    - Check how many senders have used the current topic as well as one other
    topic used by $a$ in the past

Two decay functions may be used in the calculation of the similarity score for 
each event. 



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




