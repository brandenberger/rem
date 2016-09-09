\name{similarityStat}
\alias{similarityStat}
\alias{similarity}
\title{Calculate similarity statistics}
%
\description{Calculate the endogenous network statistic \code{similarity} for relational event models. \code{similarityStat} measures the tendency for senders to adapt their behavior to that of their peers.}
\usage{
similarityStat(data, time, sender, target, 
    senderOrTarget = 'sender', 
    whichSimilarity = NULL, 
    halflifeLastEvent = NULL, 
    halflifeTimeBetweenEvents = NULL, 
    eventtypevar = NULL, 
    eventfiltervar = NULL, 
    eventfiltervalue = NULL, 
    eventvar = NULL,
    variablename = 'similarity',
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

\item{senderOrTarget}{ \code{sender} or \code{target}. Indicates on which variable (sender or target) the similarity should be calculated on. Sender similarity measures how many targets the current sender has in common with other senders who used the same targets in the past. Target similarity measures how many senders have used the current target as well as another target that the current sender used in the past.}

\item{whichSimilarity}{ \code{"total"} or \code{"average"}. Indicates how the variable should be aggregated. \code{"total"} counts the number of similar events there are in the past event history.  \code{"average"} divides the count of similar events by the number of senders or the number of  targets, depending on which mode of similarity is chosen.}

\item{halflifeLastEvent}{ A numeric value that is used in the decay function. The vector of past events is weighted by an exponential decay function using the specified halflife. The halflife parameter determines after how long a period the event weight should be halved. For sender similarity: The halflife determines the weight of the count of targets that two actors have in common. The further back the second sender was active, the less weight is given the similarity between this sender and the current sender. For target similarity: The halflife determines the weight of the count of targets that have used both been used by other senders in the past. The longer ago the current sender engaged in an event with the  other target, the less weight is given the count.}

\item{halflifeTimeBetweenEvents}{ A numeric value that is used in the decay function. Instead of counting each past event for the similarity statistic, each event is reduced depending on the time that passed between the current event and the past event. For sender similarity: Each target that two actors have in common is weighted by the time that passed between the two events. For target similarity: Each sender that two targets have in common is weighted by the time that passed between the two events. }

\item{eventtypevar}{ An optional dummy variable that represents the type of the event. If specified, only past events are considered for the count that reflect the same type as the current event (typematch).}

\item{eventfiltervar}{ An optional variable that filters past events by the \code{eventfiltervalue} specified. }

\item{eventfiltervalue}{ A string that represents an event attribute by which all past events have to be filtered by.}

\item{eventvar}{ An optional dummy variable with 0 values for null-events and 1 values for true events. If the \code{data} is in the form of counting process data, use the \code{eventvar}-option to specify which variable contains the 0/1-dummy for event occurrence. If this variable is not specified, all events in the past will be considered for the calulation of the similarity statistic, regardless if they occurred or not (= are null-events). Misspecification could result in grievous errors in the calculation of the network statistic.}

\item{variablename}{ An optional value (or values) with the name the similarity statistic variable should be given. To be used if \code{returnData = TRUE}.}

\item{returnData}{ \code{TRUE/FALSE}. Set to \code{FALSE} by default. The new  variable(s) are bound directly to the \code{data.frame} provided and the data frame is returned in full.}

\item{showprogressbar}{\code{TRUE/FALSE}. To be implemented.}

\item{inParallel}{ \code{TRUE/FALSE}. An optional boolean to specify if the loop should be run in parallel.}

\item{cluster}{ An optional numeric or character value that defines the cluster. By specifying a single number, the cluster option uses the provided number of nodes to parallellize. By specifying a cluster using the \code{makeCluster}-command in the \code{doParallel}-package, the loop can be run on multiple nodes/cores. E.g., \code{cluster = makeCluster(12, type="FORK")}.}

}
\details{

The \code{similiarityStat()}-function calculates an endogenous statistic that measures whether sender (or targets) have a tendency to cluster together. Tow distinct types of similarity measures can be calculated: sender similarity or target similarity. 

Sender similarity: How many targets does the current sender have in common with senders who used the current target in the past? How likely is it that two senders are alike?

The function proceeds as follows: 
\enumerate{
	\item First it filters out all the targets that the present sender \eqn{a} used in the past 
	\item Next it filters out all the senders that have also used the current target \eqn{b}
	\item For each of the senders found in (2) it compiles a list of targets that this sender has used in the past
	\item For each of the senders found in (2) it cross-checks the two lists generated in (1) and (3) and count how many targets the two senders have in common. 
}

% Procedure of the sender similarity calcuation:
% list - Filter out all the targets that $a$ has used in the past
%     - Filter out all the senders that also used the current target $b$
%     - Filter out events that involve senders that have used $b$ in the past
%     - Check how many targets the current sender and the past sender (that also
%     used $b$) have in common

Target similarity: How many senders have used the same two concepts that the current sender has used (in the past and is currently using)? For each target that the current sender has used in the past, how many senders have also used these past targets as well as the current target? How likely is it that two targets are used together?

The function proceeds as follows: 
\enumerate{
	\item First filter out all the targets that the current sender \eqn{a} has used in the past
	\item Next it filters out all the senders that have also used the current target \eqn{b}
	\item For each target found in (1) it compiles a list of senders that have also used this target in the past
	\item For each target found in (1) it cross-checks the list of senders that have used \eqn{b} (found under (2)) and the list of senders that also used one other target that \eqn{a} used (found under (3))
	}

% list - - Filter out all the targets that $a$ has used in the past
%     - Filter out all the senders that also used the current target $b$
%     - Filter out events that involve a target that $a$ has used in the past
%     - Check how many senders have used the current topic as well as one other
%     topic used by $a$ in the past

Two decay functions may be used in the calculation of the similarity score for each event. 

%% Additional info needed here?

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




