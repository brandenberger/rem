\name{createRemDataset}
\alias{createRemDataset}
\title{Create REM data set with dynamic risk sets}
%
\description{The function creates counting process data sets with dynamic risk sets for relational event models. For each event in the event sequence, null-events are generated and represent possible events that could have happened at that time but did not. A data set with true and null-events is returned with an event dummy for whether the event occurred or was simply possible (variable \code{eventdummy}). The returned data set also includes a variable \code{eventTime} which represents the true time of the reported event.}
\usage{
createRemDataset(data, sender, target, eventSequence, 
	eventAttribute = NULL, time = NULL, 
	start = NULL, startDate = NULL, 
	end = NULL, endDate = NULL, 
	timeformat = NULL,
	atEventTimesOnly = TRUE, untilEventOccurrs = TRUE,
	includeAllPossibleEvents = FALSE, possibleEvents = NULL, 
	returnInputData = FALSE)
}
\arguments{
\item{data}{ A data frame containing all the events.}
\item{sender}{  A string (or factor or numeric) variable that represents the sender of the event.}
\item{target}{  A string (or factor or numeric) variable that represents the target of the event.}
\item{eventSequence}{ Numeric variable that represents the event sequence. The variable has to be sorted in ascending order.}
\item{eventAttribute}{ An optional variable that represents an attribute to an event. Repeated events affect the construction of the counting process data set. Use the \code{eventAttribute}-variable to specify the uniqueness of an event. If \code{eventAttribute = NULL}, events are defines as sender-target nodes only.}
\item{time}{  An optional date variable that represents the date an event took place. The variable is used if \code{startDate} or \code{endDate} are specified. \code{timeformat} should be used to specify which format the date variable is in, in case it was not yet converted to a Date-variable.}
\item{start}{  An optional numeric variable that indicates at which point in the event sequence a specific event was at risk. The variable has to be numerical and correspond to the variable \code{eventSequence}. If this option is used, each event in the event data set will be considered at risk from the specified value onwards. If it is not specified, \code{start} is defined as the first value in the event sequence. In case of repeated events, the start-value for each duplicated event is one event-unit after the last such event.}
\item{startDate}{ An optional date variable that represents the date an event started being at risk. \code{timeformat} should be used to specify which format the date variable is in, incase it was not yet converted to a Date-variable.}
\item{end}{ An optional numeric variable that indicates at which point in the event sequence a specific event stopped being at risk. The variable has to be numerical and correspond to the variable \code{eventSequence}. If this option is used, each event in the event data set will be considered at risk until the specified value.}
\item{endDate}{ An optional date variable that represents the date an event stoped being at risk. \code{timeformat} should be used to specify which format the date variable is in, incase it was not yet converted to a Date-variable.}
\item{timeformat}{ A character string indicating the format of the \code{datevar}. see \code{\link{as.Date}}}
\item{atEventTimesOnly}{ \code{TRUE/FALSE}. Boolean option for continuous event sequences. If \code{atEventTimesOnly = TRUE}, null-events are only created at times, when an event occurred. If \code{atEventTimesOnly = FALSE}, null-events are created on each event-unit from \code{min(eventSequence):max(eventSequence)}.
For instance: Given an event sequence with three events at \code{c(1, 4, 6)}: 
If \code{atEventTimesOnly = TRUE} null events are created for events 1, 4 and 6. If \code{atEventTimesOnly = FALSE} null-events are also created for days 2, 3 and 5. }
\item{untilEventOccurrs}{ \code{TRUE/FALSE}. Boolean option to define whether null events should be an option even after an event takes place. If \code{untilEventOccurrs = TRUE} a conditional logisitc logic is applied in that events are only at risk as long as they have not taken place yet. If \code{untilEventOccurrs = FALSE} events continue to be at risk after they have occurred. Note that \code{untilEventOccurrs = TRUE} overwrites the \code{end}-Variable, if specified.}
\item{includeAllPossibleEvents}{ \code{TRUE/FALSE}. Boolean option to allow a more dynamic and specified creation of the risk set. If \code{includeAllPossibleEvents = TRUE}, a data set has to be provided to \code{possibleEvents}.}
\item{possibleEvents}{ An optional data set with the form: column 1 = sender, column 2 = target, 3 = start, 4 = end, 5 = event attribute, 6... . The data set provides all possible events for the entire event sequence and gives each possible event a start and end value to determine when each event could have been possible. This is useful if the risk set follows a complex pattern that cannot be resolved with the above options. E.g., providing a \code{startDate}-variable and setting \code{atEventTimesOnly == FALSE} will result in an error since in a continuous time setting the start variable will be matched to the closest date, rather than to the exact value of said date in the event sequence. Manually coding the possible events is neccessary.}
\item{returnInputData}{ \code{TRUE/FALSE}. Boolean option to check the original data set (handed over in \code{data}) against the created start and stop variables. If \code{returnInputData = TRUE}, a list of two data sets is returned. The first data set is the counting process data set with null-events, the second the modified \code{data}.}
}
\details{
To follow.
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
## Example 1: standard conditional logistic set-up
dt <- data.frame(
  sender = c('a', 'c', 'd', 'a', 'a', 'f', 'c'), 
  target = c('b', 'd', 'd', 'b', 'b', 'a', 'd'), 
  eventSequence = c(1, 2, 2, 3, 3, 4, 6)
)
count.data <- createRemDataset(
  data = dt, sender = dt$sender, 
  target = dt$target, eventSequence = dt$eventSequence, 
  eventAttribute = NULL, time = NULL, 
  start = NULL, startDate = NULL, 
  end = NULL, endDate = NULL, 
  timeformat = NULL,
  atEventTimesOnly = TRUE, untilEventOccurrs = TRUE,
  includeAllPossibleEvents = FALSE, possibleEvents = NULL, 
  returnInputData = FALSE)

## Example 2: add 2 attributes to the event-classification
dt <- data.frame(
  sender = c('a', 'c', 'd', 'a', 'a', 'f', 'c'), 
  target = c('b', 'd', 'd', 'b', 'b', 'a', 'd'), 
  pro.con = c('pro', 'pro', 'con', 'pro', 'con', 'pro', 'pro'),
  attack = c('yes', 'no', 'no', 'yes', 'yes', 'no', 'yes'),
  eventSequence = c(1, 2, 2, 3, 3, 4, 6)
)
count.data <- createRemDataset(
  data = dt, sender = dt$sender, 
  target = dt$target, eventSequence = dt$eventSequence, 
  eventAttribute = paste0(dt$pro.con, dt$attack), time = NULL, 
  start = NULL, startDate = NULL, 
  end = NULL, endDate = NULL, 
  timeformat = NULL,
  atEventTimesOnly = TRUE, untilEventOccurrs = TRUE,
  includeAllPossibleEvents = FALSE, possibleEvents = NULL, 
  returnInputData = FALSE)

## Example 3: adding start and end variables
# Note: the start and end variables will be overwritten 
# if there are duplicate events. If you want to 
# keep the strict start and stop values that you set, use
# includeAllPossibleEvents = TRUE and specify a 
# possibleEvents-data set.
# Note 2: if untilEventOccurrs = TRUE and an end
# variable is provided, this end variable is 
# overwritten. Set untilEventOccurrs 0 FALSE and 
# provide the end variable if you want the events 
# possibilities to stop at these exact event times.
dt <- data.frame(
  sender = c('a', 'c', 'd', 'a', 'a', 'f', 'c'), 
  target = c('b', 'd', 'd', 'b', 'b', 'a', 'd'), 
  eventSequence = c(1, 2, 2, 3, 3, 4, 6),
  start = c(0, 0, 1, 1, 1, 3, 3), 
  end = rep(6, 7)
)
count.data <- createRemDataset(
  data = dt, sender = dt$sender, 
  target = dt$target, eventSequence = dt$eventSequence, 
  eventAttribute = NULL, time = NULL, 
  start = dt$start, startDate = NULL, 
  end = dt$end, endDate = NULL, 
  timeformat = NULL,
  atEventTimesOnly = TRUE, untilEventOccurrs = TRUE,
  includeAllPossibleEvents = FALSE, possibleEvents = NULL, 
  returnInputData = FALSE)

## Example 4: using start (and stop) dates
dt <- data.frame(
  sender = c('a', 'c', 'd', 'a', 'a', 'f', 'c'), 
  target = c('b', 'd', 'd', 'b', 'b', 'a', 'd'), 
  eventSequence = c(1, 2, 2, 3, 3, 4, 6),
  date = c('01.02.1971', rep('02.02.1971', 2), 
rep('03.02.1971', 2), '04.02.1971', '06.02.1971'),
  dateAtRisk = c(rep('21.01.1971', 2), rep('01.02.1971', 5)), 
  dateRiskEnds = rep('01.03.1971', 7)
)
count.data <- createRemDataset(
  data = dt, sender = dt$sender, target = dt$target, 
  eventSequence = dt$eventSequence, 
  eventAttribute = NULL, time = dt$date, 
  start = NULL, startDate = dt$dateAtRisk, 
  end = NULL, endDate = NULL, 
  timeformat = '\%d.\%m.\%Y',
  atEventTimesOnly = TRUE, untilEventOccurrs = TRUE,
  includeAllPossibleEvents = FALSE, possibleEvents = NULL, 
  returnInputData = FALSE)
# if you want to include null-events at times when no event happened, 
# either see Example 5 or create a start-variable by yourself 
# by using the eventSequence()-command with the option 
# 'returnDateSequenceData = TRUE' in this package. With the
# generated sequence, dates from startDate can be matched
# to the event sequence values (using the match()-command).

## Example 5: using start and stop dates and including 
# possible events whenever no event occurred. 
possible.events <- data.frame(
  sender = c('a', 'c', 'd', 'f'), 
  target = c('b', 'd', 'd', 'a'), 
  start = c(0, 0, 1, 1), 
  end = c(rep(8, 4)))
count.data <- createRemDataset(
  data = dt, sender = dt$sender, target = dt$target, 
  eventSequence = dt$eventSequence, 
  eventAttribute = NULL, time = NULL, 
  start = NULL, startDate = NULL, 
  end = NULL, endDate = NULL, 
  timeformat = NULL,
  atEventTimesOnly = TRUE, untilEventOccurrs = TRUE,
  includeAllPossibleEvents = TRUE, possibleEvents = possible.events, 
  returnInputData = FALSE)
# now you can set 'atEventTimesOnly = FALSE' to include 
# null-events where none occurred until the events happened
count.data <- createRemDataset(
  data = dt, sender = dt$sender, target = dt$target, 
  eventSequence = dt$eventSequence, 
  eventAttribute = NULL, time = NULL, 
  start = NULL, startDate = NULL, 
  end = NULL, endDate = NULL, 
  timeformat = NULL,
  atEventTimesOnly = FALSE, untilEventOccurrs = TRUE,
  includeAllPossibleEvents = TRUE, possibleEvents = possible.events, 
  returnInputData = FALSE)
# plus you can set  to get the full range of the events 
# (bounded by max(possible.events$end))
count.data <- createRemDataset(
  data = dt, sender = dt$sender, target = dt$target, 
  eventSequence = dt$eventSequence, 
  eventAttribute = NULL, time = NULL, 
  start = NULL, startDate = NULL, 
  end = NULL, endDate = NULL, 
  timeformat = NULL,
  atEventTimesOnly = FALSE, untilEventOccurrs = FALSE,
  includeAllPossibleEvents = TRUE, possibleEvents = possible.events, 
  returnInputData = FALSE)
}
%\keyword{key}




