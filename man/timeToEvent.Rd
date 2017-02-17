\name{timeToEvent}
\alias{timeToEvent}
\title{Calculate the time-to-next-event or the time-since-date for a REM data set.}
%
\description{Calculate time-to-next-event or time-since-date for a REM data set.}
\usage{
timeToEvent(time, type = 'time-to-next-event', timeEventPossible = NULL)
}
\arguments{
\item{time}{ A integer or Date variable reflecting the time of the event. Note: make sure to specify event time not the event sequence in a counting process data set.}
\item{type}{ Either 'time-to-next-event' or 'time-since-date'. 
\code{type = 'time-to-next-event'} calculates the time between the current event and the event closes to the current in the past. 
\code{type = 'time-since-date'} uses the \code{time}-variable as well as the \code{timeEventPossible}-variable to calculate how much time has passed between the two variables, i.e., how long the event took to come true. }
\item{timeEventPossible}{ An optional integer or Date variable to be used if \code{type = 'time-since-date'} is specified.}
}
\details{

To come.

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
## get some random data
dt <- data.frame(
  sender = c('a', 'c', 'd', 'a', 'a', 'f', 'c'), 
  target = c('b', 'd', 'd', 'b', 'b', 'a', 'd'), 
  date = c(rep('10.01.90',2), '11.01.90', '04.01.90',
 '05.01.90', rep('10.01.90',2)),
  start = c(0, 0, 1, 1, 1, 3, 3), 
  end = rep(6, 7), 
  targetAvailableSince = c(rep(-10,6), -2),
  dateTargetAvailable = c(rep('31.12.89',6), '01.01.90')
)

## create event sequence
dt <- eventSequence(dt$date, dateformat = '\%d.\%m.\%y', data = dt,
              type = "continuous", byTime = "daily",
              excludeDate = '07.01.90',
              returnData = TRUE, sortData = TRUE,
              returnDateSequenceData = FALSE)
## also return the sequenceData
dt.seq <- eventSequence(dt$date, dateformat = '\%d.\%m.\%y', data = dt,
                    type = "continuous", byTime = "daily",
                    excludeDate = '07.01.90',
                    returnDateSequenceData = TRUE)

## create counting process data set
dts <- createRemDataset(
	data = dt, sender = dt$sender, target = dt$target, 
	eventSequence = dt$event.seq.cont, 
	eventAttribute = NULL, time = NULL, 
	start = dt$start, startDate = NULL, 
	end = dt$end, endDate = NULL, 
	timeformat = NULL, atEventTimesOnly = TRUE, 
	untilEventOccurrs = TRUE,
	includeAllPossibleEvents = FALSE, 
	possibleEvents = NULL, returnInputData = TRUE)	
## divide up the results: counting process data = 1, original data = 2
dt.rem <- dts[[1]]
dt <- dts[[2]]

## merge all necessary event attribute variables back in
dt.rem$targetAvailableSince <- dt$targetAvailableSince[match(dt.rem$eventID, 
	dt$eventID)]
dt.rem$dateTargetAvailable <- dt$dateTargetAvailable[match(dt.rem$eventID, 
	dt$eventID)]

## add dates to the eventTime
dt.rem$eventDate <- dt.seq$date.sequence[match(dt.rem$eventTime,
	dt.seq$event.sequence)]

## sort the dataframe according to eventTime
dt.rem <- dt.rem[order(dt.rem$eventTime), ]

## 1. numeric, time-to-next-event
dt.rem$timeToNextEvent <- timeToEvent(as.integer(dt.rem$eventTime))

## 2. numeric, time-since
dt.rem$timeSince <- timeToEvent(dt.rem$eventTime, 
                                type = 'time-since-date', 
                                dt.rem$targetAvailableSince)

## 3. Date, time-to-next-event
# since the event sequence excluded 06.01.90 => time to next event differs
# for the two specification with the integr (1) and the Date-variable (2).
# To be consistent, pick the eventTime instead of the Date-variable.
dt.rem$timeToNextEvent2 <- timeToEvent(as.Date(dt.rem$eventDate, '\%d.\%m.\%y'))


## 4. Date, time-since
dt.rem$timeSince2 <- timeToEvent(
	as.Date(dt.rem$eventDate, '\%d.\%m.\%y'),
	type = 'time-since-date', 
	as.Date(dt.rem$dateTargetAvailable, '\%d.\%m.\%y'))

}
%\keyword{key}




