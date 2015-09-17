\name{create.event.sequence}
\alias{create.event.sequence}
\alias{eventSequence}
\alias{event.sequence}
\alias{event sequence}
\alias{rem}
\title{Create the event sequence}
%
\description{}
\usage{
create.event.sequence(datevar, dateformat = NULL, data = NULL, 
  type = "continuous", byTime = "1 day",excludeDate = NULL, 
  excludeTypeOfDay = NULL, excludeYear = NULL, excludeFrom = NULL, 
  excludeTo = NULL, returnData = FALSE, sortData = TRUE, ...)
}
\arguments{
\item{datevar}{ The variable containing the information on the date and/or time of the event.}
\item{dateformat}{ A character string indicating the format of the \code{datevar}. see \code{\link{as.Date}}}
\item{data}{ An optional data frame containing all the variables.}
\item{type}{ "`\code{continuous}"' or "`\code{ordinal}"'. Specifies whether the event sequence is to be created as a continuous sequence or an ordinal sequence.}
\item{byTime}{ String value. Specifies at what interval the event sequence is created. Use for instance: "1 day", "2 days", "1 month", "1 year".}
\item{excludeDate}{ An optional string or string vector containing one or more dates that should be excluded from the event.sequence. The dates have to be in the same format as provided in \code{dateformat}. Only valid for continous event sequences.}
\item{excludeTypeOfDay}{ String value or vector naming the day(s) that should be excluded from the event sequence. Depending on the locale the weekdays may be named differently. Use \code{Sys.getlocale("LC_TIME")} to find which locale is installed.}
\item{excludeYear}{ A string value or vector naming the year(s) that should be excluded from the event sequence.}
\item{excludeFrom}{ A string value with the start value of the date from which on the event sequence should not be affected. The value has to be in the same format as specified in the \code{dateformat}.}
\item{excludeTo}{ A string value with the end value of the date to which time the event sequence should not be affected. The value has to be in the same format as specified in the \code{dateformat}.}
\item{returnData}{ \code{TRUE/FALSE}. Default set to \code{TRUE}. The data frame provided is returned in full, together with the new variable for the event sequence. }
\item{sortData}{ \code{TRUE/FALSE}. Default set to \code{TRUE}. Should only be used if \code{returnData = TRUE}. The entire data.frame will be ordered acording to the event sequence.}
}
\details{

In order to estimate relational event models, the events have to be ordered, either according to an ordinal or a continous event sequence. The ordinal event sequence simply orders the events and gives each event a place in the sequence. 
The continuous event sequence creates an artificial sequence ranging from \code{min(datevar)} to \code{max(datevar)} and matches each event with its place in the artificial event sequence. Dates, years or Weekdays can be excluded from the artificial event sequence. This is useful for excluding specific holidays, weekends etc..

Where two or more events occurr at the same time, they are given the same value in the event sequence. 

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




