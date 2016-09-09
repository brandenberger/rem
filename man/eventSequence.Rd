\name{eventSequence}
\alias{eventSequence}
\alias{event.sequence}
\alias{event sequence}
\title{Create event sequence}
%
\description{ Create the event sequence for relational event models. Continuous or ordinal sequences can be created. Various dates may be excluded from the sequence (e.g. special holidays, specific weekdays or longer time spans).}
\usage{
eventSequence(datevar, 
    dateformat = NULL, data = NULL, 
    type = "continuous", byTime = "daily",
    excludeDate = NULL, excludeTypeOfDay = NULL, 
    excludeYear = NULL, excludeFrom = NULL, 
    excludeTo = NULL, returnData = FALSE, 
    sortData = FALSE, 
    returnDateSequenceData = FALSE)
}
\arguments{
\item{datevar}{ The variable containing the information on the date and/or time of the event.}
\item{dateformat}{ A character string indicating the format of the \code{datevar}. see \code{\link{as.Date}}}
\item{data}{ An optional data frame containing all the variables.}
\item{type}{ "`\code{continuous}"' or "`\code{ordinal}"'. Specifies whether the event sequence is to be created as a continuous sequence or an ordinal sequence.}
\item{byTime}{ String value. Specifies at what interval the event sequence is created. Use "daily", "monthly" or "yearly".}
\item{excludeDate}{ An optional string or string vector containing one or more dates that should be excluded from the event.sequence. The dates have to be in the same format as provided in \code{dateformat}. Only valid for continuous event sequences.}
\item{excludeTypeOfDay}{ String value or vector naming the day(s) that should be excluded from the event sequence. Depending on the locale the weekdays may be named differently. Use \code{Sys.getlocale("LC_TIME")} to find which locale is installed.}
\item{excludeYear}{ A string value or vector naming the year(s) that should be excluded from the event sequence.}
\item{excludeFrom}{ A string value (or a vector of strings) with the start value of the date from (from-value included) which the event sequence should not be affected. The value has to be in the same format as specified in \code{dateformat}.}
\item{excludeTo}{ A string value (or a vector of strings) with the end value of the date to which time the event sequence should not be affected (to-value included). The value has to be in the same format as specified in \code{dateformat}.}
\item{returnData}{ \code{TRUE/FALSE}. Default set to \code{FALSE}. The data frame provided is returned in full, together with the new variable for the event sequence. }
\item{sortData}{ \code{TRUE/FALSE}. Default set to \code{FALSE}. Should only be used if \code{returnData = TRUE}. The entire data.frame will be ordered according to the event sequence.}
\item{returnDateSequenceData}{ \code{TRUE/FALSE}. Boolean option to return
the full information on which date matches to which sequence number instead
of the event sequence (and corresponding data frame).}
}
\details{

In order to estimate relational event models, the events have to be ordered, either according to an ordinal or a continuous event sequence. The ordinal event sequence simply orders the events and gives each event a place in the sequence. 
The continuous event sequence creates an artificial sequence ranging from \code{min(datevar)} to \code{max(datevar)} and matches each event with its place in the artificial event sequence. Dates, years or Weekdays can be excluded from the artificial event sequence. This is useful for excluding specific holidays, weekends etc..

Where two or more events occur at the same time, they are given the same value in the event sequence. 

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

# combine them into a data.frame
dt <- data.frame(sender, target, time)

# create continuous event sequence: return the data with the 
# event sequence and sort the data according to the event sequence.
dt <- eventSequence(datevar = dt$time, dateformat = '\%y\%m\%d', 
                    data = dt, type = 'continuous', 
                    byTime = 'daily', returnData = TRUE,
                    sortData = TRUE)

# alternative : create variable with the continuous event 
# sequence, unsorted
dt$eventSeq <- eventSequence(datevar = dt$time, 
                             dateformat = '\%y\%m\%d', 
                             data = dt, type = 'continuous',
                             byTime = 'daily', 
                             returnData = FALSE, 
                             sortData = FALSE)
# manually sort the data set
dt <- dt[order(dt$eventSeq), ]

# create the sequence by month
dt$eventSeqMonthly <- eventSequence(datevar = dt$time, 
                                    dateformat = '\%y\%m\%d', 
                                    data = dt, 
                                    type = 'continuous', 
                                    byTime = 'monthly', 
                                    returnData = FALSE, 
                                    sortData = FALSE)

# create the sequence by year
dt$eventSeqYearly <- eventSequence(datevar = dt$time, 
                                   dateformat = '\%y\%m\%d', 
                                   data = dt, 
                                   type = 'continuous', 
                                   byTime = 'yearly', 
                                   returnData = FALSE, 
                                   sortData = FALSE)

# create an ordinal event sequence
dt$eventSeqOrdinal <- eventSequence(datevar = dt$time, 
                                    dateformat = '\%y\%m\%d', 
                                    data = dt, 
                                    type = 'ordinal', 
                                    byTime = 'daily', 
                                    returnData = FALSE, 
                                    sortData = FALSE)

# exclude certain dates
dt$eventSeqEx <- eventSequence(datevar = dt$time,
                               dateformat = '\%y\%m\%d', 
                               data = dt, type = 'continuous',
                               byTime = 'daily', 
                               excludeDate = c('800108', '800112'),
                               returnData = FALSE, 
                               sortData = FALSE)

# return the sequence data set, where all values in the event sequence
# correspond to the date of the events. Useful to calculate
# start-variables for the createRemDataset-command.
seq.data <- eventSequence(datevar = dt$time,
                          dateformat = "\%y\%m\%d", 
                          data = dt, type = "continuous",
                          byTime = "daily", 
                          excludeDate = c("800108", "800112"),
                          returnData = FALSE, 
                          sortData = FALSE, 
                          returnDateSequenceData = TRUE)
}
%\keyword{key}




