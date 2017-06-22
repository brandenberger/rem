## File: Functions to generate variables that capture effects for REMs
## Author: Laurence Brandenberger
## Date: Created: 17. August, last updated: 17.August (00:27)

################################################################################
################################################################################
################################################################################

# TODO 1: add decay-function-option in the values
# TODO 2: add progress bar to functions (and set showprogressbar = TRUE)
# TODO 3: add real function names to all the stop()-outputs - if possible
# TODO 4: OPENMP - implement paralells in cpp-Functions
# TODO 5: tidy up functions - within 80char/line
# TODO 6: missing values - how to handle
# TODO 7: they all have to be of equal length!

################################################################################
##  Inertia
################################################################################

inertiaStat <- function(data, time, sender, target, halflife, 
                        weight = NULL, eventtypevar = NULL, 
                        eventtypevalue = "valuematch",
                        eventfiltervar = NULL, eventfiltervalue = NULL, 
                        eventvar = NULL, variablename = "inertia", 
                        returnData = FALSE, showprogressbar = FALSE, 
                        inParallel = FALSE, cluster = NULL){
  
  ####### check inputs
  ## check if sender and target inputs are available
  if ( is.null(sender) ) {
    stop("No 'sender' argument was provided.")
  }else{
    sender <- as.character(sender)
  }
  
  if ( is.null(target) ) {
    stop("No 'target' argument was provided.")
  }else{
    target <- as.character(target)
  }
  
  ## check if all variables have the same length.
  if (length(sender) != length(target)){
    stop("'sender' and 'target' are not of the same length.")
  }
  
  ## check if event.sequence is well defined (numeric and ever-increasing)
  if ( is.null(time) ) {
    stop("No 'time' argument was provided.")
  }else{
    #test if weight-var is in ascending order
    if ( is.unsorted(time) ) {
      stop("'", time, "' is not sorted. Sort data frame according to the event sequence.")
    }
  }
  
  ## check if time has the requested length
  if (length(sender) != length(time)){
    stop("'sender' and 'time' are not of the same length.")
  }
  
  ## check if weight-var is defined (if not -> create it)
  if ( is.null(weight) ) {
    weight <- rep(1, length(time))
  }
  if ( !is.numeric(weight) ) {
    stop("'", as.name(weight), "' variable is not numeric.") #TODO: deparse(substitute(eventattributevar)) ?
  }
  
  ## check if event-type inputs are available and correctly specified
  if ( !is.null(eventtypevar) ) {
    # check length
    if (length(sender) != length(eventtypevar)){
      stop("'eventtypevar' and 'sender' are not of the same length.")
    }
    # transform
    eventtypevar <- as.character(eventtypevar)
    if ( is.null(eventtypevalue) ){
      stop("No 'eventtypevalue' provided. Use default 'valuematch', or 'valuemix' or string value(s) to determine by which values the events should be filtered.", )
    }
    # check if eventtypevalue is part of the variable
    if ( length(eventtypevalue) > 1  ){
      for ( i in 1:length(eventtypevalue) ){
        if ( length(grep(eventtypevalue[i], eventtypevar)) == 0 ) {
          ##TODO: #deparse(substitute(eventtypevar))
          stop("Value '", eventtypevalue[i], "' is not an element of '", deparse(substitute(eventtypevar)) , "'.") ##deparse(substitute(eventtypevar))
        }
      }#closes i-loop  
      if ( length(unique(duplicated(eventtypevalue))) == 2 ) {
        stop("Duplicate values in 'eventtypevalue'.") 
      }
    }else if ( eventtypevalue != "valuematch" & eventtypevalue != "valuemix" ) {
      if ( length(grep(eventtypevalue, eventtypevar)) == 0 ) {
        ##TODO: #deparse(substitute(eventtypevar))
        stop("Value '", eventtypevalue, "' is not an element of '", deparse(substitute(eventtypevar)) , "'.") ##deparse(substitute(eventtypevar))
      }
    }	
  }
  
  ## check if event-attribute inputs are available and correctly specified
  if ( is.null(eventfiltervar) == FALSE ) {
    # check length
    if (length(sender) != length(eventfiltervar)){
      stop("'eventfiltervar' and 'sender' are not of the same length.")
    }
    # transform
    eventfiltervar <- as.character(eventfiltervar)
    if ( is.null(eventfiltervalue) ){
      stop("No 'eventfiltervalue' provided. Which value should be filtered for?", )
    }
    # check if eventattributevalue is part of the variable
    if ( length(eventfiltervalue) > 0 ){
      for ( i in 1:length(eventfiltervalue) ){
        if ( length(grep(eventfiltervalue[i], eventfiltervar)) == 0 ) {
          stop("Value '", eventfiltervalue[i], "' is not an element of '", as.name(eventfiltervar), "'.")  ##deparse(substitute(eventattributevar))
        }
      }#closes i-loop  
      if ( length(unique(duplicated(eventfiltervalue))) == 2 ) {
        stop("Duplicate values in 'eventfiltervalue'.") 
      }
    }
  }
  
  ## check event-var
  if(is.null(eventvar) == FALSE){
    if(length(unique(eventvar)) == 2){
      if( ( sort(unique(eventvar))[1] == 0 & sort(unique(eventvar))[2] == 1  ) == FALSE){
        stop('eventvar has to be a dummy variable with values 0 for non-events and
             1 for true events.')
      }
    }else{
      stop('eventvar has to be a dummy variable with values 0 for non-events and
             1 for true events.')
    }
  }
  
  ## cannot take parallel and progress bar
  if(isTRUE(inParallel) & isTRUE(showprogressbar)){
    stop('Cannot spit out progress of the function whilst running the 
         loop in parallel. Turn showprogressbar to FALSE.')
  }
  
  ## cannot have parallel without cluster
  if(isTRUE(inParallel) & is.null(cluster)){
    stop('By choosing to run the loop in parallel, you need to define a 
         cluster. For instance: makeCluster(12, type="FORK"). Alternatively, 
         hand over the number of nodes you would like to run the function on.')
  }
  
  ## check if variablename makes sense (no " " etc.)
  variablename <- gsub(" ", "", variablename, fixed = TRUE)
  
  ## create simple data set to be returned for inertia calcuations with more than 1 output-variable
  ##TODO: should there be an event-id-variable?? => that would be useful here
  data.short <- data.frame(time)
  
  ## calculate part of decay function
  xlog <- log(2)/halflife 
  
  ####### calculate stat
  
  ## use event-filter if counting process data is used
  if(is.null(eventvar)){
    countingProcessVar <- rep(1, length(sender))
  }else{
    countingProcessVar <- eventvar
  }
  
  ## 
  result <- rep(NA, length(sender))
  
  ## calculate the inertia effects for each event
  # (1): no type, no filter
  # (2): no type, with filter
  # (3): valuematch, no filter
  # (4): valuematch, with filter
  # (5): valuemix/values provided, no filter
  # (6): valuemix/values provided, with filter
  if ( is.null(eventtypevar) ) {
    if ( is.null(eventfiltervar) ) {
      ################ (1) start off with simple inertia function: no type, no filter
      # run in parallel?
      if(isTRUE(inParallel)){
        
        ##
        doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
        res <- foreach::foreach(i=1:length(sender), .combine=rbind)%dopar%{
          
          ## create vector of times, sender--target tie was made before
          timesIntertiaPast <- time[sender == sender[i] & target == target[i] & 
                                      time < time[i] & countingProcessVar == 1]
          ## get weight
          if(is.null(weight)){
            weightInertiaPast <- rep(1, length(timesIntertiaPast))
          }else{
            weightInertiaPast <- weight[sender == sender[i] & target == target[i] & 
                                          time < time[i] & countingProcessVar == 1]
          }
          
          ## run cpp-loop for all times
          result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                            weightInertiaPast) 
          ## rbind the variable:
          result 
        }# closes foreach-loop
        
        ## transform result variable
        result <- as.numeric(as.character(res))
        
      }else{ # run in standard form
        
        if(isTRUE(showprogressbar)){
          pb <- txtProgressBar(min = 1, max = length(sender), style = 3)
        }
        for(i in 1:length(sender)){
          if(isTRUE(showprogressbar)){
            setTxtProgressBar(pb, i)
          }
          
          
          ## create vector of times, sender--target tie was made before
          timesIntertiaPast <- time[sender == sender[i] & target == target[i] & 
                                      time < time[i] & countingProcessVar == 1]
          ## get weight
          if(is.null(weight)){
            weightInertiaPast <- rep(1, length(timesIntertiaPast))
          }else{
            weightInertiaPast <- weight[sender == sender[i] & target == target[i] & 
                                          time < time[i] & countingProcessVar == 1]
          }
          
          ## run cpp-loop for all times
          result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                               weightInertiaPast) 
          
        } #closes i-loop
      } # closes if no parallel
      
      ## return results
      ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
      if ( returnData == TRUE ) {
        ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 inertia variable that was generated
        return(result)
      }
      ################ done (1)
      
    }else{ # if eventfiltervar = given
      ################ (2) no type,  with filter
      # run in parallel?
      if(isTRUE(inParallel)){
        
        ##
        doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
        res <- foreach::foreach(i=1:length(sender), .combine=rbind)%dopar%{
          
          ## create vector of times, sender--target tie was made before
          timesIntertiaPast <- time[sender == sender[i] & target == target[i] & 
                                      time < time[i] & countingProcessVar == 1 &
                                      eventfiltervar == eventfiltervalue]
          ## get weight
          if(is.null(weight)){
            weightInertiaPast <- rep(1, length(timesIntertiaPast))
          }else{
            weightInertiaPast <- weight[sender == sender[i] & target == target[i] & 
                                          time < time[i] & countingProcessVar == 1 &
                                          eventfiltervar == eventfiltervalue]
          }
          
          ## run cpp-loop for all times
          result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                            weightInertiaPast) 
          ## rbind the variable:
          result 
        }# closes foreach-loop
        
        ## transform result variable
        result <- as.numeric(as.character(res))
        
      }else{ # run in standard form
        
        if(isTRUE(showprogressbar)){
          pb <- txtProgressBar(min = 1, max = length(sender), style = 3)
        }
        for(i in 1:length(sender)){
          if(isTRUE(showprogressbar)){
            setTxtProgressBar(pb, i)
          }
          
          
          ## create vector of times, sender--target tie was made before
          timesIntertiaPast <- time[sender == sender[i] & target == target[i] & 
                                      time < time[i] & countingProcessVar == 1 &
                                      eventfiltervar == eventfiltervalue]
          ## get weight
          if(is.null(weight)){
            weightInertiaPast <- rep(1, length(timesIntertiaPast))
          }else{
            weightInertiaPast <- weight[sender == sender[i] & target == target[i] & 
                                          time < time[i] & countingProcessVar == 1 &
                                          eventfiltervar == eventfiltervalue]
          }
          
          ## run cpp-loop for all times
          result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                               weightInertiaPast) 
          
        } #closes i-loop
      } # closes if no parallel
      
      ## return results
      ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
      if ( returnData == TRUE ) {
        ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
        data <- cbind(data, result)
        names(data)[length(data)] <- paste(variablename, "filtered", sep = ".")
        
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 inertia variable that was generated
        return(result)
      }
      ################ done (2)
    } # closes if-eventfiltervar == given
    
  }else{ # if eventtypevar = given
    if(eventtypevalue == 'valuematch'){
      if ( is.null(eventfiltervar) ) {
        
        ################ ## (3) valuematch, no filter
        # run in parallel?
        if(isTRUE(inParallel)){
          
          ##
          doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
          res <- foreach::foreach(i=1:length(sender), .combine=rbind)%dopar%{
            
            
            ## create vector of times, sender--target tie was made before
            timesIntertiaPast <- time[sender == sender[i] & target == target[i] & 
                                        time < time[i] & countingProcessVar == 1 &
                                        eventtypevar == eventtypevar[i]]
            ## get weight
            if(is.null(weight)){
              weightInertiaPast <- rep(1, length(timesIntertiaPast))
            }else{
              weightInertiaPast <- weight[sender == sender[i] & target == target[i] & 
                                            time < time[i] & countingProcessVar == 1 &
                                            eventtypevar == eventtypevar[i]]
            }
            
            ## run cpp-loop for all times
            result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                              weightInertiaPast) 
            ## rbind the variable:
            result 
          }# closes foreach-loop
          
          ## transform result variable
          result <- as.numeric(as.character(res))
          
        }else{ # run in standard form
          
          if(isTRUE(showprogressbar)){
            pb <- txtProgressBar(min = 1, max = length(sender), style = 3)
          }
          for(i in 1:length(sender)){
            if(isTRUE(showprogressbar)){
              setTxtProgressBar(pb, i)
            }
            
            
            ## create vector of times, sender--target tie was made before
            timesIntertiaPast <- time[sender == sender[i] & target == target[i] & 
                                        time < time[i] & countingProcessVar == 1 &
                                        eventtypevar == eventtypevar[i]]
            ## get weight
            if(is.null(weight)){
              weightInertiaPast <- rep(1, length(timesIntertiaPast))
            }else{
              weightInertiaPast <- weight[sender == sender[i] & target == target[i] & 
                                            time < time[i] & countingProcessVar == 1 &
                                            eventtypevar == eventtypevar[i]]
            }
            
            ## run cpp-loop for all times
            result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                 weightInertiaPast) 
            
          } #closes i-loop
        } # closes if no parallel
        
        ## return results
        ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
        if ( returnData == TRUE ) {
          ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "typematch", sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 inertia variable that was generated
          return(result)
        }
        ################ done (3)
        
      }else{ # if eventfiltervar = given
        
        ################ ## (4) valuematch, with filter
        # run in parallel?
        if(isTRUE(inParallel)){
          
          ##
          doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
          res <- foreach::foreach(i=1:length(sender), .combine=rbind)%dopar%{
            
            ## create vector of times, sender--target tie was made before
            timesIntertiaPast <- time[sender == sender[i] & target == target[i] & 
                                        time < time[i] & countingProcessVar == 1 &
                                        eventfiltervar == eventfiltervalue & 
                                        eventtypevar == eventtypevar[i]]
            ## get weight
            if(is.null(weight)){
              weightInertiaPast <- rep(1, length(timesIntertiaPast))
            }else{
              weightInertiaPast <- weight[sender == sender[i] & target == target[i] & 
                                            time < time[i] & countingProcessVar == 1 &
                                            eventfiltervar == eventfiltervalue & 
                                            eventtypevar == eventtypevar[i]]
            }
            
            ## run cpp-loop for all times
            result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                              weightInertiaPast) 
            ## rbind the variable:
            result 
          }# closes foreach-loop
          
          ## transform result variable
          result <- as.numeric(as.character(res))
          
        }else{ # run in standard form
          
          if(isTRUE(showprogressbar)){
            pb <- txtProgressBar(min = 1, max = length(sender), style = 3)
          }
          for(i in 1:length(sender)){
            if(isTRUE(showprogressbar)){
              setTxtProgressBar(pb, i)
            }
            
            
            ## create vector of times, sender--target tie was made before
            timesIntertiaPast <- time[sender == sender[i] & target == target[i] & 
                                        time < time[i] & countingProcessVar == 1 &
                                        eventfiltervar == eventfiltervalue & 
                                        eventtypevar == eventtypevar[i]]
            ## get weight
            if(is.null(weight)){
              weightInertiaPast <- rep(1, length(timesIntertiaPast))
            }else{
              weightInertiaPast <- weight[sender == sender[i] & target == target[i] & 
                                            time < time[i] & countingProcessVar == 1 &
                                            eventfiltervar == eventfiltervalue &
                                            eventtypevar == eventtypevar[i]]
            }
            
            ## run cpp-loop for all times
            result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                 weightInertiaPast) 
            
          } #closes i-loop
        } # closes if no parallel
        
        ## return results
        ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
        if ( returnData == TRUE ) {
          ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "typematch.filtered", sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 inertia variable that was generated
          return(result)
        }
        ################ done (4)
        
      }
      
    }else{
      ## create unique values - either for a valuemix variable or a eventtypevalues variable
      if(eventtypevalue == 'valuemix'){
        uniqueEventTypeValues <- unique(eventtypevar)
      }else{
        uniqueEventTypeValues <- eventtypevalue
      }
      
      if ( is.null(eventfiltervar) ) {
        ## (5) valuemix/values provided, no filter
        
        for (a in uniqueEventTypeValues){ #current event type
          for (b in uniqueEventTypeValues){ #past event type
            if ( a != b ){
              
              # run in parallel?
              if(isTRUE(inParallel)){
                
                ##
                doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
                res <- foreach::foreach(i=1:length(sender), .combine=rbind)%dopar%{
                  
                  ## 
                  if (eventtypevar[i] == a){
                    
                    ## create vector of times, sender--target tie was made before
                    timesIntertiaPast <- time[sender == sender[i] & target == target[i] & 
                                                time < time[i] & countingProcessVar == 1 & 
                                                eventtypevar == b]
                    ## get weight
                    if(is.null(weight)){
                      weightInertiaPast <- rep(1, length(timesIntertiaPast))
                    }else{
                      weightInertiaPast <- weight[sender == sender[i] & target == target[i] & 
                                                    time < time[i] & countingProcessVar == 1 &
                                                    eventtypevar == b]
                    }
                    
                    ## run cpp-loop for all times
                    result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                      weightInertiaPast) 
                  }else{
                    result <- 0
                  }
                  ## rbind the variable:
                  result 
                }# closes foreach-loop
                
                ## transform result variable
                result <- as.numeric(as.character(res))
                
              }else{ # run in standard form
                
                if(isTRUE(showprogressbar)){
                  pb <- txtProgressBar(min = 1, max = length(sender), style = 3)
                }
                for(i in 1:length(sender)){
                  if(isTRUE(showprogressbar)){
                    setTxtProgressBar(pb, i)
                  }
                  
                  ##
                  if(eventtypevar[i] == a){
                    
                    
                    ## create vector of times, sender--target tie was made before
                    timesIntertiaPast <- time[sender == sender[i] & target == target[i] & 
                                                time < time[i] & countingProcessVar == 1 &
                                                eventtypevar == b]
                    ## get weight
                    if(is.null(weight)){
                      weightInertiaPast <- rep(1, length(timesIntertiaPast))
                    }else{
                      weightInertiaPast <- weight[sender == sender[i] & target == target[i] & 
                                                    time < time[i] & countingProcessVar == 1 &
                                                    eventtypevar == b]
                    }
                    
                    ## run cpp-loop for all times
                    result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                         weightInertiaPast) 
                  }else{ #if eventtypevar != a
                    result[i] <- 0
                  }
                  
                } #closes i-loop
              } # closes if no parallel
              
              ##
              data.short <- cbind(data.short, result)
              
              ## change the name of the data.short's last column
              names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventattributevar))
                                                             gsub(" ", "", a, fixed = TRUE), #represents current event-actor type
                                                             gsub(" ", "", b, fixed = TRUE), #represents actor type of past actions
                                                             sep = ".")									
            }
          } # closes b-loop
          
          ## calculate effect where both events are of the same type (both a)
          # run in parallel?
          if(isTRUE(inParallel)){
            
            ##
            doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
            res <- foreach::foreach(i=1:length(sender), .combine=rbind)%dopar%{
              
              ## 
              if (eventtypevar[i] == a){
                
                ## create vector of times, sender--target tie was made before
                timesIntertiaPast <- time[sender == sender[i] & target == target[i] & 
                                            time < time[i] & countingProcessVar == 1 &
                                            eventtypevar == a]
                ## get weight
                if(is.null(weight)){
                  weightInertiaPast <- rep(1, length(timesIntertiaPast))
                }else{
                  weightInertiaPast <- weight[sender == sender[i] & target == target[i] & 
                                                time < time[i] & countingProcessVar == 1 &
                                                eventtypevar == a]
                }
                
                ## run cpp-loop for all times
                result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                  weightInertiaPast) 
              }else{
                result <- 0
              }
              ## rbind the variable:
              result 
            }# closes foreach-loop
            
            ## transform result variable
            result <- as.numeric(as.character(res))
            
          }else{ # run in standard form
            
            if(isTRUE(showprogressbar)){
              pb <- txtProgressBar(min = 1, max = length(sender), style = 3)
            }
            for(i in 1:length(sender)){
              if(isTRUE(showprogressbar)){
                setTxtProgressBar(pb, i)
              }
              
              ##
              if( eventtypevar[i] == a){
                
                
                ## create vector of times, sender--target tie was made before
                timesIntertiaPast <- time[sender == sender[i] & target == target[i] & 
                                            time < time[i] & countingProcessVar == 1 &
                                            eventtypevar == a]
                ## get weight
                if(is.null(weight)){
                  weightInertiaPast <- rep(1, length(timesIntertiaPast))
                }else{
                  weightInertiaPast <- weight[sender == sender[i] & target == target[i] & 
                                                time < time[i] & countingProcessVar == 1 &
                                                eventtypevar == a]
                }
                
                ## run cpp-loop for all times
                result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                     weightInertiaPast) 
                
              }else{
                result[i] <- 0
              }
              
            } #closes i-loop
          } # closes if no parallel
          
          ##
          data.short <- cbind(data.short, result)
          
          ## change the name of the data.short's last column
          names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventattributevar))
                                                         gsub(" ", "", a, fixed = TRUE), #represents current event-actor type
                                                         gsub(" ", "", a, fixed = TRUE), #represents actor type of past actions
                                                         sep = ".")	
          
        } # closes a-loop
        
        ## return data
        ## return data frames
        if ( returnData == TRUE ) {
          data <- cbind(data, data.short)
          return(data)
        }else{
          return(data.short)
        }
        ################ done (5)
        
      }else{
        ################ (6) valuemix/values provided, with filter 
        
        for (a in uniqueEventTypeValues){ #current event type
          for (b in uniqueEventTypeValues){ #past event type
            if ( a != b ){
              
              # run in parallel?
              if(isTRUE(inParallel)){
                
                ##
                doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
                res <- foreach::foreach(i=1:length(sender), .combine=rbind)%dopar%{
                  
                  ## 
                  if (eventtypevar[i] == a){
                    
                    ## create vector of times, sender--target tie was made before
                    timesIntertiaPast <- time[sender == sender[i] & target == target[i] & 
                                                time < time[i] & countingProcessVar == 1 & 
                                                eventtypevar == b  & 
                                                eventfiltervar == eventfiltervalue]
                    ## get weight
                    if(is.null(weight)){
                      weightInertiaPast <- rep(1, length(timesIntertiaPast))
                    }else{
                      weightInertiaPast <- weight[sender == sender[i] & target == target[i] & 
                                                    time < time[i] & countingProcessVar == 1 &
                                                    eventtypevar == b & 
                                                    eventfiltervar == eventfiltervalue]
                    }
                    
                    ## run cpp-loop for all times
                    result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                      weightInertiaPast) 
                  }else{
                    result <- 0
                  }
                  ## rbind the variable:
                  result 
                }# closes foreach-loop
                
                ## transform result variable
                result <- as.numeric(as.character(res))
                
              }else{ # run in standard form
                
                if(isTRUE(showprogressbar)){
                  pb <- txtProgressBar(min = 1, max = length(sender), style = 3)
                }
                for(i in 1:length(sender)){
                  if(isTRUE(showprogressbar)){
                    setTxtProgressBar(pb, i)
                  }
                  
                  ##
                  if(eventtypevar[i] == a){
                    
                    
                    ## create vector of times, sender--target tie was made before
                    timesIntertiaPast <- time[sender == sender[i] & target == target[i] & 
                                                time < time[i] & countingProcessVar == 1 &
                                                eventtypevar == b & 
                                                eventfiltervar == eventfiltervalue]
                    ## get weight
                    if(is.null(weight)){
                      weightInertiaPast <- rep(1, length(timesIntertiaPast))
                    }else{
                      weightInertiaPast <- weight[sender == sender[i] & target == target[i] & 
                                                    time < time[i] & countingProcessVar == 1 &
                                                    eventtypevar == b & 
                                                    eventfiltervar == eventfiltervalue]
                    }
                    
                    ## run cpp-loop for all times
                    result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                         weightInertiaPast) 
                  }else{ #if eventtypevar != a
                    result[i] <- 0
                  }
                  
                } #closes i-loop
              } # closes if no parallel
              
              ##
              data.short <- cbind(data.short, result)
              
              ## change the name of the data.short's last column
              names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventattributevar))
                                                             gsub(" ", "", a, fixed = TRUE), #represents current event-actor type
                                                             gsub(" ", "", b, fixed = TRUE), #represents actor type of past actions
                                                             sep = ".")									
            }
          } # closes b-loop
          
          ## calculate effect where both events are of the same type (both a)
          # run in parallel?
          if(isTRUE(inParallel)){
            
            ##
            doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
            res <- foreach::foreach(i=1:length(sender), .combine=rbind)%dopar%{
              
              ## 
              if (eventtypevar[i] == a){
                
                ## create vector of times, sender--target tie was made before
                timesIntertiaPast <- time[sender == sender[i] & target == target[i] & 
                                            time < time[i] & countingProcessVar == 1 &
                                            eventtypevar == a & 
                                            eventfiltervar == eventfiltervalue]
                ## get weight
                if(is.null(weight)){
                  weightInertiaPast <- rep(1, length(timesIntertiaPast))
                }else{
                  weightInertiaPast <- weight[sender == sender[i] & target == target[i] & 
                                                time < time[i] & countingProcessVar == 1 &
                                                eventtypevar == a & 
                                                eventfiltervar == eventfiltervalue]
                }
                
                ## run cpp-loop for all times
                result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                  weightInertiaPast) 
              }else{
                result <- 0
              }
              ## rbind the variable:
              result 
            }# closes foreach-loop
            
            ## transform result variable
            result <- as.numeric(as.character(res))
            
          }else{ # run in standard form
            
            if(isTRUE(showprogressbar)){
              pb <- txtProgressBar(min = 1, max = length(sender), style = 3)
            }
            for(i in 1:length(sender)){
              if(isTRUE(showprogressbar)){
                setTxtProgressBar(pb, i)
              }
              
              ##
              if( eventtypevar[i] == a){
                
                
                ## create vector of times, sender--target tie was made before
                timesIntertiaPast <- time[sender == sender[i] & target == target[i] & 
                                            time < time[i] & countingProcessVar == 1 &
                                            eventtypevar == a & 
                                            eventfiltervar == eventfiltervalue]
                ## get weight
                if(is.null(weight)){
                  weightInertiaPast <- rep(1, length(timesIntertiaPast))
                }else{
                  weightInertiaPast <- weight[sender == sender[i] & target == target[i] & 
                                                time < time[i] & countingProcessVar == 1 &
                                                eventtypevar == a & 
                                                eventfiltervar == eventfiltervalue]
                }
                
                ## run cpp-loop for all times
                result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                     weightInertiaPast) 
                
              }else{
                result[i] <- 0
              }
              
            } #closes i-loop
          } # closes if no parallel
          
          ##
          data.short <- cbind(data.short, result)
          
          ## change the name of the data.short's last column
          names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventattributevar))
                                                         gsub(" ", "", a, fixed = TRUE), #represents current event-actor type
                                                         gsub(" ", "", a, fixed = TRUE), #represents actor type of past actions
                                                         sep = ".")	
          
        } # closes a-loop
        
        ## return data
        ## return data frames
        if ( returnData == TRUE ) {
          data <- cbind(data, data.short)
          return(data)
        }else{
          return(data.short)
        }
        
        ################ done (6)
      } 
      
    } # closes eventtypevalues provided/valuemix 
  } #closes if-eventtypevar != NULL
} #closes inertiaStat-function()


################################################################################
##	Degree calculation
################################################################################

degreeStat <- function(data, time, degreevar, halflife, weight = NULL,
                       eventtypevar = NULL, eventtypevalue = "valuematch", 
                       eventfiltervar = NULL, 
                       eventfiltervalue = NULL, 
                       eventvar = NULL,
                       degreeOnOtherVar = NULL,
                       variablename = "degree", returnData = FALSE, 
                       dataPastEvents = NULL,
                       showprogressbar = FALSE, 
                       inParallel = FALSE, cluster = NULL){
  
  ####### check inputs
  ## check if degreevar input is available
  if ( is.null(degreevar) ) {
    stop("No 'degreevar' argument was provided.")
  }else{
    degreevar <- as.character(degreevar)
  }
  
  ## check if event.sequence is well defined (numeric and ever-increasing)
  if ( is.null(time) ) {
    stop("No 'time' argument was provided.")
  }else{
    #test if weight-var is in ascending order
    if ( is.unsorted(time) ) {
      stop("'", time, "' is not sorted. Sort data frame according to the event sequence.")
    }
  }
  
  ## check if degree and time are of same length
  if (length(degreevar) != length(time)){
    stop("'degreevar' and 'time' are not of the same length.")
  }
  
  ## check if weight-var is defined (if not -> create it)
  if ( is.null(weight) ) {
    weight <- rep(1, length(time))
  }
  if ( !is.numeric(weight) ) {
    stop("'", as.name(weight), "' variable is not numeric.") #TODO: deparse(substitute(eventattributevar)) ?
  }
  
  ## check if degree.on.other.var and degreevar are of same length
  if ( !is.null(degreeOnOtherVar)){
    if ( length(degreevar) != length(degreeOnOtherVar) ){
      stop("'degree.on.other.var' and 'degreevar' are not of same length.")
    }
    degreeOnOtherVar <- as.character(degreeOnOtherVar)
  }
  
  ## check if event-type inputs are available and correctly specified
  if ( !is.null(eventtypevar) ) {
    # check if degreevar and eventtypevar are of same length
    if (length(degreevar) != length(eventtypevar)){
      stop("'eventtypevar' and 'degreevar' are not of the same length.")
    }
    
    # transform eventtypevar
    eventtypevar <- as.character(eventtypevar)
    if ( is.null(eventtypevalue) ){
      stop("No 'eventtypevalue' provided. Use default 'valuematch', or 'valuemix' or string value(s) to determine by which values the events should be filtered.", )
    }
    # check if eventtypevalue is part of the variable
    if ( length(eventtypevalue) > 1  ){
      for ( i in 1:length(eventtypevalue) ){
        if ( length(grep(eventtypevalue[i], eventtypevar)) == 0 ) {
          ##TODO: #deparse(substitute(eventtypevar))
          stop("Value '", eventtypevalue[i], "' is not an element of '", deparse(substitute(eventtypevar)) , "'.") ##deparse(substitute(eventtypevar))
        }
      }#closes i-loop  
      if ( length(unique(duplicated(eventtypevalue))) == 2 ) {
        stop("Duplicate values in 'eventtypevalue'.") 
      }
    }else if ( eventtypevalue != "valuematch" & eventtypevalue != "valuemix" ) {
      if ( length(grep(eventtypevalue, eventtypevar)) == 0 ) {
        ##TODO: #deparse(substitute(eventtypevar))
        stop("Value '", eventtypevalue, "' is not an element of '", deparse(substitute(eventtypevar)) , "'.") ##deparse(substitute(eventtypevar))
      }
    }	
  }
  
  ## check if event-attribute inputs are available and correctly specified
  if ( is.null(eventfiltervar) == FALSE ) {
    # check length
    if (length(degreevar) != length(eventfiltervar)){
      stop("'eventfiltervar' and 'degreevar' are not of the same length.")
    }
    # transform
    eventfiltervar <- as.character(eventfiltervar)
    if ( is.null(eventfiltervalue) ){
      stop("No 'eventfiltervalue' provided. Which value should be filtered for?", )
    }
    # check if eventattributevalue is part of the variable
    if ( length(eventfiltervalue) > 0 ){
      for ( i in 1:length(eventfiltervalue) ){
        if ( length(grep(eventfiltervalue[i], eventfiltervar)) == 0 ) {
          stop("Value '", eventfiltervalue[i], "' is not an element of '", as.name(eventfiltervar), "'.")  ##deparse(substitute(eventattributevar))
        }
      }#closes i-loop  
      if ( length(unique(duplicated(eventfiltervalue))) == 2 ) {
        stop("Duplicate values in 'eventfiltervalue'.") 
      }
    }
  }
  
  ## check event-var
  if(is.null(eventvar) == FALSE){
    if(length(unique(eventvar)) == 2){
      if( ( sort(unique(eventvar))[1] == 0 & sort(unique(eventvar))[2] == 1  ) == FALSE){
        stop('eventvar has to be a dummy variable with values 0 for non-events and
             1 for true events.')
      }
    }else{
      stop('eventvar has to be a dummy variable with values 0 for non-events and
             1 for true events.')
    }
  }
  
  ## cannot take parallel and progress bar
  if(isTRUE(inParallel) & isTRUE(showprogressbar)){
    stop('Cannot spit out progress of the function whilst running the 
         loop in parallel. Turn showprogressbar to FALSE.')
  }
  
  ## cannot have parallel without cluster
  if(isTRUE(inParallel) & is.null(cluster)){
    stop('By choosing to run the loop in parallel, you need to define a 
         cluster. For instance: makeCluster(12, type="FORK"). Alternatively, 
         hand over the number of nodes you would like to run the function on.')
  }
  
  ## check if variablename makes sense (no " " etc.)
  variablename <- gsub(" ", "", variablename, fixed = TRUE)
  
  ## create simple data set to be returned for degree calcuations with more than 1 output-variable
  ##TODO: should there be an event-id-variable?? => that would be useful here
  data.short <- data.frame(time)
  
  ## calculate part of decay function
  xlog <- log(2)/halflife 
  
  if(is.null(dataPastEvents)){
    ## use event-filter if counting process data is used
    if(is.null(eventvar)){
      countingProcessVar <- rep(1, length(degreevar))
    }else{
      countingProcessVar <- eventvar
    }
    ## time
    pastTimeVar <- time
    
    ## weight
    pastWeightVar <- weight
    
    ## eventfiltervar
    pastEventFilterVar <- eventfiltervar
    
    ## eventtypevar
    pastEventTypeVar <- eventtypevar
    
    ## specify degreevariable
    if(is.null(degreeOnOtherVar)){
      currentEventDegreeVar <- degreevar
      pastEventDegreeVar <- degreevar
    }else{
      currentEventDegreeVar <- degreevar
      pastEventDegreeVar <- degreeOnOtherVar
    }
  }else{ # dataPastEvents is specified:
    ##
    if(is.null(eventvar)){
      countingProcessVar <- rep(1, nrow(dataPastEvents))
    }else{
      countingProcessVar <- dataPastEvents[,4]
    }
    ## degree
    if(is.null(degreeOnOtherVar)){
      currentEventDegreeVar <- degreevar
      pastEventDegreeVar <- dataPastEvents[,2]
      # nope doesn't work: dataPastEvents[,which(names(dataPastEvents) == as.name(degreevar) )]
    }else{
      currentEventDegreeVar <- degreevar
      pastEventDegreeVar <- dataPastEvents[,3]
    }
    ## time
    pastTimeVar <- dataPastEvents[,1]
    
    ## weight
    pastWeightVar <- dataPastEvents[,5]
    
    ## eventfiltervar
    pastEventFilterVar <- dataPastEvents[,7]
    
    ## eventtypevar
    pastEventTypeVar <- dataPastEvents[,6]
    
  }
  ## 
  result <- rep(NA, length(degreevar))
  
  ####### calculate stat
  ## calculate the degree effects for each event
  ## calculate the degree effects for each event
  # (1): no type, no filter
  # (2): no type, with filter
  # (3): valuematch, no filter
  # (4): valuematch, with filter
  # (5): valuemix/values provided, no filter
  # (6): valuemix/values provided, with filter
  if ( is.null(eventtypevar) ) {
    if ( is.null(eventfiltervar) ) {
      ################ (1) start off with simple inertia function: no type, no filter
      # run in parallel?
      if(isTRUE(inParallel)){
        
        ##
        doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
        res <- foreach::foreach(i=1:length(degreevar), .combine=rbind)%dopar%{
          
          ## create vector of times, degreevar--target tie was made before
          timesIntertiaPast <- pastTimeVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                             pastTimeVar < time[i] & countingProcessVar == 1]
          ## get weight
          if(is.null(weight)){
            weightInertiaPast <- rep(1, length(timesIntertiaPast))
          }else{
            weightInertiaPast <- pastWeightVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                 pastTimeVar < time[i] & countingProcessVar == 1]
          }
          
          ## run cpp-loop for all times
          result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                            weightInertiaPast) 
          ## rbind the variable:
          result 
        }# closes foreach-loop
        
        ## transform result variable
        result <- as.numeric(as.character(res))
        
      }else{ # run in standard form
        
        if(isTRUE(showprogressbar)){
          pb <- txtProgressBar(min = 1, max = length(degreevar), style = 3)
        }
        for(i in 1:length(degreevar)){
          if(isTRUE(showprogressbar)){
            setTxtProgressBar(pb, i)
          }
          
          
          ## create vector of times, degreevar--target tie was made before
          timesIntertiaPast <- pastTimeVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                             pastTimeVar < time[i] & countingProcessVar == 1]
          ## get weight
          if(is.null(weight)){
            weightInertiaPast <- rep(1, length(timesIntertiaPast))
          }else{
            weightInertiaPast <- pastWeightVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                 pastTimeVar < time[i] & countingProcessVar == 1]
          }
          
          ## run cpp-loop for all times
          result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                               weightInertiaPast) 
          
        } #closes i-loop
      } # closes if no parallel
      
      ## return results
      ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
      if ( returnData == TRUE ) {
        ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 inertia variable that was generated
        return(result)
      }
      ################ done (1)
      
    }else{ # if eventfiltervar = given
      ################ (2) no type,  with filter
      # run in parallel?
      if(isTRUE(inParallel)){
        
        ##
        doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
        res <- foreach::foreach(i=1:length(degreevar), .combine=rbind)%dopar%{
          
          ## create vector of times, sender--target tie was made before
          timesIntertiaPast <- pastTimeVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                             pastTimeVar < time[i] & countingProcessVar == 1 &
                                             pastEventFilterVar == eventfiltervalue]
          ## get weight
          if(is.null(weight)){
            weightInertiaPast <- rep(1, length(timesIntertiaPast))
          }else{
            weightInertiaPast <- pastWeightVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                 pastTimeVar < time[i] & countingProcessVar == 1 &
                                                 pastEventFilterVar == eventfiltervalue]
          }
          
          ## run cpp-loop for all times
          result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                            weightInertiaPast) 
          ## rbind the variable:
          result 
        }# closes foreach-loop
        
        ## transform result variable
        result <- as.numeric(as.character(res))
        
      }else{ # run in standard form
        
        if(isTRUE(showprogressbar)){
          pb <- txtProgressBar(min = 1, max = length(degreevar), style = 3)
        }
        for(i in 1:length(degreevar)){
          if(isTRUE(showprogressbar)){
            setTxtProgressBar(pb, i)
          }
          
          
          ## create vector of times, sender--target tie was made before
          timesIntertiaPast <- pastTimeVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                             pastTimeVar < time[i] & countingProcessVar == 1 &
                                             pastEventFilterVar == eventfiltervalue]
          ## get weight
          if(is.null(weight)){
            weightInertiaPast <- rep(1, length(timesIntertiaPast))
          }else{
            weightInertiaPast <- pastWeightVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                 pastTimeVar < time[i] & countingProcessVar == 1 &
                                                 pastEventFilterVar == eventfiltervalue]
          }
          
          ## run cpp-loop for all times
          result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                               weightInertiaPast) 
          
        } #closes i-loop
      } # closes if no parallel
      
      ## return results
      ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
      if ( returnData == TRUE ) {
        ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
        data <- cbind(data, result)
        names(data)[length(data)] <- paste(variablename, "filtered", sep = ".")
        
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 inertia variable that was generated
        return(result)
      }
      ################ done (2)
    } # closes if-eventfiltervar == given
    
  }else{ # if eventtypevar = given
    if(eventtypevalue == 'valuematch'){
      if ( is.null(eventfiltervar) ) {
        
        ################ ## (3) valuematch, no filter
        # run in parallel?
        if(isTRUE(inParallel)){
          
          ##
          doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
          res <- foreach::foreach(i=1:length(degreevar), .combine=rbind)%dopar%{
            
            
            ## create vector of times, sender--target tie was made before
            timesIntertiaPast <- pastTimeVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                               pastTimeVar < time[i] & countingProcessVar == 1 &
                                               pastEventTypeVar == eventtypevar[i]]
            ## get weight
            if(is.null(weight)){
              weightInertiaPast <- rep(1, length(timesIntertiaPast))
            }else{
              weightInertiaPast <- pastWeightVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                   pastTimeVar < time[i] & countingProcessVar == 1 &
                                                   pastEventTypeVar == eventtypevar[i]]
            }
            
            ## run cpp-loop for all times
            result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                              weightInertiaPast) 
            ## rbind the variable:
            result 
          }# closes foreach-loop
          
          ## transform result variable
          result <- as.numeric(as.character(res))
          
        }else{ # run in standard form
          
          if(isTRUE(showprogressbar)){
            pb <- txtProgressBar(min = 1, max = length(degreevar), style = 3)
          }
          for(i in 1:length(degreevar)){
            if(isTRUE(showprogressbar)){
              setTxtProgressBar(pb, i)
            }
            
            
            ## create vector of times, degreevar--target tie was made before
            timesIntertiaPast <- pastTimeVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                               pastTimeVar < time[i] & countingProcessVar == 1 &
                                               pastEventTypeVar == eventtypevar[i]]
            ## get weight
            if(is.null(weight)){
              weightInertiaPast <- rep(1, length(timesIntertiaPast))
            }else{
              weightInertiaPast <- pastWeightVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                   pastTimeVar < time[i] & countingProcessVar == 1 &
                                                   pastEventTypeVar == eventtypevar[i]]
            }
            
            ## run cpp-loop for all times
            result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                 weightInertiaPast) 
            
          } #closes i-loop
        } # closes if no parallel
        
        ## return results
        ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
        if ( returnData == TRUE ) {
          ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "typematch", sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 inertia variable that was generated
          return(result)
        }
        ################ done (3)
        
      }else{ # if eventfiltervar = given
        
        ################ ## (4) valuematch, with filter
        # run in parallel?
        if(isTRUE(inParallel)){
          
          ##
          doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
          res <- foreach::foreach(i=1:length(degreevar), .combine=rbind)%dopar%{
            
            ## create vector of times, sender--target tie was made before
            timesIntertiaPast <- pastTimeVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                               pastTimeVar < time[i] & countingProcessVar == 1 &
                                               pastEventFilterVar == eventfiltervalue & 
                                               pastEventTypeVar == eventtypevar[i]]
            ## get weight
            if(is.null(weight)){
              weightInertiaPast <- rep(1, length(timesIntertiaPast))
            }else{
              weightInertiaPast <- pastWeightVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                   pastTimeVar < time[i] & countingProcessVar == 1 &
                                                   pastEventFilterVar == eventfiltervalue & 
                                                   pastEventTypeVar == eventtypevar[i]]
            }
            
            ## run cpp-loop for all times
            result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                              weightInertiaPast) 
            ## rbind the variable:
            result 
          }# closes foreach-loop
          
          ## transform result variable
          result <- as.numeric(as.character(res))
          
        }else{ # run in standard form
          
          if(isTRUE(showprogressbar)){
            pb <- txtProgressBar(min = 1, max = length(degreevar), style = 3)
          }
          for(i in 1:length(degreevar)){
            if(isTRUE(showprogressbar)){
              setTxtProgressBar(pb, i)
            }
            
            
            ## create vector of times, sender--target tie was made before
            timesIntertiaPast <- pastTimeVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                               pastTimeVar < time[i] & countingProcessVar == 1 &
                                               pastEventFilterVar == eventfiltervalue & 
                                               pastEventTypeVar == eventtypevar[i]]
            ## get weight
            if(is.null(weight)){
              weightInertiaPast <- rep(1, length(timesIntertiaPast))
            }else{
              weightInertiaPast <- pastWeightVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                   pastTimeVar < time[i] & countingProcessVar == 1 &
                                                   pastEventFilterVar == eventfiltervalue &
                                                   pastEventTypeVar == eventtypevar[i]]
            }
            
            ## run cpp-loop for all times
            result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                 weightInertiaPast) 
            
          } #closes i-loop
        } # closes if no parallel
        
        ## return results
        ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
        if ( returnData == TRUE ) {
          ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "typematch.filtered", sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 inertia variable that was generated
          return(result)
        }
        ################ done (4)
        
      }
      
    }else{
      ## create unique values - either for a valuemix variable or a eventtypevalues variable
      if(eventtypevalue == 'valuemix'){
        uniqueEventTypeValues <- unique(eventtypevar)
      }else{
        uniqueEventTypeValues <- eventtypevalue
      }
      
      if ( is.null(eventfiltervar) ) {
        ## (5) valuemix/values provided, no filter
        
        for (a in uniqueEventTypeValues){ #current event type
          for (b in uniqueEventTypeValues){ #past event type
            if ( a != b ){
              
              # run in parallel?
              if(isTRUE(inParallel)){
                
                ##
                doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
                res <- foreach::foreach(i=1:length(degreevar), .combine=rbind)%dopar%{
                  
                  ## 
                  if (eventtypevar[i] == a){
                    
                    ## create vector of times, sender--target tie was made before
                    timesIntertiaPast <- pastTimeVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                       pastTimeVar < time[i] & countingProcessVar == 1 & 
                                                       pastEventTypeVar == b]
                    ## get weight
                    if(is.null(weight)){
                      weightInertiaPast <- rep(1, length(timesIntertiaPast))
                    }else{
                      weightInertiaPast <- pastWeightVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                           pastTimeVar < time[i] & countingProcessVar == 1 &
                                                           pastEventTypeVar == b]
                    }
                    
                    ## run cpp-loop for all times
                    result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                      weightInertiaPast) 
                  }else{
                    result <- 0
                  }
                  ## rbind the variable:
                  result 
                }# closes foreach-loop
                
                ## transform result variable
                result <- as.numeric(as.character(res))
                
              }else{ # run in standard form
                
                if(isTRUE(showprogressbar)){
                  pb <- txtProgressBar(min = 1, max = length(degreevar), style = 3)
                }
                for(i in 1:length(degreevar)){
                  if(isTRUE(showprogressbar)){
                    setTxtProgressBar(pb, i)
                  }
                  
                  ##
                  if(eventtypevar[i] == a){
                    
                    
                    ## create vector of times, sender--target tie was made before
                    timesIntertiaPast <- pastTimeVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                       pastTimeVar < time[i] & countingProcessVar == 1 &
                                                       pastEventTypeVar == b]
                    ## get weight
                    if(is.null(weight)){
                      weightInertiaPast <- rep(1, length(timesIntertiaPast))
                    }else{
                      weightInertiaPast <- pastWeightVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                           pastTimeVar < time[i] & countingProcessVar == 1 &
                                                           pastEventTypeVar == b]
                    }
                    
                    ## run cpp-loop for all times
                    result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                         weightInertiaPast) 
                  }else{ #if eventtypevar != a
                    result[i] <- 0
                  }
                  
                } #closes i-loop
              } # closes if no parallel
              
              ##
              data.short <- cbind(data.short, result)
              
              ## change the name of the data.short's last column
              names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventattributevar))
                                                             gsub(" ", "", a, fixed = TRUE), #represents current event-actor type
                                                             gsub(" ", "", b, fixed = TRUE), #represents actor type of past actions
                                                             sep = ".")									
            }
          } # closes b-loop
          
          ## calculate effect where both events are of the same type (both a)
          # run in parallel?
          if(isTRUE(inParallel)){
            
            ##
            doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
            res <- foreach::foreach(i=1:length(degreevar), .combine=rbind)%dopar%{
              
              ## 
              if (eventtypevar[i] == a){
                
                ## create vector of times, sender--target tie was made before
                timesIntertiaPast <- pastTimeVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                   pastTimeVar < time[i] & countingProcessVar == 1 &
                                                   pastEventTypeVar == a]
                ## get weight
                if(is.null(weight)){
                  weightInertiaPast <- rep(1, length(timesIntertiaPast))
                }else{
                  weightInertiaPast <- pastWeightVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                       pastTimeVar < time[i] & countingProcessVar == 1 &
                                                       pastEventTypeVar == a]
                }
                
                ## run cpp-loop for all times
                result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                  weightInertiaPast) 
              }else{
                result <- 0
              }
              ## rbind the variable:
              result 
            }# closes foreach-loop
            
            ## transform result variable
            result <- as.numeric(as.character(res))
            
          }else{ # run in standard form
            
            if(isTRUE(showprogressbar)){
              pb <- txtProgressBar(min = 1, max = length(degreevar), style = 3)
            }
            for(i in 1:length(degreevar)){
              if(isTRUE(showprogressbar)){
                setTxtProgressBar(pb, i)
              }
              
              ##
              if( eventtypevar[i] == a){
                
                
                ## create vector of times, sender--target tie was made before
                timesIntertiaPast <- pastTimeVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                   pastTimeVar < time[i] & countingProcessVar == 1 &
                                                   pastEventTypeVar == a]
                ## get weight
                if(is.null(weight)){
                  weightInertiaPast <- rep(1, length(timesIntertiaPast))
                }else{
                  weightInertiaPast <- pastWeightVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                       pastTimeVar < time[i] & countingProcessVar == 1 &
                                                       pastEventTypeVar == a]
                }
                
                ## run cpp-loop for all times
                result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                     weightInertiaPast) 
                
              }else{
                result[i] <- 0
              }
              
            } #closes i-loop
          } # closes if no parallel
          
          ##
          data.short <- cbind(data.short, result)
          
          ## change the name of the data.short's last column
          names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventattributevar))
                                                         gsub(" ", "", a, fixed = TRUE), #represents current event-actor type
                                                         gsub(" ", "", a, fixed = TRUE), #represents actor type of past actions
                                                         sep = ".")	
          
        } # closes a-loop
        
        ## return data
        ## return data frames
        if ( returnData == TRUE ) {
          data <- cbind(data, data.short)
          return(data)
        }else{
          return(data.short)
        }
        ################ done (5)
        
      }else{
        ################ (6) valuemix/values provided, with filter 
        
        for (a in uniqueEventTypeValues){ #current event type
          for (b in uniqueEventTypeValues){ #past event type
            if ( a != b ){
              
              # run in parallel?
              if(isTRUE(inParallel)){
                
                ##
                doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
                res <- foreach::foreach(i=1:length(degreevar), .combine=rbind)%dopar%{
                  
                  ## 
                  if (eventtypevar[i] == a){
                    
                    ## create vector of times, degreevar--target tie was made before
                    timesIntertiaPast <- pastTimeVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                       pastTimeVar < time[i] & countingProcessVar == 1 & 
                                                       pastEventTypeVar == b  & 
                                                       pastEventFilterVar == eventfiltervalue]
                    ## get weight
                    if(is.null(weight)){
                      weightInertiaPast <- rep(1, length(timesIntertiaPast))
                    }else{
                      weightInertiaPast <- pastWeightVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                           pastTimeVar < time[i] & countingProcessVar == 1 &
                                                           pastEventTypeVar == b & 
                                                           pastEventFilterVar == eventfiltervalue]
                    }
                    
                    ## run cpp-loop for all times
                    result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                      weightInertiaPast) 
                  }else{
                    result <- 0
                  }
                  ## rbind the variable:
                  result 
                }# closes foreach-loop
                
                ## transform result variable
                result <- as.numeric(as.character(res))
                
              }else{ # run in standard form
                
                if(isTRUE(showprogressbar)){
                  pb <- txtProgressBar(min = 1, max = length(degreevar), style = 3)
                }
                for(i in 1:length(degreevar)){
                  if(isTRUE(showprogressbar)){
                    setTxtProgressBar(pb, i)
                  }
                  
                  ##
                  if(eventtypevar[i] == a){
                    
                    
                    ## create vector of times, sender--target tie was made before
                    timesIntertiaPast <- pastTimeVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                       pastTimeVar < time[i] & countingProcessVar == 1 &
                                                       pastEventTypeVar == b & 
                                                       pastEventFilterVar == eventfiltervalue]
                    ## get weight
                    if(is.null(weight)){
                      weightInertiaPast <- rep(1, length(timesIntertiaPast))
                    }else{
                      weightInertiaPast <- pastWeightVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                           pastTimeVar < time[i] & countingProcessVar == 1 &
                                                           pastEventTypeVar == b & 
                                                           pastEventFilterVar == eventfiltervalue]
                    }
                    
                    ## run cpp-loop for all times
                    result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                         weightInertiaPast) 
                  }else{ #if eventtypevar != a
                    result[i] <- 0
                  }
                  
                } #closes i-loop
              } # closes if no parallel
              
              ##
              data.short <- cbind(data.short, result)
              
              ## change the name of the data.short's last column
              names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventattributevar))
                                                             gsub(" ", "", a, fixed = TRUE), #represents current event-actor type
                                                             gsub(" ", "", b, fixed = TRUE), #represents actor type of past actions
                                                             sep = ".")									
            }
          } # closes b-loop
          
          ## calculate effect where both events are of the same type (both a)
          # run in parallel?
          if(isTRUE(inParallel)){
            
            ##
            doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
            res <- foreach::foreach(i=1:length(degreevar), .combine=rbind)%dopar%{
              
              ## 
              if (eventtypevar[i] == a){
                
                ## create vector of times, degreevar--target tie was made before
                timesIntertiaPast <- pastTimeVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                   pastTimeVar < time[i] & countingProcessVar == 1 &
                                                   pastEventTypeVar == a & 
                                                   pastEventFilterVar == eventfiltervalue]
                ## get weight
                if(is.null(weight)){
                  weightInertiaPast <- rep(1, length(timesIntertiaPast))
                }else{
                  weightInertiaPast <- pastWeightVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                       pastTimeVar < time[i] & countingProcessVar == 1 &
                                                       pastEventTypeVar == a & 
                                                       pastEventFilterVar == eventfiltervalue]
                }
                
                ## run cpp-loop for all times
                result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                  weightInertiaPast) 
              }else{
                result <- 0
              }
              ## rbind the variable:
              result 
            }# closes foreach-loop
            
            ## transform result variable
            result <- as.numeric(as.character(res))
            
          }else{ # run in standard form
            
            if(isTRUE(showprogressbar)){
              pb <- txtProgressBar(min = 1, max = length(degreevar), style = 3)
            }
            for(i in 1:length(degreevar)){
              if(isTRUE(showprogressbar)){
                setTxtProgressBar(pb, i)
              }
              
              ##
              if( eventtypevar[i] == a){
                
                
                ## create vector of times, degreevar--target tie was made before
                timesIntertiaPast <- pastTimeVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                   pastTimeVar < time[i] & countingProcessVar == 1 &
                                                   pastEventTypeVar == a & 
                                                   pastEventFilterVar == eventfiltervalue]
                ## get weight
                if(is.null(weight)){
                  weightInertiaPast <- rep(1, length(timesIntertiaPast))
                }else{
                  weightInertiaPast <- pastWeightVar[pastEventDegreeVar == currentEventDegreeVar[i] &
                                                       pastTimeVar < time[i] & countingProcessVar == 1 &
                                                       pastEventTypeVar == a & 
                                                       pastEventFilterVar == eventfiltervalue]
                }
                
                ## run cpp-loop for all times
                result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                     weightInertiaPast) 
                
              }else{
                result[i] <- 0
              }
              
            } #closes i-loop
          } # closes if no parallel
          
          ##
          data.short <- cbind(data.short, result)
          
          ## change the name of the data.short's last column
          names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventattributevar))
                                                         gsub(" ", "", a, fixed = TRUE), #represents current event-actor type
                                                         gsub(" ", "", a, fixed = TRUE), #represents actor type of past actions
                                                         sep = ".")	
          
        } # closes a-loop
        
        ## return data
        ## return data frames
        if ( returnData == TRUE ) {
          data <- cbind(data, data.short)
          return(data)
        }else{
          return(data.short)
        }
        
        ################ done (6)
      } 
      
    } # closes eventtypevalues provided/valuemix 
  } #closes if-eventtypevar != NULL
}


################################################################################
##	FourCycle calculation
################################################################################

fourCycleStat <- function(data, time, sender, target, halflife, weight = NULL,
                          eventtypevar = NULL, eventtypevalue = 'standard', 
                          eventfiltervar = NULL, 
                          eventfilterAB = NULL, eventfilterAJ = NULL, 
                          eventfilterIB = NULL, eventfilterIJ = NULL, 
                          eventvar = NULL,
                          variablename = 'fourCycle', returnData = FALSE, 
                          dataPastEvents = NULL,
                          showprogressbar = FALSE, 
                          inParallel = FALSE, cluster = NULL){
  
  ####### check inputs
  ## check if sender input is available
  if ( is.null(sender) ) {
    stop("No 'sender' argument was provided.")
  }else{
    sender <- as.character(sender)
  }
  
  ## check if target input is available
  if ( is.null(target) ) {
    stop("No 'target' argument was provided.")
  }else{
    target <- as.character(target)
  }	
  
  ## check if event.sequence is well defined (numeric and ever-increasing)
  if ( is.null(time) ) {
    stop("No 'time' argument was provided.")
  }else{
    #test if weight-var is in ascending order
    if ( is.unsorted(time) ) {
      stop("'", time, "' is not sorted. Sort data frame according to the event sequence.")
    }
  }
  
  ## check if all variables are of same length
  if( length(sender) != length(target) ){
    stop("'sender' and 'target' are not of same length.")
  }
  if ( length(sender) != length(time) ){
    stop("'sender' and 'time' are not of same length.")
  }
  
  ## check if weight-var is defined (if not -> create it)
  if ( is.null(weight) ) {
    weight <- rep(1, length(time))
  }
  if ( !is.numeric(weight) ) {
    stop("'", as.name(weight), "' variable is not numeric.") #TODO: deparse(substitute(eventattributevar)) ?
  }
  
  ## check if event-type inputs are available and correctly specified
  if ( !is.null(eventtypevar) ) {
    # length test
    if ( length(sender) != length(eventtypevar) ){
      stop("'sender' and 'eventtypevar' are not of same length.")
    }
    # transform variable
    eventtypevar <- as.character(eventtypevar)
    if ( length(unique(eventtypevar)) != 2 ){ 
      stop("'eventtypevar' is not a dummy variable.")
    }
    if ( is.null(eventtypevalue) ){
      stop("No 'eventtypevalue' provided. Use default 'standard', or 'positive' or 'negative' to determine the overall type of the four-cycle effect.", )
    }
    if ( eventtypevalue == "standard" | eventtypevalue == "positive" | eventtypevalue == "negative"){}else{
      stop("'eventtypevalue' not specified correctly. Use default 'standard', or 'positive' or 'negative' to determine the overall type of the four-cycle effect.", )
    }	
  }else{
    eventtypevalue <- "standard"
  }
  
  ## check if event-attribute inputs are available and correctly specified
  if ( is.null(eventfiltervar) == FALSE ) {
    # length test
    if ( length(sender) != length(eventfiltervar) ){
      stop("'sender' and 'eventfiltervar' are not of same length.")
    }
    # transform variable
    eventfiltervar <- as.character(eventfiltervar)
    if ( is.null(eventfilterAB) & is.null(eventfilterAJ) & is.null(eventfilterIB) & is.null(eventfilterIJ) ){
      stop("No 'eventattribute__' provided. Provide a string value by which the events are filtered.", )
    }
    # check if eventattributevalue is part of the variable
    if ( is.null(eventfilterAB) == FALSE){
      if ( length(grep(eventfilterAB, eventfiltervar)) == 0 ) {
        ##TODO: #deparse(substitute(eventtypevar))
        stop("Value '", eventfilterAB, "' is not an element of '", deparse(substitute(eventfiltervar)) , "'.") 
      }
    }
    if ( is.null(eventfilterAJ) == FALSE){
      if ( length(grep(eventfilterAJ, eventfiltervar)) == 0 ) {
        ##TODO: #deparse(substitute(eventtypevar))
        stop("Value '", eventfilterAJ, "' is not an element of '", deparse(substitute(eventfiltervar)) , "'.") 
      }
    }
    if ( is.null(eventfilterIB) == FALSE){
      if ( length(grep(eventfilterIB, eventfiltervar)) == 0 ) {
        ##TODO: #deparse(substitute(eventtypevar))
        stop("Value '", eventfilterIB, "' is not an element of '", deparse(substitute(eventfiltervar)) , "'.") 
      }
    }
    if ( is.null(eventfilterIJ) == FALSE){
      if ( length(grep(eventfilterIJ, eventfiltervar)) == 0 ) {
        ##TODO: #deparse(substitute(eventtypevar))
        stop("Value '", eventfilterIJ, "' is not an element of '", deparse(substitute(eventfiltervar)) , "'.") 
      }
    }
  }
  
  ## check event-var
  if(is.null(eventvar) == FALSE){
    if(length(unique(eventvar)) == 2){
      if( ( sort(unique(eventvar))[1] == 0 & sort(unique(eventvar))[2] == 1  ) == FALSE){
        stop('eventvar has to be a dummy variable with values 0 for non-events and
             1 for true events.')
      }
    }else{
      stop('eventvar has to be a dummy variable with values 0 for non-events and
             1 for true events.')
    }
  }
  
  ## cannot take parallel and progress bar
  if(isTRUE(inParallel) & isTRUE(showprogressbar)){
    stop('Cannot spit out progress of the function whilst running the 
         loop in parallel. Turn showprogressbar to FALSE.')
  }
  
  ## cannot have parallel without cluster
  if(isTRUE(inParallel) & is.null(cluster)){
    stop('By choosing to run the loop in parallel, you need to define a 
         cluster. For instance: makeCluster(12, type="FORK"). Alternatively, 
         hand over the number of nodes you would like to run the function on.')
  }
  
  ## check if variablename makes sense (no " " etc.)
  variablename <- gsub(" ", "", variablename, fixed = TRUE)
  
  ## create simple data set to be returned for degree calcuations with more than 1 output-variable
  ##TODO: should there be an event-id-variable?? => that would be useful here
  data.short <- data.frame(time)
  ## 
  result <- rep(NA, length(sender))
  
  ## calculate part of decay function
  xlog <- log(2)/halflife 
  
  ## prepare the data set
  if(is.null(dataPastEvents)){
    ## prepare the data set
    if(is.null(eventvar)){
      senderLoop <- sender
      targetLoop <- target
      weightLoop <- weight
      timeLoop <- time
      if(is.null(eventtypevar)){
        eventtypevar <- rep("1", length(sender)) # if not given, define at as 1 for each event
        eventtypevarLoop <- rep("1", length(sender))
        eventtypevalueLoop <- 'standard'
      }else{
        eventtypevarLoop <- eventtypevar
        eventtypevalueLoop <- eventtypevalue
      }
      if(is.null(eventfilterAB)){
        eventfiltervarABLoop <- rep("1", length(sender))
        eventfilterABLoop <- "1"
      }else{
        eventfiltervarABLoop <- eventfiltervar
        eventfilterABLoop <- eventfilterAB
      }
      if(is.null(eventfilterAJ)){
        eventfiltervarAJLoop <- rep("1", length(sender))
        eventfilterAJLoop <- "1"
      }else{
        eventfiltervarAJLoop <- eventfiltervar
        eventfilterAJLoop <- eventfilterAJ
      }
      if(is.null(eventfilterIB)){
        eventfiltervarIBLoop <- rep("1", length(sender))
        eventfilterIBLoop <- "1"
      }else{
        eventfiltervarIBLoop <- eventfiltervar
        eventfilterIBLoop <- eventfilterIB
      }
      if(is.null(eventfilterIJ)){
        eventfiltervarIJLoop <- rep("1", length(sender))
        eventfilterIJLoop <- "1"
      }else{
        eventfiltervarIJLoop <- eventfiltervar
        eventfilterIJLoop <- eventfilterIJ
      }
    }else{ # counting process data is used
      senderLoop <- sender[eventvar == 1]
      targetLoop <- target[eventvar == 1]
      weightLoop <- weight[eventvar == 1]
      timeLoop <- time[eventvar == 1]
      if(is.null(eventtypevar)){
        eventtypevar <- rep("1", length(sender))
        eventtypevarLoop <- rep("1", length(senderLoop))
        eventtypevalueLoop <- 'standard'
      }else{
        eventtypevarLoop <- eventtypevar[eventvar == 1]
        eventtypevalueLoop <- eventtypevalue
      }
      if(is.null(eventfilterAB)){
        eventfiltervarABLoop <- rep("1", length(sender)) # same length als sender
        eventfilterABLoop <- "1"
      }else{
        eventfiltervarABLoop <- eventfiltervar # same length as sender
        eventfilterABLoop <- eventfilterAB
      }
      if(is.null(eventfilterAJ)){
        eventfiltervarAJLoop <- rep("1", length(senderLoop))
        eventfilterAJLoop <- "1"
      }else{
        eventfiltervarAJLoop <- eventfiltervar[eventvar == 1]
        eventfilterAJLoop <- eventfilterAJ
      }
      if(is.null(eventfilterIB)){
        eventfiltervarIBLoop <- rep("1", length(senderLoop))
        eventfilterIBLoop <- "1"
      }else{
        eventfiltervarIBLoop <- eventfiltervar[eventvar == 1]
        eventfilterIBLoop <- eventfilterIB
      }
      if(is.null(eventfilterIJ)){
        eventfiltervarIJLoop <- rep("1", length(senderLoop))
        eventfilterIJLoop <- "1"
      }else{
        eventfiltervarIJLoop <- eventfiltervar[eventvar == 1]
        eventfilterIJLoop <- eventfilterIJ
      }
    }
  }else{ #dataPastEvents is defined: 
    ## prepare the data set
    senderLoop <- dataPastEvents[,2]
    targetLoop <- dataPastEvents[,3]
    weightLoop <- dataPastEvents[,4]
    timeLoop <- dataPastEvents[,1]
    if(is.null(eventtypevar)){
      eventtypevar <- rep("1", length(sender)) # if not given, define at as 1 for each event
      eventtypevarLoop <- rep("1", nrow(dataPastEvents))
      eventtypevalueLoop <- 'standard'
    }else{
      eventtypevarLoop <- dataPastEvents[,5]
      eventtypevalueLoop <- eventtypevalue
    }
    if(is.null(eventfilterAB)){
      eventfiltervarABLoop <- rep("1", nrow(dataPastEvents))
      eventfilterABLoop <- "1"
    }else{
      eventfiltervarABLoop <- dataPastEvents[,6]
      eventfilterABLoop <- eventfilterAB
    }
    if(is.null(eventfilterAJ)){
      eventfiltervarAJLoop <- rep("1", nrow(dataPastEvents))
      eventfilterAJLoop <- "1"
    }else{
      eventfiltervarAJLoop <- dataPastEvents[,6]
      eventfilterAJLoop <- eventfilterAJ
    }
    if(is.null(eventfilterIB)){
      eventfiltervarIBLoop <- rep("1", nrow(dataPastEvents))
      eventfilterIBLoop <- "1"
    }else{
      eventfiltervarIBLoop <- dataPastEvents[,6]
      eventfilterIBLoop <- eventfilterIB
    }
    if(is.null(eventfilterIJ)){
      eventfiltervarIJLoop <- rep("1", nrow(dataPastEvents))
      eventfilterIJLoop <- "1"
    }else{
      eventfiltervarIJLoop <- dataPastEvents[,6]
      eventfilterIJLoop <- eventfilterIJ
    }
  }
  
  ## calculate the fourCycle effects for each event
  ##
  if(isTRUE(inParallel)){
    
    ##
    doParallel::registerDoParallel(cluster)
    
    ##
    res <- foreach::foreach(i=1:length(sender), .combine=rbind)%dopar%{
      
      ## check if eventfilterAB is TRUE
      if(eventfiltervarABLoop[i] == eventfilterABLoop){
        
        ## get list of what a said in the past
        ## TODO: add timeLoop > beginTime => when you set the nr of events you can iterate back over, you have to incorporate this in these filter-functions:
        w <- targetLoop[senderLoop == sender[i] & 
                          targetLoop != target[i] &  
                          timeLoop < time[i] &
                          eventfiltervarAJLoop == eventfilterAJLoop]
        w <- unique(w)
        
        ## get list of who else said b in the past
        ## TODO: add timeLoop > beginTime => when you set the nr of events you can iterate back over, you have to incorporate this in these filter-functions:
        if(eventtypevalue == 'standard'){
          x <- senderLoop[targetLoop == target[i] & 
                            senderLoop != sender[i] & 
                            timeLoop < time[i] & 
                            eventfiltervarIBLoop == eventfilterIBLoop]
        }else if(eventtypevalue == 'positive'){
          x <- senderLoop[targetLoop == target[i] & 
                            senderLoop != sender[i] & 
                            timeLoop < time[i] & 
                            eventfiltervarIBLoop == eventfilterIBLoop &  
                            eventtypevarLoop == eventtypevar[i]]
        }else if(eventtypevalue == 'negative'){
          x <- senderLoop[targetLoop == target[i] & 
                            senderLoop != sender[i] & 
                            timeLoop < time[i] & 
                            eventfiltervarIBLoop == eventfilterIBLoop &  
                            eventtypevarLoop != eventtypevar[i]]
        }
        x <- unique(x)
        
        ## 
        if(length(x) == 0 | length(w) == 0){
          result <- 0
        }else{
          # find i in reduced data set
          iLoop <- length(timeLoop[timeLoop < time[i]]) #+ 1 - 1 # + 1 bc in the loop it's <; however cpp starts at 0, so -1
          
          ## cpp-loop
          result <- fourCycleCpp(senderLoop, sender[i], targetLoop, target[i], 
                                 eventtypevarLoop, eventtypevar[i], timeLoop, time[i], 
                                 weightLoop, xlog, eventfiltervarAJLoop, 
                                 eventfilterAJLoop, eventfiltervarIBLoop, 
                                 eventfilterIBLoop, eventfiltervarIJLoop, 
                                 eventfilterIJLoop, eventtypevalueLoop, 
                                 w, x, iLoop, 0) #hardcode begin-optin in fourCycle to 0 => TODO: implement begin
          
        }
      }else{ # if eventfilterAB[i] != eventfilterAB
        result <- 0
      }
      
      ## rbind the variable:
      result
    } #closes i-loop
    
    ## transform result variable
    result <- as.numeric(as.character(res))
    
  }else{ # run loop without parallelization
    
    ## 
    if(isTRUE(showprogressbar)){
      pb <- txtProgressBar(min = 1, max = length(sender), style = 3)
    }
    for(i in 1:length(sender)){
      if(isTRUE(showprogressbar)){
        setTxtProgressBar(pb, i)
      }
      
      ## check if eventfilterAB is TRUE
      if(eventfiltervarABLoop[i] == eventfilterABLoop){
        
        ## get list of what a said in the past
        ## TODO: add timeLoop > beginTime => when you set the nr of events you can iterate back over, you have to incorporate this in these filter-functions:
        w <- targetLoop[senderLoop == sender[i] & 
                          targetLoop != target[i] & 
                          timeLoop < time[i] & 
                          eventfiltervarAJLoop == eventfilterAJLoop]
        w <- unique(w)
        
        ## get list of who else said b in the past
        ## TODO: add timeLoop > beginTime => when you set the nr of events you can iterate back over, you have to incorporate this in these filter-functions:
        if(eventtypevalue == 'standard'){
          x <- senderLoop[targetLoop == target[i] & 
                            senderLoop != sender[i] & 
                            timeLoop < time[i] & 
                            eventfiltervarIBLoop == eventfilterIBLoop]
        }else if(eventtypevalue == 'positive'){
          x <- senderLoop[targetLoop == target[i] & 
                            senderLoop != sender[i] & 
                            timeLoop < time[i] & 
                            eventfiltervarIBLoop == eventfilterIBLoop & 
                            eventtypevarLoop == eventtypevar[i]]
        }else if(eventtypevalue == 'negative'){
          x <- senderLoop[targetLoop == target[i] & 
                            senderLoop != sender[i] & 
                            timeLoop < time[i] & 
                            eventfiltervarIBLoop == eventfilterIBLoop & 
                            eventtypevarLoop != eventtypevar[i]]
        }
        x <- unique(x)
        
        ## 
        if(length(x) == 0 | length(w) == 0){
          result[i] <- 0
        }else{
          # find i in reduced data set
          if(is.null(eventvar)){
            iLoop <- i-1 # bc cpp-loops start at 0 not 1
          }else{
            iLoop <- length(timeLoop[timeLoop < time[i]]) #+ 1 - 1 # + 1 bc in the loop it's <; however cpp starts at 0, so -1
          }
          
          ## cpp-loop
          result[i] <- fourCycleCpp(senderLoop, sender[i], targetLoop, target[i], 
                                    eventtypevarLoop, eventtypevar[i], timeLoop, time[i], 
                                    weightLoop, xlog, eventfiltervarAJLoop, 
                                    eventfilterAJLoop,eventfiltervarIBLoop, 
                                    eventfilterIBLoop,eventfiltervarIJLoop, 
                                    eventfilterIJLoop, eventtypevalueLoop, 
                                    w, x, iLoop, 0) #hardcode begin-optin in fourCycle to 0 => TODO: implement begin
          
        }
      }else{ # if eventfilterAB[i] != eventfilterAB
        result[i] <- 0
      }
    }# closes i-loop
  } # closes if no parallel
  
  ## return results
  ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
  if ( returnData == TRUE ) {
    ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
    data <- cbind(data, result)
    names(data)[length(data)] <- variablename
    ## return the data frame with the variable bound to it
    return(data)
  }else{ 
    ## only return the 1 inertia variable that was generated
    return(result)
  }
  
} # closing

################################################################################
##	Similarity calculation
################################################################################

similarityStat <- function(data, time, sender, target, 
                           senderOrTarget = "sender",
                           whichSimilarity = NULL, 
                           halflifeLastEvent = NULL, 
                           halflifeTimeBetweenEvents = NULL,
                           eventtypevar = NULL, 
                           eventfiltervar = NULL, 
                           eventfiltervalue = NULL,
                           eventvar = NULL,
                           variablename = "similarity", 
                           returnData = FALSE, 
                           dataPastEvents = NULL,
                           showprogressbar = FALSE, 
                           inParallel = FALSE, cluster = NULL){
	
  ####### check inputs
  ## check if sender input is available
  if ( is.null(sender) ) {
    stop("No 'sender' argument was provided.")
  }else{
    sender <- as.character(sender)
  }
  
  ## check if target input is available
  if ( is.null(target) ) {
    stop("No 'target' argument was provided.")
  }else{
    target <- as.character(target)
  }	
  
  ## check if event.sequence is well defined (numeric and ever-increasing)
  if ( is.null(time) ) {
    stop("No 'time' argument was provided.")
  }else{
    #test if time-var is in ascending order
    if ( is.unsorted(time) ) {
      stop("'", time, "' is not sorted. Sort data frame according to 
           the event sequence.")
    }
  }
  
  ## check if sender and target and time are of same length
  if ( length(sender) != length(target) ){
    stop("'sender' and 'target' are not of same length.")
  }
  if ( length(time) != length(sender) ){
    stop("'sender' and 'time' are not of same length.")
  }
  
  ## check if senderOrTarget is specified
  if ( (senderOrTarget == "sender" | senderOrTarget == "target") == FALSE ){
    stop("'senderOrTarget' not correctly specified. Choose either 'sender' 
         or 'target' for the respective similarity measure.")
  }
  
  ## check if whichSimilarity is specified
  if ( is.null(whichSimilarity) == FALSE ){
    if ( (whichSimilarity == "total" | whichSimilarity == "average" ) == FALSE ){
      stop("'whichSimilarity' not correctly specified. Choose either 'total' 
           or 'average' for the respective similarity measure or set it to 'NULL'.")
    }
  }
  
  ## average/total vs. with halflife
  if ( is.null(whichSimilarity)  & is.null(halflifeLastEvent) ){
    stop("Specify type of similarity measure - either by chosing 'total' or 
         'average' in whichSimilarity, or by specifying a value for 
         'halflifeLastEvent'.")
  }
  if ( is.null(whichSimilarity) == FALSE & is.null(halflifeLastEvent) ==FALSE ){
    stop("Cannot specify 'whichSimilarity' as well as 'halflifeLastEvent'.")
  }
  if ( is.null(halflifeLastEvent) & 
       is.null(halflifeTimeBetweenEvents)==FALSE){
    stop("Please also specify a value for 'halflifeLastEvent'.")
  }
  
  ## check if event-type inputs are available and correctly specified
  if ( !is.null(eventtypevar) ) {
    # check if sender and eventfiltervar are of same length
    if ( length(eventtypevar) != length(sender) ){
      stop("'sender' and 'eventtypevar' are not of same length.")
    }
    # transform eventtypevar
    eventtypevar <- as.character(eventtypevar)
    if ( length(unique(eventtypevar)) != 2 ){ 
      stop("'eventtypevar' is not a dummy variable. Other variable types are not yet supported")
    }
  }
  
  ## check if event-attribute inputs are available and correctly specified
  if ( is.null(eventfiltervar) == FALSE ) {
    # check if sender and eventfiltervar are of same length
    if ( length(eventfiltervar) != length(sender) ){
      stop("'sender' and 'eventfiltervar' are not of same length.")
    }
    # transform eventfiltervar
    eventfiltervar <- as.character(eventfiltervar)
    if ( is.null(eventfiltervalue) ){
      stop("No 'eventfiltervalue' provided. Provide a string value by 
           which the events are filtered.", )
    }
    # check if eventfiltervalue is part of the variable
    if ( is.null(eventfiltervalue) == FALSE){
      if ( length(grep(eventfiltervalue, eventfiltervar)) == 0 ) {
        ##TODO: #deparse(substitute(eventtypevar))
        stop("Value '", eventfiltervalue, "' is not an element of '", 
             deparse(substitute(eventfiltervar)) , "'.") 
      }
    }
  }
  
  ## check event-var
  if(is.null(eventvar) == FALSE){
    if(length(unique(eventvar)) == 2){
      if( ( sort(unique(eventvar))[1] == 0 & sort(unique(eventvar))[2] == 1  ) == FALSE){
        stop('eventvar has to be a dummy variable with values 0 for non-events and
             1 for true events.')
      }
    }else{
      stop('eventvar has to be a dummy variable with values 0 for non-events and
             1 for true events.')
    }
  }
  
  ## cannot take parallel and progress bar
  if(isTRUE(inParallel) & isTRUE(showprogressbar)){
    stop('Cannot spit out progress of the function whilst running the 
         loop in parallel. Turn showprogressbar to FALSE.')
  }
  
  ## cannot have parallel without cluster
  if(isTRUE(inParallel) & is.null(cluster)){
    stop('By choosing to run the loop in parallel, you need to define a 
         cluster. For instance: makeCluster(12, type="FORK"). Alternatively, 
         hand over the number of nodes you would like to run the function on.')
  }
  
  ## check if variablename makes sense (no " " etc.)
  variablename <- gsub(" ", "", variablename, fixed = TRUE)
  
  ## create simple data set to be returned for degree calcuations with more than 1 output-variable
  ##TODO: should there be an event-id-variable?? => that would be useful here
  data.short <- data.frame(time)
  result <- rep(NA, length(sender))
  
  ## calculate part of decay function
  xlog <- log(2)/halflifeLastEvent 
  
  ## prepare data set
  if(is.null(dataPastEvents)){
    if(is.null(eventvar)){
      # sender, target, time
      senderLoop <- sender
      targetLoop <- target
      timeLoop <- time
      # eventtypevar
      if(is.null(eventtypevar)){
        eventtypevarMirrorImage <- rep("1", length(sender)) # if not given, define at as 1 for each event
        eventtypevarLoop <- rep("1", length(sender))
      }else{
        eventtypevarLoop <- eventtypevar
      }
      # eventfiltervar
      if(is.null(eventfiltervar)){
        eventfiltervarLoop <- rep("1", length(sender))
        eventfiltervalueLoop <- "1"
      }else{
        eventfiltervarLoop <- eventfiltervar
        eventfiltervalueLoop <- eventfiltervalue
      }
    }else{ # if eventvar is specified
      # sender, target, time
      senderLoop <- sender[eventvar == 1]
      targetLoop <- target[eventvar == 1]
      timeLoop <- time[eventvar == 1]
      # eventtypevar
      if(is.null(eventtypevar)){
        eventtypevarMirrorImage <- rep("1", length(sender)) # if not given, define at as 1 for each event
        eventtypevarLoop <- rep("1", length(senderLoop))
      }else{
        eventtypevarLoop <- eventtypevar[eventvar == 1]
      }
      # eventfiltervar
      if(is.null(eventfiltervar)){
        eventfiltervarLoop <- rep("1", length(senderLoop))
        eventfiltervalueLoop <- "1"
      }else{
        eventfiltervarLoop <- eventfiltervar[eventvar == 1]
        eventfiltervalueLoop <- eventfiltervalue
      }
    }
  }else{ #data Past events are provided
    # sender, target, time
    senderLoop <- dataPastEvents[,2]
    targetLoop <- dataPastEvents[,3]
    timeLoop <- dataPastEvents[,1]
    # eventtypevar
    if(is.null(eventtypevar)){
      eventtypevarMirrorImage <- rep("1", length(sender)) # if not given, define at as 1 for each event
      eventtypevarLoop <- rep("1", length(senderLoop))
    }else{
      eventtypevarLoop <- dataPastEvents[,4]
    }
    # eventfiltervar
    if(is.null(eventfiltervar)){
      eventfiltervarLoop <- rep("1", length(senderLoop))
      eventfiltervalueLoop <- "1"
    }else{
      eventfiltervarLoop <- dataPastEvents[,5]
      eventfiltervalueLoop <- eventfiltervalue
    }
  }
  
  ####### calculate the similarity effects for each event
  ## 
  ## (2) halflife last event scenario and (3) 2 halflifes specified
  if(isTRUE(inParallel)){
    
    ##
    doParallel::registerDoParallel(cluster)
    
    ## i-loop start
    res <- foreach::foreach(i=1:length(sender), .combine=rbind)%dopar%{
      
      ## construct v and w
      # v: sender-sim: v = who else used b (match= in same way)
      # w: sender-sim: w = what else has a said?; 
      if( senderOrTarget == 'sender' ){
        v <- senderLoop[ senderLoop != sender[i] & 
                           targetLoop == target[i] & 
                           timeLoop < time[i] & 
                           eventtypevarLoop == eventtypevarMirrorImage[i] &
                           eventfiltervarLoop == eventfiltervalueLoop ]
        w <- targetLoop[ senderLoop == sender[i] & 
                           targetLoop != target[i] & 
                           timeLoop < time[i] & 
                           eventtypevarLoop == eventtypevarMirrorImage[i] &
                           eventfiltervarLoop == eventfiltervalueLoop ]
        
      }else{ # target sim:
        # v: target-sim: v = what else has a said?
        # w: target-sim: w = who else said b? (??match?? here too for target-sim?)
        v <- targetLoop[ senderLoop == sender[i] & 
                           targetLoop != target[i] & 
                           timeLoop < time[i] & 
                           eventtypevarLoop == eventtypevarMirrorImage[i] &
                           eventfiltervarLoop == eventfiltervalueLoop ]
        w <- v <- senderLoop[ senderLoop != sender[i] & 
                                targetLoop == target[i] & 
                                timeLoop < time[i] & 
                                eventtypevarLoop == eventtypevarMirrorImage[i] &
                                eventfiltervarLoop == eventfiltervalueLoop ]
      }
      v <- unique(v)
      w <- unique(w)
      
      ## check if v and w are empty
      if(length(v) == 0 | length(w) == 0){
        result <- 0
      }else{
        # find i in (reduced) data set
        
        iLoop <- length(timeLoop[timeLoop < time[i]]) #+ 1 - 1 # + 1 bc in the loop it's <; however cpp starts at 0, so -1
        
        ## hand over to cpp
        if(senderOrTarget == 'sender' & is.null(eventtypevar)){
          # 1. sender sim, no match
          if(is.null(whichSimilarity) & is.null(halflifeLastEvent) == FALSE &  is.null(halflifeTimeBetweenEvents)){
            result <- similaritySimpleCpp(senderLoop, sender[i], 
                                          targetLoop, target[i], 
                                          timeLoop, time[i], 
                                          xlog, eventfiltervarLoop, 
                                          eventfiltervalueLoop, 
                                          eventtypevarLoop, eventtypevarMirrorImage[i], 
                                          'nomatch', # or nomatch if no eventtypvar is specified
                                          'sender', # or target
                                          v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                          w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                          iLoop, 0) #begin = hardcoded to zero => TODO add this as option & make sure it is also true for v and w!
          }
          if(is.null(whichSimilarity) & is.null(halflifeLastEvent) == FALSE &  is.null(halflifeTimeBetweenEvents) == FALSE){
            result <- similarityComplexCpp(senderLoop, sender[i], 
                                           targetLoop, target[i], 
                                           timeLoop, time[i], 
                                           xlog, halflifeTimeBetweenEvents,
                                           eventfiltervarLoop, 
                                           eventfiltervalueLoop, 
                                           eventtypevarLoop, eventtypevarMirrorImage[i], 
                                           'nomatch', # or nomatch if no eventtypvar is specified
                                           'sender', # or target
                                           v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                           w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                           iLoop, 0)
          }
          if ( is.null(whichSimilarity) == FALSE ){
            if(whichSimilarity == 'total'){
              result <- similarityTotalAverageCpp(senderLoop, sender[i], 
                                                  targetLoop, target[i], 
                                                  timeLoop, time[i], 
                                                  eventfiltervarLoop, 
                                                  eventfiltervalueLoop, 
                                                  eventtypevarLoop, eventtypevarMirrorImage[i], 
                                                  'total',
                                                  'nomatch', # or nomatch if no eventtypvar is specified
                                                  'sender', # or target
                                                  v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                                  w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                                  iLoop, 0)
            }
            if(whichSimilarity == 'average'){
              result <- similarityTotalAverageCpp(senderLoop, sender[i], 
                                                  targetLoop, target[i], 
                                                  timeLoop, time[i], 
                                                  eventfiltervarLoop, 
                                                  eventfiltervalueLoop, 
                                                  eventtypevarLoop, eventtypevarMirrorImage[i], 
                                                  'average',
                                                  'nomatch', # or nomatch if no eventtypvar is specified
                                                  'sender', # or target
                                                  v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                                  w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                                  iLoop, 0)
            }
          }
        }
        if(senderOrTarget == 'target' & is.null(eventtypevar)){
          # 2. target sim, no match
          if(is.null(whichSimilarity) & is.null(halflifeLastEvent) == FALSE &  is.null(halflifeTimeBetweenEvents)){
            result <- similaritySimpleCpp(senderLoop, sender[i], 
                                          targetLoop, target[i], 
                                          timeLoop, time[i], 
                                          xlog, eventfiltervarLoop, 
                                          eventfiltervalueLoop, 
                                          eventtypevarLoop, eventtypevarMirrorImage[i], 
                                          'nomatch', # or nomatch if no eventtypvar is specified
                                          'target', # or target
                                          v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                          w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                          iLoop, 0) #begin = hardcoded to zero => TODO add this as option & make sure it is also true for v and w!
          }
          if(is.null(whichSimilarity) & is.null(halflifeLastEvent) == FALSE &  is.null(halflifeTimeBetweenEvents) == FALSE){
            result <- similarityComplexCpp(senderLoop, sender[i], 
                                           targetLoop, target[i], 
                                           timeLoop, time[i], 
                                           xlog, halflifeTimeBetweenEvents,
                                           eventfiltervarLoop, 
                                           eventfiltervalueLoop, 
                                           eventtypevarLoop, eventtypevarMirrorImage[i], 
                                           'nomatch', # or nomatch if no eventtypvar is specified
                                           'target', # or target
                                           v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                           w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                           iLoop, 0)
          }
          if ( is.null(whichSimilarity) == FALSE ){
            if(whichSimilarity == 'total'){
              result <- similarityTotalAverageCpp(senderLoop, sender[i], 
                                                  targetLoop, target[i], 
                                                  timeLoop, time[i], 
                                                  eventfiltervarLoop, 
                                                  eventfiltervalueLoop, 
                                                  eventtypevarLoop, eventtypevarMirrorImage[i], 
                                                  'total',
                                                  'nomatch', # or nomatch if no eventtypvar is specified
                                                  'target', # or target
                                                  v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                                  w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                                  iLoop, 0)
            }
            if(whichSimilarity == 'average'){
              result <- similarityTotalAverageCpp(senderLoop, sender[i], 
                                                  targetLoop, target[i], 
                                                  timeLoop, time[i], 
                                                  eventfiltervarLoop, 
                                                  eventfiltervalueLoop, 
                                                  eventtypevarLoop, eventtypevarMirrorImage[i], 
                                                  'average',
                                                  'nomatch', # or nomatch if no eventtypvar is specified
                                                  'target', # or target
                                                  v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                                  w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                                  iLoop, 0)
            }
          }
        }
        if(senderOrTarget == 'sender' & is.null(eventtypevar) == FALSE){
          # 1. sender sim, match
          if(is.null(whichSimilarity) & is.null(halflifeLastEvent) == FALSE &  is.null(halflifeTimeBetweenEvents)){
            result <- similaritySimpleCpp(senderLoop, sender[i], 
                                          targetLoop, target[i], 
                                          timeLoop, time[i], 
                                          xlog, eventfiltervarLoop, 
                                          eventfiltervalueLoop, 
                                          eventtypevarLoop, eventtypevarMirrorImage[i], 
                                          'match', # or nomatch if no eventtypvar is specified
                                          'sender', # or target
                                          v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                          w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                          iLoop, 0) #begin = hardcoded to zero => TODO add this as option & make sure it is also true for v and w!
          }
          if(is.null(whichSimilarity) & is.null(halflifeLastEvent) == FALSE &  is.null(halflifeTimeBetweenEvents) == FALSE){
            result <- similarityComplexCpp(senderLoop, sender[i], 
                                           targetLoop, target[i], 
                                           timeLoop, time[i], 
                                           xlog, halflifeTimeBetweenEvents,
                                           eventfiltervarLoop, 
                                           eventfiltervalueLoop, 
                                           eventtypevarLoop, eventtypevarMirrorImage[i], 
                                           'match', # or nomatch if no eventtypvar is specified
                                           'sender', # or target
                                           v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                           w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                           iLoop, 0)
          }
          if ( is.null(whichSimilarity) == FALSE ){
            if(whichSimilarity == 'total'){
              result <- similarityTotalAverageCpp(senderLoop, sender[i], 
                                                  targetLoop, target[i], 
                                                  timeLoop, time[i], 
                                                  eventfiltervarLoop, 
                                                  eventfiltervalueLoop, 
                                                  eventtypevarLoop, eventtypevarMirrorImage[i], 
                                                  'total',
                                                  'match', # or nomatch if no eventtypvar is specified
                                                  'sender', # or target
                                                  v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                                  w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                                  iLoop, 0)
            }
            if(whichSimilarity == 'average'){
              result <- similarityTotalAverageCpp(senderLoop, sender[i], 
                                                  targetLoop, target[i], 
                                                  timeLoop, time[i], 
                                                  eventfiltervarLoop, 
                                                  eventfiltervalueLoop, 
                                                  eventtypevarLoop, eventtypevarMirrorImage[i], 
                                                  'average',
                                                  'match', # or nomatch if no eventtypvar is specified
                                                  'sender', # or target
                                                  v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                                  w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                                  iLoop, 0)
            }
          }
        }
        if(senderOrTarget == 'target' & is.null(eventtypevar) == FALSE){
          # 2. target sim, match
          if(is.null(whichSimilarity) & is.null(halflifeLastEvent) == FALSE &  is.null(halflifeTimeBetweenEvents)){
            result <- similaritySimpleCpp(senderLoop, sender[i], 
                                          targetLoop, target[i], 
                                          timeLoop, time[i], 
                                          xlog, eventfiltervarLoop, 
                                          eventfiltervalueLoop, 
                                          eventtypevarLoop, eventtypevarMirrorImage[i], 
                                          'match', # or nomatch if no eventtypvar is specified
                                          'target', # or target
                                          v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                          w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                          iLoop, 0) #begin = hardcoded to zero => TODO add this as option & make sure it is also true for v and w!
          }
          if(is.null(whichSimilarity) & is.null(halflifeLastEvent) == FALSE &  is.null(halflifeTimeBetweenEvents) == FALSE){
            result <- similarityComplexCpp(senderLoop, sender[i], 
                                           targetLoop, target[i], 
                                           timeLoop, time[i], 
                                           xlog, halflifeTimeBetweenEvents,
                                           eventfiltervarLoop, 
                                           eventfiltervalueLoop, 
                                           eventtypevarLoop, eventtypevarMirrorImage[i], 
                                           'match', # or nomatch if no eventtypvar is specified
                                           'target', # or target
                                           v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                           w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                           iLoop, 0)
          }
          if ( is.null(whichSimilarity) == FALSE ){
            if(whichSimilarity == 'total'){
              result <- similarityTotalAverageCpp(senderLoop, sender[i], 
                                                  targetLoop, target[i], 
                                                  timeLoop, time[i], 
                                                  eventfiltervarLoop, 
                                                  eventfiltervalueLoop, 
                                                  eventtypevarLoop, eventtypevarMirrorImage[i], 
                                                  'total',
                                                  'match', # or nomatch if no eventtypvar is specified
                                                  'target', # or target
                                                  v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                                  w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                                  iLoop, 0)
            }
            if(whichSimilarity == 'average'){
              result <- similarityTotalAverageCpp(senderLoop, sender[i], 
                                                  targetLoop, target[i], 
                                                  timeLoop, time[i], 
                                                  eventfiltervarLoop, 
                                                  eventfiltervalueLoop, 
                                                  eventtypevarLoop, eventtypevarMirrorImage[i], 
                                                  'average',
                                                  'match', # or nomatch if no eventtypvar is specified
                                                  'target', # or target
                                                  v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                                  w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                                  iLoop, 0)
            }
          }
        }
      } # closes if w and v are not null
      
      ## rbind the variable:
      result
    } #closes i-loop
    
    ## transform result variable
    result <- as.numeric(as.character(res))
    
  }else{ # run loop without parallelization
    
    ## 
    if(isTRUE(showprogressbar)){
      pb <- txtProgressBar(min = 1, max = length(sender), style = 3)
    }
    for(i in 1:length(sender)){
      if(isTRUE(showprogressbar)){
        setTxtProgressBar(pb, i)
      }
      
      ## construct v and w
      # v: sender-sim: v = who else used b (match= in same way)
      # w: sender-sim: w = what else has a said?; 
      if( senderOrTarget == 'sender' ){
        v <- senderLoop[ senderLoop != sender[i] & 
                           targetLoop == target[i] & 
                           timeLoop < time[i] & 
                           eventtypevarLoop == eventtypevarMirrorImage[i] &
                           eventfiltervarLoop == eventfiltervalueLoop ]
        w <- targetLoop[ senderLoop == sender[i] & 
                           targetLoop != target[i] & 
                           timeLoop < time[i] & 
                           eventtypevarLoop == eventtypevarMirrorImage[i] &
                           eventfiltervarLoop == eventfiltervalueLoop ]
        
      }else{ # target sim:
        # v: target-sim: v = what else has a said?
        # w: target-sim: w = who else said b? (??match?? here too for target-sim?)
        v <- targetLoop[ senderLoop == sender[i] & 
                           targetLoop != target[i] & 
                           timeLoop < time[i] & 
                           eventtypevarLoop == eventtypevarMirrorImage[i] &
                           eventfiltervarLoop == eventfiltervalueLoop ]
        w <- v <- senderLoop[ senderLoop != sender[i] & 
                                targetLoop == target[i] & 
                                timeLoop < time[i] & 
                                eventtypevarLoop == eventtypevarMirrorImage[i] &
                                eventfiltervarLoop == eventfiltervalueLoop ]
      }
      v <- unique(v)
      w <- unique(w)
      
      ## check if v and w are empty
      if(length(v) == 0 | length(w) == 0){
        result[i] <- 0
      }else{
        # find i in reduced data set
        if(is.null(eventvar)){
          iLoop <- i-1 # bc cpp-loops start at 0 not 1
        }else{
          iLoop <- length(timeLoop[timeLoop < time[i]]) #+ 1 - 1 # + 1 bc in the loop it's <; however cpp starts at 0, so -1
        }
        
        ## hand over to cpp
        if(senderOrTarget == 'sender' & is.null(eventtypevar)){
          # 1. sender sim, no match
          if(is.null(whichSimilarity) & is.null(halflifeLastEvent) == FALSE &  is.null(halflifeTimeBetweenEvents)){
            result[i] <- similaritySimpleCpp(senderLoop, sender[i], 
                                             targetLoop, target[i], 
                                             timeLoop, time[i], 
                                             xlog, eventfiltervarLoop, 
                                             eventfiltervalueLoop, 
                                             eventtypevarLoop, eventtypevarMirrorImage[i], 
                                             'nomatch', # or nomatch if no eventtypvar is specified
                                             'sender', # or target
                                             v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                             w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                             iLoop, 0) #begin = hardcoded to zero => TODO add this as option & make sure it is also true for v and w!
          }
          if(is.null(whichSimilarity) & is.null(halflifeLastEvent) == FALSE &  is.null(halflifeTimeBetweenEvents) == FALSE){
            result[i] <- similarityComplexCpp(senderLoop, sender[i], 
                                              targetLoop, target[i], 
                                              timeLoop, time[i], 
                                              xlog, halflifeTimeBetweenEvents,
                                              eventfiltervarLoop, 
                                              eventfiltervalueLoop, 
                                              eventtypevarLoop, eventtypevarMirrorImage[i], 
                                              'nomatch', # or nomatch if no eventtypvar is specified
                                              'sender', # or target
                                              v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                              w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                              iLoop, 0)
          }
          if ( is.null(whichSimilarity) == FALSE ){
            if(whichSimilarity == 'total'){
              result[i] <- similarityTotalAverageCpp(senderLoop, sender[i], 
                                                     targetLoop, target[i], 
                                                     timeLoop, time[i], 
                                                     eventfiltervarLoop, 
                                                     eventfiltervalueLoop, 
                                                     eventtypevarLoop, eventtypevarMirrorImage[i], 
                                                     'total',
                                                     'nomatch', # or nomatch if no eventtypvar is specified
                                                     'sender', # or target
                                                     v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                                     w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                                     iLoop, 0)
            }
            if(whichSimilarity == 'average'){
              result[i] <- similarityTotalAverageCpp(senderLoop, sender[i], 
                                                     targetLoop, target[i], 
                                                     timeLoop, time[i], 
                                                     eventfiltervarLoop, 
                                                     eventfiltervalueLoop, 
                                                     eventtypevarLoop, eventtypevarMirrorImage[i], 
                                                     'average',
                                                     'nomatch', # or nomatch if no eventtypvar is specified
                                                     'sender', # or target
                                                     v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                                     w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                                     iLoop, 0)
            }
          }
        }
        if(senderOrTarget == 'target' & is.null(eventtypevar)){
          # 2. target sim, no match
          if(is.null(whichSimilarity) & is.null(halflifeLastEvent) == FALSE &  is.null(halflifeTimeBetweenEvents)){
            result[i] <- similaritySimpleCpp(senderLoop, sender[i], 
                                             targetLoop, target[i], 
                                             timeLoop, time[i], 
                                             xlog, eventfiltervarLoop, 
                                             eventfiltervalueLoop, 
                                             eventtypevarLoop, eventtypevarMirrorImage[i], 
                                             'nomatch', # or nomatch if no eventtypvar is specified
                                             'target', # or target
                                             v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                             w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                             iLoop, 0) #begin = hardcoded to zero => TODO add this as option & make sure it is also true for v and w!
          }
          if(is.null(whichSimilarity) & is.null(halflifeLastEvent) == FALSE &  is.null(halflifeTimeBetweenEvents) == FALSE){
            result[i] <- similarityComplexCpp(senderLoop, sender[i], 
                                              targetLoop, target[i], 
                                              timeLoop, time[i], 
                                              xlog, halflifeTimeBetweenEvents,
                                              eventfiltervarLoop, 
                                              eventfiltervalueLoop, 
                                              eventtypevarLoop, eventtypevarMirrorImage[i], 
                                              'nomatch', # or nomatch if no eventtypvar is specified
                                              'target', # or target
                                              v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                              w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                              iLoop, 0)
          }
          if ( is.null(whichSimilarity) == FALSE ){
            if(whichSimilarity == 'total'){
              result[i] <- similarityTotalAverageCpp(senderLoop, sender[i], 
                                                     targetLoop, target[i], 
                                                     timeLoop, time[i], 
                                                     eventfiltervarLoop, 
                                                     eventfiltervalueLoop, 
                                                     eventtypevarLoop, eventtypevarMirrorImage[i], 
                                                     'total',
                                                     'nomatch', # or nomatch if no eventtypvar is specified
                                                     'target', # or target
                                                     v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                                     w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                                     iLoop, 0)
            }
            if(whichSimilarity == 'average'){
              result[i] <- similarityTotalAverageCpp(senderLoop, sender[i], 
                                                     targetLoop, target[i], 
                                                     timeLoop, time[i], 
                                                     eventfiltervarLoop, 
                                                     eventfiltervalueLoop, 
                                                     eventtypevarLoop, eventtypevarMirrorImage[i], 
                                                     'average',
                                                     'nomatch', # or nomatch if no eventtypvar is specified
                                                     'target', # or target
                                                     v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                                     w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                                     iLoop, 0)
            }
          }
        }
        if(senderOrTarget == 'sender' & is.null(eventtypevar) == FALSE){
          # 1. sender sim, match
          if(is.null(whichSimilarity) & is.null(halflifeLastEvent) == FALSE &  is.null(halflifeTimeBetweenEvents)){
            result[i] <- similaritySimpleCpp(senderLoop, sender[i], 
                                             targetLoop, target[i], 
                                             timeLoop, time[i], 
                                             xlog, eventfiltervarLoop, 
                                             eventfiltervalueLoop, 
                                             eventtypevarLoop, eventtypevarMirrorImage[i], 
                                             'match', # or nomatch if no eventtypvar is specified
                                             'sender', # or target
                                             v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                             w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                             iLoop, 0) #begin = hardcoded to zero => TODO add this as option & make sure it is also true for v and w!
          }
          if(is.null(whichSimilarity) & is.null(halflifeLastEvent) == FALSE &  is.null(halflifeTimeBetweenEvents) == FALSE){
            result[i] <- similarityComplexCpp(senderLoop, sender[i], 
                                              targetLoop, target[i], 
                                              timeLoop, time[i], 
                                              xlog, halflifeTimeBetweenEvents,
                                              eventfiltervarLoop, 
                                              eventfiltervalueLoop, 
                                              eventtypevarLoop, eventtypevarMirrorImage[i], 
                                              'match', # or nomatch if no eventtypvar is specified
                                              'sender', # or target
                                              v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                              w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                              iLoop, 0)
          }
          if ( is.null(whichSimilarity) == FALSE ){
            if(whichSimilarity == 'total'){
              result[i] <- similarityTotalAverageCpp(senderLoop, sender[i], 
                                                     targetLoop, target[i], 
                                                     timeLoop, time[i], 
                                                     eventfiltervarLoop, 
                                                     eventfiltervalueLoop, 
                                                     eventtypevarLoop, eventtypevarMirrorImage[i], 
                                                     'total',
                                                     'match', # or nomatch if no eventtypvar is specified
                                                     'sender', # or target
                                                     v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                                     w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                                     iLoop, 0)
            }
            if(whichSimilarity == 'average'){
              result[i] <- similarityTotalAverageCpp(senderLoop, sender[i], 
                                                     targetLoop, target[i], 
                                                     timeLoop, time[i], 
                                                     eventfiltervarLoop, 
                                                     eventfiltervalueLoop, 
                                                     eventtypevarLoop, eventtypevarMirrorImage[i], 
                                                     'average',
                                                     'match', # or nomatch if no eventtypvar is specified
                                                     'sender', # or target
                                                     v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                                     w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                                     iLoop, 0)
            }
          }
        }
        if(senderOrTarget == 'target' & is.null(eventtypevar) == FALSE){
          # 2. target sim, match
          if(is.null(whichSimilarity) & is.null(halflifeLastEvent) == FALSE &  is.null(halflifeTimeBetweenEvents)){
            result[i] <- similaritySimpleCpp(senderLoop, sender[i], 
                                             targetLoop, target[i], 
                                             timeLoop, time[i], 
                                             xlog, eventfiltervarLoop, 
                                             eventfiltervalueLoop, 
                                             eventtypevarLoop, eventtypevarMirrorImage[i], 
                                             'match', # or nomatch if no eventtypvar is specified
                                             'target', # or target
                                             v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                             w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                             iLoop, 0) #begin = hardcoded to zero => TODO add this as option & make sure it is also true for v and w!
          }
          if(is.null(whichSimilarity) & is.null(halflifeLastEvent) == FALSE &  is.null(halflifeTimeBetweenEvents) == FALSE){
            result[i] <- similarityComplexCpp(senderLoop, sender[i], 
                                              targetLoop, target[i], 
                                              timeLoop, time[i], 
                                              xlog, halflifeTimeBetweenEvents,
                                              eventfiltervarLoop, 
                                              eventfiltervalueLoop, 
                                              eventtypevarLoop, eventtypevarMirrorImage[i], 
                                              'match', # or nomatch if no eventtypvar is specified
                                              'target', # or target
                                              v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                              w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                              iLoop, 0)
          }
          if ( is.null(whichSimilarity) == FALSE ){
            if(whichSimilarity == 'total'){
              result[i] <- similarityTotalAverageCpp(senderLoop, sender[i], 
                                                     targetLoop, target[i], 
                                                     timeLoop, time[i], 
                                                     eventfiltervarLoop, 
                                                     eventfiltervalueLoop, 
                                                     eventtypevarLoop, eventtypevarMirrorImage[i], 
                                                     'total',
                                                     'match', # or nomatch if no eventtypvar is specified
                                                     'target', # or target
                                                     v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                                     w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                                     iLoop, 0)
            }
            if(whichSimilarity == 'average'){
              result[i] <- similarityTotalAverageCpp(senderLoop, sender[i], 
                                                     targetLoop, target[i], 
                                                     timeLoop, time[i], 
                                                     eventfiltervarLoop, 
                                                     eventfiltervalueLoop, 
                                                     eventtypevarLoop, eventtypevarMirrorImage[i], 
                                                     'average',
                                                     'match', # or nomatch if no eventtypvar is specified
                                                     'target', # or target
                                                     v, # sender-sim: v = who else used b (match= in same way); target-sim: v = what else has a said?
                                                     w, # sender-sim: w = what else has a said?; target-sim: w = who else said b? (??match?? here too for target-sim?)
                                                     iLoop, 0)
            }
          }
        }
      } # closes if w and v are not null
      
    }# closes i-loop
  } # closes if no parallel
  
  ## return results
  ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
  if ( returnData == TRUE ) {
    ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
    data <- cbind(data, result)
    names(data)[length(data)] <- variablename
    ## return the data frame with the variable bound to it
    return(data)
  }else{ 
    ## only return the 1 inertia variable that was generated
    return(result)
  }
  
} # closing


################################################################################
##  Create event sequence
################################################################################

eventSequence <- function(datevar, dateformat = NULL, data = NULL,
                          type = "continuous", byTime = "daily",
                          excludeDate = NULL, excludeTypeOfDay = NULL,
                          excludeYear = NULL, excludeFrom = NULL, 
                          excludeTo = NULL, returnData = FALSE, 
                          sortData = FALSE, returnDateSequenceData = FALSE){
  
  #### check if all the inputs are correct
  ## check if date and dateformat match => then create Date-object
  if (type == "continuous"){
    date <- as.Date(as.character(datevar), format = dateformat)
    if (is.na(date[1])){
      stop("'dateformat' does not match structure of 'datevar'.")
    }
  }
  
  ## cannot specify "ordinal" and then exclude variables!
  if ( type == "ordinal"){
    if ( is.null(excludeDate) == FALSE | is.null(excludeTypeOfDay) == FALSE |
         is.null(excludeYear) == FALSE | is.null(excludeFrom) == FALSE |
         is.null(excludeTo) == FALSE ){
      stop("Cannot exclude dates if type 'ordinal' is selected.")
    }
  }
  
  ## cannot specify returnData and not provide data-object
  if (returnData == TRUE & is.null(data)){
    stop("Provide 'data'-element if 'returnData = TRUE' is set.")
  }
  
  ## check if both excludeFrom and excludeTO are set
  if ( (is.null(excludeFrom) == FALSE & is.null(excludeTo)) | 
       (is.null(excludeFrom) & is.null(excludeTo)==FALSE) ){
    stop("Both 'excludeFrom' and 'excludeTo' variables need to be specified.")
  }
  
  ## check if excludeFrom and excludeTo are of equal length
  if ( is.null(excludeFrom) == FALSE & is.null(excludeTo) == FALSE & 
       length(excludeFrom) != length(excludeTo) ){
    stop("Both 'excludeFrom' and 'excludeTo' variables need to be of the same length.")
  }
  
  ## check if excludeFrom is smaller than excludeTo
  if ( is.null(excludeFrom) == FALSE & is.null(excludeTo) == FALSE) {
    for (h in length(excludeFrom)){
      if (as.Date(excludeFrom[h], format = dateformat) > 
          as.Date(excludeTo[h], format = dateformat) ){
        stop("'excludeFrom' is smaller than 'excludeTo'.")
      }	    
    }
  }
  
  ## check if type is specified correctly (either "continuous" or "ordinal")
  if ( type == "continuous" | type == "ordinal"){  
  }else{
    stop("'type' not specified correctly. Choose 'continuous' or 'ordinal' as
         event sequence options.")
  }
  
  ## send warning: if you choose "order" but not "returnData"
  if ( sortData == TRUE & returnData == FALSE){
    warning("Do not match output variable directly to your data if your data is
            not sorted. Choose 'returnData = TRUE' and 'orderData = TRUE'
            to make sure your data set has not been tarnished.")
  }
  
  ## check byTime
  if (is.na(byTime)){
    stop("'byTime' is not specified correctly. Use 'daily', 'monthly', 'yearly'.")
  }
  if (byTime == "daily" | byTime == "monthly" | byTime == "yearly"){
  }else{
    stop("'byTime' is not specified correctly. Use 'daily', 'monthly', 'yearly'.")
  }
  
  ## cannot specify returnData = TRUE and returnDateSequenceData = TRUE
  if(isTRUE(returnData) & isTRUE(returnDateSequenceData)){
    stop('cannot return both dataset and sequence set. Choose one or the other.')
  }
  
  #### create continuous event sequence
  if (type == "continuous"){   
    
    ## create artificial sequence from start to end in days
    sequence <- data.frame(seq(min(date), max(date), "1 day"))
    names(sequence) <- "date.sequence"
    
    ## erase whatever
    ## TODO: exclude only these, that do not delete entire data set! 
    if ( is.null(excludeDate) == FALSE){
      ## exlcude some of them
      for (i in excludeDate){
        sequence <- subset(sequence, sequence$date.sequence != as.Date(i, dateformat))
      }
    }
    
    if ( is.null(excludeTypeOfDay) == FALSE){
      ## add type of day
      ## depending on your locale, the weekdays will be named different, 
      ## http://stackoverflow.com/questions/17031002/get-weekdays-in-english-in-rstudio
      sequence$weekday <- weekdays(sequence$date.sequence)
      ## exlcude some of them
      for (i in excludeTypeOfDay){
        sequence <- subset(sequence, sequence$weekday != i)
      }
    }
    
    if ( is.null(excludeYear) == FALSE){
      ## add year to sequence    
      sequence$year <- as.numeric(format(sequence$date.sequence, "%Y"))
      for (i in excludeYear){
        sequence <- subset(sequence, sequence$year != i)
      }
    }
    
    if ( is.null(excludeFrom)==FALSE & is.null(excludeTo)==FALSE ){
      sequence$erase <- 0
      for (k in 1:length(excludeFrom)){
        sequence$erase <- ifelse((sequence$date.sequence >= 
                                    as.Date(excludeFrom[k], format = dateformat) & 
                                    sequence$date.sequence <= 
                                    as.Date(excludeTo[k], format = dateformat)), 1, sequence$erase)
      }
      sequence <- subset(sequence, sequence$erase != 1)
    }
    
    ## give artificial sequence 1:length()
    if (byTime == "daily"){
      sequence$event.sequence <- 1:length(sequence$date.sequence)
    }
    if (byTime == "monthly"){
      sequence$months.year <- as.character(format(sequence$date.sequence, "%m%Y"))
      sequence$event.sequence <- as.numeric(as.factor(sequence$months.year))
    }
    if (byTime == "yearly"){
      sequence$year <- as.character(format(sequence$date.sequence, "%Y"))
      sequence$event.sequence <- as.numeric(as.factor(sequence$year))
    }
    
    ## match with datevar
    result <- sequence$event.sequence[match(date, sequence$date.sequence)]
    
    ## return sequence
    if(isTRUE(returnDateSequenceData)){
      return(sequence)
      message('The output represents a data.frame with the date in the first
              column and the corresponding event sequence value in the second
              column. To return the eventSequence calculations, set 
              returnDateSequenceData = FALSE.', domain = NULL, appendLF = TRUE)
    }else{
      ## return data
      if ( returnData == FALSE){
        if ( sortData == FALSE){
          return(result)
        }else{
          result <- sort(result)
          return(result)
        }     
      }else{
        ## unsorted:
        data <- cbind(data, result)
        
        ## remove previous event-seq-variables
        if ( length(unique(grepl("event.seq.cont", names(data)))) == 2 ){
          data <- data[,-which(colnames(data) == "event.seq.cont")]
          #data <- subset(data, select = -c(event.seq.cont))
        }
        
        ## 
        names(data)[length(data)] <- 'event.seq.cont'
        if ( sortData == FALSE ){
          return(data)
        }else{ ## sorted:
          data <- data[order(data$event.seq.cont), ]
          return(data)
        }
      }
    }
    
    #### create ordinal event sequence 
  }else if (type == "ordinal"){
    ## create ordered data-sequence    
    sequence <- data.frame(datevar)
    ## if no dateformat is specified: is there no need for it? test it and report.
    if ( is.null(dateformat) ){
      temp <- as.numeric(as.character(sequence$datevar))
      if ( is.na(temp[1])){
        stop("'datevar' is not numeric. Provide numeric 'datevar' or use 
             dateformat to specify the datevar correctly.")
      }
      sequence <- sequence[order(sequence$datevar), ]
      sequence$date.sequence <- datevar
    }else{
      sequence$date.sequence <- as.Date(as.character(sequence$datevar), 
                                        format = dateformat)
      sequence <- sequence[order(sequence$date.sequence), ]
    }
    
    ## create ordinal sequence on ordered seq
    sequence$event.sequence <- as.numeric(as.factor(sequence$date.sequence))
    
    ## match with datevar
    result <- sequence$event.sequence[match(datevar, sequence$datevar)]
    
    ## return sequence
    if(isTRUE(returnDateSequenceData)){
      print('fuck this')
      return(sequence)
      message('The output represents a data.frame with the date in the first
              column and the corresponding event sequence value in the second
              column. To return the eventSequence calculations, set 
              returnDateSequenceData = FALSE.', domain = NULL, appendLF = TRUE)
    }else{
      ## return data
      if ( returnData == FALSE){
        if ( sortData == FALSE){
          return(result)
        }else{
          result <- sort(result)
          return(result)
        }
      }else{
        ## unsorted:
        data <- cbind(data, result)
        
        ## remove previous event-seq-variables
        if ( length(unique(grepl("event.seq.ord", names(data)))) == 2 ){
          data <- data[,-which(colnames(data) == "event.seq.ord")]
          #data <- subset(data, select = -c(event.seq.ord))
        }
        
        names(data)[length(data)] <- 'event.seq.ord'
        if ( sortData == FALSE ){
          return(data)
        }else{ ## sorted:
          data <- data[order(data$event.seq.ord), ]
          return(data)
        }
      }#closes if-else returnData    
    }#closes type == ordinal
  }
}

################################################################################
##	Reciprocity (one-mode statistic)
################################################################################

reciprocityStat <- function(data, time, sender, target, halflife, 
                            weight = NULL, eventtypevar = NULL, 
                            eventtypevalue = "valuematch",
                            eventfiltervar = NULL, eventfiltervalue = NULL, 
                            eventvar = NULL, variablename = "recip", 
                            returnData = FALSE, showprogressbar = FALSE, 
                            inParallel = FALSE, cluster = NULL){
  
  ####### check inputs
  ## check if sender and target inputs are available
  if ( is.null(sender) ) {
    stop("No 'sender' argument was provided.")
  }else{
    sender <- as.character(sender)
  }
  
  if ( is.null(target) ) {
    stop("No 'target' argument was provided.")
  }else{
    target <- as.character(target)
  }
  
  ## check if all variables have the same length.
  if (length(sender) != length(target)){
    stop("'sender' and 'target' are not of the same length.")
  }
  
  ## check if event.sequence is well defined (numeric and ever-increasing)
  if ( is.null(time) ) {
    stop("No 'time' argument was provided.")
  }else{
    #test if weight-var is in ascending order
    if ( is.unsorted(time) ) {
      stop("'time' is not sorted. Sort data frame according to the event sequence.")
    }
  }
  
  ## check if time has the requested length
  if (length(sender) != length(time)){
    stop("'sender' and 'time' are not of the same length.")
  }
  
  ## check if weight-var is defined (if not -> create it)
  if ( is.null(weight) ) {
    weight <- rep(1, length(time))
  }
  if ( !is.numeric(weight) ) {
    stop("'weight' variable is not numeric.") #TODO: deparse(substitute(eventattributevar)) ?
  }
  
  ## check if event-type inputs are available and correctly specified
  if ( !is.null(eventtypevar) ) {
    # check length
    if (length(sender) != length(eventtypevar)){
      stop("'eventtypevar' and 'sender' are not of the same length.")
    }
    # transform
    eventtypevar <- as.character(eventtypevar)
    if ( is.null(eventtypevalue) ){
      stop("No 'eventtypevalue' provided. Use default 'valuematch', or 'valuemix' or string value(s) to determine by which values the events should be filtered.", )
    }
    # check if eventtypevalue is part of the variable
    if ( length(eventtypevalue) > 1  ){
      for ( i in 1:length(eventtypevalue) ){
        if ( length(grep(eventtypevalue[i], eventtypevar)) == 0 ) {
          ##TODO: #deparse(substitute(eventtypevar))
          stop("Value '", eventtypevalue[i], "' is not an element of 'eventtypevar'.") ##deparse(substitute(eventtypevar))
        }
      }#closes i-loop  
      if ( length(unique(duplicated(eventtypevalue))) == 2 ) {
        stop("Duplicate values in 'eventtypevalue'.") 
      }
    }else if ( eventtypevalue != "valuematch" & eventtypevalue != "valuemix" ) {
      if ( length(grep(eventtypevalue, eventtypevar)) == 0 ) {
        ##TODO: #deparse(substitute(eventtypevar))
        stop("Value '", eventtypevalue, "' is not an element of 'eventtypevar'.") ##deparse(substitute(eventtypevar))
      }
    }	
  }
  
  ## check if event-attribute inputs are available and correctly specified
  if ( is.null(eventfiltervar) == FALSE ) {
    # check length
    if (length(sender) != length(eventfiltervar)){
      stop("'eventfiltervar' and 'sender' are not of the same length.")
    }
    # transform
    eventfiltervar <- as.character(eventfiltervar)
    if ( is.null(eventfiltervalue) ){
      stop("No 'eventfiltervalue' provided. Which value should be filtered for?", )
    }
    # check if eventattributevalue is part of the variable
    if ( length(eventfiltervalue) > 0 ){
      for ( i in 1:length(eventfiltervalue) ){
        if ( length(grep(eventfiltervalue[i], eventfiltervar)) == 0 ) {
          stop("Value '", eventfiltervalue[i], "' is not an element of 'eventfiltervar'.")  ##deparse(substitute(eventattributevar))
        }
      }#closes i-loop  
      if ( length(unique(duplicated(eventfiltervalue))) == 2 ) {
        stop("Duplicate values in 'eventfiltervalue'.") 
      }
    }
  }
  
  ## check event-var
  if(is.null(eventvar) == FALSE){
    if(length(unique(eventvar)) == 2){
      if( ( sort(unique(eventvar))[1] == 0 & sort(unique(eventvar))[2] == 1  ) == FALSE){
        stop('eventvar has to be a dummy variable with values 0 for non-events and
             1 for true events.')
      }
    }else{
      stop('eventvar has to be a dummy variable with values 0 for non-events and
             1 for true events.')
    }
  }
  
  ## cannot take parallel and progress bar
  if(isTRUE(inParallel) & isTRUE(showprogressbar)){
    stop('Cannot spit out progress of the function whilst running the 
         loop in parallel. Turn showprogressbar to FALSE.')
  }
  
  ## cannot have parallel without cluster
  if(isTRUE(inParallel) & is.null(cluster)){
    stop('By choosing to run the loop in parallel, you need to define a 
         cluster. For instance: makeCluster(12, type="FORK"). Alternatively, 
         hand over the number of nodes you would like to run the function on.')
  }
  
  ## check if variablename makes sense (no " " etc.)
  variablename <- gsub(" ", "", variablename, fixed = TRUE)
  
  ## create simple data set to be returned for inertia calcuations with more than 1 output-variable
  ##TODO: should there be an event-id-variable?? => that would be useful here
  data.short <- data.frame(time)
  
  ## calculate part of decay function
  xlog <- log(2)/halflife 
  
  ####### calculate stat
  
  ## use event-filter if counting process data is used
  if(is.null(eventvar)){
    countingProcessVar <- rep(1, length(sender))
  }else{
    countingProcessVar <- eventvar
  }
  
  ## 
  result <- rep(NA, length(sender))
  
  
  ## calculate the inertia effects for each event
  # (1): no type, no filter
  # (2): no type, with filter
  # (3): valuematch, no filter
  # (4): valuematch, with filter
  # (5): valuemix/values provided, no filter
  # (6): valuemix/values provided, with filter
  if ( is.null(eventtypevar) ) {
    if ( is.null(eventfiltervar) ) {
      ################ (1) start off with simple inertia function: no type, no filter
      # run in parallel?
      if(isTRUE(inParallel)){
        
        ##
        doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
        res <- foreach::foreach(i=1:length(sender), .combine=rbind)%dopar%{
          
          ## create vector of times, sender--target tie was made before
          timesIntertiaPast <- time[sender == target[i] & target == sender[i] & 
                                      time < time[i] & countingProcessVar == 1]
          ## get weight
          if(is.null(weight)){
            weightInertiaPast <- rep(1, length(timesIntertiaPast))
          }else{
            weightInertiaPast <- weight[sender == target[i] & target == sender[i] & 
                                          time < time[i] & countingProcessVar == 1]
          }
          
          ## run cpp-loop for all times
          result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                            weightInertiaPast) 
          ## rbind the variable:
          result 
        }# closes foreach-loop
        
        ## transform result variable
        result <- as.numeric(as.character(res))
        
      }else{ # run in standard form
        
        if(isTRUE(showprogressbar)){
          pb <- txtProgressBar(min = 1, max = length(sender), style = 3)
        }
        for(i in 1:length(sender)){
          if(isTRUE(showprogressbar)){
            setTxtProgressBar(pb, i)
          }
          
          
          ## create vector of times, sender--target tie was made before
          timesIntertiaPast <- time[sender == target[i] & target == sender[i] & 
                                      time < time[i] & countingProcessVar == 1]
          ## get weight
          if(is.null(weight)){
            weightInertiaPast <- rep(1, length(timesIntertiaPast))
          }else{
            weightInertiaPast <- weight[sender == target[i] & target == sender[i] & 
                                          time < time[i] & countingProcessVar == 1]
          }
          
          ## run cpp-loop for all times
          result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                               weightInertiaPast) 
          
        } #closes i-loop
      } # closes if no parallel
      
      ## return results
      ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
      if ( returnData == TRUE ) {
        ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 inertia variable that was generated
        return(result)
      }
      ################ done (1)
      
    }else{ # if eventfiltervar = given
      ################ (2) no type,  with filter
      # run in parallel?
      if(isTRUE(inParallel)){
        
        ##
        doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
        res <- foreach::foreach(i=1:length(sender), .combine=rbind)%dopar%{
          
          ## create vector of times, sender--target tie was made before
          timesIntertiaPast <- time[sender == target[i] & target == sender[i] & 
                                      time < time[i] & countingProcessVar == 1 &
                                      eventfiltervar == eventfiltervalue]
          ## get weight
          if(is.null(weight)){
            weightInertiaPast <- rep(1, length(timesIntertiaPast))
          }else{
            weightInertiaPast <- weight[sender == target[i] & target == sender[i] & 
                                          time < time[i] & countingProcessVar == 1 &
                                          eventfiltervar == eventfiltervalue]
          }
          
          ## run cpp-loop for all times
          result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                            weightInertiaPast) 
          ## rbind the variable:
          result 
        }# closes foreach-loop
        
        ## transform result variable
        result <- as.numeric(as.character(res))
        
      }else{ # run in standard form
        
        if(isTRUE(showprogressbar)){
          pb <- txtProgressBar(min = 1, max = length(sender), style = 3)
        }
        for(i in 1:length(sender)){
          if(isTRUE(showprogressbar)){
            setTxtProgressBar(pb, i)
          }
          
          
          ## create vector of times, sender--target tie was made before
          timesIntertiaPast <- time[sender == target[i] & target == sender[i] & 
                                      time < time[i] & countingProcessVar == 1 &
                                      eventfiltervar == eventfiltervalue]
          ## get weight
          if(is.null(weight)){
            weightInertiaPast <- rep(1, length(timesIntertiaPast))
          }else{
            weightInertiaPast <- weight[sender == target[i] & target == sender[i] & 
                                          time < time[i] & countingProcessVar == 1 &
                                          eventfiltervar == eventfiltervalue]
          }
          
          ## run cpp-loop for all times
          result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                               weightInertiaPast) 
          
        } #closes i-loop
      } # closes if no parallel
      
      ## return results
      ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
      if ( returnData == TRUE ) {
        ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
        data <- cbind(data, result)
        names(data)[length(data)] <- paste(variablename, "filtered", sep = ".")
        
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 inertia variable that was generated
        return(result)
      }
      ################ done (2)
    } # closes if-eventfiltervar == given
    
  }else{ # if eventtypevar = given
    if(eventtypevalue == 'valuematch'){
      if ( is.null(eventfiltervar) ) {
        
        ################ ## (3) valuematch, no filter
        # run in parallel?
        if(isTRUE(inParallel)){
          
          ##
          doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
          res <- foreach::foreach(i=1:length(sender), .combine=rbind)%dopar%{
            
            
            ## create vector of times, sender--target tie was made before
            timesIntertiaPast <- time[sender == target[i] & target == sender[i] & 
                                        time < time[i] & countingProcessVar == 1 &
                                        eventtypevar == eventtypevar[i]]
            ## get weight
            if(is.null(weight)){
              weightInertiaPast <- rep(1, length(timesIntertiaPast))
            }else{
              weightInertiaPast <- weight[sender == target[i] & target == sender[i] & 
                                            time < time[i] & countingProcessVar == 1 &
                                            eventtypevar == eventtypevar[i]]
            }
            
            ## run cpp-loop for all times
            result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                              weightInertiaPast) 
            ## rbind the variable:
            result 
          }# closes foreach-loop
          
          ## transform result variable
          result <- as.numeric(as.character(res))
          
        }else{ # run in standard form
          
          if(isTRUE(showprogressbar)){
            pb <- txtProgressBar(min = 1, max = length(sender), style = 3)
          }
          for(i in 1:length(sender)){
            if(isTRUE(showprogressbar)){
              setTxtProgressBar(pb, i)
            }
            
            
            ## create vector of times, sender--target tie was made before
            timesIntertiaPast <- time[sender == target[i] & target == sender[i] & 
                                        time < time[i] & countingProcessVar == 1 &
                                        eventtypevar == eventtypevar[i]]
            ## get weight
            if(is.null(weight)){
              weightInertiaPast <- rep(1, length(timesIntertiaPast))
            }else{
              weightInertiaPast <- weight[sender == target[i] & target == sender[i] & 
                                            time < time[i] & countingProcessVar == 1 &
                                            eventtypevar == eventtypevar[i]]
            }
            
            ## run cpp-loop for all times
            result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                 weightInertiaPast) 
            
          } #closes i-loop
        } # closes if no parallel
        
        ## return results
        ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
        if ( returnData == TRUE ) {
          ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "typematch", sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 inertia variable that was generated
          return(result)
        }
        ################ done (3)
        
      }else{ # if eventfiltervar = given
        
        ################ ## (4) valuematch, with filter
        # run in parallel?
        if(isTRUE(inParallel)){
          
          ##
          doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
          res <- foreach::foreach(i=1:length(sender), .combine=rbind)%dopar%{
            
            ## create vector of times, sender--target tie was made before
            timesIntertiaPast <- time[sender == target[i] & target == sender[i] & 
                                        time < time[i] & countingProcessVar == 1 &
                                        eventfiltervar == eventfiltervalue & 
                                        eventtypevar == eventtypevar[i]]
            ## get weight
            if(is.null(weight)){
              weightInertiaPast <- rep(1, length(timesIntertiaPast))
            }else{
              weightInertiaPast <- weight[sender == target[i] & target == sender[i] & 
                                            time < time[i] & countingProcessVar == 1 &
                                            eventfiltervar == eventfiltervalue & 
                                            eventtypevar == eventtypevar[i]]
            }
            
            ## run cpp-loop for all times
            result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                              weightInertiaPast) 
            ## rbind the variable:
            result 
          }# closes foreach-loop
          
          ## transform result variable
          result <- as.numeric(as.character(res))
          
        }else{ # run in standard form
          
          if(isTRUE(showprogressbar)){
            pb <- txtProgressBar(min = 1, max = length(sender), style = 3)
          }
          for(i in 1:length(sender)){
            if(isTRUE(showprogressbar)){
              setTxtProgressBar(pb, i)
            }
            
            
            ## create vector of times, sender--target tie was made before
            timesIntertiaPast <- time[sender == target[i] & target == sender[i] & 
                                        time < time[i] & countingProcessVar == 1 &
                                        eventfiltervar == eventfiltervalue & 
                                        eventtypevar == eventtypevar[i]]
            ## get weight
            if(is.null(weight)){
              weightInertiaPast <- rep(1, length(timesIntertiaPast))
            }else{
              weightInertiaPast <- weight[sender == target[i] & target == sender[i] & 
                                            time < time[i] & countingProcessVar == 1 &
                                            eventfiltervar == eventfiltervalue &
                                            eventtypevar == eventtypevar[i]]
            }
            
            ## run cpp-loop for all times
            result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                 weightInertiaPast) 
            
          } #closes i-loop
        } # closes if no parallel
        
        ## return results
        ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
        if ( returnData == TRUE ) {
          ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "typematch.filtered", sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 inertia variable that was generated
          return(result)
        }
        ################ done (4)
        
      }
      
    }else{
      ## create unique values - either for a valuemix variable or a eventtypevalues variable
      if(eventtypevalue == 'valuemix'){
        uniqueEventTypeValues <- unique(eventtypevar)
      }else{
        uniqueEventTypeValues <- eventtypevalue
      }
      
      if ( is.null(eventfiltervar) ) {
        ## (5) valuemix/values provided, no filter
        
        for (a in uniqueEventTypeValues){ #current event type
          for (b in uniqueEventTypeValues){ #past event type
            if ( a != b ){
              
              # run in parallel?
              if(isTRUE(inParallel)){
                
                ##
                doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
                res <- foreach::foreach(i=1:length(sender), .combine=rbind)%dopar%{
                  
                  ## 
                  if (eventtypevar[i] == a){
                    
                    ## create vector of times, sender--target tie was made before
                    timesIntertiaPast <- time[sender == target[i] & target == sender[i] & 
                                                time < time[i] & countingProcessVar == 1 & 
                                                eventtypevar == b]
                    ## get weight
                    if(is.null(weight)){
                      weightInertiaPast <- rep(1, length(timesIntertiaPast))
                    }else{
                      weightInertiaPast <- weight[sender == target[i] & target == sender[i] & 
                                                    time < time[i] & countingProcessVar == 1 &
                                                    eventtypevar == b]
                    }
                    
                    ## run cpp-loop for all times
                    result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                      weightInertiaPast) 
                  }else{
                    result <- 0
                  }
                  ## rbind the variable:
                  result 
                }# closes foreach-loop
                
                ## transform result variable
                result <- as.numeric(as.character(res))
                
              }else{ # run in standard form
                
                if(isTRUE(showprogressbar)){
                  pb <- txtProgressBar(min = 1, max = length(sender), style = 3)
                }
                for(i in 1:length(sender)){
                  if(isTRUE(showprogressbar)){
                    setTxtProgressBar(pb, i)
                  }
                  
                  ##
                  if(eventtypevar[i] == a){
                    
                    
                    ## create vector of times, sender--target tie was made before
                    timesIntertiaPast <- time[sender == target[i] & target == sender[i] & 
                                                time < time[i] & countingProcessVar == 1 &
                                                eventtypevar == b]
                    ## get weight
                    if(is.null(weight)){
                      weightInertiaPast <- rep(1, length(timesIntertiaPast))
                    }else{
                      weightInertiaPast <- weight[sender == target[i] & target == sender[i] & 
                                                    time < time[i] & countingProcessVar == 1 &
                                                    eventtypevar == b]
                    }
                    
                    ## run cpp-loop for all times
                    result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                         weightInertiaPast) 
                  }else{ #if eventtypevar != a
                    result[i] <- 0
                  }
                  
                } #closes i-loop
              } # closes if no parallel
              
              ##
              data.short <- cbind(data.short, result)
              
              ## change the name of the data.short's last column
              names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventattributevar))
                                                             gsub(" ", "", a, fixed = TRUE), #represents current event-actor type
                                                             gsub(" ", "", b, fixed = TRUE), #represents actor type of past actions
                                                             sep = ".")									
            }
          } # closes b-loop
          
          ## calculate effect where both events are of the same type (both a)
          # run in parallel?
          if(isTRUE(inParallel)){
            
            ##
            doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
            res <- foreach::foreach(i=1:length(sender), .combine=rbind)%dopar%{
              
              ## 
              if (eventtypevar[i] == a){
                
                ## create vector of times, sender--target tie was made before
                timesIntertiaPast <- time[sender == target[i] & target == sender[i] & 
                                            time < time[i] & countingProcessVar == 1 &
                                            eventtypevar == a]
                ## get weight
                if(is.null(weight)){
                  weightInertiaPast <- rep(1, length(timesIntertiaPast))
                }else{
                  weightInertiaPast <- weight[sender == target[i] & target == sender[i] & 
                                                time < time[i] & countingProcessVar == 1 &
                                                eventtypevar == a]
                }
                
                ## run cpp-loop for all times
                result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                  weightInertiaPast) 
              }else{
                result <- 0
              }
              ## rbind the variable:
              result 
            }# closes foreach-loop
            
            ## transform result variable
            result <- as.numeric(as.character(res))
            
          }else{ # run in standard form
            
            if(isTRUE(showprogressbar)){
              pb <- txtProgressBar(min = 1, max = length(sender), style = 3)
            }
            for(i in 1:length(sender)){
              if(isTRUE(showprogressbar)){
                setTxtProgressBar(pb, i)
              }
              
              ##
              if( eventtypevar[i] == a){
                
                
                ## create vector of times, sender--target tie was made before
                timesIntertiaPast <- time[sender == target[i] & target == sender[i] & 
                                            time < time[i] & countingProcessVar == 1 &
                                            eventtypevar == a]
                ## get weight
                if(is.null(weight)){
                  weightInertiaPast <- rep(1, length(timesIntertiaPast))
                }else{
                  weightInertiaPast <- weight[sender == target[i] & target == sender[i] & 
                                                time < time[i] & countingProcessVar == 1 &
                                                eventtypevar == a]
                }
                
                ## run cpp-loop for all times
                result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                     weightInertiaPast) 
                
              }else{
                result[i] <- 0
              }
              
            } #closes i-loop
          } # closes if no parallel
          
          ##
          data.short <- cbind(data.short, result)
          
          ## change the name of the data.short's last column
          names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventattributevar))
                                                         gsub(" ", "", a, fixed = TRUE), #represents current event-actor type
                                                         gsub(" ", "", a, fixed = TRUE), #represents actor type of past actions
                                                         sep = ".")	
          
        } # closes a-loop
        
        ## return data
        ## return data frames
        if ( returnData == TRUE ) {
          data <- cbind(data, data.short)
          return(data)
        }else{
          return(data.short)
        }
        ################ done (5)
        
      }else{
        ################ (6) valuemix/values provided, with filter 
        
        for (a in uniqueEventTypeValues){ #current event type
          for (b in uniqueEventTypeValues){ #past event type
            if ( a != b ){
              
              # run in parallel?
              if(isTRUE(inParallel)){
                
                ##
                doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
                res <- foreach::foreach(i=1:length(sender), .combine=rbind)%dopar%{
                  
                  ## 
                  if (eventtypevar[i] == a){
                    
                    ## create vector of times, sender--target tie was made before
                    timesIntertiaPast <- time[sender == target[i] & target == sender[i] & 
                                                time < time[i] & countingProcessVar == 1 & 
                                                eventtypevar == b  & 
                                                eventfiltervar == eventfiltervalue]
                    ## get weight
                    if(is.null(weight)){
                      weightInertiaPast <- rep(1, length(timesIntertiaPast))
                    }else{
                      weightInertiaPast <- weight[sender == target[i] & target == sender[i] & 
                                                    time < time[i] & countingProcessVar == 1 &
                                                    eventtypevar == b & 
                                                    eventfiltervar == eventfiltervalue]
                    }
                    
                    ## run cpp-loop for all times
                    result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                      weightInertiaPast) 
                  }else{
                    result <- 0
                  }
                  ## rbind the variable:
                  result 
                }# closes foreach-loop
                
                ## transform result variable
                result <- as.numeric(as.character(res))
                
              }else{ # run in standard form
                
                if(isTRUE(showprogressbar)){
                  pb <- txtProgressBar(min = 1, max = length(sender), style = 3)
                }
                for(i in 1:length(sender)){
                  if(isTRUE(showprogressbar)){
                    setTxtProgressBar(pb, i)
                  }
                  
                  ##
                  if(eventtypevar[i] == a){
                    
                    
                    ## create vector of times, sender--target tie was made before
                    timesIntertiaPast <- time[sender == target[i] & target == sender[i] & 
                                                time < time[i] & countingProcessVar == 1 &
                                                eventtypevar == b & 
                                                eventfiltervar == eventfiltervalue]
                    ## get weight
                    if(is.null(weight)){
                      weightInertiaPast <- rep(1, length(timesIntertiaPast))
                    }else{
                      weightInertiaPast <- weight[sender == target[i] & target == sender[i] & 
                                                    time < time[i] & countingProcessVar == 1 &
                                                    eventtypevar == b & 
                                                    eventfiltervar == eventfiltervalue]
                    }
                    
                    ## run cpp-loop for all times
                    result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                         weightInertiaPast) 
                  }else{ #if eventtypevar != a
                    result[i] <- 0
                  }
                  
                } #closes i-loop
              } # closes if no parallel
              
              ##
              data.short <- cbind(data.short, result)
              
              ## change the name of the data.short's last column
              names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventattributevar))
                                                             gsub(" ", "", a, fixed = TRUE), #represents current event-actor type
                                                             gsub(" ", "", b, fixed = TRUE), #represents actor type of past actions
                                                             sep = ".")									
            }
          } # closes b-loop
          
          ## calculate effect where both events are of the same type (both a)
          # run in parallel?
          if(isTRUE(inParallel)){
            
            ##
            doParallel::registerDoParallel(cluster) #necessary for the foreach-loop to work, cl <- makeCluster(12, type="FORK")
            res <- foreach::foreach(i=1:length(sender), .combine=rbind)%dopar%{
              
              ## 
              if (eventtypevar[i] == a){
                
                ## create vector of times, sender--target tie was made before
                timesIntertiaPast <- time[sender == target[i] & target == sender[i] & 
                                            time < time[i] & countingProcessVar == 1 &
                                            eventtypevar == a & 
                                            eventfiltervar == eventfiltervalue]
                ## get weight
                if(is.null(weight)){
                  weightInertiaPast <- rep(1, length(timesIntertiaPast))
                }else{
                  weightInertiaPast <- weight[sender == target[i] & target == sender[i] & 
                                                time < time[i] & countingProcessVar == 1 &
                                                eventtypevar == a & 
                                                eventfiltervar == eventfiltervalue]
                }
                
                ## run cpp-loop for all times
                result <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                  weightInertiaPast) 
              }else{
                result <- 0
              }
              ## rbind the variable:
              result 
            }# closes foreach-loop
            
            ## transform result variable
            result <- as.numeric(as.character(res))
            
          }else{ # run in standard form
            
            if(isTRUE(showprogressbar)){
              pb <- txtProgressBar(min = 1, max = length(sender), style = 3)
            }
            for(i in 1:length(sender)){
              if(isTRUE(showprogressbar)){
                setTxtProgressBar(pb, i)
              }
              
              ##
              if( eventtypevar[i] == a){
                
                
                ## create vector of times, sender--target tie was made before
                timesIntertiaPast <- time[sender == target[i] & target == sender[i] & 
                                            time < time[i] & countingProcessVar == 1 &
                                            eventtypevar == a & 
                                            eventfiltervar == eventfiltervalue]
                ## get weight
                if(is.null(weight)){
                  weightInertiaPast <- rep(1, length(timesIntertiaPast))
                }else{
                  weightInertiaPast <- weight[sender == target[i] & target == sender[i] & 
                                                time < time[i] & countingProcessVar == 1 &
                                                eventtypevar == a & 
                                                eventfiltervar == eventfiltervalue]
                }
                
                ## run cpp-loop for all times
                result[i] <- weightTimesSummationCpp(timesIntertiaPast, xlog, time[i], 
                                                     weightInertiaPast) 
                
              }else{
                result[i] <- 0
              }
              
            } #closes i-loop
          } # closes if no parallel
          
          ##
          data.short <- cbind(data.short, result)
          
          ## change the name of the data.short's last column
          names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventattributevar))
                                                         gsub(" ", "", a, fixed = TRUE), #represents current event-actor type
                                                         gsub(" ", "", a, fixed = TRUE), #represents actor type of past actions
                                                         sep = ".")	
          
        } # closes a-loop
        
        ## return data
        ## return data frames
        if ( returnData == TRUE ) {
          data <- cbind(data, data.short)
          return(data)
        }else{
          return(data.short)
        }
        
        ################ done (6)
      } 
      
    } # closes eventtypevalues provided/valuemix 
  } #closes if-eventtypevar != NULL
} #closes inertiaStat-function()


################################################################################
################################################################################

triadStat <- function(data, time, sender, target, halflife, weight = NULL,
                      eventtypevar = NULL, eventtypevalues = NULL, 
                      eventfiltervar = NULL, eventfilterAI = NULL,
                      eventfilterBI = NULL, eventfilterAB = NULL,
                      eventvar = NULL,
                      variablename = "triad", returnData = FALSE, 
                      showprogressbar = FALSE, 
                      inParallel = FALSE, cluster = NULL){
  
  ####### check inputs
  ## check if sender input is available
  if ( is.null(sender) ) {
    stop("No 'sender' argument was provided.")
  }else{
    sender <- as.character(sender)
  }
  
  ## check if target input is available
  if ( is.null(target) ) {
    stop("No 'target' argument was provided.")
  }else{
    target <- as.character(target)
  }  
  
  ## check if event.sequence is well defined (numeric and ever-increasing)
  if ( is.null(time) ) {
    stop("No 'time' argument was provided.")
  }else{
    #test if weight-var is in ascending order
    if ( is.unsorted(time) ) {
      stop("'time' is not sorted. Sort data frame according to the event 
           sequence.")
    }
  }
  
  ## check if vaiables are of same length
  if ( length(sender) != length(target) ){
    stop("'sender' and 'target' are not of same length.")
  }
  if ( length(sender) != length(time) ){
    stop("'sender' and 'time' are not of same length.")
  }
  
  ## check if weight-var is defined (if not -> create it)
  if ( is.null(weight) ) {
    weight <- rep(1, length(time))
  }
  if ( !is.numeric(weight) ) {
    stop("'weight' variable is not numeric.")
  }
  
  ## check if event-type inputs are available and correctly specified
  if ( !is.null(eventtypevar) ) {
    # check if variable is of same length as sender
    if ( length(sender) != length(eventtypevar) ){
      stop("'sender' and 'eventtypevar' are not of same length.")
    }
    # transform variable
    eventtypevar <- as.character(eventtypevar)
    if ( length(unique(eventtypevar)) != 2 ){ 
      stop("'eventtypevar' is not a dummy variable.")
    }
    if ( is.null(eventtypevalues) ){
      stop("No 'eventtypevalues' provided. ")
    }
    if ( length(eventtypevalues) != 2 ){
      stop("'eventtypevalues' not specified correctly. Two values need to be
           provided that will reflect either a 'friend-of-friend', a 'friend-
           of-enemy', a 'enemy-of-friend' or a 'enemy-of-enemy' triad. The two
           values indicate which value in the 'eventtypevar' relates to
           'friend' (or 'enemy') depending on the triad type.")
    }
    if ( length(grep(eventtypevalues[1], eventtypevar)) == 0 ) {
      stop("First value '", eventtypevalues[1], "' is not an element of 'eventtypevar'.") 
    }
    if ( length(grep(eventtypevalues[2], eventtypevar)) == 0 ) {
      stop("Second value '", eventtypevalues[2], "' is not an element of 'eventtypevar'.") 
    }
  }
  
  ## check if event-attribute inputs are available and correctly specified
  if ( is.null(eventfiltervar) == FALSE ) {
    # check length of variable
    if ( length(sender) != length(eventfiltervar) ){
      stop("'sender' and 'eventfiltervar' are not of same length.")
    }
    # transform variable
    eventfiltervar <- as.character(eventfiltervar)
    if ( is.null(eventfilterAB) & is.null(eventfilterAI) & 
         is.null(eventfilterBI) ){
      stop("No 'eventfilter__' provided. Provide a string value by which the 
           events are filtered.", )
    }
    # check if eventattributevalue is part of the variable
    if ( is.null(eventfilterAB) == FALSE){
      if ( length(grep(eventfilterAB, eventfiltervar)) == 0 ) {
        stop("Value '", eventfilterAB, "' is not an element of 'eventfiltervar'.") 
      }
    }
    if ( is.null(eventfilterAI) == FALSE){
      if ( length(grep(eventfilterAI, eventfiltervar)) == 0 ) {
        stop("Value '", eventfilterAI, "' is not an element of 'eventfiltervar'.") 
      }
    }
    if ( is.null(eventfilterBI) == FALSE){
      if ( length(grep(eventfilterBI, eventfiltervar)) == 0 ) {
        stop("Value '", eventfilterBI, "' is not an element of 'eventfiltervar'.") 
      }
    }
  }
  
  ## check event-var
  if(is.null(eventvar) == FALSE){
    if(length(unique(eventvar)) == 2){
      if( ( sort(unique(eventvar))[1] == 0 & sort(unique(eventvar))[2] == 1  ) == FALSE){
        stop('eventvar has to be a dummy variable with values 0 for non-events and
             1 for true events.')
      }
    }else{
      stop('eventvar has to be a dummy variable with values 0 for non-events and
             1 for true events.')
    }
  }
  
  ## cannot take parallel and progress bar
  if(isTRUE(inParallel) & isTRUE(showprogressbar)){
    stop('Cannot spit out progress of the function whilst running the 
         loop in parallel. Turn showprogressbar to FALSE.')
  }
  
  ## cannot have parallel without cluster
  if(isTRUE(inParallel) & is.null(cluster)){
    stop('By choosing to run the loop in parallel, you need to define a 
         cluster. For instance: makeCluster(12, type="FORK"). Alternatively, 
         hand over the number of nodes you would like to run the function on.')
  }
  
  ## check if variablename makes sense (no " " etc.)
  variablename <- gsub(" ", "", variablename, fixed = TRUE)
  
  ## create simple data set to be returned for degree calcuations with more than 1 output-variable
  ##TODO: should there be an event-id-variable?? => that would be useful here
  data.short <- data.frame(time)
  
  ## 
  result <- rep(NA, length(sender))
  
  ## calculate part of decay function
  xlog <- log(2)/halflife 
  
  ## calculate event-data sets
  if(is.null(eventvar)){
    senderLoop <- sender
    targetLoop <- target
    weightLoop <- weight
    timeLoop <- time
    # eventtype
    if(is.null(eventtypevar)){
      eventtypevarLoop <- rep("1", length(sender))
      typeALoop <- "1"
      typeBLoop <- "1"
    }else{ #eventtype given
      eventtypevarLoop <- eventtypevar
      typeALoop <- eventtypevalues[1]
      typeBLoop <- eventtypevalues[2]
    }
    # eventfilter: AB
    if(is.null(eventfilterAB)){
      eventfiltervarABLoop <- rep("1", length(sender))
      eventfilterABLoop <- "1"
    }else{ #eventfilterAB given:
      eventfiltervarABLoop <- eventfiltervar
      eventfilterABLoop <- eventfilterAB
    }
    # eventfilter: AI
    if(is.null(eventfilterAI)){
      eventfiltervarAILoop <- rep("1", length(sender))
      eventfilterAILoop <- "1"
    }else{ #eventfilter given:
      eventfiltervarAILoop <- eventfiltervar
      eventfilterAILoop <- eventfilterAI
    }
    # eventfilter BI
    if(is.null(eventfilterBI)){
      eventfiltervarBILoop <- rep("1", length(sender))
      eventfilterBILoop <- "1"
    }else{ #eventfilter given:
      eventfiltervarBILoop <- eventfiltervar
      eventfilterBILoop <- eventfilterBI
    }
  }else{ # if counting process data is given:
    senderLoop <- sender[eventvar == 1]
    targetLoop <- target[eventvar == 1]
    weightLoop <- weight[eventvar == 1]
    timeLoop <- time[eventvar == 1]
    if(is.null(eventtypevar)){
      eventtypevarLoop <- rep("1", length(senderLoop))
      typeALoop <- "1"
      typeBLoop <- "1"
    }else{ #eventtype given
      eventtypevarLoop <- eventtypevar[eventvar == 1]
      typeALoop <- eventtypevalues[1]
      typeBLoop <- eventtypevalues[2]
    }
    # eventfilter: AB
    if(is.null(eventfilterAB)){
      eventfiltervarABLoop <- rep("1", length(sender)) ## eventfiltervarABLoop keeps its original length => does not get reduced!
      eventfilterABLoop <- "1"
    }else{ #eventfilterAB given:
      eventfiltervarABLoop <- eventfiltervar # not subsetted, bc it's used in the outter i-loop, not the cpp-loop with the reduced data set
      eventfilterABLoop <- eventfilterAB
    }
    # eventfilter: AI
    if(is.null(eventfilterAI)){
      eventfiltervarAILoop <- rep("1", length(senderLoop))
      eventfilterAILoop <- "1"
    }else{ #eventfilter given:
      eventfiltervarAILoop <- eventfiltervar[eventvar == 1]
      eventfilterAILoop <- eventfilterAI
    }
    # eventfilter BI
    if(is.null(eventfilterBI)){
      eventfiltervarBILoop <- rep("1", length(senderLoop))
      eventfilterBILoop <- "1"
    }else{ #eventfilter given:
      eventfiltervarBILoop <- eventfiltervar[eventvar == 1]
      eventfilterBILoop <- eventfilterBI
    }
  }
  
  ####### calculate stat
  
  ##
  if(isTRUE(inParallel)){
    
    ##
    doParallel::registerDoParallel(cluster)
    
    ##
    res <- foreach::foreach(i=1:length(sender), .combine=rbind)%dopar%{
      
      ## check if eventfilterAB is TRUE
      if(eventfiltervarABLoop[i] == eventfilterABLoop){
        
        ## get list of parters sender and target have interacted with
        # With whom (other than target B) has sender interacted in the past? => AI
        xa <- targetLoop[senderLoop == sender[i] & 
                           targetLoop != target[i] & 
                           timeLoop < time[i] & 
                           eventtypevarLoop == typeALoop & 
                           eventfiltervarAILoop == eventfilterAILoop]
        xb <- senderLoop[targetLoop == sender[i] & 
                           senderLoop != target[i] & 
                           timeLoop < time[i] & 
                           eventtypevarLoop == typeALoop & 
                           eventfiltervarAILoop == eventfilterAILoop]
        x <- unique(c(xa, xb))
        # With whom has target interacted in past? (other than current sender)
        ya <- targetLoop[senderLoop == target[i] & 
                           targetLoop != sender[i] & 
                           timeLoop < time[i] & 
                           eventtypevarLoop == typeBLoop & 
                           eventfiltervarBILoop == eventfilterBILoop]
        yb <- senderLoop[targetLoop == target[i] & 
                           senderLoop != sender[i] & 
                           timeLoop < time[i] & 
                           eventtypevarLoop == typeBLoop & 
                           eventfiltervarBILoop == eventfilterBILoop]
        y <- unique(c(ya, yb))
        
        ## Do x and y overlap?
        v <- intersect(x, y)
        
        ## if they overlap, run the cpp-loop, otherwise return result = 0
        if(length(v) == 0){
          result <- 0
        }else{ #if the sender and target have one in common => run cpp-loop
          
          # find i in reduced data set
          if(is.null(eventvar)){
            iLoop <- i-1 # bc cpp-loops start at 0 not 1
          }else{
            iLoop <- length(timeLoop[timeLoop < time[i]]) #+ 1 - 1 # + 1 bc in the loop it's <; however cpp starts at 0, so -1
          }
          
          # run the cpp-loop
          result <- triadCpp(v, senderLoop, targetLoop, timeLoop, weightLoop,
                             eventtypevarLoop, typeALoop, typeBLoop, 
                             eventfiltervarAILoop, eventfilterAILoop, 
                             eventfiltervarBILoop, eventfilterBILoop, 
                             xlog, iLoop, 
                             sender[i], target[i], time[i])
        }
      }else{ # if eventfilterAB[i] != eventfilterAB
        result <- 0
      }
      
      ## rbind the variable:
      result
    } #closes i-loop
    
    ## transform result variable
    result <- as.numeric(as.character(res))
  }else{ # run loop without parallelization
    
    ## 
    if(isTRUE(showprogressbar)){
      pb <- txtProgressBar(min = 1, max = length(sender), style = 3)
    }
    for(i in 1:length(sender)){
      if(isTRUE(showprogressbar)){
        setTxtProgressBar(pb, i)
      }
      
      ## check if eventfilterAB is TRUE
      if(eventfiltervarABLoop[i] == eventfilterABLoop){
        
        ## get list of parters sender and target have interacted with
        # With whom (other than target B) has sender interacted in the past? => AI
        xa <- targetLoop[senderLoop == sender[i] & 
                           targetLoop != target[i] & 
                           timeLoop < time[i] & 
                           eventtypevarLoop == typeALoop & 
                           eventfiltervarAILoop == eventfilterAILoop]
        xb <- senderLoop[targetLoop == sender[i] & 
                           senderLoop != target[i] & 
                           timeLoop < time[i] & 
                           eventtypevarLoop == typeALoop & 
                           eventfiltervarAILoop == eventfilterAILoop]
        x <- unique(c(xa, xb))
        # With whom has target interacted in past? (other than current sender)
        ya <- targetLoop[senderLoop == target[i] & 
                           targetLoop != sender[i] & 
                           timeLoop < time[i] & 
                           eventtypevarLoop == typeBLoop & 
                           eventfiltervarBILoop == eventfilterBILoop]
        yb <- senderLoop[targetLoop == target[i] & 
                           senderLoop != sender[i] & 
                           timeLoop < time[i] & 
                           eventtypevarLoop == typeBLoop & 
                           eventfiltervarBILoop == eventfilterBILoop]
        y <- unique(c(ya, yb))
        
        ## Do x and y overlap?
        v <- intersect(x, y)
        
        ## if they overlap, run the cpp-loop, otherwise return result = 0
        if(length(v) == 0){
          result[i] <- 0
        }else{ #if the sender and target have one in common => run cpp-loop
          
          # find i in reduced data set
          if(is.null(eventvar)){
            iLoop <- i-1 # bc cpp-loops start at 0 not 1
          }else{
            iLoop <- length(timeLoop[timeLoop < time[i]]) #+ 1 - 1 # + 1 bc in the loop it's <; however cpp starts at 0, so -1
          }
          # run the cpp-loop
          result[i] <- triadCpp(v, senderLoop, targetLoop, timeLoop, weightLoop,
                                eventtypevarLoop, typeALoop, typeBLoop, 
                                eventfiltervarAILoop, eventfilterAILoop, 
                                eventfiltervarBILoop, eventfilterBILoop, 
                                xlog, iLoop, 
                                sender[i], target[i], time[i])
        }
      }else{ # if eventfilterAB[i] != eventfilterAB
        result[i] <- 0
      }
    }# closes i-loop
  } # closes if no parallel
  
  ## return results
  ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
  if ( returnData == TRUE ) {
    ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
    data <- cbind(data, result)
    names(data)[length(data)] <- variablename
    ## return the data frame with the variable bound to it
    return(data)
  }else{ 
    ## only return the 1 inertia variable that was generated
    return(result)
  }
  
}#closing



################################################################################
## create REM Dataset with null events - various data strategy options allowed
################################################################################

createRemDataset <- function(data, sender, target, eventSequence, 
                             eventAttribute = NULL, time = NULL, 
                             start = NULL, startDate = NULL, 
                             end = NULL, endDate = NULL, 
                             timeformat = NULL,
                             atEventTimesOnly = TRUE, untilEventOccurrs = TRUE,
                             includeAllPossibleEvents = FALSE, 
                             possibleEvents = NULL, 
                             returnInputData = FALSE){
  
  ## possibleEvents has to be formatted as follows: 
  ## 1=sender, 2=target, 3=start, 4=end, 5=attribute, 6...=egal
  
  ## check if sender and target inputs are available
  if ( is.null(sender) ) {
    stop("No 'sender' argument was provided.")
  }else{
    sender <- as.character(sender)
  }
  
  if ( is.null(target) ) {
    stop("No 'target' argument was provided.")
  }else{
    target <- as.character(target)
  }
  
  ## check if all variables have the same length.
  if (length(sender) != length(target)){
    stop("'sender' and 'target' are not of the same length.")
  }
  
  ## check if event.sequence is well defined (numeric and ever-increasing)
  if(is.null(eventSequence)){
    stop("No 'eventSequence' argument was provided.")
  }else{
    #test if weight-var is in ascending order
    if(is.unsorted(eventSequence)){
      stop("'", eventSequence, "' is not sorted. Sort data frame according to the event sequence.")
    }
  }
  
  ## time needs to be specified if no end/start!
  if(is.null(time) & is.null(start) & !is.null(startDate)){
    stop('if no start-variable is provied, time-variable has to be specified.')
  }
  if(is.null(time) & is.null(end) & !is.null(endDate)){
    stop('if no end-variable is provied, time-variable has to be specified.')
  }
  
  ## create return data set
  dataNullEvents <- data.frame()
  
  ## 
  if(isTRUE(includeAllPossibleEvents) & is.null(possibleEvents)){
    stop('including all possible events (includeAllPossibleEvents = TRUE) has
         to be accompanied by a data frame specified in possibleEvents. The 
         data frame has to include a sender, target, start and end variable.')
  }
  
  ## create variable with unique event times (for atEventTimesOnly == TRUE)
  if(isTRUE(atEventTimesOnly)){
    allEventTimes <- unique(eventSequence)
  }
  
  ## create start-variable (1day/event before)
  if(includeAllPossibleEvents == FALSE){
    if(is.null(start)){
      ## if events should occurr at any possible time - the start variable has to be calculated beforehand
      if(atEventTimesOnly == FALSE){
        stop('start-variable has to be provided if null events should be created
             at every possible event-timepoint. This is due to the fact that in
             a continuous event sequence, the start variable using startDate is
             calculated to the closest true event, rather than the exact days/etc. .')
      }
      if(is.null(startDate)){
        ## repeated events start 1eventday after the last identical event occurred
        if(is.null(eventAttribute)){
          duplicateSenderTarget <- duplicated(paste0(sender, target))
        }else{
          duplicateSenderTarget <- duplicated(paste0(sender, target, eventAttribute))
        }
        ##
        for(a in 1:length(sender)){
          if(duplicateSenderTarget[a] == FALSE){
            start[a] <- 0
          }else{
            # start begins again from the next possible event
            start[a] <- ifelse(is.null(eventAttribute),
                               suppressWarnings(max(eventSequence[sender == sender[a] & 
                                                                    target == target[a] &
                                                                    eventSequence < eventSequence[a]])), 
                               suppressWarnings(max(eventSequence[sender == sender[a] & 
                                                                    target == target[a] &
                                                                    eventAttribute == eventAttribute[a] & 
                                                                    eventSequence < eventSequence[a]])))
            if(start[a] == -Inf){
              ## special case where two events occurr at the same time - the second
              ## event is duplicated, but there is no event earlier than itself, so it gets start = 0
              start[a] <- 0
            }
            # now give it the next possible value
            start[a] <- which(allEventTimes > start[a])[1]
          }
        }#closes a-loop
      }else{ # start-Date is not NULL
        if(is.null(timeformat)){
          if(class(time) != 'Date'){
            stop('time-Variable is not specified correctly. time-variable has 
                 to be a Date-class variable. Alternatively a timeformat can be provided
                 to transform the time variable into a Date-variable.')
          }
          if(class(startDate) != 'Date'){
            stop('startDate-Variable is not specified correctly. startDate-variable has 
                 to be a Date-class variable. Alternatively a timeformat can be provided
                 to transform the startDate variable into a Date-variable.')
          }
        }else{
          time <- as.Date(time, timeformat)
          startDate <- as.Date(startDate, timeformat)
        }
        ## repeated events start 1eventday after the last identical event occurred
        if(is.null(eventAttribute)){
          duplicateSenderTarget <- duplicated(paste0(sender, target))
        }else{
          duplicateSenderTarget <- duplicated(paste0(sender, target, eventAttribute))
        }
        ## for each variable:
        for(a in 1:length(sender)){
          if(duplicateSenderTarget[a] == FALSE){
            if(startDate[a] < min(time)){
              start[a] <- 0
            }else{
              start[a] <- min(eventSequence[time >= startDate[a]]) 
            }
          }else{ # there are duplicates! 
            ## two things need to be considered: 
            ## first: a start date is given
            ## and has to be matched to the time-variable for the closest
            ## possible event to that data
            ## second: if there are duplicates within the start-time interval, 
            ## the duplicate time +1 needs to be taken
            
            # (1) smallest timepoint that is the same as the start-date or a unit more
            whichStartPointB <- min(eventSequence[time >= startDate[a]]) 
            
            # (2) start begins again from the next possible event
            whichStartPointA <- ifelse(is.null(eventAttribute),
                                       suppressWarnings(max(eventSequence[sender == sender[a] & 
                                                                            target == target[a] &
                                                                            eventSequence < eventSequence[a]])), 
                                       suppressWarnings(max(eventSequence[sender == sender[a] & 
                                                                            target == target[a] &
                                                                            eventAttribute == eventAttribute[a] & 
                                                                            eventSequence < eventSequence[a]])))
            # special case where two events occurr at the same time - the second
            # event is duplicated, but there is no event earlier than itself, so it gets start = 0
            if(whichStartPointA == -Inf){
              whichStartPointA <- 0
            }
            # now give it the next possible value
            whichStartPointA <- which(allEventTimes > whichStartPointA)[1]
            ## assign it the max value of the two possibilities
            whichStartPoint <- c(whichStartPointA, whichStartPointB)
            start[a] <- max(whichStartPoint)
          }
        }#closes a-loop
      }
    }else{ #start is defined:
      ## here too: correct start-variable for duplicate events
      ## repeated events start 1eventday after the last identical event occurred
      if(is.null(eventAttribute)){
        duplicateSenderTarget <- duplicated(paste0(sender, target))
      }else{
        duplicateSenderTarget <- duplicated(paste0(sender, target, eventAttribute))
      }
      
      for(a in 1:length(sender)){
        if(duplicateSenderTarget[a] == FALSE){ # no duplicates
          # if the start value is smaller than the lowest value of the event sequence
          if(start[a] < min(eventSequence)){
            start[a] <- 0
          }else{
            # otherwise take the start value that is provided
            start[a] <- start[a]
          }
        }else{ # if there are duplicates! Attention!
          if(isTRUE(atEventTimesOnly)){
            # start begins from next possible event/
            # (1) - start-value that is given 
            whichStartPointB <- start[a]
            # (2) - start-value from a duplicate event within that timespan
            whichStartPointA <- ifelse(is.null(eventAttribute),
                                       suppressWarnings(max(eventSequence[sender == sender[a] & 
                                                                            target == target[a] &
                                                                            eventSequence < eventSequence[a]])), 
                                       suppressWarnings(max(eventSequence[sender == sender[a] & 
                                                                            target == target[a] &
                                                                            eventAttribute == eventAttribute[a] & 
                                                                            eventSequence < eventSequence[a]])))
            # special case where two events occurr at the same time - the second
            # event is duplicated, but there is no event earlier than itself, so it gets start = 0
            if(whichStartPointA == -Inf){
              whichStartPointA <- 0
            }
            # now give it the next possible value
            whichStartPointA <- which(allEventTimes > whichStartPointA)[1]
            ## Finally: assign it the max value of the two possibilities
            whichStartPoint <- c(whichStartPointA, whichStartPointB)
            start[a] <- max(whichStartPoint)
            
          }else{ #not at event times
            # (1) take the start value provided
            whichStartPointB <- start[a]
            # (2) where is the last duplicate + 1
            whichStartPointA <- ifelse(is.null(eventAttribute),
                                       suppressWarnings(max(eventSequence[sender == sender[a] & 
                                                                            target == target[a] &
                                                                            eventSequence < eventSequence[a]])), 
                                       suppressWarnings(max(eventSequence[sender == sender[a] & 
                                                                            target == target[a] &
                                                                            eventAttribute == eventAttribute[a] & 
                                                                            eventSequence < eventSequence[a]])))
            # special case where two events occurr at the same time - the second
            # event is duplicated, but there is no event earlier than itself, so it gets start = 0
            if(whichStartPointA == -Inf){
              whichStartPointA <- 0
            }
            whichStartPointA <- whichStartPointA+1 # 1 eventday after the previous event took place
            
            ## Finally: 
            whichStartPoint <- c(whichStartPointA, whichStartPointB)
            start[a] <- max(whichStartPoint)
          }
        }
        
      }#closes a-loop
    }
  }
  
  ## create end-variable
  if(includeAllPossibleEvents == FALSE){
    if(is.null(end)){
      if(is.null(endDate)){
        ## events that are repeated in future end after they occurr
        if(is.null(eventAttribute)){
          duplicateSenderTargetFromLast <- duplicated(paste0(sender, target), fromLast = TRUE)
        }else{
          duplicateSenderTargetFromLast <- duplicated(paste0(sender, target, eventAttribute), fromLast = TRUE)
        }
        ##
        for(a in 1:length(sender)){
          if(isTRUE(duplicateSenderTargetFromLast[a])){
            end[a] <- eventSequence[a]
          }else{
            if(isTRUE(untilEventOccurrs)){
              end[a] <- eventSequence[a]
            }else{
              end[a] <- max(eventSequence)
            }
          }
        }#closes a-loop
      }else{ # if endDate is specified (but 'end' is not)
        # if no time-format is set
        if(is.null(timeformat)){
          if(class(time) != 'Date'){
            stop('time-Variable is not specified correctly. time-variable has 
                 to be a Date-class variable. Alternatively a timeformat can be provided
                 to transform the time variable into a Date-variable.')
          }
          if(class(endDate) != 'Date'){
            stop('endDate-Variable is not specified correctly. endDate-variable has 
                 to be a Date-class variable. Alternatively a timeformat can be provided
                 to transform the endDate variable into a Date-variable.')
          }
        }else{ # if timeformat is set, transform date variables
          time <- as.Date(time, timeformat)
          endDate <- as.Date(endDate, timeformat)
        }
        # find duplicates
        if(is.null(eventAttribute)){
          duplicateSenderTargetFromLast <- duplicated(paste0(sender, target), fromLast = TRUE)
        }else{
          duplicateSenderTargetFromLast <- duplicated(paste0(sender, target, eventAttribute), fromLast = TRUE)
        }
        
        # loop over all events
        for(a in 1:length(sender)){
          # duplicates
          if(isTRUE(duplicateSenderTargetFromLast[a])){
            end[a] <- eventSequence[a]
          }else{
            if(isTRUE(untilEventOccurrs)){
              end[a] <- eventSequence[a]
            }else{
              # if endDate is later than max(time)
              if(endDate[a] > max(time)){
                end[a] <- max(eventSequence)
              }else{
                end[a] <- max(eventSequence[endDate[a] < time ])
              }
            }
          }
        }#closes a-loop
      }
    }else{
      if(isTRUE(untilEventOccurrs)){
        end <- eventSequence
      }else{
        end <- end
      }
    }
  }
  
  ## create an event-ID variable
  eventID <- paste0(rep('eventID', length(sender)), 1:length(sender))
  
  ## create a 
  if(is.null(eventAttribute)){
    eventAttr <- rep("1", length(sender))
    eventAttr <- as.character(eventAttr)
  }else{
    eventAttr <- as.character(eventAttribute)
  }
  
  ############## 1. and 2.
  # 1. fill risk set until 'end'-variable
  # 2. fill risk set until end-variable, use any time point (regardless if an 
  # event occurred or not) risk set contains all events that eventually took place
  if(includeAllPossibleEvents == FALSE){
    if(isTRUE(atEventTimesOnly)){
      # 1. includeAllPossibleEvents = FALSE and atEventTimesOnly = TRUE
      # above: allEventTimes <- unique(eventSequence)
      dataNullEvents <- createNullEvents(eventID, 
                                         sender, target, eventAttr,
                                         eventSequence, 
                                         start, end, allEventTimes)
    }else{
      # 2. includeAllPossibleEvents = FALSE and atEventTimesOnly = FALSE
      allEventTimes <- min(start):max(end)
      dataNullEvents <- createNullEvents(eventID, 
                                         sender, target, eventAttr, 
                                         eventSequence, 
                                         start, end, allEventTimes)
    }
    ## remove useless eventattribute-Variable if set to null
    if(is.null(eventAttribute)){
      dataNullEvents <- dataNullEvents[,1:7]
    }
  }
  ############## 
  
  ##TODO: Nr 3 and 4 are still slow - recode them in cpp
  ############## 3. 
  if(atEventTimesOnly == TRUE & includeAllPossibleEvents == TRUE){
    # for each event in the data set
    for(i in 1:length(allEventTimes)){
      # include all null-events that range within the current time
      temp <- possibleEvents[allEventTimes[i] >= possibleEvents[,3] &
                               allEventTimes[i] <= possibleEvents[,4],]
      # specify whether or not an event occurred
      temp$eventTime <- allEventTimes[i]
      temp$eventdummy <- 0
      if(i != 1){
        dataNullEvents <- rbind(dataNullEvents, temp)
      }else{
        dataNullEvents <- temp
      }
    }#closes i-loop
    
    for(j in 1:nrow(data)){
      if(is.null(eventAttribute)){
        dataNullEvents$eventdummy[dataNullEvents[,1] == sender[j] & dataNullEvents[,2] == target[j] & 
                                    dataNullEvents$eventTime == eventSequence[j]] <- 1
      }else{
        dataNullEvents$eventdummy[dataNullEvents[,1] == sender[j] & dataNullEvents[,2] == target[j] & 
                                    dataNullEvents$eventTime == eventSequence[j] & dataNullEvents[,5] == eventAttribute[j] ] <- 1
      }
    }#closes j-loop
  }
  ############## 
  
  ############## 4. 
  if(atEventTimesOnly == FALSE & includeAllPossibleEvents == TRUE){
    if(isTRUE(untilEventOccurrs)){
      for(i in 1:max(eventSequence)){
        # include all null-events that range within the current time
        temp <- possibleEvents[i >= possibleEvents[,3] &
                                 i <= possibleEvents[,4],]
        # specify whether or not an event occurred
        temp$eventTime <- i
        temp$eventdummy <- 0
        if(i != 1){
          dataNullEvents <- rbind(dataNullEvents, temp)
        }else{
          dataNullEvents <- temp
        }
      }#closes i-loop
      for(j in 1:nrow(data)){
        if(is.null(eventAttribute)){
          dataNullEvents$eventdummy[dataNullEvents[,1] == sender[j] & dataNullEvents[,2] == target[j] & 
                                      dataNullEvents$eventTime == eventSequence[j]] <- 1
        }else{
          dataNullEvents$eventdummy[dataNullEvents[,1] == sender[j] & dataNullEvents[,2] == target[j] & 
                                      dataNullEvents$eventTime == eventSequence[j] & dataNullEvents[,5] == eventAttribute[j] ] <- 1
        }
      }#closes j-loop
    }else{
      # for each event in the data set
      maxSeq <- max(c(eventSequence, possibleEvents[,4]))
      for(i in 1:max(maxSeq)){
        # include all null-events that range within the current time
        temp <- possibleEvents[i >= possibleEvents[,3] &
                                 i <= possibleEvents[,4],]
        # specify whether or not an event occurred
        temp$eventTime <- i
        temp$eventdummy <- 0
        if(i != 1){
          dataNullEvents <- rbind(dataNullEvents, temp)
        }else{
          dataNullEvents <- temp
        }
      }#closes i-loop
      for(j in 1:nrow(data)){
        if(is.null(eventAttribute)){
          dataNullEvents$eventdummy[dataNullEvents[,1] == sender[j] & dataNullEvents[,2] == target[j] & 
                                      dataNullEvents$eventTime == eventSequence[j]] <- 1
        }else{
          dataNullEvents$eventdummy[dataNullEvents[,1] == sender[j] & dataNullEvents[,2] == target[j] & 
                                      dataNullEvents$eventTime == eventSequence[j] & dataNullEvents[,5] == eventAttribute[j] ] <- 1
        }
      }#closes j-loop 
    }
  }
  ############## 
  
  ## DONE:
  if(isTRUE(returnInputData)){
    data <- cbind(eventID, data, start, end)
    return(list(dataNullEvents, data))
  }else{
    return(dataNullEvents)
  }
}

################################################################################
## Calculate time to next event or time since date
################################################################################

timeToEvent <- function(time, type = 'time-to-next-event', 
                        timeEventPossible = NULL){
  
  ############ 0. data checks
  ## is time sorted?
  if ( is.null(time) ) {
    stop("No 'time' argument was provided.")
  }else{
    #test if weight-var is in ascending order
    if ( is.unsorted(time) ) {
      stop("'", time, "' is not sorted. Sort data according to the event time")
    }
  }
  
  ## time and timeEventPossible specified correctly
  if(is.null(timeEventPossible)){
    if(class(time) == 'integer' | class(time) == 'Date'){
      
    }else{
      stop('time variable not specified correctly. Has to be either 
           integer or Date class.')
    }
  }else{
    if(class(time) == 'integer' & class(timeEventPossible) == 'Date'){
      stop('tiem and timeEventPossible not specified correclty.
             Both variables have to be either integer or Date class objects.')
    }else if(class(time) == 'Date' & class(timeEventPossible) == 'integer' ){
      stop('tiem and timeEventPossible not specified correclty.
             Both variables have to be either integer or Date class objects.')
    }
  }
  
  ##
  timeToEventVar <- NULL
  
  ############ 1. time-to-next-event
  if(type == 'time-to-next-event'){
    for(i in 1:length(time)){
      if(time[i] == min(time)){
        timeToEventVar[i] <- 1
      }else{
        timeToEventVar[i] <- as.numeric(time[i] - max(time[time < time[i]]))
      }
    }
  }
  ############# 2. time since 
  if(type == 'time-since-date'){
    ## 
    timeToEventVar <- as.numeric(time - timeEventPossible)
  }
  
  ## 
  return(timeToEventVar)
}

