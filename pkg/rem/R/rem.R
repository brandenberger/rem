## File: Functions to generate variables that capture effects for REMs
## Author: Laurence Brandenberger
## Date: Created: 17. August, last updated: 17.August (00:27)

####################################################################
####################################################################
####################################################################

# TODO 1: add decay-function-option in the values
# TODO 2: add progress bar to functions (and set showprogressbar = TRUE)
# TODO 3: add real function names to all the stop()-outputs - if possible
# TODO 4: OPENMP - implement paralells in cpp-Functions
# TODO 5: tidy up functions - within 80char/line

####################################################################
##  Inertia
####################################################################

get.inertia.stat <- function(data, time, sender, target, halflife, weight = NULL, eventtypevar = NULL, eventtypevalue = "valuematch", 
                             eventattributevar = NULL, eventattributevalue = "valuematch", variablename = "inertia", returnData = TRUE, 
                             showprogressbar = FALSE){
  
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
  
  ## check if event.sequence is well defined (numeric and ever-increasing)
  if ( is.null(time) ) {
    stop("No 'time' argument was provided.")
  }else{
    #test if weight-var is in ascending order
    if ( is.unsorted(time) ) {
      stop("'", time, "' is not sorted. Sort data frame according to the event sequence.")
    }
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
  if ( is.null(eventattributevar) == FALSE ) {
    eventattributevar <- as.character(eventattributevar)
    if ( is.null(eventattributevalue) ){
      stop("No 'eventattributevalue' provided. Use default 'valuematch', or 'valuemix' or string value(s) to determine by which values the events should be filtered.", )
    }
    # check if eventattributevalue is part of the variable
    if ( length(eventattributevalue) > 1 ){
      for ( i in 1:length(eventattributevalue) ){
        if ( length(grep(eventattributevalue[i], eventattributevar)) == 0 ) {
          stop("Value '", eventattributevalue[i], "' is not an element of '", as.name(eventattributevar), "'.")  ##deparse(substitute(eventattributevar))
        }
      }#closes i-loop  
      if ( length(unique(duplicated(eventattributevalue))) == 2 ) {
        stop("Duplicate values in 'eventattributevalue'.") 
      }
    }else if ( eventattributevalue != "valuematch" &  eventattributevalue != "valuemix") {
      if ( length(grep(eventattributevalue, eventattributevar)) == 0 ) {
        ##TODO: #deparse(substitute(eventtypevar))
        stop("Value '", eventattributevalue, "' is not an element of '", deparse(substitute(eventattributevar)) , "'.") ##deparse(substitute(eventtypevar))
      }
    }
  }
  
  ## check if variablename makes sense (no " " etc.)
  variablename <- gsub(" ", "", variablename, fixed = TRUE)
  
  ## create simple data set to be returned for inertia calcuations with more than 1 output-variable
  ##TODO: should there be an event-id-variable?? => that would be useful here
  data.short <- data.frame(time)
  
  ## calculate part of decay function
  xlog <- log(2)/halflife 
  
  ####### calculate stat
  ## create placeholder-variables to be used in the cpp-Function
  placeholder <- rep("1", length(time))
  
  ## calculate the inertia effects for each event
  if ( is.null(eventtypevar) ) {
    if ( is.null(eventattributevar) ) {
      ## (1) start off with simple inertia function: no type, no attribute
      result <- inertiaCpp(time, weight, sender, target, placeholder, "1", "1", placeholder, "1", "1", xlog, "s-t-only")		
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
    }else{
      ## all the functions with event attribute variable but no type variable:
      if ( length(eventattributevalue) == 1 ){
        if ( eventattributevalue == "valuematch" ){
          ## (2) with eventattributevalue set to "valuematch"
          result <- inertiaCpp(time, weight, sender, target, placeholder, "1", "1", eventattributevar, "1", "1", xlog, "s-t-attributematch")		
          ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
          if ( returnData == TRUE ) {
            ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
            ##TODO: also figure out the deparse(substitute(eventattributevar)) problem so that var-name can be used in varname
            data <- cbind(data, result)
            names(data)[length(data)] <- paste(variablename, "attr", "match", sep = ".") #deparse(substitute(eventattributevar))
            
            ## return the data frame with the variable bound to it
            return(data)
          }else{ 
            ## only return the 1 inertia variable that was generated
            return(result)
          }	
        } else if ( eventattributevalue != "valuemix" ) {
          ## (3) with one eventattirbutevalue selected (and used as filter)
          result <- inertiaCpp(time, weight, sender, target, placeholder, "1", "1", eventattributevar, eventattributevalue, "1", xlog, "s-t-attributefilter")		
          ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
          if ( returnData == TRUE ) {
            ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
            ##TODO: also figure out the deparse(substitute(eventattributevar)) problem so that var-name can be used in varname
            data <- cbind(data, result)
            names(data)[length(data)] <- paste(variablename, "attr", gsub(" ", "", eventattributevalue, fixed = TRUE), sep = ".") #deparse(substitute(eventattributevar))
            
            ## return the data frame with the variable bound to it
            return(data)
          }else{ 
            ## only return the 1 inertia variable that was generated
            return(result)
          }	
        } else if ( eventattributevalue == "valuemix"){
          for (i in unique(eventattributevar)){
            for (j in unique(eventattributevar)){
              if ( i != j ){
                ## calculate inertia for the two distinct attribute-values
                temp <- inertiaCpp(time, weight, sender, target, placeholder, "1", "1", eventattributevar, i, j, xlog, "s-t-attributemix")		
                data.short <- cbind(data.short, temp)
                names(data.short)[length(data.short)] <- paste(variablename, "attr", #deparse(substitute(eventattributevar))
                                                               gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                               sep = ".")									
              }
            }#closes j-loop
            ## calculate inertia for the attribute-values where i and j are the same (both i used)
            temp <- inertiaCpp(time, weight, sender, target, placeholder, "1", "1", eventattributevar, i, i, xlog, "s-t-attributemix")		
            ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
            data.short <- cbind(data.short, temp)
            
            names(data.short)[length(data.short)] <- paste(variablename, "attr", #deparse(substitute(eventattributevar))
                                                           gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                           gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                           sep = ".")
            
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }	
        }
      } else if ( length(eventattributevalue) > 1 ) {
        ## (4) with specific eventattributevalues selected
        for (i in eventattributevalue){
          for (j in eventattributevalue){
            if ( i != j ){
              ## calculate inertia for the two distinct attribute-values
              temp <- inertiaCpp(time, weight, sender, target, placeholder, "1", "1", eventattributevar, i, j, xlog, "s-t-attributemix")		
              data.short <- cbind(data.short, temp)
              names(data.short)[length(data.short)] <- paste(variablename, "attr", #deparse(substitute(eventattributevar))
                                                             gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                             gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                             sep = ".")
            }
          }#closes j-loop
          ## calculate inertia for the attribute-values where i and j are the same (both i used)
          temp <- inertiaCpp(time, weight, sender, target, placeholder, "1", "1", eventattributevar, i, i, xlog, "s-t-attributemix")		
          ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
          data.short <- cbind(data.short, temp)
          names(data.short)[length(data.short)] <- paste(variablename, "attr", #deparse(substitute(eventattributevar))
                                                         gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                         gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                         sep = ".")
          
        }#closes i-loop
        ## return data frames
        if ( returnData == TRUE ) {
          data <- cbind(data, data.short)
          return(data)
        }else{
          return(data.short)
        }
      } #closes if ( length(eventattributevalue) > 1 ) 
    } #closes if-else command "if ( is.null(eventattributevar) ) {}else{}"
    
  }else if ( is.null(eventattributevar) ) { #closes if-is.null(eventtypevar)-command 
    ## all the functions that include a type variable (and no attribute): 
    if ( length(eventtypevalue) == 1 ){
      ## (5) with eventtypevar set to "valuematch"
      if ( eventtypevalue == "valuematch" ){
        result <- inertiaCpp(time, weight, sender, target, eventtypevar, "1", "1", placeholder, "1", "1", xlog, "s-t-typematch")  	
        ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
        if ( returnData == TRUE ) {
          ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
          ##TODO: also figure out the deparse(substitute(eventattributevar)) problem so that var-name can be used in varname
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "type", "match", sep = ".") #deparse(substitute(eventtypevar))
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 inertia variable that was generated
          return(result)
        }	
      }  else if ( eventtypevalue != "valuemix" ) {
        ## (6) with only 1 eventtypevar selected (used as filter)
        result <- inertiaCpp(time, weight, sender, target, eventtypevar, eventtypevalue, "1", placeholder, "1", "1", xlog, "s-t-typefilter")  	
        ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
        if ( returnData == TRUE ) {
          ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
          ##TODO: also figure out the deparse(substitute(eventattributevar)) problem so that var-name can be used in varname
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "type", gsub(" ", "", eventtypevalue, fixed = TRUE), sep = ".") #deparse(substitute(eventtypevar))
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 inertia variable that was generated
          return(result)
        }	
      } else if ( eventtypevalue == "valuemix"){
        ## (7) with more than one eventtypevalue selected (or all of them = nodemix)
        for (i in unique(eventtypevar)){
          for (j in unique(eventtypevar)){
            if ( i != j ){
              ## calculate inertia for the two distinct attribute-values
              temp <- inertiaCpp(time, weight, sender, target, eventtypevar, i, j, placeholder, "1", "1", xlog, "s-t-typemix")  	
              data.short <- cbind(data.short, temp)
              names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                             gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                             gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                             sep = ".")
            }										
          }#closes j-loop
          ## calculate inertia for the attribute-values where i and j are the same (both i used)
          temp <- inertiaCpp(time, weight, sender, target, eventtypevar, i, i, placeholder, "1", "1", xlog, "s-t-typemix")		
          ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
          data.short <- cbind(data.short, temp)
          names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                         gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                         gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                         sep = ".")
        }#closes i-loop
        ## return data frames
        if ( returnData == TRUE ) {
          data <- cbind(data, data.short)
          return(data)
        }else{
          return(data.short)
        }	
      }
    }else if ( length(eventtypevalue) > 1 ) {
      ## (4) with specific eventattributevalues selected
      for (i in eventtypevalue ){
        for (j in eventtypevalue ){
          if ( i != j ){
            ## calculate inertia for the two distinct attribute-values
            temp <- inertiaCpp(time, weight, sender, target, eventtypevar, i, j , placeholder, "1", "1", xlog, "s-t-attributemix")  	
            data.short <- cbind(data.short, temp)
            names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventattributevar))
                                                           gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                           gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                           sep = ".")
          }
        }#closes j-loop
        ## calculate inertia for the attribute-values where i and j are the same (both i used)
        temp <- inertiaCpp(time, weight, sender, target, eventtypevar, i, i , placeholder, "1", "1", xlog, "s-t-attributemix")  	
        ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
        data.short <- cbind(data.short, temp)
        names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventattributevar))
                                                       gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                       gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                       sep = ".")
      }#closes i-loop
      ## return data frames
      if ( returnData == TRUE ) {
        data <- cbind(data, data.short)
        return(data)
      }else{
        return(data.short)
      }
    } #closes if ( length(eventtypevalue) > 1 ) {}
  } #closes if ( is.null(eventattributevar) ) {}
  
  ## if both eventtypevar and eventattributevar are selected:
  if ( is.null(eventtypevar) == FALSE & is.null(eventattributevar) == FALSE ){
    ## all the functions that include both a type variable as well as an attribute variable
    if ( length(eventtypevalue) == 1 ){
      if ( length(eventattributevalue) == 1){
        ## (8) with typevar = valuematch and attributevar = valuematch
        if ( eventtypevalue == "valuematch" & eventattributevalue == "valuematch") {
          result <- inertiaCpp(time, weight, sender, target, eventtypevar, "1", "1", eventattributevar, "1", "1", xlog, "s-t-typematch-attributematch")		
          ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
          if ( returnData == TRUE ) {
            ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
            ##TODO: also figure out the deparse(substitute(eventattributevar)) problem so that var-name can be used in varname
            data <- cbind(data, result)
            names(data)[length(data)] <- paste(variablename, "type", "match", "attr", "match", sep = ".") #deparse(substitute(eventattributevar))
            ## return the data frame with the variable bound to it
            return(data)
          }else{ 
            ## only return the 1 inertia variable that was generated
            return(result)
          }	
        }
        
        ## (9) with typevar = 1 value selected and attributevar = valuematch
        if ( eventtypevalue != "valuematch" & eventtypevalue != "valuemix" & eventattributevalue == "valuematch"){
          result <- inertiaCpp(time, weight, sender, target, eventtypevar, eventtypevalue, "1", eventattributevar, "1", "1", xlog, "s-t-typefilter-attributematch")		
          ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
          if ( returnData == TRUE ) {
            ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
            ##TODO: also figure out the deparse(substitute(eventattributevar)) problem so that var-name can be used in varname
            data <- cbind(data, result)
            names(data)[length(data)] <- paste(variablename, "type", gsub(" ", "", eventtypevalue, fixed = TRUE), "attr", "match", sep = ".") #deparse(substitute(eventattributevar))
            ## return the data frame with the variable bound to it
            return(data)
          }else{ 
            ## only return the 1 inertia variable that was generated
            return(result)
          }
        }
        
        ## (10) with typevar = valuematch and attributevar = 1 value selected
        if ( eventtypevalue == "valuematch" & eventattributevalue != "valuemix" & eventattributevalue != "valuematch"){
          result <- inertiaCpp(time, weight, sender, target, eventtypevar, "1", "1", eventattributevar, eventattributevalue, "1", xlog, "s-t-typematch-attributefilter")		
          ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
          if ( returnData == TRUE ) {
            ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
            ##TODO: also figure out the deparse(substitute(eventattributevar)) problem so that var-name can be used in varname
            data <- cbind(data, result)
            names(data)[length(data)] <- paste(variablename, "type", "match", "attr", gsub(" ", "", eventattributevalue, fixed = TRUE), sep = ".") #deparse(substitute(eventattributevar))
            ## return the data frame with the variable bound to it
            return(data)
          }else{ 
            ## only return the 1 inertia variable that was generated
            return(result)
          }
        }
        
        ## (11) with typevar = 1 value selected and attributevar = 1 value selected
        if ( eventtypevalue != "valuematch" & eventtypevalue != "valuemix" & eventattributevalue != "valuematch" & eventattributevalue != "valuemix"){
          result <- inertiaCpp(time, weight, sender, target, eventtypevar, eventtypevalue, "1", eventattributevar, eventattributevalue, "1", xlog, "s-t-typefilter-attributefilter")		
          ## if returnData = TRUE => return the entire data frame as well as the 1 additional inertia-variable
          if ( returnData == TRUE ) {
            ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
            ##TODO: also figure out the deparse(substitute(eventattributevar)) problem so that var-name can be used in varname
            data <- cbind(data, result)
            names(data)[length(data)] <- paste(variablename, "type", gsub(" ", "", eventtypevalue, fixed = TRUE), "attr", gsub(" ", "", eventattributevalue, fixed = TRUE), sep = ".") #deparse(substitute(eventattributevar))
            ## return the data frame with the variable bound to it
            return(data)
          }else{ 
            ## only return the 1 inertia variable that was generated
            return(result)
          }
        }
        
        ## (12-1) with typevar = valuematch and attributevar = valuemix
        if ( eventtypevalue == "valuematch" & eventattributevalue == "valuemix"){
          for (i in unique(eventattributevar)){
            for (j in unique(eventattributevar)){
              if ( i != j ){
                ## calculate inertia for the two distinct attribute-values
                temp <- inertiaCpp(time, weight, sender, target, eventtypevar, "1", "1", eventattributevar, i, j, xlog, "s-t-typematch-attributemix")		
                data.short <- cbind(data.short, temp)
                names(data.short)[length(data.short)] <- paste(variablename, "type", "match", "attr", #deparse(substitute(eventattributevar))
                                                               gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                               sep = ".")
              }
            }#closes j-loop
            ## calculate inertia for the attribute-values where i and j are the same (both i used)
            temp <- inertiaCpp(time, weight, sender, target, eventtypevar, "1", "1", eventattributevar, i, i, xlog, "s-t-typematch-attributemix")		
            ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
            data.short <- cbind(data.short, temp)
            names(data.short)[length(data.short)] <- paste(variablename, "type", "match", "attr", #deparse(substitute(eventattributevar))
                                                           gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                           gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                           sep = ".")
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
        ## (13-1) with typevar = 1 value selected and attributevar = valuemix
        if ( eventtypevalue != "valuematch" & eventtypevalue != "valuemix" & eventattributevalue == "valuemix"){
          for (i in unique(eventattributevar)){
            for (j in unique(eventattributevar)){
              if ( i != j ){
                ## calculate inertia for the two distinct attribute-values
                temp <- inertiaCpp(time, weight, sender, target, eventtypevar, eventtypevalue, "1", eventattributevar, i, j, xlog, "s-t-typefilter-attributemix")		
                data.short <- cbind(data.short, temp)
                names(data.short)[length(data.short)] <- paste(variablename, "type", eventtypevalue, "attr", #deparse(substitute(eventattributevar))
                                                               gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                               sep = ".")
              }
            }#closes j-loop
            ## calculate inertia for the attribute-values where i and j are the same (both i used)
            temp <- inertiaCpp(time, weight, sender, target, eventtypevar, eventtypevalue, "1", eventattributevar, i, i, xlog, "s-t-typefilter-attributemix")		
            ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
            data.short <- cbind(data.short, temp)
            names(data.short)[length(data.short)] <- paste(variablename, "type", eventtypevalue, "attr", #deparse(substitute(eventattributevar))
                                                           gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                           gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                           sep = ".")
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
        ## (14-1) with typevar = valuemix and attributevar = valuematch   
        if ( eventtypevalue == "valuemix" & eventattributevalue == "valuematch"){
          for (i in unique(eventtypevar)){
            for (j in unique(eventtypevar)){
              if ( i != j ){
                ## calculate inertia for the two distinct attribute-values
                temp <- inertiaCpp(time, weight, sender, target, eventtypevar, i, j, eventattributevar, "1", "1", xlog, "s-t-typemix-attributematch")  	
                data.short <- cbind(data.short, temp)
                names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                               gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                               "attr", "match",
                                                               sep = ".")
              }
            }#closes j-loop
            ## calculate inertia for the attribute-values where i and j are the same (both i used)
            temp <- inertiaCpp(time, weight, sender, target, eventtypevar, i, i, eventattributevar, "1", "1", xlog, "s-t-typemix-attributematch")		
            ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
            data.short <- cbind(data.short, temp)
            names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                           gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                           gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                           "attr", "match",
                                                           sep = ".")
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
        ## (15-1) with typevar = valuemix and attributevar = 1 value selected
        if ( eventtypevalue == "valuemix" & eventattributevalue != "valuematch" & eventattributevalue != "valuemix"){
          for (i in unique(eventtypevar)){
            for (j in unique(eventtypevar)){
              if ( i != j ){
                ## calculate inertia for the two distinct attribute-values
                temp <- inertiaCpp(time, weight, sender, target, eventtypevar, i, j, eventattributevar, eventattributevalue, "1", xlog, "s-t-typemix-attributefilter")  	
                data.short <- cbind(data.short, temp)
                names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                               gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                               "attr", eventattributevalue,
                                                               sep = ".")
              }
            }#closes j-loop
            ## calculate inertia for the attribute-values where i and j are the same (both i used)
            temp <- inertiaCpp(time, weight, sender, target, eventtypevar, i, i, eventattributevar, eventattributevalue, "1", xlog, "s-t-typemix-attributefilter")		
            ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
            data.short <- cbind(data.short, temp)
            names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                           gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                           gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                           "attr", eventattributevalue,
                                                           sep = ".")
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
        ## (16-1) with typevar = valuemix and attributevar = valuemix
        if ( eventtypevalue == "valuemix" & eventattributevalue == "valuemix"){
          for (i in unique(eventtypevar)){
            for (j in unique(eventtypevar)){
              for (k in unique(eventattributevar)){
                for (l in unique(eventattributevar)){
                  temp <- inertiaCpp(time, weight, sender, target, eventtypevar, i, j, eventattributevar, k, l, xlog, "s-t-typemix-attributemix")  	
                  data.short <- cbind(data.short, temp)
                  names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                                 gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                                 gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                                 "attr", #deparse(substitute(eventtypevar))
                                                                 gsub(" ", "", k, fixed = TRUE), #represents current event-actor type
                                                                 gsub(" ", "", l, fixed = TRUE), #represents actor type of past actions
                                                                 sep = ".")
                  #TODO: do not allow duplicate entries
                }#closes l-loop	
              }#closes k-loop
            }#closes j-loop
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
      }else{ #else length(eventattributevalue) > 2
        
        ## (12-2) with typevar = valuematch and attributevar = valuemix
        if ( eventtypevalue == "valuematch" ){
          for ( i in eventattributevalue ){
            for ( j in eventattributevalue ){
              if ( i != j ){
                ## calculate inertia for the two distinct attribute-values
                temp <- inertiaCpp(time, weight, sender, target, eventtypevar, "1", "1", eventattributevar, i, j, xlog, "s-t-typematch-attributemix")		
                data.short <- cbind(data.short, temp)
                names(data.short)[length(data.short)] <- paste(variablename, "type", "match", "attr", #deparse(substitute(eventattributevar))
                                                               gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                               sep = ".")
              }
            }#closes j-loop
            ## calculate inertia for the attribute-values where i and j are the same (both i used)
            temp <- inertiaCpp(time, weight, sender, target, eventtypevar, "1", "1", eventattributevar, i, i, xlog, "s-t-typematch-attributemix")		
            ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
            data.short <- cbind(data.short, temp)
            names(data.short)[length(data.short)] <- paste(variablename, "type", "match", "attr", #deparse(substitute(eventattributevar))
                                                           gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                           gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                           sep = ".")
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
        ## (13-2) with typevar = 1 value selected and attributevar = valuemix
        if ( eventtypevalue != "valuematch" & eventtypevalue != "valuemix" ){
          for (i in eventattributevalue ){
            for (j in eventattributevalue ){
              if ( i != j ){
                ## calculate inertia for the two distinct attribute-values
                temp <- inertiaCpp(time, weight, sender, target, eventtypevar, eventtypevalue, "1", eventattributevar, i, j, xlog, "s-t-typefilter-attributemix")		
                data.short <- cbind(data.short, temp)
                names(data.short)[length(data.short)] <- paste(variablename, "type", eventtypevalue, "attr", #deparse(substitute(eventattributevar))
                                                               gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                               sep = ".")
              }
            }#closes j-loop
            ## calculate inertia for the attribute-values where i and j are the same (both i used)
            temp <- inertiaCpp(time, weight, sender, target, eventtypevar, eventtypevalue, "1", eventattributevar, i, i, xlog, "s-t-typefilter-attributemix")		
            ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
            data.short <- cbind(data.short, temp)
            names(data.short)[length(data.short)] <- paste(variablename, "type", eventtypevalue, "attr", #deparse(substitute(eventattributevar))
                                                           gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                           gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                           sep = ".")
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
        ## (16-2) with typevar = valuemix and attributevar = valuemix (values)
        if ( eventtypevalue == "valuemix"){
          for ( i in unique(eventtypevar) ){
            for ( j in unique(eventtypevar) ){
              for ( k in eventattributevalue ){
                for ( l in eventattributevalue ){
                  temp <- inertiaCpp(time, weight, sender, target, eventtypevar, i, j, eventattributevar, k, l, xlog, "s-t-typemix-attributemix")  	
                  data.short <- cbind(data.short, temp)
                  names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                                 gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                                 gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                                 "attr", #deparse(substitute(eventtypevar))
                                                                 gsub(" ", "", k, fixed = TRUE), #represents current event-actor type
                                                                 gsub(" ", "", l, fixed = TRUE), #represents actor type of past actions
                                                                 sep = ".")
                  #TODO: do not allow duplicate entries
                }#closes l-loop	
              }#closes k-loop
            }#closes j-loop
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
      } #closes else{length(eventattributevalue) >2 }
    }else{ #else: length(eventtypevalue) > 1:
      if ( length(eventattributevalue) == 1) {
        
        ## (14-2) with typevar = valuemix (values) and attributevar = valuematch   
        if ( eventattributevalue == "valuematch"){
          for ( i in eventtypevalue ){
            for ( j in eventtypevalue ){
              if ( i != j ){
                ## calculate inertia for the two distinct attribute-values
                temp <- inertiaCpp(time, weight, sender, target, eventtypevar, i, j, eventattributevar, "1", "1", xlog, "s-t-typemix-attributematch")  	
                data.short <- cbind(data.short, temp)
                names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                               gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                               "attr", "match",
                                                               sep = ".")
              }
            }#closes j-loop
            ## calculate inertia for the attribute-values where i and j are the same (both i used)
            temp <- inertiaCpp(time, weight, sender, target, eventtypevar, i, i, eventattributevar, "1", "1", xlog, "s-t-typemix-attributematch")		
            ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
            data.short <- cbind(data.short, temp)
            names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                           gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                           gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                           "attr", "match",
                                                           sep = ".")
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
        ## (15-2) with typevar = valuemix (values) and attributevar = 1 value selected
        if ( eventattributevalue != "valuematch" & eventattributevalue != "valuemix"){
          for ( i in eventtypevalue ){
            for ( j in eventtypevalue ){
              if ( i != j ){
                ## calculate inertia for the two distinct attribute-values
                temp <- inertiaCpp(time, weight, sender, target, eventtypevar, i, j, eventattributevar, eventattributevalue, "1", xlog, "s-t-typemix-attributefilter")  	
                data.short <- cbind(data.short, temp)
                names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                               gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                               "attr", eventattributevalue,
                                                               sep = ".")
              }
            }#closes j-loop
            ## calculate inertia for the attribute-values where i and j are the same (both i used)
            temp <- inertiaCpp(time, weight, sender, target, eventtypevar, i, i, eventattributevar, eventattributevalue, "1", xlog, "s-t-typemix-attributefilter")		
            ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
            data.short <- cbind(data.short, temp)
            names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                           gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                           gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                           "attr", eventattributevalue,
                                                           sep = ".")
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
        ## (16-2) with typevar = valuemix (values) and attributevar = valuemix
        if ( eventattributevalue == "valuemix"){
          for (i in eventtypevalue ){
            for (j in eventtypevalue ){
              for (k in unique(eventattributevar)){
                for (l in unique(eventattributevar)){
                  temp <- inertiaCpp(time, weight, sender, target, eventtypevar, i, j, eventattributevar, k, l, xlog, "s-t-typemix-attributemix")  	
                  data.short <- cbind(data.short, temp)
                  names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                                 gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                                 gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                                 "attr", #deparse(substitute(eventtypevar))
                                                                 gsub(" ", "", k, fixed = TRUE), #represents current event-actor type
                                                                 gsub(" ", "", l, fixed = TRUE), #represents actor type of past actions
                                                                 sep = ".")
                  #TODO: do not allow duplicate entries
                }#closes l-loop	
              }#closes k-loop
            }#closes j-loop
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
      }else{ #else: length(eventattributevalue) > 1
        
        ## (16-2) with typevar = valuemix (values) and attributevar = valuemix (values)
        for ( i in eventtypevalue ){
          for ( j in eventtypevalue ){
            for ( k in eventattributevalue ){
              for ( l in eventattributevalue ){
                temp <- inertiaCpp(time, weight, sender, target, eventtypevar, i, j, eventattributevar, k, l, xlog, "s-t-typemix-attributemix")  	
                data.short <- cbind(data.short, temp)
                names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                               gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                               "attr", #deparse(substitute(eventtypevar))
                                                               gsub(" ", "", k, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", l, fixed = TRUE), #represents actor type of past actions
                                                               sep = ".")
                #TODO: do not allow duplicate entries
              }#closes l-loop	
            }#closes k-loop
          }#closes j-loop
        }#closes i-loop
        ## return data frames
        if ( returnData == TRUE ) {
          data <- cbind(data, data.short)
          return(data)
        }else{
          return(data.short)
        }
      }#closes else{length(eventattributevalue) > 1}		
    }#closes if-else "( length(eventtypevalue) == 1 ){}"
  }#closes -if both eventtypevar and eventattributevar are selected
}


####################################################################
##	Degree calculation
####################################################################

get.degree.stat <- function(data, time, degreevar, halflife, weight = NULL, eventtypevar = NULL, eventtypevalue = "valuematch", 
                            eventattributevar = NULL, eventattributevalue = "valuematch", variablename = "degree", returnData = TRUE, 
                            showprogressbar = FALSE){
  
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
  
  ## check if weight-var is defined (if not -> create it)
  if ( is.null(weight) ) {
    weight <- rep(1, length(time))
  }
  if ( !is.numeric(weight) ) {
    stop("'", as.name(weight), "' variable is not numeric.") #TODO: deparse(substitute(eventattributevar)) ?
  }
  
  ## check if event-type inputs are available and correctly specified
  if ( !is.null(eventtypevar) ) {
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
  if ( is.null(eventattributevar) == FALSE ) {
    eventattributevar <- as.character(eventattributevar)
    if ( is.null(eventattributevalue) ){
      stop("No 'eventattributevalue' provided. Use default 'valuematch', or 'valuemix' or string value(s) to determine by which values the events should be filtered.", )
    }
    # check if eventattributevalue is part of the variable
    if ( length(eventattributevalue) > 1 ){
      for ( i in 1:length(eventattributevalue) ){
        if ( length(grep(eventattributevalue[i], eventattributevar)) == 0 ) {
          stop("Value '", eventattributevalue[i], "' is not an element of '", as.name(eventattributevar), "'.")  ##deparse(substitute(eventattributevar))
        }
      }#closes i-loop  
      if ( length(unique(duplicated(eventattributevalue))) == 2 ) {
        stop("Duplicate values in 'eventattributevalue'.") 
      }
    }else if ( eventattributevalue != "valuematch" &  eventattributevalue != "valuemix") {
      if ( length(grep(eventattributevalue, eventattributevar)) == 0 ) {
        ##TODO: #deparse(substitute(eventtypevar))
        stop("Value '", eventattributevalue, "' is not an element of '", deparse(substitute(eventattributevar)) , "'.") ##deparse(substitute(eventtypevar))
      }
    }
  }
  
  ## check if variablename makes sense (no " " etc.)
  variablename <- gsub(" ", "", variablename, fixed = TRUE)
  
  ## create simple data set to be returned for degree calcuations with more than 1 output-variable
  ##TODO: should there be an event-id-variable?? => that would be useful here
  data.short <- data.frame(time)
  
  ## calculate part of decay function
  xlog <- log(2)/halflife 
  
  ####### calculate stat
  ## create placeholder-variables to be used in the cpp-Function
  placeholder <- rep("1", length(time))
  
  ## calculate the degree effects for each event
  if ( is.null(eventtypevar) ) {
    if ( is.null(eventattributevar) ) {
      ## (1) start off with simple degree function: no type, no attribute
      result <- degreeCpp(time, weight, degreevar, placeholder, "1", "1", placeholder, "1", "1", xlog, "d-only")		
      ## if returnData = TRUE => return the entire data frame as well as the 1 additional degree-variable
      if ( returnData == TRUE ) {
        ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }else{
      ## all the functions with event attribute variable but no type variable:
      if ( length(eventattributevalue) == 1 ){
        if ( eventattributevalue == "valuematch" ){
          ## (2) with eventattributevalue set to "valuematch"
          result <- degreeCpp(time, weight, degreevar, placeholder, "1", "1", eventattributevar, "1", "1", xlog, "d-attributematch")		
          ## if returnData = TRUE => return the entire data frame as well as the 1 additional degree-variable
          if ( returnData == TRUE ) {
            ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
            ##TODO: also figure out the deparse(substitute(eventattributevar)) problem so that var-name can be used in varname
            data <- cbind(data, result)
            names(data)[length(data)] <- paste(variablename, "attr", "match", sep = ".") #deparse(substitute(eventattributevar))
            
            ## return the data frame with the variable bound to it
            return(data)
          }else{ 
            ## only return the 1 degree variable that was generated
            return(result)
          }	
        } else if ( eventattributevalue != "valuemix" ) {
          ## (3) with one eventattirbutevalue selected (and used as filter)
          result <- degreeCpp(time, weight, degreevar, placeholder, "1", "1", eventattributevar, eventattributevalue, "1", xlog, "d-attributefilter")		
          ## if returnData = TRUE => return the entire data frame as well as the 1 additional degree-variable
          if ( returnData == TRUE ) {
            ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
            ##TODO: also figure out the deparse(substitute(eventattributevar)) problem so that var-name can be used in varname
            data <- cbind(data, result)
            names(data)[length(data)] <- paste(variablename, "attr", gsub(" ", "", eventattributevalue, fixed = TRUE), sep = ".") #deparse(substitute(eventattributevar))
            
            ## return the data frame with the variable bound to it
            return(data)
          }else{ 
            ## only return the 1 degree variable that was generated
            return(result)
          }	
        } else if ( eventattributevalue == "valuemix"){
          for (i in unique(eventattributevar)){
            for (j in unique(eventattributevar)){
              if ( i != j ){
                ## calculate degree for the two distinct attribute-values
                temp <- degreeCpp(time, weight, degreevar, placeholder, "1", "1", eventattributevar, i, j, xlog, "d-attributemix")		
                data.short <- cbind(data.short, temp)
                names(data.short)[length(data.short)] <- paste(variablename, "attr", #deparse(substitute(eventattributevar))
                                                               gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                               sep = ".")									
              }
            }#closes j-loop
            ## calculate degree for the attribute-values where i and j are the same (both i used)
            temp <- degreeCpp(time, weight, degreevar, placeholder, "1", "1", eventattributevar, i, i, xlog, "d-attributemix")		
            ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter? now it only filter!)
            data.short <- cbind(data.short, temp)
            names(data.short)[length(data.short)] <- paste(variablename, "attr", #deparse(substitute(eventattributevar))
                                                           gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                           gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                           sep = ".")
            
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }	
        }
      } else if ( length(eventattributevalue) > 1 ) {
        ## (4) with specific eventattributevalues selected
        for (i in eventattributevalue){
          for (j in eventattributevalue){
            if ( i != j ){
              ## calculate degree for the two distinct attribute-values
              temp <- degreeCpp(time, weight, degreevar, placeholder, "1", "1", eventattributevar, i, j, xlog, "d-attributemix")		
              data.short <- cbind(data.short, temp)
              names(data.short)[length(data.short)] <- paste(variablename, "attr", #deparse(substitute(eventattributevar))
                                                             gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                             gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                             sep = ".")
            }
          }#closes j-loop
          ## calculate degree for the attribute-values where i and j are the same (both i used)
          temp <- degreeCpp(time, weight, degreevar, placeholder, "1", "1", eventattributevar, i, i, xlog, "d-attributemix")		
          ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
          data.short <- cbind(data.short, temp)
          names(data.short)[length(data.short)] <- paste(variablename, "attr", #deparse(substitute(eventattributevar))
                                                         gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                         gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                         sep = ".")
          
        }#closes i-loop
        ## return data frames
        if ( returnData == TRUE ) {
          data <- cbind(data, data.short)
          return(data)
        }else{
          return(data.short)
        }
      } #closes if ( length(eventattributevalue) > 1 ) 
    } #closes if-else command "if ( is.null(eventattributevar) ) {}else{}"
    
  }else if ( is.null(eventattributevar) ) { #closes if-is.null(eventtypevar)-command 
    ## all the functions that include a type variable (and no attribute): 
    if ( length(eventtypevalue) == 1 ){
      ## (5) with eventtypevar set to "valuematch"
      if ( eventtypevalue == "valuematch" ){
        result <- degreeCpp(time, weight, degreevar, eventtypevar, "1", "1", placeholder, "1", "1", xlog, "d-typematch")  	
        ## if returnData = TRUE => return the entire data frame as well as the 1 additional degree-variable
        if ( returnData == TRUE ) {
          ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
          ##TODO: also figure out the deparse(substitute(eventattributevar)) problem so that var-name can be used in varname
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "type", "match", sep = ".") #deparse(substitute(eventtypevar))
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }	
      }  else if ( eventtypevalue != "valuemix" ) {
        ## (6) with only 1 eventtypevar selected (used as filter)
        result <- degreeCpp(time, weight, degreevar, eventtypevar, eventtypevalue, "1", placeholder, "1", "1", xlog, "d-typefilter")  	
        ## if returnData = TRUE => return the entire data frame as well as the 1 additional degree-variable
        if ( returnData == TRUE ) {
          ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
          ##TODO: also figure out the deparse(substitute(eventattributevar)) problem so that var-name can be used in varname
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "type", gsub(" ", "", eventtypevalue, fixed = TRUE), sep = ".") #deparse(substitute(eventtypevar))
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }	
      } else if ( eventtypevalue == "valuemix"){
        ## (7) with more than one eventtypevalue selected (or all of them = nodemix)
        for (i in unique(eventtypevar)){
          for (j in unique(eventtypevar)){
            if ( i != j ){
              ## calculate degree for the two distinct attribute-values
              temp <- degreeCpp(time, weight, degreevar, eventtypevar, i, j, placeholder, "1", "1", xlog, "d-typemix")  	
              data.short <- cbind(data.short, temp)
              names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                             gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                             gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                             sep = ".")
            }										
          }#closes j-loop
          ## calculate degree for the attribute-values where i and j are the same (both i used)
          temp <- degreeCpp(time, weight, degreevar, eventtypevar, i, i, placeholder, "1", "1", xlog, "d-typemix")		
          ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
          data.short <- cbind(data.short, temp)
          names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                         gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                         gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                         sep = ".")
        }#closes i-loop
        ## return data frames
        if ( returnData == TRUE ) {
          data <- cbind(data, data.short)
          return(data)
        }else{
          return(data.short)
        }	
      }
    }else if ( length(eventtypevalue) > 1 ) {
      ## (4) with specific eventattributevalues selected
      for (i in eventtypevalue ){
        for (j in eventtypevalue ){
          if ( i != j ){
            ## calculate degree for the two distinct attribute-values
            temp <- degreeCpp(time, weight, degreevar, eventtypevar, i, j , placeholder, "1", "1", xlog, "d-attributemix")  	
            data.short <- cbind(data.short, temp)
            names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventattributevar))
                                                           gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                           gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                           sep = ".")
          }
        }#closes j-loop
        ## calculate degree for the attribute-values where i and j are the same (both i used)
        temp <- degreeCpp(time, weight, degreevar, eventtypevar, i, i , placeholder, "1", "1", xlog, "d-attributemix")  	
        ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
        data.short <- cbind(data.short, temp)
        names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventattributevar))
                                                       gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                       gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                       sep = ".")
      }#closes i-loop
      ## return data frames
      if ( returnData == TRUE ) {
        data <- cbind(data, data.short)
        return(data)
      }else{
        return(data.short)
      }
    } #closes if ( length(eventtypevalue) > 1 ) {}
  } #closes if ( is.null(eventattributevar) ) {}
  
  ## if both eventtypevar and eventattributevar are selected:
  if ( is.null(eventtypevar) == FALSE & is.null(eventattributevar) == FALSE ){
    ## all the functions that include both a type variable as well as an attribute variable
    if ( length(eventtypevalue) == 1 ){
      if ( length(eventattributevalue) == 1){
        ## (8) with typevar = valuematch and attributevar = valuematch
        if ( eventtypevalue == "valuematch" & eventattributevalue == "valuematch") {
          result <- degreeCpp(time, weight, degreevar, eventtypevar, "1", "1", eventattributevar, "1", "1", xlog, "d-typematch-attributematch")		
          ## if returnData = TRUE => return the entire data frame as well as the 1 additional degree-variable
          if ( returnData == TRUE ) {
            ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
            ##TODO: also figure out the deparse(substitute(eventattributevar)) problem so that var-name can be used in varname
            data <- cbind(data, result)
            names(data)[length(data)] <- paste(variablename, "type", "match", "attr", "match", sep = ".") #deparse(substitute(eventattributevar))
            ## return the data frame with the variable bound to it
            return(data)
          }else{ 
            ## only return the 1 degree variable that was generated
            return(result)
          }	
        }
        
        ## (9) with typevar = 1 value selected and attributevar = valuematch
        if ( eventtypevalue != "valuematch" & eventtypevalue != "valuemix" & eventattributevalue == "valuematch"){
          result <- degreeCpp(time, weight, degreevar, eventtypevar, eventtypevalue, "1", eventattributevar, "1", "1", xlog, "d-typefilter-attributematch")		
          ## if returnData = TRUE => return the entire data frame as well as the 1 additional degree-variable
          if ( returnData == TRUE ) {
            ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
            ##TODO: also figure out the deparse(substitute(eventattributevar)) problem so that var-name can be used in varname
            data <- cbind(data, result)
            names(data)[length(data)] <- paste(variablename, "type", gsub(" ", "", eventtypevalue, fixed = TRUE), "attr", "match", sep = ".") #deparse(substitute(eventattributevar))
            ## return the data frame with the variable bound to it
            return(data)
          }else{ 
            ## only return the 1 degree variable that was generated
            return(result)
          }
        }
        
        ## (10) with typevar = valuematch and attributevar = 1 value selected
        if ( eventtypevalue == "valuematch" & eventattributevalue != "valuemix" & eventattributevalue != "valuematch"){
          result <- degreeCpp(time, weight, degreevar, eventtypevar, "1", "1", eventattributevar, eventattributevalue, "1", xlog, "d-typematch-attributefilter")		
          ## if returnData = TRUE => return the entire data frame as well as the 1 additional degree-variable
          if ( returnData == TRUE ) {
            ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
            ##TODO: also figure out the deparse(substitute(eventattributevar)) problem so that var-name can be used in varname
            data <- cbind(data, result)
            names(data)[length(data)] <- paste(variablename, "type", "match", "attr", gsub(" ", "", eventattributevalue, fixed = TRUE), sep = ".") #deparse(substitute(eventattributevar))
            ## return the data frame with the variable bound to it
            return(data)
          }else{ 
            ## only return the 1 degree variable that was generated
            return(result)
          }
        }
        
        ## (11) with typevar = 1 value selected and attributevar = 1 value selected
        if ( eventtypevalue != "valuematch" & eventtypevalue != "valuemix" & eventattributevalue != "valuematch" & eventattributevalue != "valuemix"){
          result <- degreeCpp(time, weight, degreevar, eventtypevar, eventtypevalue, "1", eventattributevar, eventattributevalue, "1", xlog, "d-typefilter-attributefilter")		
          ## if returnData = TRUE => return the entire data frame as well as the 1 additional degree-variable
          if ( returnData == TRUE ) {
            ##TODO: not simply add new variable - but check if a variable with this name already exists and replace it?
            ##TODO: also figure out the deparse(substitute(eventattributevar)) problem so that var-name can be used in varname
            data <- cbind(data, result)
            names(data)[length(data)] <- paste(variablename, "type", gsub(" ", "", eventtypevalue, fixed = TRUE), "attr", gsub(" ", "", eventattributevalue, fixed = TRUE), sep = ".") #deparse(substitute(eventattributevar))
            ## return the data frame with the variable bound to it
            return(data)
          }else{ 
            ## only return the 1 degree variable that was generated
            return(result)
          }
        }
        
        ## (12-1) with typevar = valuematch and attributevar = valuemix
        if ( eventtypevalue == "valuematch" & eventattributevalue == "valuemix"){
          for (i in unique(eventattributevar)){
            for (j in unique(eventattributevar)){
              if ( i != j ){
                ## calculate degree for the two distinct attribute-values
                temp <- degreeCpp(time, weight, degreevar, eventtypevar, "1", "1", eventattributevar, i, j, xlog, "d-typematch-attributemix")		
                data.short <- cbind(data.short, temp)
                names(data.short)[length(data.short)] <- paste(variablename, "type", "match", "attr", #deparse(substitute(eventattributevar))
                                                               gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                               sep = ".")
              }
            }#closes j-loop
            ## calculate degree for the attribute-values where i and j are the same (both i used)
            temp <- degreeCpp(time, weight, degreevar, eventtypevar, "1", "1", eventattributevar, i, i, xlog, "d-typematch-attributemix")		
            ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
            data.short <- cbind(data.short, temp)
            names(data.short)[length(data.short)] <- paste(variablename, "type", "match", "attr", #deparse(substitute(eventattributevar))
                                                           gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                           gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                           sep = ".")
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
        ## (13-1) with typevar = 1 value selected and attributevar = valuemix
        if ( eventtypevalue != "valuematch" & eventtypevalue != "valuemix" & eventattributevalue == "valuemix"){
          for (i in unique(eventattributevar)){
            for (j in unique(eventattributevar)){
              if ( i != j ){
                ## calculate degree for the two distinct attribute-values
                temp <- degreeCpp(time, weight, degreevar, eventtypevar, eventtypevalue, "1", eventattributevar, i, j, xlog, "d-typefilter-attributemix")		
                data.short <- cbind(data.short, temp)
                names(data.short)[length(data.short)] <- paste(variablename, "type", eventtypevalue, "attr", #deparse(substitute(eventattributevar))
                                                               gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                               sep = ".")
              }
            }#closes j-loop
            ## calculate degree for the attribute-values where i and j are the same (both i used)
            temp <- degreeCpp(time, weight, degreevar, eventtypevar, eventtypevalue, "1", eventattributevar, i, i, xlog, "d-typefilter-attributemix")		
            ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
            data.short <- cbind(data.short, temp)
            names(data.short)[length(data.short)] <- paste(variablename, "type", eventtypevalue, "attr", #deparse(substitute(eventattributevar))
                                                           gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                           gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                           sep = ".")
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
        ## (14-1) with typevar = valuemix and attributevar = valuematch   
        if ( eventtypevalue == "valuemix" & eventattributevalue == "valuematch"){
          for (i in unique(eventtypevar)){
            for (j in unique(eventtypevar)){
              if ( i != j ){
                ## calculate degree for the two distinct attribute-values
                temp <- degreeCpp(time, weight, degreevar, eventtypevar, i, j, eventattributevar, "1", "1", xlog, "d-typemix-attributematch")  	
                data.short <- cbind(data.short, temp)
                names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                               gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                               "attr", "match",
                                                               sep = ".")
              }
            }#closes j-loop
            ## calculate degree for the attribute-values where i and j are the same (both i used)
            temp <- degreeCpp(time, weight, degreevar, eventtypevar, i, i, eventattributevar, "1", "1", xlog, "d-typemix-attributematch")		
            ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
            data.short <- cbind(data.short, temp)
            names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                           gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                           gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                           "attr", "match",
                                                           sep = ".")
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
        ## (15-1) with typevar = valuemix and attributevar = 1 value selected
        if ( eventtypevalue == "valuemix" & eventattributevalue != "valuematch" & eventattributevalue != "valuemix"){
          for (i in unique(eventtypevar)){
            for (j in unique(eventtypevar)){
              if ( i != j ){
                ## calculate degree for the two distinct attribute-values
                temp <- degreeCpp(time, weight, degreevar, eventtypevar, i, j, eventattributevar, eventattributevalue, "1", xlog, "d-typemix-attributefilter")  	
                data.short <- cbind(data.short, temp)
                names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                               gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                               "attr", eventattributevalue,
                                                               sep = ".")
              }
            }#closes j-loop
            ## calculate degree for the attribute-values where i and j are the same (both i used)
            temp <- degreeCpp(time, weight, degreevar, eventtypevar, i, i, eventattributevar, eventattributevalue, "1", xlog, "d-typemix-attributefilter")		
            ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
            data.short <- cbind(data.short, temp)
            names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                           gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                           gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                           "attr", eventattributevalue,
                                                           sep = ".")
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
        ## (16-1) with typevar = valuemix and attributevar = valuemix
        if ( eventtypevalue == "valuemix" & eventattributevalue == "valuemix"){
          for (i in unique(eventtypevar)){
            for (j in unique(eventtypevar)){
              for (k in unique(eventattributevar)){
                for (l in unique(eventattributevar)){
                  temp <- degreeCpp(time, weight, degreevar, eventtypevar, i, j, eventattributevar, k, l, xlog, "d-typemix-attributemix")  	
                  data.short <- cbind(data.short, temp)
                  names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                                 gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                                 gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                                 "attr", #deparse(substitute(eventtypevar))
                                                                 gsub(" ", "", k, fixed = TRUE), #represents current event-actor type
                                                                 gsub(" ", "", l, fixed = TRUE), #represents actor type of past actions
                                                                 sep = ".")
                  #TODO: do not allow duplicate entries
                }#closes l-loop	
              }#closes k-loop
            }#closes j-loop
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
      }else{ #else length(eventattributevalue) > 2
        
        ## (12-2) with typevar = valuematch and attributevar = valuemix
        if ( eventtypevalue == "valuematch" ){
          for ( i in eventattributevalue ){
            for ( j in eventattributevalue ){
              if ( i != j ){
                ## calculate degree for the two distinct attribute-values
                temp <- degreeCpp(time, weight, degreevar, eventtypevar, "1", "1", eventattributevar, i, j, xlog, "d-typematch-attributemix")		
                data.short <- cbind(data.short, temp)
                names(data.short)[length(data.short)] <- paste(variablename, "type", "match", "attr", #deparse(substitute(eventattributevar))
                                                               gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                               sep = ".")
              }
            }#closes j-loop
            ## calculate degree for the attribute-values where i and j are the same (both i used)
            temp <- degreeCpp(time, weight, degreevar, eventtypevar, "1", "1", eventattributevar, i, i, xlog, "d-typematch-attributemix")		
            ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
            data.short <- cbind(data.short, temp)
            names(data.short)[length(data.short)] <- paste(variablename, "type", "match", "attr", #deparse(substitute(eventattributevar))
                                                           gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                           gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                           sep = ".")
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
        ## (13-2) with typevar = 1 value selected and attributevar = valuemix
        if ( eventtypevalue != "valuematch" & eventtypevalue != "valuemix" ){
          for (i in eventattributevalue ){
            for (j in eventattributevalue ){
              if ( i != j ){
                ## calculate degree for the two distinct attribute-values
                temp <- degreeCpp(time, weight, degreevar, eventtypevar, eventtypevalue, "1", eventattributevar, i, j, xlog, "d-typefilter-attributemix")		
                data.short <- cbind(data.short, temp)
                names(data.short)[length(data.short)] <- paste(variablename, "type", eventtypevalue, "attr", #deparse(substitute(eventattributevar))
                                                               gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                               sep = ".")
              }
            }#closes j-loop
            ## calculate degree for the attribute-values where i and j are the same (both i used)
            temp <- degreeCpp(time, weight, degreevar, eventtypevar, eventtypevalue, "1", eventattributevar, i, i, xlog, "d-typefilter-attributemix")		
            ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
            data.short <- cbind(data.short, temp)
            names(data.short)[length(data.short)] <- paste(variablename, "type", eventtypevalue, "attr", #deparse(substitute(eventattributevar))
                                                           gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                           gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                           sep = ".")
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
        ## (16-2) with typevar = valuemix and attributevar = valuemix (values)
        if ( eventtypevalue == "valuemix"){
          for ( i in unique(eventtypevar) ){
            for ( j in unique(eventtypevar) ){
              for ( k in eventattributevalue ){
                for ( l in eventattributevalue ){
                  temp <- degreeCpp(time, weight, degreevar, eventtypevar, i, j, eventattributevar, k, l, xlog, "d-typemix-attributemix")  	
                  data.short <- cbind(data.short, temp)
                  names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                                 gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                                 gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                                 "attr", #deparse(substitute(eventtypevar))
                                                                 gsub(" ", "", k, fixed = TRUE), #represents current event-actor type
                                                                 gsub(" ", "", l, fixed = TRUE), #represents actor type of past actions
                                                                 sep = ".")
                  #TODO: do not allow duplicate entries
                }#closes l-loop	
              }#closes k-loop
            }#closes j-loop
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
      } #closes else{length(eventattributevalue) >2 }
    }else{ #else: length(eventtypevalue) > 1:
      if ( length(eventattributevalue) == 1) {
        
        ## (14-2) with typevar = valuemix (values) and attributevar = valuematch   
        if ( eventattributevalue == "valuematch"){
          for ( i in eventtypevalue ){
            for ( j in eventtypevalue ){
              if ( i != j ){
                ## calculate degree for the two distinct attribute-values
                temp <- degreeCpp(time, weight, degreevar, eventtypevar, i, j, eventattributevar, "1", "1", xlog, "d-typemix-attributematch")  	
                data.short <- cbind(data.short, temp)
                names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                               gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                               "attr", "match",
                                                               sep = ".")
              }
            }#closes j-loop
            ## calculate degree for the attribute-values where i and j are the same (both i used)
            temp <- degreeCpp(time, weight, degreevar, eventtypevar, i, i, eventattributevar, "1", "1", xlog, "d-typemix-attributematch")		
            ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
            data.short <- cbind(data.short, temp)
            names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                           gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                           gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                           "attr", "match",
                                                           sep = ".")
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
        ## (15-2) with typevar = valuemix (values) and attributevar = 1 value selected
        if ( eventattributevalue != "valuematch" & eventattributevalue != "valuemix"){
          for ( i in eventtypevalue ){
            for ( j in eventtypevalue ){
              if ( i != j ){
                ## calculate degree for the two distinct attribute-values
                temp <- degreeCpp(time, weight, degreevar, eventtypevar, i, j, eventattributevar, eventattributevalue, "1", xlog, "d-typemix-attributefilter")  	
                data.short <- cbind(data.short, temp)
                names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                               gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                               "attr", eventattributevalue,
                                                               sep = ".")
              }
            }#closes j-loop
            ## calculate degree for the attribute-values where i and j are the same (both i used)
            temp <- degreeCpp(time, weight, degreevar, eventtypevar, i, i, eventattributevar, eventattributevalue, "1", xlog, "d-typemix-attributefilter")		
            ##TODO: calculate one effect each for a filtered-variable (not just match, but also filter?)
            data.short <- cbind(data.short, temp)
            names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                           gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                           gsub(" ", "", i, fixed = TRUE), #represents actor type of past actions
                                                           "attr", eventattributevalue,
                                                           sep = ".")
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
        ## (16-2) with typevar = valuemix (values) and attributevar = valuemix
        if ( eventattributevalue == "valuemix"){
          for (i in eventtypevalue ){
            for (j in eventtypevalue ){
              for (k in unique(eventattributevar)){
                for (l in unique(eventattributevar)){
                  temp <- degreeCpp(time, weight, degreevar, eventtypevar, i, j, eventattributevar, k, l, xlog, "d-typemix-attributemix")  	
                  data.short <- cbind(data.short, temp)
                  names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                                 gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                                 gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                                 "attr", #deparse(substitute(eventtypevar))
                                                                 gsub(" ", "", k, fixed = TRUE), #represents current event-actor type
                                                                 gsub(" ", "", l, fixed = TRUE), #represents actor type of past actions
                                                                 sep = ".")
                  #TODO: do not allow duplicate entries
                }#closes l-loop	
              }#closes k-loop
            }#closes j-loop
          }#closes i-loop
          ## return data frames
          if ( returnData == TRUE ) {
            data <- cbind(data, data.short)
            return(data)
          }else{
            return(data.short)
          }
        }
        
      }else{ #else: length(eventattributevalue) > 1
        
        ## (16-2) with typevar = valuemix (values) and attributevar = valuemix (values)
        for ( i in eventtypevalue ){
          for ( j in eventtypevalue ){
            for ( k in eventattributevalue ){
              for ( l in eventattributevalue ){
                temp <- degreeCpp(time, weight, degreevar, eventtypevar, i, j, eventattributevar, k, l, xlog, "d-typemix-attributemix")  	
                data.short <- cbind(data.short, temp)
                names(data.short)[length(data.short)] <- paste(variablename, "type", #deparse(substitute(eventtypevar))
                                                               gsub(" ", "", i, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", j, fixed = TRUE), #represents actor type of past actions
                                                               "attr", #deparse(substitute(eventtypevar))
                                                               gsub(" ", "", k, fixed = TRUE), #represents current event-actor type
                                                               gsub(" ", "", l, fixed = TRUE), #represents actor type of past actions
                                                               sep = ".")
                #TODO: do not allow duplicate entries
              }#closes l-loop	
            }#closes k-loop
          }#closes j-loop
        }#closes i-loop
        ## return data frames
        if ( returnData == TRUE ) {
          data <- cbind(data, data.short)
          return(data)
        }else{
          return(data.short)
        }
      }#closes else{length(eventattributevalue) > 1}		
    }#closes if-else "( length(eventtypevalue) == 1 ){}"
  }#closes -if both eventtypevar and eventattributevar are selected
}

####################################################################
##	FourCycle calculation
####################################################################

get.fourCycle.stat <- function(data, time, sender, target, halflife, weight = NULL, eventtypevar = NULL, eventtypevalue = "standard", 
                               eventattributevar = NULL, eventattributeAB = NULL, eventattributeAJ = NULL, 
                               eventattributeIB = NULL,eventattributeIJ = NULL, variablename = "fourCycle", returnData = TRUE, 
                               showprogressbar = FALSE){
  
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
  
  ## check if weight-var is defined (if not -> create it)
  if ( is.null(weight) ) {
    weight <- rep(1, length(time))
  }
  if ( !is.numeric(weight) ) {
    stop("'", as.name(weight), "' variable is not numeric.") #TODO: deparse(substitute(eventattributevar)) ?
  }
  
  ## check if event-type inputs are available and correctly specified
  if ( !is.null(eventtypevar) ) {
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
  if ( is.null(eventattributevar) == FALSE ) {
    eventattributevar <- as.character(eventattributevar)
    if ( is.null(eventattributeAB) & is.null(eventattributeAJ) & is.null(eventattributeIB) & is.null(eventattributeIJ) ){
      stop("No 'eventattribute__' provided. Provide a string value by which the events are filtered.", )
    }
    # check if eventattributevalue is part of the variable
    if ( is.null(eventattributeAB) == FALSE){
      if ( length(grep(eventattributeAB, eventattributevar)) == 0 ) {
        ##TODO: #deparse(substitute(eventtypevar))
        stop("Value '", eventattributeAB, "' is not an element of '", deparse(substitute(eventattributevar)) , "'.") 
      }
    }
    if ( is.null(eventattributeAJ) == FALSE){
      if ( length(grep(eventattributeAJ, eventattributevar)) == 0 ) {
        ##TODO: #deparse(substitute(eventtypevar))
        stop("Value '", eventattributeAJ, "' is not an element of '", deparse(substitute(eventattributevar)) , "'.") 
      }
    }
    if ( is.null(eventattributeIB) == FALSE){
      if ( length(grep(eventattributeIB, eventattributevar)) == 0 ) {
        ##TODO: #deparse(substitute(eventtypevar))
        stop("Value '", eventattributeIB, "' is not an element of '", deparse(substitute(eventattributevar)) , "'.") 
      }
    }
    if ( is.null(eventattributeIJ) == FALSE){
      if ( length(grep(eventattributeIJ, eventattributevar)) == 0 ) {
        ##TODO: #deparse(substitute(eventtypevar))
        stop("Value '", eventattributeIJ, "' is not an element of '", deparse(substitute(eventattributevar)) , "'.") 
      }
    }
  }
  
  ## check if variablename makes sense (no " " etc.)
  variablename <- gsub(" ", "", variablename, fixed = TRUE)
  
  ## create simple data set to be returned for degree calcuations with more than 1 output-variable
  ##TODO: should there be an event-id-variable?? => that would be useful here
  data.short <- data.frame(time)
  
  ## calculate part of decay function
  xlog <- log(2)/halflife 
  
  ####### calculate stat
  ## create placeholder-variables to be used in the cpp-Function
  placeholder <- rep("1", length(time))
  
  ## calculate the fourCycle effects for each event
  
  if (eventtypevalue == "standard"){
    
    ## (1) fourCycle without filter
    if ( is.null(eventattributeAB) & is.null(eventattributeAJ) & is.null(eventattributeIB) & is.null(eventattributeIJ)){
      result <- fourCycleCpp(sender, target, placeholder, time, weight, xlog, placeholder, "1", placeholder , "1", placeholder, "1", placeholder, "1", eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (1)
    
    ## (2) fourCycle with AB-filter
    if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAJ) & is.null(eventattributeIB) & is.null(eventattributeIJ)){		
      result <- fourCycleCpp(sender, target, placeholder, time, weight, xlog, eventattributevar, eventattributeAB, placeholder , "1", placeholder, "1", placeholder, "1", eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (2)
    
    ## (3) fourCycle with AJ-filter
    if ( is.null(eventattributeAB) & is.null(eventattributeAJ) == FALSE & is.null(eventattributeIB) & is.null(eventattributeIJ)){		
      result <- fourCycleCpp(sender, target, placeholder, time, weight, xlog, placeholder, "1", eventattributevar , eventattributeAJ, placeholder, "1", placeholder, "1", eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (3)
    
    ## (4) fourCycle with IB-filter
    if ( is.null(eventattributeAB) & is.null(eventattributeAJ) & is.null(eventattributeIB) == FALSE & is.null(eventattributeIJ)){		
      result <- fourCycleCpp(sender, target, placeholder, time, weight, xlog, placeholder, "1", placeholder , "1", eventattributevar, eventattributeIB, placeholder, "1", eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (4)
    
    ## (5) fourCycle with IJ-filter
    if ( is.null(eventattributeAB) & is.null(eventattributeAJ) & is.null(eventattributeIB) & is.null(eventattributeIJ)== FALSE ){		
      result <- fourCycleCpp(sender, target, placeholder, time, weight, xlog, placeholder, "1", placeholder , "1", placeholder, "1", eventattributevar, eventattributeIJ, eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (5)
    
    ## (6) fourCycle with AB + AJ filter
    if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAJ) == FALSE & is.null(eventattributeIB) & is.null(eventattributeIJ)){		
      result <- fourCycleCpp(sender, target, placeholder, time, weight, xlog, eventattributevar, eventattributeAB, eventattributevar , eventattributeIB, placeholder, "1", placeholder, "1", eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (6)
    
    ## (7) fourCycle with AB + IB filter
    if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAJ) & is.null(eventattributeIB) == FALSE & is.null(eventattributeIJ)){		
      result <- fourCycleCpp(sender, target, placeholder, time, weight, xlog, eventattributevar, eventattributeAB, placeholder , "1", eventattributevar, eventattributeIB, placeholder, "1", eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (7)
    
    ## (8) fourCycle with AB + IJ filter
    if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAJ) & is.null(eventattributeIB) & is.null(eventattributeIJ) == FALSE){		
      result <- fourCycleCpp(sender, target, placeholder, time, weight, xlog, eventattributevar, eventattributeAB, placeholder , "1", placeholder, "1", eventattributevar, eventattributeIJ, eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (8)
    
    ## (9) fourCycle with AJ + IB filter
    if ( is.null(eventattributeAB) & is.null(eventattributeAJ) == FALSE & is.null(eventattributeIB) == FALSE & is.null(eventattributeIJ)  ){		
      result <- fourCycleCpp(sender, target, placeholder, time, weight, xlog, placeholder, "1", eventattributevar , eventattributeAJ, eventattributevar, eventattributeIB, placeholder, "1", eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (9)
    
    ## (10) fourCycle with AJ + IJ filter
    if ( is.null(eventattributeAB)  & is.null(eventattributeAJ) == FALSE & is.null(eventattributeIB)  & is.null(eventattributeIJ) == FALSE ){		
      result <- fourCycleCpp(sender, target, placeholder, time, weight, xlog, placeholder, "1", eventattributevar , eventattributeAJ, placeholder, "1", eventattributevar, eventattributeIJ, eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (10)
    
    ## (11) fourCycle with IB + IJ filter
    if ( is.null(eventattributeAB) & is.null(eventattributeAJ) & is.null(eventattributeIB) == FALSE & is.null(eventattributeIJ) == FALSE ){		
      result <- fourCycleCpp(sender, target, placeholder, time, weight, xlog, placeholder, "1", placeholder , "1", eventattributevar, eventattributeIB, eventattributevar, eventattributeIJ, eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (11)
    
    ## (12) fourCycle with AB, AJ, IB filter
    if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAJ) == FALSE & is.null(eventattributeIB) == FALSE & is.null(eventattributeIJ) ){		
      result <- fourCycleCpp(sender, target, placeholder, time, weight, xlog, eventattributevar, eventattributeAB, eventattributevar , eventattributeAJ, eventattributevar, eventattributeIB, placeholder, "1", eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (12)
    
    ## (13) fourCycle with AB, AJ, IJ filter
    if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAJ) == FALSE & is.null(eventattributeIB)  & is.null(eventattributeIJ) == FALSE ){		
      result <- fourCycleCpp(sender, target, placeholder, time, weight, xlog, eventattributevar, eventattributeAB, eventattributevar , eventattributeAJ, placeholder, "1", eventattributevar, eventattributeIJ, eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (13)
    
    ## (14) fourCycle with AB, IB, IJ filter
    if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAJ)  & is.null(eventattributeIB) == FALSE & is.null(eventattributeIJ) == FALSE ){		
      result <- fourCycleCpp(sender, target, placeholder, time, weight, xlog, eventattributevar, eventattributeAB, placeholder , "1", eventattributevar, eventattributeIB, eventattributevar, eventattributeIJ, eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (16)
    
    ## (15) fourCycle with AJ, IB, IJ filter
    if ( is.null(eventattributeAB)  & is.null(eventattributeAJ) == FALSE & is.null(eventattributeIB) == FALSE & is.null(eventattributeIJ) == FALSE ){		
      result <- fourCycleCpp(sender, target, placeholder, time, weight, xlog, placeholder, "1", eventattributevar , eventattributeAJ, eventattributevar, eventattributeIB, eventattributevar, eventattributeIJ, eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (16)
    
    ## (16) fourCycle with AB, AJ, IB, IJ filter
    if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAJ) == FALSE & is.null(eventattributeIB) == FALSE & is.null(eventattributeIJ) == FALSE ){		
      result <- fourCycleCpp(sender, target, placeholder, time, weight, xlog, eventattributevar, eventattributeAB, eventattributevar , eventattributeAJ, eventattributevar, eventattributeIB, eventattributevar, eventattributeIJ, eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (16)
    
  }else{ #else: if eventtypevar is set
    
    ## (1) fourCycle without filter
    if ( is.null(eventattributeAB) & is.null(eventattributeAJ) & is.null(eventattributeIB) & is.null(eventattributeIJ)){		
      result <- fourCycleCpp(sender, target, eventtypevar, time, weight, xlog, placeholder, "1", placeholder , "1", placeholder, "1", placeholder, "1", eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (1)
    
    ## (2) fourCycle with AB-filter
    if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAJ) & is.null(eventattributeIB) & is.null(eventattributeIJ)){		
      result <- fourCycleCpp(sender, target, eventtypevar, time, weight, xlog, eventattributevar, eventattributeAB, placeholder , "1", placeholder, "1", placeholder, "1", eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (2)
    
    ## (3) fourCycle with AJ-filter
    if ( is.null(eventattributeAB) & is.null(eventattributeAJ) == FALSE & is.null(eventattributeIB) & is.null(eventattributeIJ)){		
      result <- fourCycleCpp(sender, target, eventtypevar, time, weight, xlog, placeholder, "1", eventattributevar , eventattributeAJ, placeholder, "1", placeholder, "1", eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (3)
    
    ## (4) fourCycle with IB-filter
    if ( is.null(eventattributeAB) & is.null(eventattributeAJ) & is.null(eventattributeIB) == FALSE & is.null(eventattributeIJ)){		
      result <- fourCycleCpp(sender, target, eventtypevar, time, weight, xlog, placeholder, "1", placeholder , "1", eventattributevar, eventattributeIB, placeholder, "1", eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (4)
    
    ## (5) fourCycle with IJ-filter
    if ( is.null(eventattributeAB) & is.null(eventattributeAJ) & is.null(eventattributeIB) & is.null(eventattributeIJ)== FALSE ){		
      result <- fourCycleCpp(sender, target, eventtypevar, time, weight, xlog, placeholder, "1", placeholder , "1", placeholder, "1", eventattributevar, eventattributeIJ, eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (5)
    
    ## (6) fourCycle with AB + AJ filter
    if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAJ) == FALSE & is.null(eventattributeIB) & is.null(eventattributeIJ)){		
      result <- fourCycleCpp(sender, target, eventtypevar, time, weight, xlog, eventattributevar, eventattributeAB, eventattributevar , eventattributeIB, placeholder, "1", placeholder, "1", eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (6)
    
    ## (7) fourCycle with AB + IB filter
    if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAJ) & is.null(eventattributeIB) == FALSE & is.null(eventattributeIJ)){		
      result <- fourCycleCpp(sender, target, eventtypevar, time, weight, xlog, eventattributevar, eventattributeAB, placeholder , "1", eventattributevar, eventattributeIB, placeholder, "1", eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (7)
    
    ## (8) fourCycle with AB + IJ filter
    if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAJ) & is.null(eventattributeIB) & is.null(eventattributeIJ) == FALSE){		
      result <- fourCycleCpp(sender, target, eventtypevar, time, weight, xlog, eventattributevar, eventattributeAB, placeholder , "1", placeholder, "1", eventattributevar, eventattributeIJ, eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (8)
    
    ## (9) fourCycle with AJ + IB filter
    if ( is.null(eventattributeAB) & is.null(eventattributeAJ) == FALSE & is.null(eventattributeIB) == FALSE & is.null(eventattributeIJ)  ){		
      result <- fourCycleCpp(sender, target, eventtypevar, time, weight, xlog, placeholder, "1", eventattributevar , eventattributeAJ, eventattributevar, eventattributeIB, placeholder, "1", eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (9)
    
    ## (10) fourCycle with AJ + IJ filter
    if ( is.null(eventattributeAB)  & is.null(eventattributeAJ) == FALSE & is.null(eventattributeIB)  & is.null(eventattributeIJ) == FALSE ){		
      result <- fourCycleCpp(sender, target, eventtypevar, time, weight, xlog, placeholder, "1", eventattributevar , eventattributeAJ, placeholder, "1", eventattributevar, eventattributeIJ, eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (10)
    
    ## (11) fourCycle with IB + IJ filter
    if ( is.null(eventattributeAB) & is.null(eventattributeAJ) & is.null(eventattributeIB) == FALSE & is.null(eventattributeIJ) == FALSE ){		
      result <- fourCycleCpp(sender, target, eventtypevar, time, weight, xlog, placeholder, "1", placeholder , "1", eventattributevar, eventattributeIB, eventattributevar, eventattributeIJ, eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (11)
    
    ## (12) fourCycle with AB, AJ, IB filter
    if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAJ) == FALSE & is.null(eventattributeIB) == FALSE & is.null(eventattributeIJ) ){		
      result <- fourCycleCpp(sender, target, eventtypevar, time, weight, xlog, eventattributevar, eventattributeAB, eventattributevar , eventattributeAJ, eventattributevar, eventattributeIB, placeholder, "1", eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (12)
    
    ## (13) fourCycle with AB, AJ, IJ filter
    if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAJ) == FALSE & is.null(eventattributeIB)  & is.null(eventattributeIJ) == FALSE ){		
      result <- fourCycleCpp(sender, target, eventtypevar, time, weight, xlog, eventattributevar, eventattributeAB, eventattributevar , eventattributeAJ, placeholder, "1", eventattributevar, eventattributeIJ, eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (13)
    
    ## (14) fourCycle with AB, IB, IJ filter
    if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAJ)  & is.null(eventattributeIB) == FALSE & is.null(eventattributeIJ) == FALSE ){		
      result <- fourCycleCpp(sender, target, eventtypevar, time, weight, xlog, eventattributevar, eventattributeAB, placeholder , "1", eventattributevar, eventattributeIB, eventattributevar, eventattributeIJ, eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (16)
    
    ## (15) fourCycle with AJ, IB, IJ filter
    if ( is.null(eventattributeAB)  & is.null(eventattributeAJ) == FALSE & is.null(eventattributeIB) == FALSE & is.null(eventattributeIJ) == FALSE ){		
      result <- fourCycleCpp(sender, target, eventtypevar, time, weight, xlog, placeholder, "1", eventattributevar , eventattributeAJ, eventattributevar, eventattributeIB, eventattributevar, eventattributeIJ, eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (16)
    
    ## (16) fourCycle with AB, AJ, IB, IJ filter
    if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAJ) == FALSE & is.null(eventattributeIB) == FALSE & is.null(eventattributeIJ) == FALSE ){		
      result <- fourCycleCpp(sender, target, eventtypevar, time, weight, xlog, eventattributevar, eventattributeAB, eventattributevar , eventattributeAJ, eventattributevar, eventattributeIB, eventattributevar, eventattributeIJ, eventtypevalue )		
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 degree variable that was generated
        return(result)
      }
    }#closes (16)
    
  }#closes else eventtypevar != null
}

####################################################################
##	Similarity calculation
####################################################################

get.similarity.stat <- function(data, time, sender, target, 
                                senderOrTarget = "sender",
                                whichSimilarity = NULL, 
                                halflife.last.event = NULL, 
                                halflife.time.between.events = NULL,
                                eventtypevar = NULL, 
                                eventattributevar = NULL, eventattributevalue = NULL,
                                variablename = "similarity", returnData = TRUE, 
                                showprogressbar = FALSE){
  
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
      stop("'", time, "' is not sorted. Sort data frame according to 
           the event sequence.")
    }
    }
  
  ## check if senderOrTarget is specified
  if ( (senderOrTarget == "sender" | senderOrTarget == "target") == FALSE ){
    stop("'senderOrTarget' not correctly specified. Choose either 'sender' 
         or 'target' for the respective similarity measure.")
  }
  
  ## check if whichSimilarity is specified
  if ( is.null(whichSimilarity) == FALSE ){
    if ( (whichSimilarity == "total" | whichSimilarity == "average" ) == FALSE ){
      stop("'whichSimilarity' not correctly specified. Choose either 'NULL', 
           'total' or 'average' for the respective similarity measure.")
    }
    }
  
  ## average/total vs. with halflife
  if ( is.null(whichSimilarity)  & is.null(halflife.last.event) ){
    stop("Specify type of similarity measure - either by chosing 'total' or 
         'average' in whichSimilarity, or by specifying a value for 
         'halflife.last.event'.")
  }
  if ( is.null(whichSimilarity) == FALSE & is.null(halflife.last.event) ==FALSE ){
    stop("Cannot specify 'whichSimilarity' as well as 'halflife.last.event'.")
  }
  if ( is.null(halflife.last.event) & 
         is.null(halflife.time.between.events)==FALSE){
    stop("Please also specify a value for 'halflife.last.event'.")
  }
  
  ## check if event-type inputs are available and correctly specified
  if ( !is.null(eventtypevar) ) {
    eventtypevar <- as.character(eventtypevar)
    if ( length(unique(eventtypevar)) != 2 ){ 
      stop("'eventtypevar' is not a dummy variable.")
    }
  }
  
  ## check if event-attribute inputs are available and correctly specified
  if ( is.null(eventattributevar) == FALSE ) {
    eventattributevar <- as.character(eventattributevar)
    if ( is.null(eventattributevalue) ){
      stop("No 'eventattributevalue' provided. Provide a string value by 
           which the events are filtered.", )
    }
    # check if eventattributevalue is part of the variable
    if ( is.null(eventattributevalue) == FALSE){
      if ( length(grep(eventattributevalue, eventattributevar)) == 0 ) {
        ##TODO: #deparse(substitute(eventtypevar))
        stop("Value '", eventattributevalue, "' is not an element of '", 
             deparse(substitute(eventattributevar)) , "'.") 
      }
    }
    }
  
  ## check if variablename makes sense (no " " etc.)
  variablename <- gsub(" ", "", variablename, fixed = TRUE)
  
  ## create simple data set to be returned for degree calcuations with more than 1 output-variable
  ##TODO: should there be an event-id-variable?? => that would be useful here
  data.short <- data.frame(time)
  
  ## calculate part of decay function
  xlog.last.event <- log(2)/halflife.last.event 
  
  ####### calculate stat
  ## create placeholder-variables to be used in the cpp-Function
  placeholder <- rep("1", length(time))
  
  ## calculate the similarity effects for each event
  # sender similarity
  if ( senderOrTarget == "sender" ){
    
    if ( is.null(whichSimilarity) == FALSE ){
      if ( whichSimilarity == "total" ){
        ##########
        if ( is.null(eventtypevar) & is.null(eventattributevar)){
          ## (1a) sender, total				
          result <- similarityTotalAverageCpp(sender, target, time, 
                                              placeholder, "1", placeholder, 
                                              "total", "nomatch", "sender" )		
          if ( returnData == TRUE ) {
            data <- cbind(data, result)
            names(data)[length(data)] <- paste(variablename, "total", sep = ".")
            ## return the data frame with the variable bound to it
            return(data)
          }else{ 
            ## only return the 1 degree variable that was generated
            return(result)
          }		
        }
        if ( is.null(eventtypevar) & is.null(eventattributevar)==FALSE ){
          ## (2a) sender, total, filter
          result <- similarityTotalAverageCpp(sender, target, time, 
                                              eventattributevar, 
                                              eventattributevalue, placeholder, 
                                              "total", "nomatch", "sender" )		
          if ( returnData == TRUE ) {
            data <- cbind(data, result)
            names(data)[length(data)] <- paste(variablename, "total", 
                                               eventattributevalue, sep = ".")
            ## return the data frame with the variable bound to it
            return(data)
          }else{ 
            ## only return the 1 degree variable that was generated
            return(result)
          }
          
        }
        if ( is.null(eventtypevar) == FALSE & is.null(eventattributevar)==FALSE ){
          ## (3a) sender, total, match, filter
          result <- similarityTotalAverageCpp(sender, target, time, 
                                              eventattributevar, 
                                              eventattributevalue, eventtypevar, 
                                              "total", "match", "sender" )		
          if ( returnData == TRUE ) {
            data <- cbind(data, result)
            names(data)[length(data)] <- paste(variablename, "total", 
                                               eventattributevalue, 
                                               "sameType", sep = ".")
            ## return the data frame with the variable bound to it
            return(data)
          }else{ 
            ## only return the 1 degree variable that was generated
            return(result)
          }
        }
        if ( is.null(eventtypevar) == FALSE & is.null(eventattributevar) ){
          ## (4a) sender, total, match
          result <- similarityTotalAverageCpp(sender, target, time, placeholder,
                                              "1", eventtypevar, "total", 
                                              "match", "sender" )		
          if ( returnData == TRUE ) {
            data <- cbind(data, result)
            names(data)[length(data)] <- paste(variablename, "total", 
                                               "sameType", sep = ".")
            ## return the data frame with the variable bound to it
            return(data)
          }else{ 
            ## only return the 1 degree variable that was generated
            return(result)
          }
        }	
        #closes whichSimilarity == "total"
      } else if (  whichSimilarity == "average"){ 
        ##########
        if ( is.null(eventtypevar) & is.null(eventattributevar)){
          ## (1b) sender, average				
          result <- similarityTotalAverageCpp(sender, target, time, placeholder,
                                              "1", placeholder, "average", 
                                              "nomatch", "sender" )		
          if ( returnData == TRUE ) {
            data <- cbind(data, result)
            names(data)[length(data)] <- paste(variablename, "average", sep = ".")
            ## return the data frame with the variable bound to it
            return(data)
          }else{ 
            ## only return the 1 degree variable that was generated
            return(result)
          }		
        }
        if ( is.null(eventtypevar) & is.null(eventattributevar)==FALSE ){
          ## (2b) sender, average, filter
          result <- similarityTotalAverageCpp(sender, target, time, 
                                              eventattributevar, 
                                              eventattributevalue, placeholder, 
                                              "average", "nomatch", "sender" )		
          if ( returnData == TRUE ) {
            data <- cbind(data, result)
            names(data)[length(data)] <- paste(variablename, "average", 
                                               eventattributevalue, sep = ".")
            ## return the data frame with the variable bound to it
            return(data)
          }else{ 
            ## only return the 1 degree variable that was generated
            return(result)
          }
          
        }
        if ( is.null(eventtypevar) == FALSE & is.null(eventattributevar)==FALSE ){
          ## (3b) sender, average, match, filter
          result <- similarityTotalAverageCpp(sender, target, time, 
                                              eventattributevar, 
                                              eventattributevalue, eventtypevar,
                                              "average", "match", "sender" )		
          if ( returnData == TRUE ) {
            data <- cbind(data, result)
            names(data)[length(data)] <- paste(variablename, "average", 
                                               eventattributevalue, "sameType",
                                               sep = ".")
            ## return the data frame with the variable bound to it
            return(data)
          }else{ 
            ## only return the 1 degree variable that was generated
            return(result)
          }
        }
        if ( is.null(eventtypevar) == FALSE & is.null(eventattributevar) ){
          ## (4b) sender, average, match
          result <- similarityTotalAverageCpp(sender, target, time, placeholder, 
                                              "1", eventtypevar, "average", 
                                              "match", "sender" )		
          if ( returnData == TRUE ) {
            data <- cbind(data, result)
            names(data)[length(data)] <- paste(variablename, "average", 
                                               "sameType", sep = ".")
            ## return the data frame with the variable bound to it
            return(data)
          }else{ 
            ## only return the 1 degree variable that was generated
            return(result)
          }
        }
      }#closes whichSimilarity == "average"	
    }#closes is.null(whichSimilarity) == FALSE 
    ##########
    # with 1 hallife parameter set
    if ( is.null(whichSimilarity) & is.null(halflife.last.event)==FALSE & 
           is.null(halflife.time.between.events)){
      if ( is.null(eventtypevar) & is.null(eventattributevar)){
        ## (1c) sender, 1time				
        result <- similaritySimpleCpp(sender, target, time, xlog.last.event, 
                                      placeholder, "1", placeholder, 
                                      "nomatch", "sender" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "overTime", sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }		
      }
      if ( is.null(eventtypevar) & is.null(eventattributevar)==FALSE ){
        ## (2c) sender, 1time, filter
        result <- similaritySimpleCpp(sender, target, time, xlog.last.event, 
                                      eventattributevar, eventattributevalue, 
                                      placeholder,  "nomatch", "sender" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "overTime", 
                                             eventattributevalue, sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }
        
      }
      if ( is.null(eventtypevar) == FALSE & is.null(eventattributevar)==FALSE ){
        ## (3c) sender, 1time, match, filter
        result <- similaritySimpleCpp(sender, target, time, xlog.last.event, 
                                      eventattributevar, eventattributevalue, 
                                      eventtypevar,  "match", "sender" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "overTime", 
                                             eventattributevalue, "sameType",
                                             sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }
      }
      if ( is.null(eventtypevar) == FALSE & is.null(eventattributevar) ){
        ## (4c) sender, 1time, match
        result <- similaritySimpleCpp(sender, target, time, xlog.last.event,
                                      placeholder, "1", eventtypevar, "match",
                                      "sender" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "overTime", 
                                             "sameType", sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }
      }
    }#closes if --1halflife parameter set--command
    
    ##########
    if ( is.null(whichSimilarity) & is.null(halflife.last.event)==FALSE & 
           is.null(halflife.time.between.events)==FALSE){
      if ( is.null(eventtypevar) & is.null(eventattributevar)){
        ## (1d) sender, 2times				
        result <- similarityComplexCpp(sender, target, time, xlog.last.event, 
                                       halflife.time.between.events, 
                                       placeholder, "1", placeholder, "nomatch",
                                       "sender" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "overTimes", 
                                             sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }		
      }
      if ( is.null(eventtypevar) & is.null(eventattributevar)==FALSE ){
        ## (2d) sender, 2times, filter
        result <- similarityComplexCpp(sender, target, time, xlog.last.event, 
                                       halflife.time.between.events, 
                                       eventattributevar, eventattributevalue,
                                       placeholder,  "nomatch", "sender" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "overTimes", 
                                             eventattributevalue, sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }
        
      }
      if ( is.null(eventtypevar) == FALSE & is.null(eventattributevar)==FALSE ){
        ## (3d) sender, 2times, match, filter
        result <- similarityComplexCpp(sender, target, time, xlog.last.event, 
                                       halflife.time.between.events, 
                                       eventattributevar, eventattributevalue, 
                                       eventtypevar,  "match", "sender" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "overTimes", 
                                             eventattributevalue, "sameType",
                                             sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }
      }
      if ( is.null(eventtypevar) == FALSE & is.null(eventattributevar) ){
        ## (4d) sender, 2times, match
        result <- similarityComplexCpp(sender, target, time, xlog.last.event, 
                                       halflife.time.between.events, 
                                       placeholder, "1", eventtypevar, "match", 
                                       "sender" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "overTimes", 
                                             "sameType", sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }
      }
    }#closes if --two-halflife-parameters-set--command
    
    # target similarity
  } else if ( senderOrTarget == "target"){
    if ( whichSimilarity == "total"){
      ##########
      if ( is.null(eventtypevar) & is.null(eventattributevar)){
        ## (1a) sender, total				
        result <- similarityTotalAverageCpp(sender, target, time, placeholder,
                                            "1", placeholder, "total", 
                                            "nomatch", "target" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "total", sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }		
      }
      if ( is.null(eventtypevar) & is.null(eventattributevar)==FALSE ){
        ## (2a) sender, total, filter
        result <- similarityTotalAverageCpp(sender, target, time, 
                                            eventattributevar,
                                            eventattributevalue, placeholder,
                                            "total", "nomatch", "target" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "total",
                                             eventattributevalue, sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }
        
      }
      if ( is.null(eventtypevar) == FALSE & is.null(eventattributevar)==FALSE ){
        ## (3a) sender, total, match, filter
        result <- similarityTotalAverageCpp(sender, target, time, 
                                            eventattributevar, 
                                            eventattributevalue, eventtypevar, 
                                            "total", "match", "target" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "total", 
                                             eventattributevalue, "sameType",
                                             sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }
      }
      if ( is.null(eventtypevar) == FALSE & is.null(eventattributevar) ){
        ## (4a) sender, total, match
        result <- similarityTotalAverageCpp(sender, target, time, placeholder, 
                                            "1", eventtypevar, "total", "match",
                                            "target" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "total", "sameType",
                                             sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }
      }	
    }#closes whichSimilarity == "total"
    if (  whichSimilarity == "average"){
      ##########
      if ( is.null(eventtypevar) & is.null(eventattributevar)){
        ## (1b) sender, average				
        result <- similarityTotalAverageCpp(sender, target, time, placeholder, 
                                            "1", placeholder, "average", 
                                            "nomatch", "target" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "average", sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }		
      }
      if ( is.null(eventtypevar) & is.null(eventattributevar)==FALSE ){
        ## (2b) sender, average, filter
        result <- similarityTotalAverageCpp(sender, target, time, 
                                            eventattributevar, 
                                            eventattributevalue, placeholder,
                                            "average", "nomatch", "target" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "average", 
                                             eventattributevalue, sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }
        
      }
      if ( is.null(eventtypevar) == FALSE & is.null(eventattributevar)==FALSE ){
        ## (3b) sender, average, match, filter
        result <- similarityTotalAverageCpp(sender, target, time, 
                                            eventattributevar, 
                                            eventattributevalue, eventtypevar, 
                                            "average", "match", "target" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "average", 
                                             eventattributevalue, "sameType",
                                             sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }
      }
      if ( is.null(eventtypevar) == FALSE & is.null(eventattributevar) ){
        ## (4b) sender, average, match
        result <- similarityTotalAverageCpp(sender, target, time, placeholder, 
                                            "1", eventtypevar, "average",
                                            "match", "target" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "average",
                                             "sameType", sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }
      }
    }#closes whichSimilarity == "average"		
    ##########
    # with 1 hallife parameter set
    if ( is.null(whichSimilarity) & is.null(halflife.last.event)==FALSE & 
           is.null(halflife.time.between.events)){
      if ( is.null(eventtypevar) & is.null(eventattributevar)){
        ## (1c) sender, 1time				
        result <- similaritySimpleCpp(sender, target, time, xlog.last.event, 
                                      placeholder, "1", placeholder, 
                                      "nomatch", "target" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "overTime", 
                                             sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }		
      }
      if ( is.null(eventtypevar) & is.null(eventattributevar)==FALSE ){
        ## (2c) sender, 1time, filter
        result <- similaritySimpleCpp(sender, target, time, xlog.last.event, 
                                      eventattributevar, eventattributevalue, 
                                      placeholder,  "nomatch", "target" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "overTime", 
                                             eventattributevalue, sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }
        
      }
      if ( is.null(eventtypevar) == FALSE & is.null(eventattributevar)==FALSE ){
        ## (3c) sender, 1time, match, filter
        result <- similaritySimpleCpp(sender, target, time, xlog.last.event, 
                                      eventattributevar, eventattributevalue, 
                                      eventtypevar,  "match", "target" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "overTime", 
                                             eventattributevalue, "sameType",
                                             sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }
      }
      if ( is.null(eventtypevar) == FALSE & is.null(eventattributevar) ){
        ## (4c) sender, 1time, match
        result <- similaritySimpleCpp(sender, target, time, xlog.last.event, 
                                      placeholder, "1", eventtypevar, "match", 
                                      "target" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "overTime", 
                                             "sameType", sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }
      }
    }#closes if --1halflife parameter set--command
    
    ##########
    if ( is.null(whichSimilarity) & is.null(halflife.last.event)==FALSE & 
           is.null(halflife.time.between.events)==FALSE){
      if ( is.null(eventtypevar) & is.null(eventattributevar)){
        ## (1d) sender, 2times				
        result <- similarityComplexCpp(sender, target, time, xlog.last.event, 
                                       halflife.time.between.events, 
                                       placeholder, "1", placeholder, "nomatch",
                                       "target" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "overTimes", 
                                             sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }		
      }
      if ( is.null(eventtypevar) & is.null(eventattributevar)==FALSE ){
        ## (2d) sender, 2times, filter
        result <- similarityComplexCpp(sender, target, time, xlog.last.event,
                                       halflife.time.between.events, 
                                       eventattributevar, eventattributevalue,
                                       placeholder,  "nomatch", "target" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "overTimes", 
                                             eventattributevalue, sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }
        
      }
      if ( is.null(eventtypevar) == FALSE & is.null(eventattributevar)==FALSE ){
        ## (3d) sender, 2times, match, filter
        result <- similarityComplexCpp(sender, target, time, xlog.last.event,
                                       halflife.time.between.events, 
                                       eventattributevar, eventattributevalue,
                                       eventtypevar,  "match", "target" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "overTimes", 
                                             eventattributevalue, "sameType",
                                             sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }
      }
      if ( is.null(eventtypevar) == FALSE & is.null(eventattributevar) ){
        ## (4d) sender, 2times, match
        result <- similarityComplexCpp(sender, target, time, xlog.last.event, 
                                       halflife.time.between.events, 
                                       placeholder, "1", eventtypevar, "match", 
                                       "target" )		
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- paste(variablename, "overTimes", 
                                             "sameType", sep = ".")
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 degree variable that was generated
          return(result)
        }
      }
    }#closes if --two-halflife-parameters-set--command
    
  }##closes if senderOrTarget == "target"
  }

####################################################################
##  Create event sequence
####################################################################

create.event.sequence <- function(datevar, dateformat = NULL, data = NULL,
                                  type = "continuous", byTime = "1 day",
                                  excludeDate = NULL, excludeTypeOfDay = NULL,
                                  excludeYear = NULL, excludeFrom = NULL, 
                                  excludeTo = NULL, returnData = FALSE, 
                                  sortData = TRUE, ...){
  
  #### check if all the inputs are correct
  ## check if date and dateformat match => then create Date-object
  if (type == "continuous"){
    date <- as.Date(datevar, format = dateformat)
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
  
  ## check if excludeFrom is smaller than excludeTo
  if ( is.null(excludeFrom) == FALSE & is.null(excludeTo) == FALSE) {
    if ( as.Date(excludeFrom, 
                 format = dateformat) > as.Date(excludeTo,format = dateformat) ){
      stop("'excludeFrom' is smaller than 'excludeTo'.")
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
  
  #### create continuous event sequence
  if (type == "continuous"){   
    
    ## create artificial sequence from start to end
    sequence <- data.frame(seq(min(date), max(date), byTime))
    names(sequence) <- "date.sequence"
    
    ## erase whatever
    ## TODO: exclude only these, that do not delete entire data set! 
    if ( is.null(excludeDate) == FALSE){
      ## exlcude some of them
      for (i in excludeDate){
        sequence <- subset(sequence, sequence$date.sequence != i)
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
      sequence$erase <- ifelse((sequence$date.sequence < 
                                  as.Date(excludeFrom, format = dateformat) | 
                                  sequence$date.sequence > 
                                  as.Date(excludeTo, format = dateformat)), 1, 0)
      sequence <- subset(sequence, sequence$erase == 1)
    }
    
    ## give artificial sequence 1:length()
    sequence$event.sequence <- 1:length(sequence$date.sequence)
    
    ## match with datevar
    result <- sequence$event.sequence[match(date, sequence$date.sequence)]
    
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
      names(data)[length(data)] <- "event.seq.cont"
      if ( sortData == FALSE ){
        return(data)
      }else{ ## sorted:
        data <- data[order(data$event.seq.cont), ]
        return(data)
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
      names(data)[length(data)] <- "event.seq.ord"
      if ( sortData == FALSE ){
        return(data)
      }else{ ## sorted:
        data <- data[order(data$event.seq.ord), ]
        return(data)
      }
    }#closes if-else returnData    
  }#closes type == ordinal
}


