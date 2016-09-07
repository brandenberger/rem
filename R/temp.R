################################################################################
##	Triads (one-mode statistic)
################################################################################

triadStatold <- function(data, time, sender, target, halflife, weight = NULL,
                         eventtypevar = NULL, eventtypevalues = NULL, 
                         eventattributevar = NULL, eventattributeAI = NULL,
                         eventattributeBI = NULL, eventattributeAB = NULL,
                         variablename = "triad", returnData = FALSE, 
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
      stop("'", time, "' is not sorted. Sort data frame according to the event 
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
    stop("'", as.name(weight), "' variable is not numeric.")
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
      stop("First value '", eventattributeAB, "' is not an element of '", 
           deparse(substitute(eventattributevar)) , "'.") 
    }
    if ( length(grep(eventtypevalues[2], eventtypevar)) == 0 ) {
      stop("Second value '", eventattributeAB, "' is not an element of '", 
           deparse(substitute(eventattributevar)) , "'.") 
    }
    }
  
  ## check if event-attribute inputs are available and correctly specified
  if ( is.null(eventattributevar) == FALSE ) {
    # check length of variable
    if ( length(sender) != length(eventattributevar) ){
      stop("'sender' and 'eventattributevar' are not of same length.")
    }
    # transform variable
    eventattributevar <- as.character(eventattributevar)
    if ( is.null(eventattributeAB) & is.null(eventattributeAI) & 
         is.null(eventattributeBI) ){
      stop("No 'eventattribute__' provided. Provide a string value by which the 
           events are filtered.", )
    }
    # check if eventattributevalue is part of the variable
    if ( is.null(eventattributeAB) == FALSE){
      if ( length(grep(eventattributeAB, eventattributevar)) == 0 ) {
        stop("Value '", eventattributeAB, "' is not an element of '", 
             deparse(substitute(eventattributevar)) , "'.") 
      }
    }
    if ( is.null(eventattributeAI) == FALSE){
      if ( length(grep(eventattributeAI, eventattributevar)) == 0 ) {
        stop("Value '", eventattributeAI, "' is not an element of '", 
             deparse(substitute(eventattributevar)) , "'.") 
      }
    }
    if ( is.null(eventattributeBI) == FALSE){
      if ( length(grep(eventattributeBI, eventattributevar)) == 0 ) {
        stop("Value '", eventattributeBI, "' is not an element of '", 
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
  xlog <- log(2)/halflife 
  
  ####### calculate stat
  ## create placeholder-variables to be used in the cpp-Function
  placeholder <- rep("1", length(time))
  
  ## calculate the triad effects for each event
  
  ## all the statistics without an event type
  if ( is.null(eventtypevar) ){
    ## all stats without an event type and an event attribute
    if ( is.null(eventattributevar) ){
      ## (1) no type, no attribute. Simple triad-effect
      result <- triadCpp(sender, target, time, weight, placeholder, "1", "1", 
                         placeholder, "1", placeholder, "1", placeholder, "1",
                         xlog )  	
      if ( returnData == TRUE ) {
        data <- cbind(data, result)
        names(data)[length(data)] <- variablename
        ## return the data frame with the variable bound to it
        return(data)
      }else{ 
        ## only return the 1 triad variable that was generated
        return(result)
      }      
    }else{
      ## all stats without event type but with event attribute
      ## (2) no type, attributeAB
      if ( is.null(eventattributeAI) & is.null(eventattributeBI) & is.null(eventattributeAB) == FALSE ){
        result <- triadCpp(sender, target, time, weight, placeholder, "1", "1", 
                           eventattributevar, eventattributeAB, placeholder, "1", 
                           placeholder, "1", xlog )    
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- variablename
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 triad variable that was generated
          return(result)
        }
      }
      ## (3) no type, attributeAI
      if ( is.null(eventattributeAB) & is.null(eventattributeBI) & is.null(eventattributeAI) == FALSE ){
        result <- triadCpp(sender, target, time, weight, placeholder, "1", "1", 
                           placeholder, "1", eventattributevar, eventattributeAI,
                           placeholder, "1", xlog )    
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- variablename
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 triad variable that was generated
          return(result)
        } 
      }
      ## (4) no type, attributeBI
      if ( is.null(eventattributeAB) & is.null(eventattributeAI) & is.null(eventattributeBI) == FALSE ){
        result <- triadCpp(sender, target, time, weight, placeholder, "1", "1", 
                           placeholder, "1", placeholder, "1", eventattributevar,
                           eventattributeBI, xlog )    
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- variablename
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 triad variable that was generated
          return(result)
        } 
      }
      ## (5) no type, attributeAB & attributeAI
      if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAI) == FALSE & is.null(eventattributeBI) ){
        result <- triadCpp(sender, target, time, weight, placeholder, "1", "1", 
                           eventattributevar, eventattributeAB, eventattributevar, 
                           eventattributeAI, placeholder, "1", xlog )    
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- variablename
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 triad variable that was generated
          return(result)
        }      
      }
      ## (6) no type, attribute AB & attributeBI
      if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAI)  & is.null(eventattributeBI) == FALSE ){
        result <- triadCpp(sender, target, time, weight, placeholder, "1", "1", 
                           eventattributevar, eventattributeAB, placeholder, 
                           "1", eventattributevar, eventattributeBI, xlog )    
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- variablename
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 triad variable that was generated
          return(result)
        }
      }
      ## (7) no type, attribute AI & attributeBI
      if ( is.null(eventattributeAB) & is.null(eventattributeAI) == FALSE & is.null(eventattributeBI) == FALSE ){
        result <- triadCpp(sender, target, time, weight, placeholder, "1", "1", 
                           placeholder, "1", eventattributevar, 
                           eventattributeAI, eventattributevar, eventattributeBI, xlog )    
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- variablename
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 triad variable that was generated
          return(result)
        }
      }
      ## (8) no type, attribute AB & attributeAI & attributeBI
      if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAI) == FALSE & is.null(eventattributeBI) == FALSE ){
        result <- triadCpp(sender, target, time, weight, placeholder, "1", "1", 
                           eventattributevar, eventattributeAB, eventattributevar, 
                           eventattributeAI, eventattributevar, eventattributeBI, xlog )    
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- variablename
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 triad variable that was generated
          return(result)
        }
      }
      
    }#closes else attributevar != null   
  }else{
    ## all the statistics with an event type
    if ( is.null(eventattributevar) ){
      ## with type, but no attribute
      ## (9) type, no attribute
      if ( is.null(eventattributeAB)  & is.null(eventattributeAI)  & is.null(eventattributeBI)  ){
        result <- triadCpp(sender, target, time, weight, eventtypevar, eventtypevalues[1], eventtypevalues[2], 
                           placeholder, "1", placeholder, 
                           "1", placeholder, "1", xlog )    
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- variablename
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 triad variable that was generated
          return(result)
        }
      }
      
    }else{
      ## all stats with type and attribute
      ## (10) type, attributeAB
      if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAI)  & is.null(eventattributeBI)  ){
        result <- triadCpp(sender, target, time, weight, eventtypevar, eventtypevalues[1], eventtypevalues[2], 
                           eventattributevar, eventattributeAB, placeholder, 
                           "1", placeholder, "1", xlog )    
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- variablename
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 triad variable that was generated
          return(result)
        }
      }
      
      ## (11) type, attributeAI
      if ( is.null(eventattributeAB)  & is.null(eventattributeAI) == FALSE & is.null(eventattributeBI)  ){
        result <- triadCpp(sender, target, time, weight, eventtypevar, eventtypevalues[1], eventtypevalues[2], 
                           placeholder, "1", eventattributevar, 
                           eventattributeAI, placeholder, "1", xlog )    
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- variablename
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 triad variable that was generated
          return(result)
        }
      }
      
      ## (12) type, attributeBI
      if ( is.null(eventattributeAB)  & is.null(eventattributeAI)  & is.null(eventattributeBI) == FALSE ){
        result <- triadCpp(sender, target, time, weight, eventtypevar, eventtypevalues[1], eventtypevalues[2], 
                           placeholder, "1", placeholder, 
                           "1", eventattributevar, eventattributeBI, xlog )    
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- variablename
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 triad variable that was generated
          return(result)
        }
      }
      
      ## (13) type, attributeAB & attributeAI
      if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAI) == FALSE & is.null(eventattributeBI) ){
        result <- triadCpp(sender, target, time, weight, eventtypevar, eventtypevalues[1], eventtypevalues[2], 
                           eventattributevar, eventattributeAB, eventattributevar, 
                           eventattributeAI, placeholder, "1", xlog )    
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- variablename
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 triad variable that was generated
          return(result)
        }
      }
      
      ## (14) type, attribute AB & attributeBI
      if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAI)  & is.null(eventattributeBI) == FALSE ){
        result <- triadCpp(sender, target, time, weight, eventtypevar, eventtypevalues[1], eventtypevalues[2], 
                           eventattributevar, eventattributeAB, placeholder, 
                           "1", eventattributevar, eventattributeBI, xlog )    
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- variablename
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 triad variable that was generated
          return(result)
        }
      }
      
      ## (15) type, attribute AI & attributeBI
      if ( is.null(eventattributeAB)  & is.null(eventattributeAI) == FALSE & is.null(eventattributeBI) == FALSE ){
        result <- triadCpp(sender, target, time, weight, eventtypevar, eventtypevalues[1], eventtypevalues[2], 
                           placeholder, "1", eventattributevar, 
                           eventattributeAI, eventattributevar, eventattributeBI, xlog )    
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- variablename
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 triad variable that was generated
          return(result)
        }
      }
      
      ## (16) type, attribute AB & attributeAI & attributeBI
      if ( is.null(eventattributeAB) == FALSE & is.null(eventattributeAI) == FALSE & is.null(eventattributeBI) == FALSE ){
        result <- triadCpp(sender, target, time, weight, eventtypevar, eventtypevalues[1], eventtypevalues[2], 
                           eventattributevar, eventattributeAB, eventattributevar, 
                           eventattributeAI, eventattributevar, eventattributeBI, xlog )    
        if ( returnData == TRUE ) {
          data <- cbind(data, result)
          names(data)[length(data)] <- variablename
          ## return the data frame with the variable bound to it
          return(data)
        }else{ 
          ## only return the 1 triad variable that was generated
          return(result)
        }
      }
      
    }##closes else attr-var != null
  }## closes else-type-var != null
  
  }#closing

