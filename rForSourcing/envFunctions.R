# functinos that haven't made it into getWBData yet...
# file location: projects/westbrook/code/westbrookKata/getWBData/R in osensei

#'Add environmental data to core westbrook data for intervals for survey == 'shock'. Use addEnvironmentalDaily to add daily environmental data to reader observations
#'@return Environmental data added to a coreData
#'@param coreData a data.frame created using createCoreData
#'@param sampleFlow Logical indicating whether to add the discharge on detectionDate to the output
#'@param fun A quoted name of a function used to characterize the environmental data between samples (e.g., mean, median, min, max). Custom functions can be used but must have two arguments: x and na.rm.
#'@details View this function as a template to be altered to get different env variables or different summaries. / This function add columns of mean daily temperature and mean daily flow to coreData. For each individual, it returns the mean env values between capture dates or between median capture dates of all individuals of a fish was not captured on an occasion.
#' The internal function getIntervalMean returns means for the observed river or for all rivers if an individual was unobserved.
#' Could alter getIntervalMean to return data for the mean daily min or max env data or for data from the last observed river, instead of all rivers.
#'
#'@export 

addEnvironmentalInterval <-function( coreData, sampleFlow=F , funName="mean"){
  func<-get(funName)
  # Only want to calculate stats on intervals for 'shock' samples
  # Env data for closest env observation and the surrounding 24 hrs are calculated in addEnvrionmentalDaily
  coreDataNotShock <- coreData %>% filter(survey != "shock")
  coreData<- coreData %>% filter(survey == "shock")
  # get env data from database
  whichDrainage<-"west"
  if(all(!unique(coreData$river) %in% c("west brook","wb jimmy","wb mitchell","wb obear"))){
    whichDrainage<-"stanley"
  }
  if(whichDrainage=="west"){
    envData<-tbl(conDplyr,"data_daily_temperature") %>%
      collect(n=Inf) %>%
      full_join(tbl(conDplyr,"data_flow_extension") %>% collect(n=Inf),
                by=c("river","date")) %>%
      select(-source) %>%
      dplyr::filter(date <= max(coreData$detectionDate),
                    date >= min(coreData$detectionDate)) %>%
      rename(temperature = daily_mean_temp,
             flow = qPredicted) %>%
      data.frame()
  } else {
    envData<-tbl(conDplyr,"stanley_environmental") %>%
      filter(section==11) %>%
      select(datetime,temperature,depth) %>%
      collect(n=Inf) %>%
      rename(flow = depth,date = datetime) %>%
      data.frame()
    warning("Depth was inserted into flow column because that is what is available in Stanley")
  }
  ###########################################################################
  # set up table of intervals based on date for each capture interval by fish
  coreData <- coreData %>%
    group_by( tag ) %>%
    mutate( lagDetectionDate = lead( detectionDate ) ) %>%
    ungroup()
  if(whichDrainage=="west"){ #splitting by drainage here because river isn't required for Stanley, not a good way to achieve this :/
    # function to get mean environmental data for each row of coreData
    getIntervalMean <- function( start,end,r,e ,fun=func){
      d <- envData$date 
      if( e == "Temperature" ) {
        envCol <- "temperature"
        if ( is.na( r ) ) meanTemp <- fun( envData[ d >= start & d <= end , envCol ], na.rm = T ) #use data from all rivers
        if ( !is.na( r ) ) meanTemp <- fun( envData[ d >= start & d <= end & envData$river == r, envCol ], na.rm = T )
      }
      if( e == "Flow" ) {
        envCol <- "flow"
        meanTemp <- fun( envData[ d >= start & d <= end , envCol ], na.rm = T )
      }
      return( meanTemp )
    }
    # get unique start and end dates from coreData and calc env means for the intervals
    coreDataUniqueDates <- coreData %>%
      select( river, detectionDate, lagDetectionDate ) %>%
      unique() %>%
      group_by( river, detectionDate, lagDetectionDate ) %>% #just a loop, probably a better way to do this.
      mutate( meanTemperature = getIntervalMean( detectionDate, lagDetectionDate, river, "Temperature" ),
              meanFlow = getIntervalMean( detectionDate, lagDetectionDate, river, "Flow" )) %>%
      ungroup()
    coreData <- left_join( coreData,coreDataUniqueDates,
                           by=c("detectionDate","river","lagDetectionDate"))
  } else {
    # function to get mean environmental data for each row of coreData
    getIntervalMean <- function( start,end,e ,fun=func){
      d <- envData$date
      meanEnv <- fun( envData[ d >= start & d <= end , tolower(e) ], na.rm = T )
      return( meanEnv )
    }
    # get unique start and end dates from coreData and calc env means for the intervals
    coreDataUniqueDates <- coreData %>%
      select( detectionDate, lagDetectionDate ) %>%
      unique() %>%
      group_by( detectionDate, lagDetectionDate ) %>% #just a loop, probably a better way to do this.
      mutate( meanTemperature = getIntervalMean( detectionDate, lagDetectionDate, "Temperature" ),
              meanFlow = getIntervalMean( detectionDate, lagDetectionDate, "Flow" )) %>%
      ungroup()
    coreData <- left_join( coreData,coreDataUniqueDates,
                           by=c("detectionDate","lagDetectionDate"))
  }
  if(sampleFlow){
    coreData<-coreData %>%
      mutate(date=as.Date(detectionDate)) %>%
      filter(enc==1) %>%
      select(sampleName,date) %>%
      group_by(sampleName,date) %>%
      summarize(n=n()) %>%
      ungroup() %>%
      left_join(envData %>%
                  filter(!is.na(flow)) %>%
                  mutate(date=as.Date(date)) %>%
                  select(date,flow) %>%
                  rename(flowForP=flow) %>%
                  unique(),
                by=c("date")) %>%
      group_by(sampleName) %>%
      summarize(flowForP=sum(flowForP*n)/(sum(n))) %>%
      ungroup() %>%
      right_join(coreData,by="sampleName")
    # coreData <- envData %>%
    # filter(!is.na(flow)) %>%
    # select(date,flow) %>%
    # rename(flowForP=flow) %>%
    # unique() %>%
    # right_join (coreData,by=c("date"="detectionDate")) %>%
    # rename(detectionDate=date)
  }
  names(coreData)[which(names(coreData)=="meanTemperature")]<-paste0(funName,"Temperature")
  names(coreData)[which(names(coreData)=="meanFlow")]<-paste0(funName,"Flow")
  coreData <- bind_rows( coreData, coreDataNotShock ) %>% arrange( tag, detectionDate )
  return( coreData )
}

#'Add environmental data to core westbrook data. Just for reader observations. Use addEnvironmentalInterval to add env data between sample occasions.
#'@return Environmental data added to a coreData
#'@param coreData a data.frame created using createCoreData
#'@param sampleFlow Logical indicating whether to add the discharge on detectionDate to the output
#'@param fun A quoted name of a function used to characterize the environmental data between samples (e.g., mean, median, min, max). Custom functions can be used but must have two arguments: x and na.rm.
#'@details View this function as a template to be altered to get different env variables or different summaries. / This function add columns of mean daily temperature and mean daily flow to coreData. For each individual, it returns the mean env values between capture dates or between median capture dates of all individuals of a fish was not captured on an occasion.
#' The internal function getIntervalMean returns means for the observed river or for all rivers if an individual was unobserved.
#' Could alter getIntervalMean to return data for the mean daily min or max env data or for data from the last observed river, instead of all rivers.
#'
#'@export
addEnvironmentalDaily <-function( coreData, sampleFlow=F , funName="mean"){
  func<-get(funName)
  # Only want to calculate stats on intervals for 'reader' samples
  coreDataShock <- coreData %>% filter(survey == "shock")
  coreData<- coreData %>% filter(survey != "shock")
  # get env data from database
  whichDrainage<-"west"
  if(all(!unique(coreData$river) %in% c("west brook","wb jimmy","wb mitchell","wb obear"))){
    whichDrainage<-"stanley"
  }
  if(whichDrainage=="west"){
    envData<-tbl(conDplyr,"data_daily_temperature") %>%
      collect(n=Inf) %>%
      full_join(tbl(conDplyr,"data_flow_extension") %>% collect(n=Inf),
                by=c("river","date")) %>%
      select(-source) %>%
      dplyr::filter(date <= max(coreData$detectionDate),
                    date >= min(coreData$detectionDate)) %>%
      rename(temperature = daily_mean_temp,
             flow = qPredicted) %>%
      data.frame()
  } else {
    envData<-tbl(conDplyr,"stanley_environmental") %>%
      filter(section==11) %>%
      select(datetime,temperature,depth) %>%
      collect(n=Inf) %>%
      rename(flow = depth,date = datetime) %>%
      data.frame()
    warning("Depth was inserted into flow column because that is what is available in Stanley")
  }
  if( whichDrainage == 'stanley') {
    # get daily mean for each day and merge into each antenna observation
    e <- envData %>%
      mutate( date = date(date) ) %>%
      group_by( date ) %>%
      summarize( meanDailyTemp = mean( temperature, rm.na = T ),
                 meanDailyFlow = mean( flow, rm.na = T))
    coreData$date <- date(coreData$detectionDate)
    coreData <- left_join( coreData, e, by = 'date' )
    # rounds detectionDate to nearest 15 or 60 min to join tidal stage and temperature data based on time
    # coreData$time15 <- (as.POSIXct(round(as.numeric(coreData$detectionDate)/(15*60))*(15*60), origin='1970-01-01'))
    # coreData$time60 <- (as.POSIXct(round(as.numeric(coreData$detectionDate)/(60*60))*(60*60), origin='1970-01-01'))
    # get daily closest obs and daily mean for each reader observation
    #coreData <- left_join( coreData, envData, by = c('time15' = "date"))
    # rounds detectionDate to nearest 15 or 60 min to join tidal stage and temperature data based on time
    coreData$time15 <- (as.POSIXct(round(as.numeric(coreData$detectionDate)/(15*60))*(15*60), origin='1970-01-01'))
    coreData$time30 <- (as.POSIXct(round(as.numeric(coreData$detectionDate)/(30*60))*(30*60), origin='1970-01-01'))
    coreData$time60 <- (as.POSIXct(round(as.numeric(coreData$detectionDate)/(60*60))*(60*60), origin='1970-01-01'))
    # get daily closest obs and daily mean for each reader observation
    coreData <- left_join( coreData, envData, by = c('time30' = "date"))
  }
  # put the data frames back together
  coreData <- bind_rows( coreData, coreDataShock ) %>% arrange( tag, detectionDate )
  return( coreData )
}

