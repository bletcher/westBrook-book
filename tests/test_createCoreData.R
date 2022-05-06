
d0 <- 
  createCoreData(
    sampleType = "electrofishing", #"stationaryAntenna","portableAntenna"),
    whichDrainage = "west",
    columnsToAdd =
      c("sampleNumber",
        "river",
        "riverMeter",
        "survey",
        "pass",
        'observedLength',
        'observedWeight')
  ) 

d1 <- d0 %>%
  addTagProperties(
    columnsToAdd = 
      c("cohort",
        "species",
        "dateEmigrated",
        "sex",
        "species")
  )

d2 <- d1 %>%
  dplyr::filter(!is.na(tag), 
                area %in% c("trib","inside","below","above"), 
                !is.na(sampleNumber) 
  ) 

d3 <- d2 %>%
  createCmrData(maxAgeInSamples = maxAgeInSamples, 
                inside = F, 
                censorDead = F, 
                censorEmigrated = F)  # may want to change censorEmigrated = T to = F
  # sample 83 is the last tagging sample

d4 <- d3 %>% 
  filter(sampleNumber <= 83)

d5 <- d4 %>%
  addSampleProperties() 

d6 <- d5 %>%
  addEnvironmental() 

  # these functions do not work for CMR data - they separate out shock and non-shock samples
  #addEnvironmentalDaily() %>%
  #addEnvironmentalInterval() %>%
 
d7 <- d6 %>% 
 addKnownZ2() 

d8 <- d7 %>%
  addFirstLast() 

d9 <- d8 %>%
  fillRiver()

d10 <- d9 %>%
  scaleEnvData()

d8 %>% 
  select(tag, cohort, detectionDate, lagDetectionDate, season, year, sampleNumber, survey, ageInSamples, river, enc, meanFlow, 
         meanTemperature, knownZ) %>% 
  data.frame() %>% 
  head(500)

d8 %>% 
  select(tag, cohort, detectionDate, season, year, sampleNumber, survey, ageInSamples, river, enc, meanFlow, 
         meanTemperature, knownZ, firstObserved, lastObserved, isFirstObserved, isLastObserved) %>% 
  data.frame() %>% 
  head(500)
