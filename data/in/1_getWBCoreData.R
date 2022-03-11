library(getWBData)
coreData<-createCoreData(sampleType="electrofishing",
                         columnsToAdd=c("sampleNumber","river",
                                        'observedLength','observedWeight')) %>% 
  addTagProperties() %>%
  dplyr::filter( species %in% c( "bkt","bnt","ats" ) ) %>%
#  createCmrData( maxAgeInSamples=20 ) %>%
  addSampleProperties() %>%
  addEnvironmental()

#  addKnownZ()


