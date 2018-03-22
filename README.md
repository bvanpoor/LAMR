# LAMR
Length-based Age-structured Mark-recapture Estimation (also known as LAME)

This estimation model is designed to use seasonal (or annual) growth of animals (we use fish) to help remove bias that occurs when capture probability varies with size. If capture probability increases with size and size-at-age varies among individuals (which is almost universal), fast growing animals will be marked first and mean capture probability of marked fish will be higher than unmarked fish: a violation of mark-recapture. This method estimates variation in length at age and predicts the differences in capture probability between marked and unmarked fish, thereby removing this bias. 

This paper is being reviewed in Ecological Modelling. The attached code is that used in the manuscript and the datafiles are based on the rainbow trout populations used as a case study there. 
