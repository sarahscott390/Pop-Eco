#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                   ##  
## Here we import cleaned data under a robust design for occurrence  ##
#  observations of Piute ground squirrels at the NCA. We run a      ##
## dynamic occupancy analysis. See Mackenzie et al. (2003) Ecology   ##
## for details of the model. The dynamic occupancy model is hierarchical #
# with two components: (1) an ecological submodel that describes how ##
## site colonization and extinction processes explain site occupancy  ##
## with the opportunity to link these to environmental predictors at ##
## the site. (2) describes the observation submodel linking detection ##            
## probability to relevant predictors.                                ##
##                                                                   ##
#######################################################################
##### Set up your workspace and load relevant packages -----------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()

# Install new packages from "CRAN" repository. # 

# load packages:
library( tidyverse )#includes dplyr, tidyr and ggplot2
library( unmarked ) #run occupancy and abundance models
library( AICcmodavg ) #model evaluation 
## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are:
datadir <- paste( getwd(), "/Data/", sep = "" )
# load cleaned data with robust format
robdf <- read.csv( file = paste( datadir, "opendf.csv", sep = ""),
                   header = TRUE )
#view
head( robdf ); dim( robdf ) 
#### End of data load -------------
####################################################################
##### Ready data for analysis ----------------------------
# What predictors do we think drive colonization, extinction and # 
# detection of Piute ground squirrels at the NCA? #
# We think that cheatgrass and sagebrush are both potentially        #
# important to occupancy. Breaking those relationships down into the #
# components of our dynamic process we predict that:
# Extinction:
# Cheatgrass is an invasive species that increases the likelihood of # 
# more frequent fires in the system so it may increase the probability #
# of a site becoming extinct. We also assess effects of min temperatures # 
# in Feb, assuming that years with colder Feb temperatures lower survival of  #
# adults and reproduction if it limits food availability when ground squirrels #
# come out of hibernation. Particularly sites with low abundances, may become #
# extinct as a result. 
# Colonization:
# Sites with more sagebrush are more likely to become colonized as sagebrush #
# provides food, and cover from weather and predators. Years with cooler summers#
# may increase survival of young, and therefore the amount of individuals, #
# spreading to new areas.

# Detection:
# We expect observer effects influence detection. We also test the #
# effects of sagebrush again. #

#collapse wide presences and observer ids into long format:
longdf <- robdf %>% pivot_longer( cols = contains("j"), 
                                  names_to = c(".value", "survey"),
                                  names_pattern = "(.+).(.+)" )

#check
head( longdf );dim(longdf )
#do we have the right number of rows?
#answer: yes, we went from 1 to 3 rows

# Arrange columns in the correct order for unmarked fomatMult() function:
longdf <- longdf %>% select( year, o.sites, survey, pres =  pres., 
                             cheatgrass, sagebrush, Feb.minT, AprMay.maxT,
                             observer = observer. ) 
# Use formatMulti to get data in correct format for robust analysis:
umf <- formatMult( longdf )
#view
summary( umf )

#now scale yeasrXsite level predictors:
sc <- apply( yearlySiteCovs(umf), MARGIN = 2, FUN = scale )
# We replace the predictors in our unmarked dataframe with the scaled values:
yearlySiteCovs( umf ) <- sc

### end data prep -----------
### Analyze data ------------------------------------------
# We are now ready to perform our analysis. Since the number of predictors #
# is reasonable for the sample size, and there were no issues with #
# correlation, we focus on a single full, additive model:
fm.dyn <- colext( #define detection submodel:
  pformula = ~ 1 + observer + sagebrush,
  #define occupancy submodel for year 1:
  psiformula = ~ 1,
  #define extinction submodel for years 2:T:
  epsilonformula = ~ 1 +  cheatgrass + Feb.minT, 
  #define colonization submodel for years 2:T:
  gammaformula = ~ 1 + sagebrush + AprMay.maxT, 
  #data to use:
  data = umf )
# Why didn't we put predictors on psi?
# Answer: So we can determine baseline occupancy probability? Need an intercept-only model.
#
# View model results:
fm.dyn

# Estimate confidence intervals for coefficients in ecological submodels:
plogis(confint( fm.dyn, type = "psi" ))
# what does that plogis function do? What does the result tell us?
# Answer:  Plogis applies the logistic transformation to a set of values.  Confint() estimates confidence intervals for the occupancy coefficients (psi) in the model and plogis() converts the estimated log-odds back to probabilities between 0 and 1.  
#
# Why don't we do the same for colonization and extinction?
# Answer: Colonization (gamma) and extinction (epsilon) represent rates of change, not occupancy probabilities. Plogis() would not return meaningful results.
#
confint( fm.dyn, type = "col" )
confint( fm.dyn, type = "ext" )
# How do we interpret these results?
## Answer: Sagebrush has a strong positive effect on the probability of colonization (CI: 0.8704710  1.9020003; entirely above zero).  Cheatgrass has a strong positive effect on the probability of extinction (CI: 2.4703049  4.60450386; entirely above zero), so a negative impact on ground squirrel survival.  Temperature does not seem to have much of an effect on colonization or extinction (both CIs cross zero).
#
# Estimate confidence intervals for coefficients for detection submodel:
confint( fm.dyn, type = 'det' )
# How do we interpret these results?
# Answer: Sagebrush's CI crossed zero, indicating that sagebrush does not have an effect on the probability of detection.  Observer 2's CI is entirely below zero, indicating that their presence decreases the probability of detection compared to Observer 1.  Observer 3 and Observer 4's CI is entirely below zero, indicating that their presence decreases the probability of detection compared to Observer 1.  They also seem to have an even lower detection probability than Observer 2, as their sub-zero values are even more negative than Observer 2.  
#
################ How are they different to our single-season detection model?
###### Answer:  The great-horned owl one?  the output was in Bayesian in the eval script I believe, not sure how to compare the two
#

##########################################################################
# Model fit and evaluation -----------------------------------------------

# Now that we looked at the initial output we can evaluate our model to #
# decide if we are happy to proceed or need to modify our analysis #
# somehow.

# We start with goodness of fit (GoF) on detection frequencies, which relies on a #
# Pearson chi-square to assess fit as suggested by Mackenzie and Bailey (2004) #
# J. Agr., Bio. & Env. Stats. 9: 300-318. 
# This test is extended in AICmodavg package to dynamic occupancy models of #
# MacKenzie et al. (2003) by using the occupancy estimates for each season obtained #
# from the model. These estimates are then used to compute the predicted and #
# observed frequencies separately within each season. The chi-squares are then #
# summed to be used as the test statistic for the dynamic occupancy model.
gof.boot <- AICcmodavg::mb.gof.test( fm.dyn, nsim = 1000, print.table = TRUE )
#view
gof.boot
# What does the output tell us about our model fit?
## Answer:  There is variation between seasons in goodness of fit (ex. season 10 fit is a bit worse than other seasons at 11.9 and Season 4 is has the best fit of all years at 1.0).  Our p-value of 0.517 suggests that there is no significant evidence to reject the null hypothesis of good fit.  The bootstrapped histogram looks solid (I think?  the median line falls in the visible center of the values).  Our c-hat is close to 1, so I dont think our model is over or under dispersed in a problematic way.  
#
#
# If we want to look at each season to see if any of them had particularly bad fit:
gof.boot$chisq.table$tables
# Is there a season that was particularly bad? Which?
## Answer: Yes, season 10 
# Remember that higher chi-squared values represent worse fit

# We also evaluate how well our full model did against the null model # 
# by estimating pseudo-R^2, based on Nagelkerke, N.J.D. (2004) A Note #
# on a General Definition of the Coefficient of Determination. Biometrika 78,#
# pp. 691-692.#
# We run the null model
fm.null <- colext( #define detection submodel:
  pformula = ~ 1,
  #define occupancy submodel for year 1:
  psiformula = ~ 1,
  #define extinction submodel for years 2:T:
  epsilonformula = ~ 1,
  #define colonization submodel for years 2:T:
  gammaformula = ~ 1,
  #data to use:
  data = umf )
#view
fm.null
# Now build the list with the two models of interest:
rms <- fitList( 'full' = fm.dyn,
                'null' = fm.null )
# Then use model selection function from unmarked, defining which is the null:
unmarked::modSel(rms, nullmod = "null" )

#########################################################################
##### Summarizing model output ##############
# Now we plot the results we are interested in.
# We can see how mean occupancy changed through time by extracting values from
# the projected.mean data table:
fm.dyn@projected.mean
# select occupied row and combine it with year
data.frame( year = sort( unique( robdf$year) ),
            occupied = as.vector( fm.dyn@projected.mean["occupied",] ) ) %>%
  ggplot(., aes( x = year, y = occupied ) ) +
  theme_bw(base_size = 15 ) + 
  geom_line( size = 2 )

# What is happening to Piute ground squirrels at the NCA?
## Answer: Occupancy is declining consistently and dramatically
# 
# We now see the effects that our predictors are having on this trend. #
# by plotting partial prediction plots for our ecological submodels #
# Here I focus only on those with 95% CIs not overlapping zero:
# We start by creating our datasets to predict over
# how many values do we use:
n <- 100
# we use the observed values to define our range:
sagebrush <- seq( min( robdf[,"sagebrush"]),max( robdf[,"sagebrush"]),
                  length.out = n )
cheatgrass <- seq( min( robdf[,"cheatgrass"]),max( robdf[,"cheatgrass"]),
                   length.out = n )
#standardize them
sage.std <- scale( sagebrush )
cheat.std <- scale( cheatgrass )
#combine standardized predictor into a new dataframe to predict partial relationship
# with sagebrush. We replace value of other predictor with its mean
colData <- data.frame( sagebrush = sage.std, AprMay.maxT = 0 )

#predict partial relationship between sagebrush and occupancy
pred.col.sage <- predict( fm.dyn, type = "col", newdata = colData, 
                          appendData = TRUE )
#view
head( pred.col.sage ); dim( pred.col.sage )

#combine standardized predictor into a new dataframe to predict partial relationship
# with cheatgrass. We replace value of other predictor with its mean
extData <- data.frame( Feb.minT = 0, cheatgrass = cheat.std )
#predict partial relationship between sagebrush and occupancy
pred.ext.cheat <- predict( fm.dyn, type = "ext", newdata = extData, 
                           appendData = TRUE )

# create plots for ecological submodels
#Starting with sagebrush and colonization:
# select the predicted values we want to plot and combine with unscaled predictor
sagep <- cbind( pred.col.sage[,c("Predicted", "lower", "upper") ], sagebrush ) %>%
  # define x and y values
  ggplot(., aes( x = sagebrush, y = Predicted ) ) + 
  #choose preset look
  theme_bw( base_size = 15 ) +
  # add labels
  labs( x = "Sagebrush (%)", y = "Estimated colonization" ) +
  # add band of confidence intervals
  geom_smooth( aes(ymin = lower, ymax = upper ), 
               stat = "identity",
               size = 1.5, alpha = 0.5, color = "grey" ) +
  # add mean line on top
  geom_line( size = 2 ) 
#view
sagep
# How do you interpret this relationship?
# Is there a potential threshold beyond which colonization becomes unlikely?
## Answer:  Sagebrush has a positive effect on the probability of colonization.  A potential threshold beyond which colonization becomes unlikely looks to to be ~15%.
#
#Now plot cheatgrass and extinction:
# select the predicted values we want to plot and combine with unscaled predictor
cheatp <- cbind( pred.ext.cheat[,c("Predicted", "lower", "upper") ], cheatgrass ) %>%
  # define x and y values
  ggplot(., aes( x = cheatgrass, y = Predicted ) ) + 
  #choose preset look
  theme_bw( base_size = 15 ) +
  # add labels
  labs( x = "Cheatgrass (%)", y = "Estimnate extinction" ) +
  # add band of confidence intervals
  geom_smooth( aes(ymin = lower, ymax = upper ), 
               stat = "identity",
               size = 1.5, alpha = 0.5, color = "grey" ) +
  # add mean line on top
  geom_line( size = 2 ) 
#view
cheatp
# How would you interpret this relationship?
## Answer: Cheatgrass has a strong positive effect on the probability of extinction.  The slope of this relationship is steep and the effect begins at a low % cheatgrass.  
# 

#### On your own, plot partial predictions for your detection submodel. 

# Prepare the data for predictions
# Create a new data frame for each observer
observer_levels <- unique(longdf$observer)
n <- 100  # Number of values for the prediction

# Create an empty list to store predictions
predictions <- list()

# Loop through each observer to get predictions
for (obs in observer_levels) {
  # Create a data frame for each observer with mean values for other predictors
  new_data <- data.frame(
    observer = obs,
    sagebrush = seq(min(robdf$sagebrush), max(robdf$sagebrush), length.out = n)
  )
  
  # Predict detection probabilities
  pred <- predict(fm.dyn, type = "det", newdata = new_data, appendData = TRUE)
  predictions[[obs]] <- data.frame(
    sagebrush = new_data$sagebrush,
    Predicted = pred$Predicted,
    lower = pred$lower,
    upper = pred$upper,
    observer = obs
  )
}

# Combine all predictions into a single data frame
all_predictions <- bind_rows(predictions)


# Which partial prediction plots are you focusing on and why?
## Answer: Sagebrush percentage vs detection probability and differences in detection probability between observers

# Plot:

ggplot(all_predictions, aes(x = sagebrush, y = Predicted, color = observer)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  labs(x = "Sagebrush (%)", y = "Predicted Detection Probability") +
  theme_minimal() +
  ggtitle("Partial Predictions of Detection Probability by Observer")

#
######
############################################################################
################## Save your data and workspace ###################

# Save workspace:
save.image( "RobOccResults.RData" )

#save the plot objects you need for your presentation
#start by calling the file where you will save it
tiff( 'Data/SageXCol.tiff',
      height = 10, width = 12, units = 'cm', 
      compression = "lzw", res = 400 )
#call the plot
sagep
#end connection
dev.off()
#Now the cheatgrass x occupancy plot:
tiff( 'Data/CheatXExt.tiff',
      height = 10, width = 12, units = 'cm', 
      compression = "lzw", res = 400 )
#call the plot
cheatp
#end connection
dev.off()

########## End of saving section ##################################

############# END OF SCRIPT #####################################