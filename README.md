# -Facilitation-buffers-a-Pinus-edulis-population
data and code to estimate the effect of facilitation via nurse plants on a Pinus edulis (piñon pine) population using experimental seedling survival data and an integral projection model from a 10-year census.

All data needed is in /data.

The idea was to take this IPM for piñon, and feed my seedling data into it. Because of models and observations indicating increasing drought frequency for the southwest US, I want to evaluate if survival differences we see based on WHERE seedlings are planted (under tree canopies or not) buffer the whole population from decline, but ask whether this matters more or less as drought frequency increases. 

First, I used SPEI-36 with a 1 year lag to determine if a year was in "drought". This is because, in my mind, the long-run effects of prolonged reductions in precipitation and high potential evapotranspiration is really what drives piñon mortality. An SPEI-36 lag1 < or equal to -1 indicated a "drought year". Doing it this way shows that 6/10 years in the census period were considered drought. I obtained SPEI values from SPEIbase.

Next, since we cannot infer recruitment from this data, and do not have recruitment data in general, i attempted to define a multiplier that cone counts will be multiplied by to quantify recruitment (in this case, seedlings entering the first size class of the model, and thus population). I defined this recruitment multiplier as R = seeds per cone * germination rate of seeds * first two-year survival probability (from my data) * probability of progressing to the minimum size class. I found some literature references for seeds per cone and germination rate, but had to make some assumptions, MUCH MORE INFO outlined in the 01_facilitation_lambda.R script. 

Then, using a sequential projection, I attempted to determine 1) the lamdba of the population based on whether seedling come up underneath tree canopies or out in the open, and 2) how that lambda differs based on differing frequencies of drought years. 

This essentially flips a weighted coin, where the probability of selecting a drought year is equal to specified drought frequency (if 20% drought frequency, 20% chance of encountering a drought year). What this is "pulling" at random based on that frequency is the vital rates of the population from drought or non-drought years based on SPEI-36 lag1. The R multiplier stays the same, but the vital rates change based on what years are randomly pulled. 
