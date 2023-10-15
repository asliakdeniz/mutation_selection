# Evolution and the Ultimatum Game by Akdeniz and van Veelen (2023)
Simulation code and data used in Akdeniz and van Veelen (2023, Games and Economic Behavior). 

The abstract of the paper can be found below. The full paper can be found in https://www.sciencedirect.com/science/article/pii/S0899825623001173#ab0010. 

Abstract: In this paper we review, upgrade, and synthesize existing models from evolutionary game theory, all of which aim at explaining human behaviour in the ultimatum game. Our new and improved versions of Gale et al. (1995), Nowak et al. (2000), and Rand et al. (2013) avoid shortcomings that the original versions have, one of which is that the results in the first and the last are driven by bias in the mutations. We also compare the predictions of these three models with the existing experimental evidence by looking at properties of the distributions of minimal acceptable offers. We find that the observed distributions do not conform to the predictions from Gale et al. (1995), Rand et al. (2013), or any other model in which there is no fitness benefit to rejecting. This does not rule out commitment-based explanations, such as Nowak et al. (2000).

Simulation code:
Matlab code used for the ultimatum game simulations in "Evolution and the Ultimatum Game" by Akdeniz and van Veelen (2023). The model settings are described in Section 2 of the paper. The four Matlab files in this directory are listed and explained below.

Global_Co: Ultimatum game with global mutations

Local_Co: Ultimatum game with local and co-occurring mutations

Local_NotCo: Ultimatum game with local and independent mutations

Local_NotCo_Partial: Ultimatum game with local and independent mutations, and with (partial) observability

The code can be run directly on Matlab. 

Data:
Data.zip file contains the data used in "Evolution and the Ultimatum Game" by Akdeniz and van Veelen (2023). 

For the direct-response method, the data is from the ultimatum game experiments by: Andersen et al. (2011); Barmettler et al. (2012); Bornstein and Yaniv (1998); Cameron (1999); Carpenter et al., 2005a, Carpenter et al., 2005b; Croson (1996); Forsythe et al. (1994); Lightner et al. (2017); Ruffle (1998); Slonim and Roth (1998); Tomlin (2015).

For the strategy method, the data is from the ultimatum game experiments by: Bader et al. (2021); Bahry and Wilson (2006); Benndorf et al. (2017); Soo Hong Chew et al. (2013); Demiral and Mollerstrom (2020); Inaba et al. (2018); Keuschnigg et al. (2016); Peysakhovich et al. (2014). 
