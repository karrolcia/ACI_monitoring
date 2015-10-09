# ACI_monitoring
Aerosol - Cloud Interactions Monitoring Software
ACI Monitoring Scheme Data Processing program
   
This program filters data in accordance with the data selesction criteria of the ACI Monitoring Scheme. 
The criteria are as follows:

 1. ACI monitoring scheme applies only to liquid water clouds on top of 
    the boundary layer in well-mixed conditions. This is doneto ensures 
    that the cloud is not decoupled from the boundary layer and the 
    aerosol background below the cloud (Feingold et al., 2006).
 2. Any form of precipitation, including drizzle, needs to be eliminated 
    from the profile as it can obscure the formative stage of a cloud 
    (Feingold et al., 2003). We use the Cloudnet categorisation data for
    the classification of the observed targets. 
 3. The assumption of ACI is that the variation in the aerosol concentration
    affects the variation in the cloud properties. Thus, both aerosol 
    and cloud parameters need to vary to observe ACI. 
 4. The ACI monitoring scheme relies on measurements from three separate 
    instruments. Only profiles where all three instruments provide good 
    quality data are analysed.
 5. We apply a constraint on LWP. In most cases we divide the data into 
    bins of LWP of 10 gm-2.  LWP % should be above 30 gm-2 and below 
    150 gm-2. Values below 30 gm-2 are disregarded because of the 
    uncertainty of LWP calculated from MWR, which is around 15 gm-2 
    (Turner et al., 2007). The values above 150 gm-2 170 are excluded to 
    avoid precipitating clouds.

After applying filtering data is further prepared for the analysis. This preparation includes converting Radar 
Reflectivity Factor to a linear values, integrating the Attenuated Backscatter Coefficient and Linear 
Radar Reflectivity Factor. If the depolarisation data is available, the Attenuated Backscatter Coefficient is 
divided between the parallel and perpendicular channel to detect the cloud base.

Additional part of this program is the ACI_scheme_plots.m script which creates time series and scatter plots 
and performs the statistical analysis of the data.
 
version 1.2
*AUTHOR*
Karolina Sarna
*DATE*
2015-07-16
 
k.sarna@tudelft.nl
karolina.sarna@gmail.com
