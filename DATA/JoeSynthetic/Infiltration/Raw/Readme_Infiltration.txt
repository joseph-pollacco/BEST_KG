Dear Joseph,
 
Please, find the data refined for short times. Now you have 40 points by runs. These curves were obtained by sampling the numerical result. It corresponds to curves with more than thousands points. Then I have sampled these data with a constant increment to avoid too many points, which explained the very high values of the first point. In this new version, I have combined previous samples with a better sampling of the first part of the curve.
 
This question of sampling is crucial because, in your inverting procedure, the number of points will represent the weight you put into any part of the infiltration curve. The weight between the beginning of the curve (transient state) and the final part of the curve (steady state) will play a role in your inversion quality and on the estimated parameters. You must not miss the transient state and put enough points in the first part of the curve.
 
The second dataset (a_joseph_2new)  corresponds to the first (a_josephh_new) except that I removed the very first point (obtained for t = 0.001 min). In some cases, the infiltrations rates at t = 0.001 min were tending towards infinity which produced a step for cumulative infiltration. That is not physically sound, so I prefer to remove the first point.
 
You can try with one or the other of the dataset and see if there is any difference. The second one could be better.
 
Best
 
Laurent
