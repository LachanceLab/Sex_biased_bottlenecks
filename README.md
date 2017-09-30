# Sex_biased_bottlenecks

This directory includes MATLAB code to generate Watterson's Theta and Pi statistics for bottlenecked populations.

Each bottleneck can involve a different number of males and females.
Note that a serial founder effect model is used (e.g. CHB simulations assume that a GIH bottleneck occured prior to the CHB bottleck, and JPT simulations assume that GIH and CHB bottlenecks has occured prior to the JPT bottleneck).

bnldev.m is used speed up binomial sampling during the genetic drift phase of each simulated generation 
http://www.mathworks.com/matlabcentral/fileexchange/42464-binomial-random-number-generator?focused=3789359&tab=function
