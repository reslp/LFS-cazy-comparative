#!~/bin/cafe
version
date
tree %%tree%%


#analysis 1:
load -i %%datafile%% -l %%prefix%%_cafe.log -filter -t 8 -p 0.05 -r 1000
lambdamu -s -t %%cafetree%%

report %%prefix%%_cafe_results






