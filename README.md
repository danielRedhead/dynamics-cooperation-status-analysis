Analysis scripts and data for "The dynamics of cooperation and social status in a small-scale society" 
============

Requirements: 

- R (3.5.2 or greater): https://cran.r-project.org
- RSiena (1.2-4 or greater): https://www.stats.ox.ac.uk/~snijders/siena/
- Statnet (2018.10 or greater):https://statnet.org


Instructions:

In R, set the working directory to that containing this readme file. On Mac, for example, you could say

```
setwd('~/Desktop/dynamics_cooperation_status_analysis')
```

Check to see if you're in the right place by typing dir() and see whether this readme file is present. 


The longitudinal analysis uses nine data files as input: 

```
'village1_cooperation_t1.csv' - The cooperation network at time point 1 \
'village1_cooperation_t2.csv' - The cooperation network at time point 2 \
'village1_cooperation_t3.csv' - The cooperation network at time point 3 \
'village1_kinship_t1.csv' - The kinship network at time point 1 \
'village1_kinship_t2.csv' - The kinship network at time point 2 \
'village1_status.csv' - The status ratings of all actors at each time point \
'village1_strength_size.csv' - The physical strength and size of all actors at the first two time points \
'village1_income.csv' - The income of all actors at the first two time points \
'village1_log_age.csv' - The log age of all actors at the first two time points \
```

Please note that only data from the first two time points are used for varying covariates, and dyadic varying covariates (kinship), in the analysis. Another important note is that the rows in the cooperation and kinship matrices represent the individuals nominating. At each time point rows do not change and represent the same individuals. 

When the project folder is the working directory, the longitudinal analysis may run itself (assuming that you have installed all of the dependencies) by calling

```
source('./village1_SAOM.R')
```

However, I would advise to run the analysis script in blocks. 

The cross-sectional analysis uses six data files as input: 

```
'village2_cooperation.csv' - The cooperation network \
'village2_kinship.csv' - The kinship network \
'village2_attributes.csv' - The attributes of actors (i.e. status, strength and size, income, log age)
```

There are three data files that only include participants with full data:

```
'village2_cooperation_deleted.csv' - The cooperation network \
'village2_kinship_deleted.csv' - The kinship network \
'village2_attributes_deleted.csv' - The attributes of actors (i.e. status, strength and size, income, log age)
```

When the project folder is the working directory, the cross-sectional analysis may run itself (assuming that you have installed all of the dependencies) by calling

```
source('village2_ERGM.R')
```

However, I would advise to run the analysis script in blocks. 

The anaylses may take some time. The total time until completion may vary by machine. 

The project is maintained by Daniel Redhead (daniel_redhead@eva.mpg.de) and is hosted at https://github.com/danielRedhead
