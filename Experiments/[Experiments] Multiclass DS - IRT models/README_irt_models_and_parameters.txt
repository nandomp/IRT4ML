R VARIABLES

************************************************

all_param (in file "irt_parameters.RData"): 

- List of the IRT parameters for each dataset. The element all_param[[i]] is realted to a dataset and 
stores a matrix of size "ninstances x 3", with the guessing, difficulty and discrimination parameters for each instance.  

Example: The parameters of the 100th instance of the 2nd dataset

> all_param[[2]][100,] 
  Gussng      Dffclt      Dscrmn 
 0.07305088 -0.94614917  1.18836896 


************************************************

all_abilities (in file "algor_abilities.RData"):

- Matrix of size "ndatasets x nalgorithms" with the abilities estimated by the IRT model. We have in this 
matrix ndatasets = 23 and nalgorithms = 15.

Example: algorithms' abilities in the 3rd dataset 

> all_abilities[3,]
         J48       J48Unp          LMT       logist          SMO           NB 
 0.684764891  0.684529552 -0.008338582 -0.062735081 -3.247900414  2.828140205 
       Stump          IBK         JRip         PART          SMV          Ada 
 0.017366766  0.771090101  0.651328714  0.034509551 -3.247900414  0.018821796 
     Bagging   LogitBoost     Stacking 
 0.020994640  0.557145705 -3.822916005 

************************************************

all_accuracies (in file "algor_accuracies.RData"):

- Stores the accuracies of the 15 algorithms in the 23 datasets

Example: algorithms' accuracies in the 2nd dataset 

> all_accuracies[2,]

       J48     J48Unp        LMT     logist        SMO         NB      Stump 
 0.9513591  0.9470672  0.9599428  0.9642346  0.9656652  0.9642346  0.9098712 
       IBK       JRip       PART        SMV        Ada    Bagging LogitBoost 
 0.9699571  0.9527897  0.9599428  0.9656652  0.9513591  0.9542203  0.9613734 
  Stacking 
 0.6552217 

************************************************

results (in file "results_responses.RData"):

- List of binary responses of the algorithms for each dataset. The element results[[i]] is related to a dataset and stores a matrix of size "ninstances x nalgorithms" with 1's (right responses) and 0's (wrong responses). It 
is useful to generate other measures like accuracy of the algorithms and average 0\1 per instance.

Example: responses in the 1st instance of the 3rd dataset

> results[[3]][1,]
       J48     J48Unp        LMT     logist        SMO         NB      Stump 
         1          1          1          1          0          0          0 
       IBK       JRip       PART        SMV        Ada    Bagging LogitBoost 
         1          1          0          0          1          1          1 
  Stacking 
         0 


************************************************

all_models (in file "all_3P_IRT_models.RData"):

- List of the IRT models (for while 23 models - a 3P model for each dataset). It can used to produce the ICC 
plots. 

Example: plot the ICC of the 31st instance of the 11st dataset
    
> plot_ICC(all_models,results,all_abilities,11,31)


*************************************************

List of datasets and indices for the above variables

[[1]] = 06_balance-scale
[[2]] = 07_energy    
[[3]] = 08_seeds
[[4]] = 09_teaching        
[[5]] = 10-vertebral-column
[[6]] = 11-balance-scale-v2  
[[7]] = 12-ecoli
[[8]] = 13-flags         







