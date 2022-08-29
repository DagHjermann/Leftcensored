# leftcensored  

An R package for linear or non-linear regression analysis of left-censored data (and later, right-censored as well?) with or without measurement error in the y value. 
  
Specifications:  

- The package currently does only regression with one response variable depending on one predictor variable (both continuous), and where the error is normally distributed. I.e., models on the form   
```
y = F(x) + error  
```

- The function F can be either linear (i.e. F(x) = a + b*x) or non-linear (using thin-plate splines)   
- The threshold for left-censoring can vary among observations. For instance, one can easily accomodate chemical data where for instance the level of quantification (LOQ) is varying among years or even among samples   
- The user can specify a measurement error in y; this error can also either be constanct or vary between observations   
- Parameter estimation is done using Bayesian methodology (Markov Chain Monte Carlo) with non-informative prior distributions  
- The method of Qi et al. (2022) is the basis for parameter estimation, and the method also returns the DIC (Deviance information criterion) for models to aid model selection between e.g. linear and non-linear models   

Future developments:  

- Regression of several related response variables with similar (or opposite) trends in relation to x. This may be modelled as having an "undelying" but unobserved (latent) variable (y_unobserved) which dictates the observed variables y1-y.k:      
```
y1 = a1 + b1*y_unobserved + error  
y2 = a2 + b2*y_unobserved + error  
...
y.k = a.k + b.k*y_unobserved + error  

y_unobserved = F(x) + error  

```

- 

Files and functions:  

- `leftcensored_lm.R`     
    - lc_linear   
    
- `leftcensored_lm_qi_measerror.R`    
    - lc_linear_qi_measerror  
    - lc_linear_qi_measerror_uncens  
    
- `leftcensored_lm_dinterval.R`    
    - lc_linear_dinterval  
    
- `leftcensored_lm_qi.R`    
    - lc_linear_qi  
    - lc_linear_qi_uncens   
    
- `leftcensored_fixedsplines_tp.R`     
    - lc_fixedsplines_tp  
    - get_jags_model_code      
    
- `leftcensored_fixedsplines_qi.R`     
    - lc_fixedsplines_tp  

    


    
