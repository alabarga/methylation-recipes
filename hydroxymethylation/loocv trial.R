

glm.fit=glm(group ~ cpg1 + cpg2,
            +             data = data ,family = binomial)
glm.fit$coefficients
(Intercept)        cpg1        cpg2 
142.93970   -31.22971   -22.43628 
glm.fit$residuals

install.packages("boot")

cv.err <- cv.glm(data, glm.fit)
Warning messages:
  1: glm.fit: algorithm did not converge 
2: glm.fit: fitted probabilities numerically 0 or 1 occurred 
3: glm.fit: algorithm did not converge 
4: glm.fit: fitted probabilities numerically 0 or 1 occurred 
5: glm.fit: algorithm did not converge 
6: glm.fit: fitted probabilities numerically 0 or 1 occurred 


# LOGIT REGRESSION

> cpglist <- uniqueSelectedVarsCpG.mval

> data <- data.frame(names = colnames(mval.dmp.fdr10), group = as.factor(pd$Sample_Group), 
                   cpg1 = mval.dmp.fdr10[cpglist[1],], cpg2 = mval.dmp.fdr10[cpglist[2],], 
                   stringsAsFactors = FALSE)

> glm.fit=glm(group ~ cpg1 + cpg2,
            data = data ,family = binomial)

> library(boot)
> cv.err <- cv.glm(data, glm.fit)

> summary(glm.fit)

Call:
  glm(formula = group ~ cpg1 + cpg2, family = binomial, data = data)

Deviance Residuals: 
  Min        1Q    Median        3Q       Max  
-1.11019  -0.00092  -0.00003   0.00035   1.61243  

Coefficients:
  Estimate Std. Error z value Pr(>|z|)
(Intercept)   142.94     146.24   0.977    0.328
cpg1          -31.23      35.47  -0.880    0.379
cpg2          -22.44      22.64  -0.991    0.322

(Dispersion parameter for binomial family taken to be 1)

Null deviance: 38.6733  on 27  degrees of freedom
Residual deviance:  4.5912  on 25  degrees of freedom
AIC: 10.591

Number of Fisher Scoring iterations: 10

> cv.err$delta
[1] 0.1091066 0.1083301

> cv.err$K
[1] 28

> cv.err$seed
[1]         403           1 -1584343791   322326848   513556766  1049212393 -2000701669 -1344700118
[9] -2018174388 -1436760561    50358597   490042108  2000562514  -779632947   311154311  -196591250
