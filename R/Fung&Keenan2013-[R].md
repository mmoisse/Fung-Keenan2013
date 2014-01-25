<center>
Supporting Webpage 2 for Fung & Keenan (2014): A Description of R functions used to calculate  confidence intervals
========================================================

Tak Fung^1,2 and Kevin Keenan^3
----------------------------------------------------

<h6>
<sup>1</sup> National University of Singapore, Department of Biological Sciences, 14 Science Drive 4, Singapore 117543

<sup>2</sup> Queens University Belfast, School of Biological Sciences, Belfast BT9 7BL, UK

<sup>3</sup> Queens University Belfast, Institute for Global Food Security, School of Biological Sciences, Belfast BT9 7BL, UK

</center>
</h6>

## Introduction
This document describes the functionality of the `R` code converted from the _Mathematica_ code used in [Fung & Keenan (2014)](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0085925). The code was written and tested using _Mathematica_ v5.0[1] and subsequently converted and tested in `R`. This document describes four separate functions named `pmfSamplingDistYiN`, `AcceptanceRegion`, `CIforpiCasePiiUnknown` and `CIforpiCasePiiKnown`. In addition, code is presented at the end to calculate a CI for Jost's $D$, for the butterfly example examined in the main text of Fung & Keenan (2014). The original _Mathematica_ version of this web document can be found <a href="http://rpubs.com/kkeenan02/Fung-Keenan-Mathematica/" target="_blank">here</a>

## pmfSamplingDistYiN
This function returns $P(Y_{i,N} = y_{i,N})$ as specified by equation (9) in the main text, given $M$, $N$, $p_{i}$, $P_{ii}$ and $y_{i,n}$. Here, $Y_{i,N}$ is the random variable specifying the number of copies of allele $A_{i}$ in a sample of size $N$ taken from a finite diploid population of size $M$, with the frequency of allele $A_{i}$ in the population being $p_{i}$ and the frequency of homozygotes of allele $A_{i}$ in the population being $P_{ii}$.




### pmfSamplingDistYiN source code (`R`)

```r
pmfSamplingDistYiN <- function(M, NN, p_i, P_ii, yiN){
  MaxFunc <- max((yiN/2) - (M*p_i) + (M*P_ii), yiN - NN, 0)
  MinFunc <- min(M*P_ii, yiN/2, M-NN+yiN-(2*M*p_i)+(M*P_ii))
  LowerBound <- ceiling(MaxFunc)
  UpperBound <- floor(MinFunc)
  
  # Calculate P(YiN = yiN) according to eqn 9
  Numerator1 <- 0
  
  for(i in LowerBound:UpperBound){
    Summand <- choose(M*P_ii, i)*choose(2*M*(p_i-P_ii), yiN-(2*i))*
      choose(M+(M*P_ii)-(2*M*p_i), NN+i-yiN)
    Numerator1 <- Numerator1 + Summand
  }
  
  probOut <- Numerator1/choose(M, NN)
  return(probOut)
}
```


### pmfSamplingDistYiN example

```r
# test timing
system.time(res <- pmfSamplingDistYiN(1000, 10, 0.1, 0.04, 1))
```

```
   user  system elapsed 
      0       0       0 
```

```r
# print results
print(res)
```

```
[1] 0.250394
```


## AcceptanceRegion
This function tests the null hypothesis $H_{o}:p_{i}=p_{i,0}, P_{ii}=P_{ii,0}$ for an observed value of $y_{i,N}$, $latex \hat{y}_{i,N}$, given the sampling scenario considered. It does this by calculating the acceptance region for a specified significance level, $latex \alpha$, and then determining whether $latex \hat{y}_{i,N}$ lies within this region. The outputs are bounds of the acceptance region and an indication of whether $latex \hat{y}_{i,N}$ falls within this region or not ('1' == TRUE, '0' == FALSE, respectively).


### AcceptanceRegion source code (`R`)

```r
AcceptanceRegion <- function(M, NN, p_i0, P_ii0, yiNobs, alpha){
  
  # find the lower bound of the acceptance region
  SumProb <- 0
  yiNlow <- 0
  dummy1 <- 0
  while(SumProb <= (alpha/2)){
    dummy1 <- pmfSamplingDistYiN(M, NN, p_i0, P_ii0, yiNlow)
    SumProb <- SumProb + dummy1
    yiNlow <- yiNlow + 1
  }
  yiNlow <- yiNlow - 1
  dummy1 <- pmfSamplingDistYiN(M, NN, p_i0, P_ii0, yiNlow)
  SumProb <- SumProb - dummy1
  # find the upper bound of the acceptance region
  SumProb2 <- 1
  yiNup <- 2*NN
  while(SumProb2 >= (1-(alpha/2))){
    dummy1 <- pmfSamplingDistYiN(M, NN, p_i0, P_ii0, yiNup)
    SumProb2 <- SumProb2 - dummy1
    yiNup <- yiNup - 1
  }
  yiNup <- yiNup + 1
  dummy1 <- pmfSamplingDistYiN(M, NN, p_i0, P_ii0, yiNup)
  SumProb2 <- SumProb2 + dummy1
  # Test if yiNobs is with in the accptance region
  test <- yiNobs >= yiNlow && yiNobs <= yiNup
  out <- list(lowerbound = yiNlow, 
              upperbound = yiNup,
              result = test)
  return(unlist(out))
}
```


### AcceptanceRegion example

```r
# test timing
system.time(res <- AcceptanceRegion(1000, 30, 0.625, 0.25, 50, 0.05))
```

```
   user  system elapsed 
   0.02    0.00    0.02 
```

```r
# print results
print(res)
```

```
lowerbound upperbound     result 
        33         42          0 
```


## CIforpiCasePiiUnknown
This function calculates $latex \geq 100(1-\alpha)$ % confidence intervals (CI's) for $p_i$ and $P_{ii}$ given $M$, $N$ and $latex \hat{y}_{i,N}$. These CI's are computed using equations (14a) and (14b) in the main text, and uses `pmfSamplingDistYiN` and `AcceptanceRegion`. The outputs are the lower and upper limits of the CI for $p_i$ followed by those of the CI for $P_{ii}$.

### CIforpiCasePiiUnknown source code (`R`)

```r
CIforpiCasePiiUnknown <- function(M, NN, yiNobs, alpha){
  pi0List <- vector()
  Pii0List <- vector()
  out <- matrix()
  i = 0
  while(i  <= (2*M)){
    p_i0 <- i/(2*M)
    if((p_i0 >= (yiNobs/(2*M))) && (p_i0 <= (1 - (((2*NN)-yiNobs)/(2*M))))){
      k=0
      while(k <= i){
        P_ii0 <- k/(2*M)
        if((P_ii0 >= max(0, (2*p_i0)-1)) && (P_ii0 <= p_i0)){
          out1 <- AcceptanceRegion(M, NN, p_i0, P_ii0, yiNobs, alpha)
          if(out1[3] == 1L){
            pi0List <- c(pi0List, p_i0)
            Pii0List <- c(Pii0List, P_ii0)
          }
          out2 <- c(out1,p_i0, P_ii0)
          if(all(is.na(out))){
            out <- out2
          } else{
            out <- rbind(out, out2) 
          }
        }
        k <- k + 1
      }
    }
    i <- i + 1
  }
  if(length(pi0List) > 0){
    piCIlow <- min(pi0List)
    piCIup <- max(pi0List)
    PiiCIlow <- min(Pii0List)
    PiiCIup <- max(Pii0List)
  } else if(length(pi0List) == 0){
    largestAlpha <- 0
    outCounter <- 1
    while(i <= (2*M)){
      p_i0 <- i/(2*M)
      if((p_i0 >= (yiNobs/(2*M))) && (p_i0 <= (1 - (((2*NN)-yiNobs)/(2*M))))){
        while(k <= i){
          P_ii0 <- k/(2*M)
          if((P_ii0 = max[0, (2*p_i0) - 1]) && (P_ii0 = p_i0)){
            if(yiNobs < out[outCounter, 1]){
              Sumprob <- 0
              m = 0
              while(m <= (yiNobs - 1)){
                Sumprob <- Sumprob + pmfSamplingDistYiN(M, NN, p_i0List[i],
                                                        P_ii0List[i], m)
                m <- m + 1
              }
              largestAlpha <- c(largestAlpha, Sumprob*2) 
            }
            if(yiNobs > out[outCounter, 2]){
              Sumprob2 <- 0
              m = (2*NN)
              while(m >= (yiNobs +1)){
                Sumprob2 <- Sumprob2 + pmfSamplingDistYiN(M, NN, p_i0List[i],
                                                          P_ii0List[i], m)
                m <- m - 1
              }
              largestAlpha <- c(largestAlpha, Sumprob2*2)
            }            
          }
          outCounter <- outCounter + 1
        }
      }
    }
    alphaMax <- -1
    alphaMaxIdx <- -1
    j = 1
    while(j <= length(largestAlpha)){
      if(largestAlpha[j] > alphaMax){
        alphaMax <- largestAlpha[j]
        alphaMaxIdx <- j
      }
      j <- j+1
    }
    piCIlow <- out[alphaMaxIdx, 4]
    piCIup <- piCIlow
    PiiCIlow <- out[alphaMaxIdx, 5]
    PiiCIup <- out[alphaMaxIdx, 5]
  }
  list(piCI = c(piCIlow, piCIup), 
       PiiCI = c(PiiCIlow, PiiCIup))
}
```


### CIforpiCasePiiUnknown example

```r
# test timing
system.time(res <- CIforpiCasePiiUnknown(100, 30, 5, 0.05))
```

```
   user  system elapsed 
  58.32    0.02   59.61 
```

```r
# print results
print(res)
```

```
$piCI
[1] 0.025 0.200

$PiiCI
[1] 0.000 0.195
```


## CIforpiCasePiiKnown

This function calculates a $latex \geq 100(1-\alpha)$ % CI for $p_{i}$, given $M$, $N$ and $latex \hat{y}_{i,N}$, for a population with maximum homozygosity ($p_{i}=P_{ii}$). The function uses equation (13) in the main text to calculate the CI, and uses `pmfSamplingDistYiN` and `AcceptanceRegion`. The outputs are the lower and upper bounds for the CI. The code of `CIforpiCasePiiKnown` can be easily adapted to calculate CI's under the scenario of HWE or the Scenario of minimum homozygosity, by altering two lines indicated within the code.

### CIforpiCasePiiKnown source code (`R`)

```r
CIforpiCasePiiKnown <- function(M, NN, yiNobs, alpha){
  out <- matrix()
  pi0List <- vector()
  i <- 0
  while(i <= (2*M)){
    p_i0 <- i/(2*M)
    if((p_i0 >= (yiNobs/(2*M))) && (p_i0 <= (1-(((2*NN)-yiNobs)/(2*M))))){
      # Max homoxygosity (Pii0 <- pi0)
      # HWE (Pii0 <- pi0^2)
      # Min homozygosity (Pii0 <- max(0, (2^pi0)-1))
      P_ii0 <- p_i0
      # make sure Pii * M is an integer
      if((P_ii0 >= max(0, (2*p_i0)-1)) && (P_ii0 <= p_i0) && (P_ii0*M == as.integer(P_ii0*M))){
        out1 <- AcceptanceRegion(M, NN, p_i0, P_ii0, yiNobs, alpha)
        if(out1[3] == 1L){
          pi0List <- c(pi0List, p_i0)
          out2 <- c(out1 , p_i0, P_ii0) 
          if(all(is.na(out))){
            out <- out2
          } else{
            out <- rbind(out, out2) 
          }
        }
      }
    }
    i <- i+1
  }
  if(length(pi0List) > 0){
    piCIlow <- min(pi0List)
    piCIup <- max(pi0List)
  } else if(length(pi0List) == 0){
    largestAlpha <- 0
    outCounter <- 1
    i <- 1
    while(i <= (2*M)){
      p_i0 <- i/(2*M)
      if((p_i0 >= (yiNobs/(2*M))) && (p_i0 <= (1- (((2*NN)-yiNobs)/(2*M))))){
        # Max homoxygosity (Pii0 == pi0)
        # HWE (Pii0 == pi0^2)
        # Min homozygosity (Pii0 == max(0, (2^pi0)-1))
        P_ii0 <- p_i0
        # make sure Pii * M is an integer
        if((P_ii0 >= max(0, (2*p_i0)-1)) && (P_ii0 <= p_i0) && (P_ii0*M == as.integer(P_ii0*M))){
          if(yiNobs <= out[outCounter, 1]){
            Sumprob <- 0
            m <- 0
            while(m <= (yiNobs - 1)){
              Sumprob <- Sumprob + pmfSamplingDistYiN(M, NN, p_i0List[i],
                                                      P_ii0List[i], m)
              m <- m+1
            }
            largestAlpha <- c(largestAlpha, Sumprob*2) 
          }
          if(yiNobs > out[largestAlpha, 2]){
            Sumprob2 <- 0
            m <- (2*NN)
            while(m >= (yiNobs + 1)){
              Sumprob2 <- Sumprob2 + pmfSamplingDistYiN(M, NN, p_i0List[i],
                                                        P_ii0List[i], m)
              m <- m - 1
            }
            largestAlpha <- c(largestAlpha, Sumprob2*2)
          }
          outCounter <- outCounter + 1
        }
      }
      i <- i + 1
    }
    alphaMax <- -1
    alphaMaxIdx <- -1
    j <- 1
    while(j <= length(largestAlpha)){
      if(largestAlpha[j] > alphaMax){
        alphaMax <- largestAlpha[j]
        alphaMaxIdx <- j
      }
      j <- j + 1
    }
    piCIlow <- out[alphaMaxIdx, 4]
    piCIup <- piCIlow
  }
  res <- list(piCIlower = piCIlow,
              piCIupper = piCIup)
  return(unlist(res))
}  
```


### CIforpiCasePiiKnown example

```r
# test timing
system.time({res <- CIforpiCasePiiKnown(438, 53, 1, 0.05/3)})
```

```
   user  system elapsed 
   1.95    0.00    2.01 
```

```r
# print results
print(res)
```

```
 piCIlower  piCIupper 
0.00228311 0.07990868 
```


## Calculate $latex \geq$ 95% CI for Jost's $D$
In this example, the $latex \geq$ 95% CI for Jost's $D$ is calculated for the butterfly example given in the main text. This code uses `CIforpiCasePiiKnown` to calculate the CI's for the population frequencies of each of the three alleles in each of the two populations. Then it uses these CI's to minimize and maximize Jost's $D$, which corresponds to the lower and upper bounds for its CI, respectively.


```r
# Find CI's for p1, p2 and p3, for Prasto population
CIforp1Prasto <- CIforpiCasePiiKnown(2*219, 53, 1, 0.05/3)
CIforp2Prasto <- CIforpiCasePiiKnown(2*219, 53, 80, 0.05/3)
CIforp3Prasto <- CIforpiCasePiiKnown(2*219, 53, 24, 0.05/3)

# Find CI's for q1, q2 and q3, for Finstrom population
CIforq1Finstrom <- CIforpiCasePiiKnown(7*219, 74, 4, 0.05/3)
CIforq2Finstrom <- CIforpiCasePiiKnown(7*219, 74, 123, 0.05/3)
CIforq3Finstrom <- CIforpiCasePiiKnown(7*219, 74, 4, 0.05/3)

# Define Jost's D
JostD <- function(freq){
  plist <- c(freq[1:3], 1-sum(freq[1:3]))
  qlist <- c(freq[4:6], 1-sum(freq[4:6]))
  JostT <- sum(((plist + qlist)/2)^2)
  JostS <- (sum(plist^2) + sum(qlist^2))/2
  return(((JostT/JostS) - 1)/((1/2) - 1))
}

# JostDMinus is -JostD, required because optimization program
# used below minimizes functions, and minimizing -JostD
# is equivalent to maximizing JostD
JostDMinus=function(x){
  -JostD(x)
}

# Check if the packages 'Rsolnp' needs to be installed before loading it. 
# Rsolnp is used for the optimisation step.
if("Rsolnp" %in% rownames(installed.packages()) == FALSE){ 
  install.packages("Rsolnp", repo = "http://cran.rstudio.com",
                   dep = TRUE)
}
library("Rsolnp")


# define inequalities
ineqfun1 <- function(freq){
  return(c(freq[1:6], sum(freq[1:3]), sum(freq[4:6])))
}

# Define lower and upper bounds on inequalities to be used in 
# optimization program
ineqLB1=c(CIforp1Prasto[1], CIforp2Prasto[1], CIforp3Prasto[1],
          CIforq1Finstrom[1], CIforq2Finstrom[1], CIforq3Finstrom[1],
          0, 0)
ineqUB1=c(CIforp1Prasto[2], CIforp2Prasto[2], CIforp3Prasto[2],
          CIforq1Finstrom[2], CIforq2Finstrom[2], CIforq3Finstrom[2],
          1, 1)

# Define initial values of population allele frequencies to be used in optimization program,
# which are mid-points of CI's
x0=c(mean(CIforp1Prasto), mean(CIforp2Prasto), mean(CIforp3Prasto),
     mean(CIforq1Finstrom), mean(CIforq2Finstrom), mean(CIforq3Finstrom))

# Find lower bound of CI for Jost's D
FindMinJostD <- solnp(x0, fun = JostD, ineqfun = ineqfun1, 
                      ineqLB = ineqLB1, ineqUB = ineqUB1)
```

```

Iter: 1 fn: 0.00002344	 Pars:  0.01278 0.87206 0.11416 0.01451 0.87382 0.10894
Iter: 2 fn: 0.00002344	 Pars:  0.01277 0.87208 0.11416 0.01450 0.87383 0.10894
solnp--> Completed in 2 iterations
```

```r
MinJostD <- FindMinJostD$values[length(FindMinJostD$values)]

# Find upper bound of CI for Jost's D
FindMaxJostD <- solnp(x0, fun = JostDMinus, ineqfun = ineqfun1, 
                      ineqLB = ineqLB1, ineqUB = ineqUB1)
```

```

Iter: 1 fn: -0.1858	 Pars:  0.020548 0.598174 0.381279 0.002609 0.913894 0.002609
Iter: 2 fn: -0.1858	 Pars:  0.020548 0.598174 0.381279 0.002609 0.913894 0.002609
solnp--> Completed in 2 iterations
```

```r
MaxJostD <- -FindMaxJostD$values[length(FindMaxJostD$values)]

# Lower bound
MinJostD
```

```
[1] 2.34359e-05
```

```r
# Upper bound
MaxJostD
```

```
[1] 0.185774
```


## Source code
In the interest of reproducibility, all source code, both for the <a href="http://rpubs.com/kkeenan02/Fung-Keenan-Mathematica/" target="_blank">_Mathematica_</a> and the `R` versions of this document can be freely accessed at <a href="https://github.com/kkeenan02/Fung-Keenan2013/" target="_blank">github</a>.

### References for Supporting Webpage 2
  1. Wolfram Research Inc. (2003) Mathematica Edition: Version 5.0. Champaign, Illinois, USA, Wolfram Research Inc.
  
  2. R Core Team (2013). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL: http://www.R-project.org/.
  
  3. Alexios Ghalanos and Stefan Theussl (2012). Rsolnp: General
Non-linear Optimization Using Augmented Lagrange Multiplier Method.
R package version 1.14.
  
  4. Fung T, Keenan K (2014) Confidence Intervals for Population Allele Frequencies: The General Case of Sampling from a Finite Diploid Population of Any Size. PLoS ONE 9(1): e85925. doi:10.1371/journal.pone.0085925

## Reproducibility

```
R Under development (unstable) (2014-01-19 r64835)
Platform: x86_64-w64-mingw32/x64 (64-bit)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] knitr_1.5.1

loaded via a namespace (and not attached):
[1] digest_0.6.4   evaluate_0.5.1 formatR_0.10   stringr_0.6.2 
[5] tools_3.1.0   
```

