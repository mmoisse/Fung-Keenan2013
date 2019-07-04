# This script contains all of the R code necessary to 
# calculate CI's from Fung & Keenan (2013)."

# Tak Fung & Kevin Keenan 2013
library("gmp")

BinomCoeff <- function(nn,rr)
{ 
  if(rr<0){return(0)}
  if(rr==0 || rr==nn){return(1)}
  
  # Recur
  return(BinomCoeff(nn-1, rr-1)*(nn/rr))
}

pmfSamplingDistYiN <- function(M, NN, p_i, P_ii, yiN){
   res <- pmfSamplingDistYiN_1(M, NN, p_i, P_ii, yiN)
   if(is.na(res)) {
      res <- pmfSamplingDistYiN_2(M, NN, p_i, P_ii, yiN)
   }
   return(res)
}

# New version of pmfSamplingDistYiN that works for large integers
pmfSamplingDistYiN_2 <- function(M, NN, p_i, P_ii, yiN){
  MaxFunc <- max((yiN/2) - (M*p_i) + (M*P_ii), yiN - NN, 0)
  MinFunc <- min(M*P_ii, yiN/2, M-NN+yiN-(2*M*p_i)+(M*P_ii))
  LowerBound <- ceiling(MaxFunc)
  UpperBound <- floor(MinFunc)
  
  # Calculate P(YiN = yiN) according to eqn 9
  Numerator1 <- 0

  for(i in LowerBound:UpperBound){
    Summand <- BinomCoeff(as.bigz(M*P_ii), as.bigz(i))*BinomCoeff(as.bigz(2*M*(p_i-P_ii)), as.bigz(yiN-(2*i)))*
      BinomCoeff(as.bigz(M+(M*P_ii)-(2*M*p_i)), as.bigz(NN+i-yiN))
    Numerator1 <- Numerator1 + Summand
  }
  
  probOut <- Numerator1/BinomCoeff(as.bigz(M), as.bigz(NN))
  return(asNumeric(probOut))
}

pmfSamplingDistYiN_1 <- function(M, NN, p_i, P_ii, yiN){
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

# load Rsolnp package for optimization
library(Rsolnp)

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
MinJostD <- FindMinJostD$values[length(FindMinJostD$values)]

# Find upper bound of CI for Jost's D
FindMaxJostD <- solnp(x0, fun = JostDMinus, ineqfun = ineqfun1, 
                      ineqLB = ineqLB1, ineqUB = ineqUB1)
MaxJostD <- -FindMaxJostD$values[length(FindMaxJostD$values)]

# Lower bound
MinJostD
# Upper bound
MaxJostD

## END
