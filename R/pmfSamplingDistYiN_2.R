
# Load gmp package
library("gmp")


# Recursive function to calculate binomial coefficient
# Works with variables of type "bigz", unlike "choose" function
BinomCoeff <- function(nn,rr)
{ 
  if(rr==0 || rr==nn){return(1)}
  
  # Recur
  return(BinomCoeff(nn-1, rr-1)*(nn/rr))
}


# Check that BinomCoeff gives same results as choose
choose(10,3)
BinomCoeff(10,3)
choose(145,31)
BinomCoeff(145,31)


# Website code for pmfSamplingDistYiN
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


# Check that pmfSamplingDistYiN_2 gives same results as pmfSamplingDistYiN
# for relatively small integers
pmfSamplingDistYiN(1000, 10, 0.1, 0.04, 1)
pmfSamplingDistYiN_2(1000, 10, 0.1, 0.04, 1)


# Check that pmfSamplingDistYiN_2 works for large integers, unlike
# pmfSamplingDistYiN
pmfSamplingDistYiN(5000, 600, 0.02, 0.02^2, 0)
pmfSamplingDistYiN_2(5000, 600, 0.02, 0.02^2, 0)
pmfSamplingDistYiN(10000, 1000, 0.01, 0.01^2, 0)
pmfSamplingDistYiN_2(10000, 1000, 0.01, 0.01^2, 0)
# Output from pmfSamplingDistYiN_2 matches that from Mathematica




