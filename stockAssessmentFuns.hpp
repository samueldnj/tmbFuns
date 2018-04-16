// <><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>
// stockAssessmentFuns.hpp
// 
// A library of stock assessment related C++ functions, defined originally
// in separate models but compiled here as a library for future
// use.
// 
// Author: Samuel D. N. Johnson
// Date: 16 April, 2018
// 
// Last updated: 16 April, 2018
// <><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><><>><><>><><>><>


// posfun()
// Compares a quantity x to a threshold epsilon and if below the threshold, 
// increases x above the threshold and increments a penalty
// inputs:    x = value to be checked
//            eps = threshold
//            & pen = external penalty variable
// output:    y conditional on value of x-eps:
//            y = x if x > eps
//            y = eps / (2 - eps/x)
// side-effs: pen += 0.01 * (x - eps)^2
// Usage:     required for state space models when catch may
//            become larger than biomass under stochastic conditions.
// Source:    kaskr.github.com/adcomp wiki, modified by SDNJ
//            to actually make output > eps
template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x-eps,2), Type(0.));
  return CppAD::CondExpGe(x, eps, x, eps/(Type(2.0)-eps/x));
}


// square()
// Shortcut wrapper for squaring quantities
// inputs:    x = quantity to be squared
// ouputs:    y = x^2
// Usage:     when squaring.
// Source:    Almost surely stolen from somewhere.
template <class Type> 
Type square(Type x){return pow(x,2);}
VECTORIZE1_t(square)


// addCompNoise() 
// Adds noise to compositional data assuming a
// multivariate logistic distribution
// inputs:    inputComp = vector of input compositional poportions
//            noise = vector of random errors to add to compositional data
// ouputs:    outputComp = new compositional data with error
// Usage:     Used when making random walks in compositional data
//            or for simulating noise (e.g. ageing error)
// Author:    S. D. N. Johnson
template<class Type>
vector<Type> addCompNoise(  vector<Type>  inputComp, 
                            vector<Type>  noise )
{
  int nComps = inputComp.size();
  // Take logs 
  vector<Type> logInputComp(nComps);
  logInputComp.fill(0.0);
  for( int c = 0; c < nComps; c++)
   logInputComp(c) = log( inputComp(c) );
  // Now add noise
  vector<Type> outputComp(nComps);
  outputComp = exp(logInputComp + noise);
  Type tmpTotal = outputComp.sum();
  outputComp /= tmpTotal;

  return outputComp;
  
} // End addCompNoise()


// CRmort()
// A Chapman-Robson total mortality (Z) estimator
// inputs:    ageComp = vector of age composition data (proportions)
//            kage = age (k) of full selectivity/recruitment
//            Aplus = plus group age A
//            minObs = minimum number of observations to truncate
//                     age composition data
//            & Ztmp = external total mortality value
// ouputs:    NA, void function
// side-effs: total mortality estimate is saved in variable
//            passed in as &Ztmp
// Author:    S. D. N. Johnson
// Source:    Chapman and Robson, 1960;
//            Dunn et al, 2002
template<class Type>
void CRmort(  vector<Type> ageComp, 
              int kage, 
              int Aplus, 
              int minObs,
              Type& Ztmp ){
  // Restrict to entries between kage and Aplus
  int maxAges = Aplus - kage + 1;
  // Create a vector to hold observations, fill with zeroes
  vector<Type> ageObs( maxAges );
  ageObs.fill(0);
  // Average age
  Type abar = 0.0;

  // Loop and fill vector
  for( int a = kage-1; a < Aplus; a ++ )
  {
    if( ageComp(a) >= minObs )
    {
      ageObs(a - kage + 1) = ageComp(a);
      abar += (a - kage + 1) * ageComp(a);
    } else break;
  }
  
  // Now compute total observations and complete abar calc
  Type     N = ageObs.sum();
  abar      /= N;


  // Return Z estimate if there are any age observations
  if( abar == 0 ) Ztmp =  -1;
  else  Ztmp =  log( ( 1 + abar - 1/N ) / abar );
} // End CRmort()

// solveBaranovDD()
// Newton-Rhapson solver for Baranov catch equation for a population
// modeled with no age classes (e.g. Delay Difference formulation) at
// a given time step
// inputs:    nIter = number of NR iterations
//            Bstep = fraction of NR step (Jacobian) to take at each iteration
//            C = Catch
//            M = natural mortality
//            B = Biomass
//            & Z = total mortality (external variable)
//            & F = Fishing mortality (external variable)
// returns:   NA, void function
// Side-effs: variables passed as Z, F overwritten with total, fishing mortality
// Author:    Modified by S. D. N. Johnson from S. Rossi and S. P. Cox
template<class Type>
void solveBaranovDD(  int   nIter,
                      Type  Bstep,
                      Type  C,
                      Type  M,
                      Type  B,
                      Type& Z,
                      Type& F)
{
  Type f    = 0.;   // Function value
  Type J    = 0.;   // Jacobian
  Type newZ = 0.;   // Updated Z
  Type tmp  = 0.;   // predicted catch given F

  // Initial approximation of F
  F = C / (C+B);
  
  newZ = M + F;
  Z    = M + F;

  // Refine F
  for( int i=0; i<nIter; i++ )
  {
    // Total mortality
    Z     = newZ;
    newZ  = M;
    // Predicted catch given F
    tmp   = B*(1.-exp(-Z))*F/Z;

    // Function value: difference of pred - obs catch
    f   = C - tmp;
    // Jacobian
    J   = -B * ((1. - exp(-Z)) * M / pow(Z,2) + exp( -Z ) * F / Z);

   
    // Updated fishing mortality
    F -= Bstep * f / J;

    // Updated total mortality
    newZ += F;

  }  // end i

}  // end solveBaranovDD()

// negLogLogisticNormal()
// Calculates the negative log density for a logistic normal
// distribution.
// inputs:    y = vector of observed proportions
//            p = vector of parameters (true class proportions)
//            var = variance of logistic normal distribution
// outputs:   nld = negative log density of logistic normal distribution
// Usage:     For computing the likelihood of observed compositional data
// Source:    S. D. N. Johnson
// Reference: Schnute and Haigh, 2007
template<class Type>
Type negLogLogisticNormal(  vector<Type> y, 
                            vector<Type> p, 
                            Type var )

{
  // Get dimension of MV distribution
  int N = y.size();

  // Take log var
  Type lnvar = log( var );

  // Variable to returning
  Type nld = 0;

  // Take geometric means
  Type ytilde = pow(y.prod(),1/N);
  Type ptilde = pow(p.prod(),1/N);

  // Add variance term
  nld += (N - 1) * lnvar / 2;

  // Now loop over dimensions, add residuals
  for( int i = 0; i < N; i++ )
    nld += square(log(y(i)/ytilde) - log(p(i)/ptilde))/2/var;

  return nld;
} // end negLogLogisticNormal()


