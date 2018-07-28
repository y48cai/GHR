// -----------------------------------------------------------------------------
//
// Rejection sampler for the Generalized Inverse Gaussian (GIG) distribution
//
// Algorithm: Devroy, L. (2014), "Random variate generation for the generalized
// inverse Gaussian distribution". Statistics and Computing 24(2), pp 239–246.
//
// Implementation: Copyright Martin Lysy, 2018. <mlysy@uwaterloo.ca>
//
// If X ~ GIG(omega, eta, lambda), then its PDF is
//
// f(x) = K * x^(lambda - 1) * exp(-omega/2 * (x/eta + eta/x),  for   x > 0.
//
// -----------------------------------------------------------------------------

#include <Rcpp.h>
using namespace Rcpp;

class randGIG {
private:
  //static const int MAX_ITER = 1e3;
  // variables
  bool neg_lambda;
  double lambda, omega;
  double p, q, r;
  double s, t;
  double s1, t1;
  double eta, zeta, theta, xi;
  double fq, fqr;
  double U, V, W, X;
  double alpha, chi;
  // required functions
  double psi(double x) {
    return -alpha * (cosh(x) - 1.0) - lambda * (exp(x) - x - 1);
  }
  double psi1(double x) {
    return -alpha * sinh(x) - lambda * (exp(x) - 1.0);
  }
public:
  randGIG() {}; // uninitialized constructor
  randGIG(double _omega, double _lambda); // initialized constructor
  void init(double _omega, double _lambda); // do precomputations
  double rand(); // random number generation
};

// pre-compute all internal constants
inline void randGIG::init(double _omega, double _lambda) {
  double tmp, tmp2;
  omega = _omega;
  lambda = _lambda;
  neg_lambda = (lambda < 0);
  if(neg_lambda) {
    lambda = -lambda;
  }
  alpha = sqrt(omega*omega + lambda*lambda) - lambda;
  // set t and s
  tmp = -psi(1.0);
  if(tmp < .5) {
    t = log(4.0/(alpha + 2.0*lambda));
  } else if(tmp <= 2.0) {
    t = 1.0;
  } else {
    t = sqrt(2.0/(alpha + lambda));
  }
  tmp = -psi(-1.0);
  if(tmp < .5) {
    s = 1.0/alpha;
    s = log(1.0 + s + sqrt(s*(s + 2.0)));
    tmp2 = 1.0/lambda;
    s = tmp2 < s ? tmp2 : s;
  } else if(tmp <= 2.0) {
    s = 1.0;
  } else {
    s = sqrt(4.0/(alpha * cosh(1) + lambda));
  }
  // set (eta, zeta, theta, xi)
  eta = -psi(t);
  zeta = -psi1(t);
  theta = -psi(-s);
  xi = psi1(-s);
  // set p, r
  p = 1/xi;
  r = 1/zeta;
  // set s1, t1, q
  t1 = t - r*eta;
  s1 = s - p*theta;
  q = t1 + s1;
  // set q/(p+q+r) and (q+r)/(p+q+r)
  fq = 1.0/(p+q+r);
  fqr = (q+r) * fq;
  fq *= q;
  return;
}

// constructor
inline randGIG::randGIG(double _omega, double _lambda) {
  init(_omega, _lambda);
}

// random number generation
inline double randGIG::rand() {
  //int count = 0;
  do {
    U = unif_rand();
    V = unif_rand();
    W = unif_rand();
    // set X
    if(U < fq) {
      X = -s1 + q*V;
    } else if(U < fqr) {
      X = t1 - r * log(V);
    } else {
      X = -s1 + p * log(V);
    }
    // set chi(X)
    if(X < -s1) {
      chi = exp(-theta + xi * (X+s));
    } else if(X <= t1) {
      chi = 1.0;
    } else {
      chi = exp(-eta - zeta * (X-t));
    }
  // } while(count++ < MAX_ITER && W * chi > exp(psi(X)));
  } while(W * chi > exp(psi(X)));
  chi = lambda/omega;
  X = (chi + sqrt(1 + chi*chi)) * exp(X);
  if(neg_lambda) {
    X = 1.0/X;
  }
  return X;
}

// try to prevent R users from interacting with C++ code directly,
// as it's much easier and cleaner to do argument checking at R level.
//[[Rcpp::export(".GenerateGIG")]]
NumericVector GenerateGIG(int n, NumericVector omega, NumericVector eta,
			  NumericVector lambda) {
  // check which parameters are single-valued
  bool isSingleOmega = omega.length() == 1;
  bool isSingleEta = eta.length() == 1;
  bool isSingleLambda = lambda.length() == 1;
  NumericVector x(n);
  randGIG gig; // GIG random number generator
  if(isSingleOmega && isSingleLambda) {
    // precompute internal constants once
    gig.init(omega[0], lambda[0]);
    for(int ii=0; ii<n; ii++) {
      x[ii] = eta[(!isSingleEta)*ii] * gig.rand();
    }
  } else {
    // compute internal constants at each step
    for(int ii=0; ii<n; ii++) {
      gig.init(omega[(!isSingleOmega)*ii], lambda[(!isSingleLambda)*ii]);
      x[ii] = eta[(!isSingleEta)*ii] * gig.rand();
    }
  }
  return x;
}
