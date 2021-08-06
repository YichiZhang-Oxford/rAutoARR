# include <RcppArmadillo.h>
# include <algorithm>
# include <string>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

              /* --------------------------------- */    
/*------------ Tuning-Free Huber Regression Funtions ---------------*/
              /* --------------------------------- */ 


// [[Rcpp::export]]
int sgn(const double x)
{
  return (x > 0) - (x < 0);
}

// [[Rcpp::export]]
double f1(const double x, const arma::vec &resSq, const int n, const double rhs)
{
  return arma::mean(arma::min(resSq / x, arma::ones(n))) - rhs;
}

// [[Rcpp::export]]
double rootf1(const arma::vec &resSq, const int n, const double rhs, double low, double up, const double tol = 0.001, const int maxIte = 500)
{
  int ite = 1;
  while (ite <= maxIte && up - low > tol)
  {
    double mid = 0.5 * (up + low);
    double val = f1(mid, resSq, n, rhs);
    if (val < 0)
    {
      up = mid;
    }
    else
    {
      low = mid;
    }
    ite++;
  }
  return 0.5 * (low + up);
}

// [[Rcpp::export]]
double f2(const double x, const arma::vec &resSq, const int N, const double rhs)
{
  return arma::mean(arma::min(resSq / x, arma::ones(N))) - rhs;
}

// [[Rcpp::export]]
double rootf2(const arma::vec &resSq, const int n, const int d, const int N, const double rhs, double low, double up, const double tol = 0.001,
              const int maxIte = 500)
{
  int ite = 0;
  while (ite <= maxIte && up - low > tol)
  {
    double mid = 0.5 * (up + low);
    double val = f2(mid, resSq, N, rhs);
    if (val < 0)
    {
      up = mid;
    }
    else
    {
      low = mid;
    }
    ite++;
  }
  return 0.5 * (low + up);
}

// [[Rcpp::export]]
double g1(const double x, const arma::vec &resSq, const int n, const double rhs)
{
  return arma::mean(arma::min(resSq / x, arma::ones(n))) - rhs;
}

// [[Rcpp::export]]
double rootg1(const arma::vec &resSq, const int n, const double rhs, double low, double up, const double tol = 0.001, const int maxIte = 500)
{
  int ite = 0;
  while (ite <= maxIte && up - low > tol)
  {
    double mid = 0.5 * (up + low);
    double val = g1(mid, resSq, n, rhs);
    if (val < 0)
    {
      up = mid;
    }
    else
    {
      low = mid;
    }
    ite++;
  }
  return 0.5 * (low + up);
}

// [[Rcpp::export]]
double huberDer(const arma::vec &res, const double tau, const int n)
{
  double rst = 0.0;
  for (int i = 0; i < n; i++)
  {
    double cur = res(i);
    rst -= std::abs(cur) <= tau ? cur : tau * sgn(cur);
  }
  return rst / n;
}

// [[Rcpp::export]]
double huberMean(arma::vec X, const double tol = 0.001, const int iteMax = 500)
{
  int n = X.n_rows;
  double rhs = std::log(n) / n;
  double mx = arma::mean(X);
  X -= mx;
  double tau = arma::stddev(X) * std::sqrt((long double)n / std::log(n));
  double derOld = huberDer(X, tau, n);
  double mu = -derOld, muDiff = -derOld;
  arma::vec res = X - mu;
  arma::vec resSq = arma::square(res);
  tau = std::sqrt((long double)rootf1(resSq, n, rhs, arma::min(resSq), arma::accu(resSq)));
  double derNew = huberDer(res, tau, n);
  double derDiff = derNew - derOld;
  int ite = 1;
  while (std::abs(derNew) > tol && ite <= iteMax)
  {
    double alpha = 1.0;
    double cross = muDiff * derDiff;
    if (cross > 0)
    {
      double a1 = cross / derDiff * derDiff;
      double a2 = muDiff * muDiff / cross;
      alpha = std::min(std::min(a1, a2), 100.0);
    }
    derOld = derNew;
    muDiff = -alpha * derNew;
    mu += muDiff;
    res = X - mu;
    resSq = arma::square(res);
    tau = std::sqrt((long double)rootf1(resSq, n, rhs, arma::min(resSq), arma::accu(resSq)));
    derNew = huberDer(res, tau, n);
    derDiff = derNew - derOld;
    ite++;
  }
  return mu + mx;
}

// [[Rcpp::export]]
double hMeanCov(const arma::vec &Z, const int n, const int d, const int N, double rhs, const double epsilon = 0.0001, const int iteMax = 500)
{
  double muOld = 0;
  double muNew = arma::mean(Z);
  double tau = arma::stddev(Z) * std::sqrt((long double)n / (2 * std::log(d) + std::log(n)));
  int iteNum = 0;
  arma::vec res(n), resSq(n), w(n);
  while ((std::abs(muNew - muOld) > epsilon) && iteNum < iteMax)
  {
    muOld = muNew;
    res = Z - muOld;
    resSq = arma::square(res);
    tau = std::sqrt((long double)rootf2(resSq, n, d, N, rhs, arma::min(resSq), arma::accu(resSq)));
    w = arma::min(tau / arma::abs(res), arma::ones(N));
    muNew = arma::as_scalar(Z.t() * w) / arma::accu(w);
    iteNum++;
  }
  return muNew;
}

// [[Rcpp::export]]
Rcpp::List huberCov(const arma::mat& X)
{
  const int n = X.n_rows;
  const int p = X.n_cols;
  
  double rhs2 = (2 * std::log(p) + std::log(n)) / n;
  arma::vec mu(p);
  arma::mat sigmaHat(p, p);
  for (int j = 0; j < p; j++)
  {
    mu(j) = huberMean(X.col(j));
    double theta = huberMean(arma::square(X.col(j)));
    double temp = mu(j) * mu(j);
    if (theta > temp)
    {
      theta -= temp;
    }
    sigmaHat(j, j) = theta;
  }
  int N = n * (n - 1) >> 1;
  arma::mat Y(N, p);
  for (int i = 0, k = 0; i < n - 1; i++)
  {
    for (int j = i + 1; j < n; j++)
    {
      Y.row(k++) = X.row(i) - X.row(j);
    }
  }
  for (int i = 0; i < p - 1; i++)
  {
    for (int j = i + 1; j < p; j++)
    {
      sigmaHat(i, j) = sigmaHat(j, i) = hMeanCov(0.5 * Y.col(i) % Y.col(j), n, p, N, rhs2);
    }
  }
  return Rcpp::List::create(Rcpp::Named("means") = mu, Rcpp::Named("cov") = sigmaHat);
}

// [[Rcpp::export]]
double mad(const arma::vec &x)
{
  return 1.482602 * arma::median(arma::abs(x - arma::median(x)));
}

// [[Rcpp::export]]
arma::mat standardize(arma::mat X, const arma::rowvec &mx, const arma::vec &sx, const int p)
{
  for (int i = 0; i < p; i++)
  {
    X.col(i) = (X.col(i) - mx(i)) / sx(i);
  }
  return X;
}

// [[Rcpp::export]]
void updateHuber(const arma::mat &Z, const arma::vec &res, arma::vec &der, arma::vec &grad, const int n, const double tau, const double n1)
{
  for (int i = 0; i < n; i++)
  {
    double cur = res(i);
    if (std::abs(cur) <= tau)
    {
      der(i) = -cur;
    }
    else
    {
      der(i) = -tau * sgn(cur);
    }
  }
  grad = n1 * Z.t() * der;
}

// [[Rcpp::export]]
arma::vec adaHuberReg(const arma::mat& X, arma::vec Y, const double tol = 0.0001, const int iteMax = 5000)
{
  int n = X.n_rows;
  int p = X.n_cols;
  const double n1 = 1.0 / n;
  double rhs = n1 * (p + std::log(n * p));
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx = arma::stddev(X, 0, 0).t();
  double my = arma::mean(Y);
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx, p));
  Y -= my;
  double tau = 1.345 * mad(Y);
  arma::vec der(n);
  arma::vec gradOld(p + 1), gradNew(p + 1);
  updateHuber(Z, Y, der, gradOld, n, tau, n1);
  arma::vec beta = -gradOld, betaDiff = -gradOld;
  arma::vec res = Y - Z * beta;
  arma::vec resSq = arma::square(res);
  tau = std::sqrt((long double)rootg1(resSq, n, rhs, arma::min(resSq), arma::accu(resSq)));
  updateHuber(Z, res, der, gradNew, n, tau, n1);
  arma::vec gradDiff = gradNew - gradOld;
  int ite = 1;
  while (arma::norm(gradNew, "inf") > tol && ite <= iteMax)
  {
    double alpha = 1.0;
    double cross = arma::as_scalar(betaDiff.t() * gradDiff);
    if (cross > 0)
    {
      double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
      double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
      alpha = std::min(std::min(a1, a2), 100.0);
    }
    gradOld = gradNew;
    betaDiff = -alpha * gradNew;
    beta += betaDiff;
    res -= Z * betaDiff;
    resSq = arma::square(res);
    tau = std::sqrt((long double)rootg1(resSq, n, rhs, arma::min(resSq), arma::accu(resSq)));
    updateHuber(Z, res, der, gradNew, n, tau, n1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  beta.rows(1, p) /= sx;
  beta(0) = huberMean(Y + my - X * beta.rows(1, p));
  return beta;
}

// [[Rcpp::export]]
arma::vec huberReg(const arma::mat& X, arma::vec Y, const double tol = 0.0001, const double constTau = 1.345, const int iteMax = 5000)
{
  int n = X.n_rows;
  int p = X.n_cols;
  const double n1 = 1.0 / n;
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx = arma::stddev(X, 0, 0).t();
  double my = arma::mean(Y);
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx, p));
  Y -= my;
  double tau = constTau * mad(Y);
  arma::vec der(n);
  arma::vec gradOld(p + 1), gradNew(p + 1);
  updateHuber(Z, Y, der, gradOld, n, tau, n1);
  arma::vec beta = -gradOld, betaDiff = -gradOld;
  arma::vec res = Y - Z * beta;
  tau = constTau * mad(res);
  updateHuber(Z, res, der, gradNew, n, tau, n1);
  arma::vec gradDiff = gradNew - gradOld;
  int ite = 1;
  while (arma::norm(gradNew, "inf") > tol && ite <= iteMax)
  {
    double alpha = 1.0;
    double cross = arma::as_scalar(betaDiff.t() * gradDiff);
    if (cross > 0)
    {
      double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
      double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
      alpha = std::min(std::min(a1, a2), 100.0);
    }
    gradOld = gradNew;
    betaDiff = -alpha * gradNew;
    beta += betaDiff;
    res -= Z * betaDiff;
    tau = constTau * mad(res);
    updateHuber(Z, res, der, gradNew, n, tau, n1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  beta.rows(1, p) /= sx;
  beta(0) = huberMean(Y + my - X * beta.rows(1, p));
  
  return beta;
}

              /* ------------------------------------ */    
/*------------ Alternating Gradient Descent Functions ------------*/
              /* ------------------------------------ */    

// [[Rcpp::export]]
double LnVal(const arma::vec &Y, double mu, double tau, const int n, double z)
{
  return (arma::accu(arma::sqrt(tau * tau + arma::square(Y - mu)) - tau) / (z * std::sqrt(n)) + z * tau / std::sqrt(n)) / n;
}

// [[Rcpp::export]]
double gradientMuVal(const arma::vec &Y, double mu, double tau, const int n, double z)
{
  return -(1 / std::sqrt(n)) * arma::accu((Y - mu) / (z * arma::sqrt(tau * tau + arma::square(Y - mu)))) / n;
}

// [[Rcpp::export]]
double gradientTauVal(const arma::vec &Y, double mu, double tau, const int n, double z)
{
  return (1 / std::sqrt(n)) * arma::accu(tau / (z * arma::sqrt(tau * tau + arma::square(Y - mu)))) / n - (std::sqrt(n) / z - z / std::sqrt(n)) / n;
}

// [[Rcpp::export]]
double agd(const arma::vec &Y, double epsilon = 1e-5, int iteMax = 5000)
{
  int n = Y.n_elem;
  double eta1 = 1.0;
  double eta2 = std::sqrt(n);
  double muOld = std::numeric_limits<double>::infinity();
  double muNew = arma::mean(Y);
  double tauOld = std::numeric_limits<double>::infinity();
  double delta = 0.05;
  double z = std::min(std::sqrt(std::log(n)), std::sqrt(std::log(1 / delta)));
  double tauNew = arma::stddev(Y) * std::sqrt(n) / (std::sqrt(2) * z);
  int iteNum = 0;
  double LMuOld = std::numeric_limits<double>::infinity();
  double LMuNew = LnVal(Y, muNew, tauNew, n, z);
  
  while (((LMuOld - LMuNew > 1e-10) && std::abs(muOld - muNew) > epsilon) && std::abs(tauOld - tauNew) > epsilon && iteNum < iteMax)
  {
    muOld = muNew;
    tauOld = tauNew;
    LMuOld = LMuNew;
    
    double gradientLnTauVal = gradientTauVal(Y, muOld, tauOld, n, z);
    tauNew = tauOld - eta2 * gradientLnTauVal;
    
    double gradientLnMuVal = gradientMuVal(Y, muOld, tauNew, n, z);
    muNew = muOld - eta1 * gradientLnMuVal;
    LMuNew = LnVal(Y, muNew, tauNew, n, z);
    iteNum += 1;
  }
  return muNew;
}

// [[Rcpp::export]]
double agdBB(const arma::vec &Y, double epsilon = 1e-5, int iteMax = 5000)
{
  int n = Y.n_elem;
  double eta1 = 1.0;
  double eta2 = std::sqrt(n);
  double muOld = std::numeric_limits<double>::infinity();
  double muNew = arma::mean(Y);
  double tauOld = std::numeric_limits<double>::infinity();
  double delta = 0.05;
  double z = std::min(std::sqrt(std::log(n)), std::sqrt(std::log(1 / delta)));
  double tauNew = arma::stddev(Y) * std::sqrt(n) / (std::sqrt(2) * z);
  int iteNum = 0;
  
  double stepSizeTau = eta2;
  double stepSizeMu = eta1;
  
  double LMuOld = std::numeric_limits<double>::infinity();
  double LMuNew = LnVal(Y, muNew, tauNew, n, z);
  
  double lastTauTilde;
  double lastGradTau;
  double lastMuTilde;
  double lastGradMu;
  
  while (((LMuOld - LMuNew) > epsilon * epsilon && std::abs(muOld - muNew) > epsilon) && std::abs(tauOld - tauNew) > epsilon && iteNum < iteMax)
  {
    LMuOld = LMuNew;
    muOld = muNew;
    tauOld = tauNew;
    
    double gradTau = gradientTauVal(Y, muOld, tauOld, n, z);
    if (iteNum > 0)
    {
      double sTau = tauOld - lastTauTilde;
      double yTau = gradTau - lastGradTau;
      stepSizeTau = std::max(eta2, std::min(std::min((sTau * yTau) / (yTau * yTau), (sTau * sTau) / std::max(0.000001, sTau * yTau)), 1e3));
    }
    lastGradTau = gradTau;
    lastTauTilde = tauOld;
    
    tauNew = std::max(1.35, tauOld - stepSizeTau * (gradientTauVal(Y, muOld, tauOld, n, z)));
    
    double gradMu = gradientMuVal(Y, muOld, tauNew, n, z);
    if (iteNum > 0)
    {
      double sMu = muOld - lastMuTilde;
      double yMu = gradMu - lastGradMu;
      stepSizeMu = std::max(eta1, std::min(std::min((sMu * yMu) / (yMu * yMu), (sMu * sMu) / std::max(0.000001, sMu * yMu)), 1e3));
    }
    lastGradMu = gradMu;
    lastMuTilde = muOld;
    
    muNew = muOld - stepSizeMu * (gradientMuVal(Y, muOld, tauNew, n, z));
    LMuNew = LnVal(Y, muNew, tauNew, n, z);
    iteNum += 1;
  }
  return muNew;
}

// [[Rcpp::export]]
double agdBacktracking(const arma::vec &Y, double s1 = 1.0, double gamma1 = 0.5, double beta1 = 0.8, double s2 = 1.0, 
                                 double gamma2 = 0.5, double beta2 = 0.8, double epsilon = 1e-5, int iteMax = 5000)
{    
  int n = Y.n_elem;
  double muOld = std::numeric_limits<double>::infinity();
  double muNew = arma::mean(Y);
  double tauOld = std::numeric_limits<double>::infinity();
  double delta = 0.05;
  double z = std::min(std::sqrt(std::log(n)), std::sqrt(std::log(1 / delta)));
  double tauNew = arma::stddev(Y) * std::sqrt(n) / (std::sqrt(2) * z);
  int iteNum = 0;
  
  while (std::abs(muOld - muNew) > epsilon && std::abs(tauOld - tauNew) > epsilon && iteNum < iteMax)
  {
    muOld = muNew;
    tauOld = tauNew;
    double LTauOld = LnVal(Y, muOld, tauOld, n, z);
    double gradientLnTauVal = gradientTauVal(Y, muOld, tauOld, n, z);
    
    double eta2 = s2;
    tauNew = tauOld - eta2*gradientLnTauVal;
    double LTauNew = LnVal(Y, muOld, tauNew, n ,z);
    
    while (LTauNew > LTauOld + gamma2*eta2*gradientLnTauVal)
    {
      eta2 = beta2*eta2;
      tauNew = tauOld - eta2*gradientLnTauVal;
      LTauNew = LnVal(Y, muOld, tauNew, n ,z);
      gradientLnTauVal = gradientTauVal(Y, muOld, tauNew, n, z);
    }
    
    double LMuOld = LTauNew;
    double gradientLnMuVal = gradientMuVal(Y, muOld, tauNew, n, z);
    
    double eta1 = s1;
    muNew = muOld - eta1*gradientLnMuVal;
    double LMuNew = LnVal(Y, muNew, tauNew, n ,z);
    
    while (LMuNew > LMuOld + gamma1*eta1*gradientLnMuVal)
    {
      eta1 = beta1*eta1;
      muNew = muOld - eta1*gradientLnMuVal;
      LMuNew = LnVal(Y, muNew, tauNew, n ,z);
      gradientLnMuVal = gradientMuVal(Y, muNew, tauNew, n, z);
    }
    iteNum += 1;
  }
  return muNew;
}
