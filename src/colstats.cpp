/******************************************************************************/

#include <bigstatsr/BMCodeAcc.h>

/******************************************************************************/

// [[Rcpp::export]]
ListOf<NumericVector> snp_colstats(Environment BM,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   int ncores) {

  XPtr<FBM> xpBM = BM["address"];
  if (BM.exists("code256")) {
    SubBMCode256Acc macc(xpBM, rowInd, colInd, BM["code256"], 1);
    return snp_colstats0(macc, rowInd, colInd, ncores);
  } else {
    switch(xpBM->matrix_type()) {
    case 6:
    {
      SubBMAcc<float> macc(xpBM, rowInd, colInd, 1);
      return snp_colstats0(macc, rowInd, colInd, ncores);
    }
    default:
      throw Rcpp::exception(ERROR_TYPE);
    }
  }
}


ListOf<NumericVector> snp_colstats0(C macc,
                                   const IntegerVector& rowInd,
                                   const IntegerVector& colInd,
                                   int ncores) {
  size_t n = macc.nrow();
  size_t m = macc.ncol();

  NumericVector sumX(m), denoX(m);

  #pragma omp parallel for num_threads(ncores)
  for (size_t j = 0; j < m; j++) {
    double xSum = 0, xxSum = 0;
    for (size_t i = 0; i < n; i++) {
      double x = macc(i, j);
      xSum += x;
      xxSum += x*x;
    }
    sumX[j] = xSum;
    denoX[j] = xxSum - xSum * xSum / n;
  }

  return List::create(_["sumX"]  = sumX,
                      _["denoX"] = denoX);
}

/******************************************************************************/
