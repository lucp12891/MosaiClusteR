// MosaiClusteR C++ kernels for the ABC / M-ABC family.
// Compiled on demand by enable_cpp_acceleration() via Rcpp::sourceCpp().
// Original sampler: mcpp_sampling.cpp (IntegraClust); co-clustering kernel added.
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <numeric>

extern "C" {
#include <R.h>
#include <Rmath.h>
}

using namespace Rcpp;

static inline int findInterval01(double u, const double* boundaries, int nBoundaries) {
    int lo = 0;
    int hi = nBoundaries - 1;
    while (lo + 1 < hi) {
        int mid = (lo + hi) >> 1;
        if (u >= boundaries[mid]) lo = mid;
        else                      hi = mid;
    }
    return lo;
}

static void calculateBoundaries(const double* weights, double* boundaries, int populationSize) {
    double acc = 0.0;
    boundaries[0] = 0.0;
    for (int i = 1; i <= populationSize; ++i) { acc += weights[i - 1]; boundaries[i] = acc; }
}

static void removeWeightAndNormalize(double* weights, int indexToRemove, int populationSize) {
    double wrem = weights[indexToRemove];
    weights[indexToRemove] = 0.0;
    double remain = 1.0 - wrem;
    if (remain > 0) for (int i = 0; i < populationSize; ++i) weights[i] /= remain;
}

static void sampleWithReplacementWithWeights(int sampleSize, int populationSize,
    double* weights, int* sampledIndices) {
    RNGScope scope;
    const int nB = populationSize + 1;
    std::vector<double> boundaries(nB);
    calculateBoundaries(weights, boundaries.data(), populationSize);
    for (int i = 0; i < sampleSize; ++i) {
        double u = unif_rand();
        sampledIndices[i] = findInterval01(u, boundaries.data(), nB);
    }
}

static void sampleWithoutReplacementWithWeights(int sampleSize, int populationSize,
    double* weights, int* sampledIndices) {
    RNGScope scope;
    const int nB = populationSize + 1;
    std::vector<double> w(weights, weights + populationSize);
    std::vector<double> boundaries(nB);
    for (int i = 0; i < sampleSize; ++i) {
        calculateBoundaries(w.data(), boundaries.data(), populationSize);
        double u = unif_rand();
        int idx = findInterval01(u, boundaries.data(), nB);
        sampledIndices[i] = idx;
        removeWeightAndNormalize(w.data(), idx, populationSize);
    }
}

static void sampleWithReplacement(int sampleSize, int populationSize, int* sampledIndices) {
    RNGScope scope;
    for (int i = 0; i < sampleSize; ++i) sampledIndices[i] = (int)(unif_rand() * populationSize);
}

static void sampleWithoutReplacement(int sampleSize, int populationSize, int* sampledIndices) {
    RNGScope scope;
    std::vector<int> idx(populationSize);
    std::iota(idx.begin(), idx.end(), 0);
    int last = populationSize - 1;
    for (int i = 0; i < sampleSize; ++i) {
        int k = (int)(unif_rand() * (last + 1));
        sampledIndices[i] = idx[k];
        std::swap(idx[k], idx[last]);
        --last;
    }
}

// [[Rcpp::export]]
IntegerVector mcpp_rf_sample(int populationSize, int sampleSize,
    bool useWeights = false, bool withReplacement = true,
    Nullable<NumericVector> weights = R_NilValue) {
    IntegerVector out(sampleSize);
    if (withReplacement) {
        if (useWeights) {
            NumericVector w = weights.isNull()
                ? NumericVector(populationSize, 1.0 / populationSize) : NumericVector(weights);
            sampleWithReplacementWithWeights(sampleSize, populationSize, w.begin(), out.begin());
        } else sampleWithReplacement(sampleSize, populationSize, out.begin());
    } else {
        if (useWeights) {
            NumericVector w = weights.isNull()
                ? NumericVector(populationSize, 1.0 / populationSize) : NumericVector(weights);
            sampleWithoutReplacementWithWeights(sampleSize, populationSize, w.begin(), out.begin());
        } else sampleWithoutReplacement(sampleSize, populationSize, out.begin());
    }
    return out; // 0-based indices
}

// [[Rcpp::export]]
IntegerVector mcpp_sample_weighted_no_replace(NumericVector z, int ng) {
    int p = z.size();
    IntegerVector ord = Rcpp::seq(0, p - 1);
    std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) { return z[a] > z[b]; });
    NumericVector w(p);
    for (int r = 0; r < p; ++r) w[ord[r]] = 1.0 / ((r + 1) + 100.0);
    double sumw = std::accumulate(w.begin(), w.end(), 0.0);
    if (sumw <= 0) sumw = 1.0;
    for (int i = 0; i < p; ++i) w[i] /= sumw;
    ng = std::min(ng, p);
    IntegerVector out(ng);
    std::vector<double> ww(w.begin(), w.end());
    for (int k = 0; k < ng; ++k) {
        RNGScope scope;
        double u = unif_rand();
        double acc = 0.0; int pick = p - 1;
        for (int i = 0; i < p; ++i) { acc += ww[i]; if (u <= acc) { pick = i; break; } }
        out[k] = pick;
        double rem = ww[pick];
        ww[pick] = 0.0;
        double remain = 1.0 - rem;
        if (remain > 0) for (int i = 0; i < p; ++i) ww[i] /= remain;
    }
    return out; // 0-based
}

// Consensus co-clustering counts over a (K*R) x N label matrix (0 = unselected).
// [[Rcpp::export]]
List count_coclustering_cpp(IntegerMatrix res) {
    int R = res.nrow();
    int N = res.ncol();
    NumericMatrix co_sel(N, N);
    NumericMatrix not_co_clust(N, N);
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            int cs = 0, ncc = 0;
            for (int r = 0; r < R; r++) {
                int li = res(r, i);
                int lj = res(r, j);
                if (li != 0 && lj != 0) { cs++; if (li != lj) ncc++; }
            }
            co_sel(i, j) = cs;  co_sel(j, i) = cs;
            not_co_clust(i, j) = ncc; not_co_clust(j, i) = ncc;
        }
    }
    return List::create(Named("co_sel") = co_sel, Named("not_co_clust") = not_co_clust);
}
