// KDTreeBallQuery.cpp
// [idx, dist] = KDTreeBallQuery(inPts, queryPts, radii)
// Self-contained R2018a C++-API MEX (no external KD-tree dependency)
//
// Inputs:
//   inPts    : MxK double
//   queryPts : NxK double
//   radii    : scalar or Nx1/1xN double
//
// Outputs:
//   idx  : Nx1 cell of index vectors (1-based indices into inPts)
//   dist : Nx1 cell of corresponding Euclidean distances
//
// Compile (R2024b example):
//   clear mex
//   mex -v -R2018a ...
//       CXXFLAGS="$CXXFLAGS -std=c++17 -fno-omit-frame-pointer" ...
//       CXXOPTIMFLAGS="-O3 -DNDEBUG" ...
//       KDTreeBallQuery.cpp
//
// Author: (ChatGPT) for Dr. Han

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"

#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits>

using matlab::data::Array;
using matlab::data::ArrayFactory;
using matlab::data::TypedArray;
using matlab::data::TypedArrayRef;
using matlab::data::CellArray;
using matlab::mex::ArgumentList;
using matlab::mex::Function;

namespace {

// Column-major accessor for MxK matrix stored as MATLAB double
inline double getMK(const double* A, size_t M, size_t /*K*/, size_t r, size_t c) {
    // MATLAB column-major: element (r,c) sits at c*M + r
    return A[c * M + r];
}

// Validate a 2-D double matrix
void ensureDouble2D(const Array& a, const char* name) {
    if (a.getType() != matlab::data::ArrayType::DOUBLE || a.getNumberOfElements() == 0) {
        throw std::invalid_argument(std::string(name) + " must be non-empty double.");
    }
    const auto dims = a.getDimensions();
    if (dims.size() != 2) {
        throw std::invalid_argument(std::string(name) + " must be 2-D.");
    }
}

// Extract radii for each query (N). Accepts scalar, Nx1, or 1xN.
std::vector<double> makeRadiiPerQuery(const TypedArray<double>& radiiArr, size_t N) {
    std::vector<double> radii(N, 0.0);
    const auto dims = radiiArr.getDimensions();

    if (radiiArr.getNumberOfElements() == 1) {
        // scalar radius
        double r = *radiiArr.begin();
        if (!(r >= 0.0)) { // NaN or negative
            throw std::invalid_argument("radii must be nonnegative.");
        }
        std::fill(radii.begin(), radii.end(), r);
        return radii;
    }

    if (dims.size() != 2) {
        throw std::invalid_argument("radii must be scalar or a 2-D vector of size Nx1 or 1xN.");
    }

    size_t R = dims[0];
    size_t C = dims[1];

    if (R == N && C == 1) {
        size_t i = 0;
        for (double v : radiiArr) {
            if (!(v >= 0.0)) throw std::invalid_argument("radii must be nonnegative.");
            radii[i++] = v;
        }
    } else if (R == 1 && C == N) {
        // row vector
        size_t i = 0;
        for (double v : radiiArr) {
            if (!(v >= 0.0)) throw std::invalid_argument("radii must be nonnegative.");
            radii[i++] = v;
        }
    } else {
        throw std::invalid_argument("radii must be scalar or match the number of query points (N).");
    }
    return radii;
}

} // namespace

class MexFunction : public Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) override {
        ArrayFactory factory;

        // --- Validate number of args
        if (inputs.size() != 3) {
            error_("KDTreeBallQuery requires 3 inputs: inPts, queryPts, radii.");
        }
        if (outputs.size() > 2) {
            error_("KDTreeBallQuery returns at most two outputs: [idx, dist].");
        }

        // --- Validate and extract inputs
        // inPts
        if (inputs[0].getType() != matlab::data::ArrayType::DOUBLE) {
            error_("inPts must be double.");
        }
        if (inputs[1].getType() != matlab::data::ArrayType::DOUBLE) {
            error_("queryPts must be double.");
        }
        if (inputs[2].getType() != matlab::data::ArrayType::DOUBLE) {
            error_("radii must be double.");
        }

        TypedArray<double> inPts   = inputs[0];
        TypedArray<double> queryPts= inputs[1];
        TypedArray<double> radiiIn = inputs[2];

        ensureDouble2D(inPts,   "inPts");
        ensureDouble2D(queryPts,"queryPts");

        auto dimsIn  = inPts.getDimensions();    // M x K
        auto dimsQ   = queryPts.getDimensions(); // N x K

        size_t M = dimsIn[0];
        size_t K = dimsIn[1];
        size_t N = dimsQ[0];
        size_t Kq= dimsQ[1];

        if (K == 0 || Kq == 0) {
            error_("inPts/queryPts must have nonzero number of columns.");
        }
        if (K != Kq) {
            error_("Dimensionality mismatch: size(inPts,2) must equal size(queryPts,2).");
        }
        if (K > 3) {
            error_("Only 1D/2D/3D point sets are supported (K must be 1..3).");
        }

        // Copy raw data buffers (column-major)
        std::vector<double> inBuf;
        inBuf.reserve(M * K);
        for (double v : inPts) inBuf.push_back(v);

        std::vector<double> qBuf;
        qBuf.reserve(N * K);
        for (double v : queryPts) qBuf.push_back(v);

        // Radii per query
        std::vector<double> radii = makeRadiiPerQuery(radiiIn, N);

        // --- Prepare outputs (N x 1 cell arrays)
        CellArray idxCell  = factory.createCellArray({ N, 1 });
        CellArray distCell = factory.createCellArray({ N, 1 });

        // --- Brute-force radius search
        // For each query i, compute distances to all M input points, keep those <= radii[i]
        for (size_t i = 0; i < N; ++i) {
            double r  = radii[i];
            double r2 = r * r;

            // First pass: count matches to allocate exactly
            size_t count = 0;
            for (size_t m = 0; m < M; ++m) {
                double d2 = 0.0;
                // accumulate squared distance across K dims
                for (size_t c = 0; c < K; ++c) {
                    double q  = getMK(qBuf.data(),  N, K, i, c);
                    double ip = getMK(inBuf.data(), M, K, m, c);
                    double diff = ip - q;
                    d2 += diff * diff;
                }
                if (d2 <= r2) ++count;
            }

            // Allocate MATLAB vectors (column) for indices and distances
            TypedArray<double> idxVec  = factory.createArray<double>({count, 1});
            TypedArray<double> dstVec  = factory.createArray<double>({count, 1});

            // Second pass: fill
            size_t w = 0;
            for (size_t m = 0; m < M; ++m) {
                double d2 = 0.0;
                for (size_t c = 0; c < K; ++c) {
                    double q  = getMK(qBuf.data(),  N, K, i, c);
                    double ip = getMK(inBuf.data(), M, K, m, c);
                    double diff = ip - q;
                    d2 += diff * diff;
                }
                if (d2 <= r2) {
                    // MATLAB uses 1-based indices
                    idxVec[w] = static_cast<double>(m + 1);
                    dstVec[w] = std::sqrt(d2);
                    ++w;
                }
            }

            idxCell[i]  = std::move(idxVec);
            distCell[i] = std::move(dstVec);
        }

        // --- Set outputs
        outputs[0] = std::move(idxCell);
        if (outputs.size() > 1) {
            outputs[1] = std::move(distCell);
        }
    }

private:
    // convenience to throw MATLAB error
    [[noreturn]] void error_(const char* msg) {
        getEngine()->feval(u"error", 1, std::vector<Array>{ ArrayFactory().createScalar(msg) });
        // feval(error) above throws into MATLAB; noreturn to satisfy compiler
        throw std::runtime_error(msg);
    }
};
// MATLAB_MEX_FUNCTION_CLASS_EXPORT(MexFunction);