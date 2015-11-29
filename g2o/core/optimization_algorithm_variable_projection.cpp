// g2o - General Graph Optimization
// Copyright (C) 2011 R. Kuemmerle, G. Grisetti, W. Burgard
// All rights reserved.
//
// Copyright (C) 2015 K. Khosoussi, S. Huang and G. Dissanayake
// Centre for Autonomous Systems (CAS)
// University of Technology Sydney
// (https://kasra.github.io)
// See "Exploiting the Separable Structure of SLAM" @ RSS 2015.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "optimization_algorithm_variable_projection.h"

#include <iostream>

#include "g2o/stuff/timeutil.h"
#include "g2o/stuff/macros.h"

#include "solver.h"
#include "batch_stats.h"
#include "sparse_optimizer.h"

using namespace std;

// Projection threshold - should be in (0,1)
// (>1) means projection happens in every iteration
#define GAMMA_THRESHOLD 2
// is the position-position block of Omega spherical?
#define SPHERICAL false

namespace g2o {

  OptimizationAlgorithmVP::OptimizationAlgorithmVP(Solver* solver) :
    OptimizationAlgorithmWithHessian(solver),
    _projectNextIter(true),
    _gammaThreshold(GAMMA_THRESHOLD)
  {
  }

  OptimizationAlgorithmVP::~OptimizationAlgorithmVP()
  {
  }

  OptimizationAlgorithm::SolverResult OptimizationAlgorithmVP::solve(int iteration, bool online)
  {
    assert(_optimizer && "_optimizer not set");
    assert(_solver->optimizer() == _optimizer && "underlying linear solver operates on different graph");
    bool ok = true;
    
    //here so that correct component for max-mixtures can be computed before the build structure
    double t=get_monotonic_time();
    _optimizer->computeActiveErrors();
    G2OBatchStatistics* globalStats = G2OBatchStatistics::globalStats();
    if (globalStats) {
      globalStats->timeResiduals = get_monotonic_time()-t;
    }

    double preProjCost = 0;

    if (_projectNextIter) preProjCost = _optimizer->activeRobustChi2();

    if (iteration == 0 && !online) { // built up the CCS structure, here due to easy time measure
      if (_projectNextIter) {
        ok = _solver->buildStructureLinear();
        assert(ok);
      }
      ok = _solver->buildStructure();
      if (! ok) {
        cerr << __PRETTY_FUNCTION__ << ": Failure while building CCS structure" << endl;
        return OptimizationAlgorithm::Fail;
      }
    }

    // is the position-position block of omega spherical?
    bool isSpherical = SPHERICAL;

    if (iteration == 0)
      cerr << "[VP]: isSpherical = " << isSpherical << " - Gamma_t = " << _gammaThreshold << "\n";

    if (_projectNextIter) {
      if (iteration > 0 && isSpherical == true) {
        _solver->buildSystemSpherical();
        _solver->project(isSpherical);
        _optimizer->updateLinear(_solver->xLinear());
        _optimizer->computeActiveErrors();
      } else { // not spherical or iteration = 0
        _solver->buildSystemLinear();
        _solver->project();
        _optimizer->updateLinear(_solver->xLinear());
        _optimizer->computeActiveErrors();
        assert(preProjCost);

        double postProjCost = _optimizer->activeRobustChi2();
        cerr << "[VP]: Cost Before Proj = " << preProjCost << " -> Cost After Proj = " << postProjCost << "\n";

        /* gamma the gain */
        double gamma = postProjCost/preProjCost;

        cerr << "[VP]: Gamma = " << FIXED(gamma) << endl;

        if (gamma > _gammaThreshold && isSpherical == false) { // no more projection step
          _projectNextIter = false;
          cerr << "[VP]: _projNextIter : true -> false" << endl;
        }
      }
    }

    t=get_monotonic_time();
    _solver->buildSystem();
    if (globalStats) {
      globalStats->timeQuadraticForm = get_monotonic_time()-t;
      t=get_monotonic_time();
    }

    ok = _solver->solve();


    if (globalStats) {
      globalStats->timeLinearSolution = get_monotonic_time()-t;
      t=get_monotonic_time();
    }

    _optimizer->update(_solver->x());


    if (globalStats) {
      globalStats->timeUpdate = get_monotonic_time()-t;
    }
    if (ok)
      return OK;
    else
      return Fail;
  }

  void OptimizationAlgorithmVP::printVerbose(std::ostream& os) const
  {
    os
      << "\t schur= " << _solver->schur();
  }

} // end namespace
