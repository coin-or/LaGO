// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id: LaGOAlgorithm.hpp 135 2007-10-09 19:37:41Z stefan $

#include "LaGOQuadraticRelaxTest.hpp"

#include "LaGOQuadraticOrConvexApproximation.hpp"

#include "BonBonminSetup.hpp"
#include "BonCbc.hpp"

namespace LaGO {

QuadraticRelaxTest::QuadraticRelaxTest(MINLPData& data_)
: data(data_), decomp(data_), curvcheck(data_), conprob(data_), guess(data_), convexify(data_)
{
#ifdef COIN_HAS_FILIB
	filib::fp_traits<double>::setup();
#endif
}

void QuadraticRelaxTest::preprocessing() {
	decomp.decompose();
	
//	conprob.print_level=5;
	conprob.initDependencyGraph();
	BoxReductionStatistics statistics=conprob.reduceBox();
	cout << statistics;
	if (statistics.empty_box) {
		cout << "Problem infeasible." << endl;
		return;
	}
//	BoxReductionStatistics::printBox(cout, data);
	
	guess.guessBounds();
	cout << "Variables with guessed bounds: " << guess.numGuessedBounds() << endl; 
	
	curvcheck.computeCurvature();
	
	quadest.computeEstimators(data);
	
	if (guess.numGuessedBounds())
		convexify.convexify(guess.bound_is_guessed);
	else 
		convexify.convexify();
	

//	cout << data;
}

void QuadraticRelaxTest::solve_relax(SmartPtr<QuadraticOrConvexApproximation> quad) {
	Bonmin::BonminSetup bonmin_setup;
	bonmin_setup.initialize(GetRawPtr(quad));

	Bonmin::Bab bb;
	bb(bonmin_setup);	
}

	
void QuadraticRelaxTest::run() {
	preprocessing();
	
	SmartPtr<QuadraticOrConvexApproximation> quad;
	
	int nr_new_estimators;
	int iter=0;
	do {
		quad=new QuadraticOrConvexApproximation(data, true);
//		cout << *quad;
		
		solve_relax(quad);
		
		Bonmin::TMINLP::SolverReturn status=quad->getSolutionStatus();
		if (status!=Bonmin::TMINLP::SUCCESS) return;

		DenseVector x(data.numVariables());
		quad->getSolution(x);
		cout << "Infeasibilities: " << endl;
		for (int c=0; c<data.numConstraints(); ++c) {
			double infeas=data.getConstraint(c).getInfeasibility(x);
			if (infeas)
				cout << data.getConstraint(c).name << ":\t " << infeas << endl;
		}
		cout << "Objective function value: " << data.getObjective().eval(x) << endl;
		
		nr_new_estimators=quadest.computeImprovingEstimators(data, x);

		if (guess.numGuessedBounds())
			convexify.convexifyEstimators(guess.bound_is_guessed);
		else 
			convexify.convexifyEstimators();

	} while (iter++<1);
		
}

} // namespace LaGO
