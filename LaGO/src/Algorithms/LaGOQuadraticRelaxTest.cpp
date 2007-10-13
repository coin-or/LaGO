// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id: LaGOAlgorithm.hpp 135 2007-10-09 19:37:41Z stefan $

#include "LaGOQuadraticRelaxTest.hpp"

#include "LaGOQuadraticOrConvexApproximation.hpp"

#include "BonBonminSetup.hpp"
#include "BonCbc.hpp"
#include "BonTMINLP2TNLP.hpp"

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
	
	int count=quadest.computeEstimators(data);
	cout << "Number of quadratic estimators computed: " << count << endl; 
	
	if (guess.numGuessedBounds())
		convexify.convexify(guess.bound_is_guessed);
	else 
		convexify.convexify();
	

//	cout << data;
}

bool QuadraticRelaxTest::solveMINLPRelax(SmartPtr<QuadraticOrConvexApproximation> quad) {
	Bonmin::BonminSetup bonmin_setup;
	bonmin_setup.initialize(GetRawPtr(quad));

	Bonmin::Bab bb;
	bb(bonmin_setup);	

	Bonmin::TMINLP::SolverReturn status=quad->getSolutionStatus();
	switch (status) {
		case Bonmin::TMINLP::LIMIT_EXCEEDED:
			clog << "Bonmin: some limit exceeded." << endl;
			if (!bb.bestSolution())
				return false;
		case Bonmin::TMINLP::SUCCESS:
			clog << "Bonmin best bound: " << bb.bestBound() << endl;
			return true;
		default:
			clog << "Bonmin failed. Status: " << status << endl;
			return false;
	}
}

bool QuadraticRelaxTest::solveNLPRelax(SmartPtr<QuadraticOrConvexApproximation> quad) {
	SmartPtr<Ipopt::TNLP> nlp=new Bonmin::TMINLP2TNLP(GetRawPtr(quad));
	
	Ipopt::IpoptApplication ipopt;
//	ipopt.Options()->SetIntegerValue("print_level", 4);
//	ipopt.Options()->SetIntegerValue("max_iter", 10000);
//	ipopt.Initialize("");
	ipopt.Initialize(); // this reads ipopt.opt
	
	Ipopt::ApplicationReturnStatus status=ipopt.OptimizeTNLP(GetRawPtr(nlp));
	if (status!=Ipopt::Solve_Succeeded)
		clog << "Ipopt return: " << status << endl;
	return (status==Ipopt::Solve_Succeeded);
}

	
void QuadraticRelaxTest::run() {
	preprocessing();
	
	quadest.testEstimators(data);
	
	SmartPtr<QuadraticOrConvexApproximation> quad;
	
	int nr_new_estimators;
	int iter=0;
	do {
		quadest.testEstimators(data);
		
		quad=new QuadraticOrConvexApproximation(data, true);
		
		bool success=solveMINLPRelax(quad);
		cout << "Solved: " << success << endl;
		if (!success) {
//			cout << *quad;
			break;
		} 

		DenseVector x(data.numVariables());
		quad->getSolution(x);
		cout << "Infeasibilities w.r.t. original problem: " << endl;
		for (int c=0; c<data.numConstraints(); ++c) {
			double infeas=data.getConstraint(c).getInfeasibility(x);
			if (infeas>1e-4)
				cout << data.getConstraint(c).name << ":\t " << infeas << endl;
		}
		cout << "Iteration " << iter << "\t Relaxation optimal value: " << quad->getSolutionValue() << endl;
		cout << "Iteration " << iter << "\t Original objective function value: " << data.getObjective().eval(x) << endl;
		
		if (++iter>3) {
			cout << "Stop after " << iter << " iterations." << endl;
			break;
		}

		nr_new_estimators=quadest.computeImprovingEstimators(data, x);
		cout << "Iteration " << iter << "\t Number of generated estimators: " << nr_new_estimators << endl; 

		if (nr_new_estimators) {
			if (guess.numGuessedBounds())
				convexify.convexifyEstimators(guess.bound_is_guessed);
			else 
				convexify.convexifyEstimators();
		}
		
	} while (nr_new_estimators);
		
}

} // namespace LaGO
