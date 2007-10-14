// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id: LaGOAlgorithm.hpp 135 2007-10-09 19:37:41Z stefan $

#ifndef LAGOGAMSSOLVER_HPP_
#define LAGOGAMSSOLVER_HPP_

#include "LaGObase.hpp"
#include "LaGOGamsReader.hpp"

namespace LaGO {

/** Class that calls a GAMS (NLP) solver on the original problem with fixed discrete variables.
 */  
class GamsSolver : public ReferencedObject {
private:
	const GamsReader& gams;
	const vector<int>& discrete_var;

	DenseVector lower;
	DenseVector upper;
	DenseVector startpoint;
	DenseVector discr_fix;

	void* iolibsave;

	bool first_call;
	char* tmpsolfn; // temporary space to store solutionfilename
	char** args; // arguments for subsolver-call

	int subsolver;
	char* solvername;
	int optfile;
	
	int model_status;
	int solver_status;
	int iteration_number;
	double time;
	double solution_value;
	DenseVector solution_primal;
	DenseVector solution_conval;
	DenseVector solution_dualcon;
	DenseVector solution_dualvar;
	
public:
	GamsSolver(const GamsReader& gams_, const vector<int>& discrete_var_);
	
	~GamsSolver();
	
	void setStartPoint(const DenseVector& x);
	
	void setBounds(const DenseVector& lower, const DenseVector& upper);
	
	/** Calls a GAMS solver.
	 * @return True on success, and false on a failure.
	 */ 
	bool callSolver();
	
	int getModelStatus() const { return model_status; }
	int getSolverStatus() const { return solver_status; }
	int getNumIterations() const { return iteration_number; }
	double getTime() const { return time; }
	
	double getOptimalValue() const { return solution_value; }  
	const DenseVector& getOptimalPoint() const { return solution_primal; }
	
}; // class GamsSolver

} // namespace LaGO

#endif /*LAGOGAMSSOLVER_HPP_*/
