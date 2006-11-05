// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef _MINLPDATA_H_
#define _MINLPDATA_H_

#include "standard.h"
#include "func.h"
#include "problem.h"

#include <string>

class MINLP;

class MinlpOpt;
class Reformulation;
class LinearizedConCutGenerator;
class PolynomialUnderestimator2;
class QuadraticUnderestimator;

/** Storage for the data of a MINLP.
 */
class MINLPData {
	friend class MINLP;
	friend class MinlpOpt;
	friend class Reformulation;
	friend class LinearizedConCutGenerator;
	friend class PolynomialUnderestimator2;
	friend class QuadraticUnderestimator;
private:
	/** Storage for the data of an objective function or constraint.
	 */
	class ObjCon {
	public:
		/** Name of the objective or constraint.
		 */
		string name;
		/** The function.
		 */
		Pointer<SepQcFunc> func;
		/** Type of the function.
		 */
		SepQcFunc::ftype functype;
		/** Curvature of the function.
		 */
		Func::CurvatureType curvtype;
		
		/** The quadratic underestimator of func, if different from func.
		 */
		Pointer<SepQcFunc> polynomial_underestimator;
		/** Blockwise constants in the polynomial approximation functions.
		 */
		map<int, double> polynomial_approx_constants_lower;
		
		/** Characteristica of convex alpha-underestimator.
		 */
		map<int, double> convexification_characteristica_lower;
		
		/** The convex underestimator of func, if different from func.
		 */
		Pointer<SepQcFunc> convex_underestimator;
		
		/** Tells which block was moved to which constraint in the reformulation.
		 */
		map<int, int> reformulation_constraints_lower;
		
		ObjCon(const Pointer<SepQcFunc>& func_, const string& name_=string());
		
		ObjCon()
		: functype(SepQcFunc::CONSTANT), curvtype(Func::LINEAR)
		{ }
		
		friend ostream& operator<<(ostream& out, const ObjCon& objcon);

		// TODO: information about usage of blocks, sparsity, relaxations, ...
	};
	friend ostream& operator<<(ostream& out, const ObjCon& objcon);
	
public:
	/** Storage for the data of a variable.
	 */
	class Variable {
	public:
		/** Index of variable in problem.
		 */
		int index;
		/** Number of the block where the variable belongs to.
		 */
		int block_nr;
		/** Index of the variable in the block.
		 */
		int index_in_block;

		/** Name of the variable.
		 */		
		string name;
		
		/** Lower bound.
		 */
		double lower;
		/** Upper bound.
		 */
		double upper;
		/** Discrete or not.
		 * TODO: use an enumerate to distinguish between continuous, discrete (lower or upper), and integer.
		 */
		bool discrete;

		Variable(int index_, double lower_, double upper_, bool discrete_=false, const string& name_=string())
		: index(index_), block_nr(-1), index_in_block(-1), name(name_), lower(lower_), upper(upper_), discrete(discrete_)
		{ }

	};

	/** Storage for the data of a constraint.
	 * Constraints are of the form g(x)<=0 or g(x)==0.
	 */	
	class Constraint : public ObjCon {
	public:
		/** Index of constraint in problem.
		 */
		int index;
		/** Equality constraint or not.
		 */
		bool equality;

		/** The quadratic underestimator of -func, if different from func.
		 */
		Pointer<SepQcFunc> polynomial_overestimator;
		/** Blockwise constants in the polynomial approximation functions.
		 */
		map<int, double> polynomial_approx_constants_upper;

		/** Characteristica of concave alpha-overestimator.
		 */
		map<int, double> convexification_characteristica_upper;
		/** The convex underestimator of -func, if needed.
		 */
		Pointer<SepQcFunc> concave_overestimator;

		map<int, int> reformulation_constraints_upper;

		/** If this constraint is an equality and was split in two inequalities, then this is the index of the second inequality constraint.
		 */
//		int ineq_index;
		
		Constraint(int index_, const Pointer<SepQcFunc>& func_, bool equality_, const string& name_=string())
		: ObjCon(func_, name_), index(index_), equality(equality_)
		{ }

		friend ostream& operator<<(ostream& out, const Constraint& con);
	};
	friend ostream& operator<<(ostream& out, const Constraint& con);
	
	/** Storage for the data of an objective function.
	 */
	class Objective : public ObjCon {
	public:
		Objective(const Pointer<SepQcFunc>& func_, const string& name_=string())
		: ObjCon(func_, name_)
		{ }
		
		Objective()
		: ObjCon()
		{ }
	};
	
	typedef vector<int> Block;
	typedef vector<Block> BlockStructure;
	
	
private:
	/** Give a name to the problem.
	 */
	string name;
	/** Variables.
	 */
	vector<Variable> var;
	
	/** Constraints.
	 */
	vector<Constraint> con;
	
	/** Objective.
	 */
	Objective obj;
	
	/** Information on block structure.
	 */
	BlockStructure block;
	
	/** Indices of discrete variables.
	 */
	vector<int> discrete_var;


public:
	MINLPData();
	MINLPData(const MinlpProblem& prob);
	virtual ~MINLPData();
};

#endif //_MINLPDATA_H_
