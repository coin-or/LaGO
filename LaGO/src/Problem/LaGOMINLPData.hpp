// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOMINLPDATA_H_
#define LAGOMINLPDATA_H_

#include "LaGObase.hpp"
#include "LaGOBlockFunction.hpp"
#include "LaGOSparsity.hpp"
#include "LaGOCurvature.hpp"

namespace LaGO {

class GamsReader;
class Decomposition;
class CurvatureCheck;

/** Storage for the data of a MINLP.
 */
class MINLPData : public ReferencedObject {
	friend class GamsReader;
	friend class Decomposition;
	friend class CurvatureCheck;
public:
	/** Storage for the data of a variable.
	 */
	class Variable : public ReferencedObject {
		friend class GamsReader;
	private:
		/** Index of variable in problem.
		 */
		int index;

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
		 */
		bool discrete;

	public:
		Variable(int index_, double lower_, double upper_, bool discrete_=false, const string& name_=string())
		: index(index_), name(name_), lower(lower_), upper(upper_), discrete(discrete_)
		{ }
		
		bool isDiscrete() const { return discrete; }
		
		double getLower() const { return lower; }
		double getUpper() const { return upper; }
		
		const string& getName() const { return name; }
		
		friend ostream& operator<<(ostream& out, const Variable& var);
	};
	
	/** Storage for the data of an objective function or constraint.
	 */
	class ObjCon : public ReferencedObject {
		friend class GamsReader;
		friend class Decomposition;
		friend class CurvatureCheck;
	public:
		/** Name of the objective or constraint.
		 */
		string name;
		
		/** The original function (before decomposition).
		 */
		SmartPtr<Function> origfuncNL;
		SmartPtr<SparseVector> origfuncLin;
		double origfuncConstant;
		
		/** Stores the sparsity graph of origfuncNL.
		 */
		SmartPtr<SparsityGraph> sparsitygraph;
		
		/** The function in decomposed form.
		 */
		vector<SmartPtr<BlockFunction> > decompfuncNL;
		SmartPtr<SparseVector> decompfuncLin;
		double decompfuncConstant;
		/** Tells for each variable in which blocks at which position it appears.
		 * for (k,j) in decompmapping[i]: decompfuncNL[k]->indices[j]==i
		 */ 
		vector<vector<pair<int,int> > > decompmapping;

	public:
		ObjCon(const SmartPtr<Function>& origfuncNL_=NULL, const SmartPtr<SparseVector>& origfuncLin_=NULL, double origfuncConstant_=0., const string& name_=string());
		
		virtual ~ObjCon();
		
		bool isLinear() const { return IsNull(origfuncNL); }
		bool isConstant() const { return IsNull(origfuncLin) && IsNull(origfuncNL); }
		Curvature getCurvature() const;
		
		void print(ostream& out, const vector<MINLPData::Variable>& var) const;

		friend ostream& operator<<(ostream& out, const ObjCon& objcon);
	};

	/** Storage for the data of a constraint.
	 * Constraints are of the form g(x)<=0 or g(x)==0.
	 */	
	class Constraint : public ObjCon {
	public:
		/** Index of constraint in problem.
		 */
		int index;
		/** Lower and upper bounds on the constraint.
		 */
		double lower, upper;
		
		Constraint(int index_, double lower_=0., double upper_=0., 
		  const SmartPtr<Function>& origfuncNL_=NULL, const SmartPtr<SparseVector>& origfuncLin_=NULL, double origfuncConstant_=0.,		
		  const string& name_=string())
		: ObjCon(origfuncNL_, origfuncLin_, origfuncConstant_, name_), index(index_), lower(lower_), upper(upper_)
		{ }

		void print(ostream& out, const vector<MINLPData::Variable>& var) const;

		friend ostream& operator<<(ostream& out, const Constraint& con);
	};
	
	/** Storage for the data of an objective function.
	 */
	class Objective : public ObjCon {
	public:
		Objective(const SmartPtr<Function>& origfuncNL_=NULL, const SmartPtr<SparseVector>& origfuncLin_=NULL, double origfuncConstant_=0.,		
		  const string& name_=string())
		: ObjCon(origfuncNL_, origfuncLin_, origfuncConstant_, name_)
		{ }
	};
	
	
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

	/** Indices of discrete variables.
	 */
	vector<int> discrete_var;

	/** A collection of start points.
	 */
	vector<DenseVector> start_points;

public:
	MINLPData();
	virtual ~MINLPData();
	
	int numVariables() const { return var.size(); }
	int numDiscrVariables() const { return discrete_var.size(); }
	int numConstraints() const { return con.size(); }
	
	const Variable& getVariable(int index) const { return var.at(index); }
	const Constraint& getConstraint(int index) const { return con.at(index); }
	const Objective& getObjective() const { return obj; }
	const vector<int>& getDiscrVariables() const { return discrete_var; }
	const vector<DenseVector>& getStartPoints() const { return start_points; }

	/** Creates vectors with lower and upper bounds of those variables that are listed in indices.
	 */ 	
	void getBox(DenseVector& lower, DenseVector& upper, const vector<int>& indices) const;
	void getBox(DenseVector& lower, DenseVector& upper) const;
	/** Whether the MINLP is convex.
	 * The MINLP is convex, if each constraint by its own is convex.
	 * If the curvature of some constraints is not known, false is returned.  
	 */
	bool isConvex() const;
	
	friend ostream& operator<<(ostream& out, const MINLPData& data);
	
	//TODO: some methods for fast evaluation of jacobian and hessian of lagrangian
}; // class MINLPData

} //namespace LaGO

#endif // LAGOMINLPDATA_H_
