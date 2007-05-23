// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOMINLPDATA_H_
#define LAGOMINLPDATA_H_

#include "LaGObase.hpp"
#include "LaGOBlockFunction.hpp"

namespace LaGO {

class GamsReader;

/** Storage for the data of a MINLP.
 */
class MINLPData : public ReferencedObject {
	friend class GamsReader;
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
		
		friend ostream& operator<<(ostream& out, const Variable& var);
	};
	
	/** Storage for the data of an objective function or constraint.
	 */
	class ObjCon : public ReferencedObject {
		friend class GamsReader;
	protected:
		/** Name of the objective or constraint.
		 */
		string name;
		
		/** The original function (before decomposition).
		 */
		SmartPtr<Function> origfuncNL;
		SmartPtr<SparseVector> origfuncLin;
		double origfuncConstant;
		
		/** The function in decomposed form.
		 */
		vector<SmartPtr<BlockFunction> > decompfuncNL;
		SmartPtr<SparseVector> decompfuncLin;
		double decompfuncConstant;

	public:
		ObjCon(const SmartPtr<Function>& origfuncNL_=NULL, const SmartPtr<SparseVector>& origfuncLin_=NULL, double origfuncConstant_=0., const string& name_=string());
		
		virtual ~ObjCon();
		
		bool isLinear() { return IsNull(origfuncNL); }
		bool isConstant() { return IsNull(origfuncLin) && IsNull(origfuncNL); }
				
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
	
	friend ostream& operator<<(ostream& out, const MINLPData& data);
	
	//TODO: some methods for fast evaluation of jacobian and hessian of lagrangian
}; // class MINLPData

} //namespace LaGO

#endif // LAGOMINLPDATA_H_
