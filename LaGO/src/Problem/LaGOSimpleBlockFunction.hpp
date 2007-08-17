#ifndef LAGOSIMPLEBLOCKFUNCTION_HPP_
#define LAGOSIMPLEBLOCKFUNCTION_HPP_

#include "LaGObase.hpp"

namespace LaGO {
	
/** A function that is only defined on a subset of the variables.
 * A kind of reverse for of an RestrictedFunction.
 */ 
class SimpleBlockFunction : public Function {
private:
	SmartPtr<const Function> f;
	vector<int> indices;
	
	mutable DenseVector my_x;

	vector<int> sparsity;
public:
	SimpleBlockFunction(const SmartPtr<const Function>& f_, const vector<int>& indices_)
	: f(f_), indices(indices_), my_x(indices.size())
	{ assert(IsValid(f));
		if (f->haveSparsity()) {
			const vector<int>& f_sparsity(f->getSparsity());
			sparsity.resize(f_sparsity.size());
			for (int i=0; i<(int)f_sparsity.size(); ++i)
				sparsity[i]=indices[f_sparsity[i]];
		} else {
			sparsity=indices;
		}
	} 

	double eval(const DenseVector& x) const {
		my_x.setToBlock(x, indices);
		return f->eval(my_x);
	}
	
	void gradient(DenseVector& grad, const DenseVector& x) const;
	
	void evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const;

	void hessianVectorProduct(DenseVector& product, const DenseVector& x, const DenseVector& factor) const;

#ifdef COIN_HAS_FILIB
	bool canIntervalEvaluation() const { return f->canIntervalEvaluation(); }
	
	interval<double> eval(const IntervalVector& x) const;

	void evalAndGradient(interval<double>& value, IntervalVector& grad, const IntervalVector& x) const;
#endif
	
	bool haveSparsity() const { return true; } 

	virtual const vector<int>& getSparsity() const { return sparsity; }
	
	void print(ostream& out) const { f->print(out); }	
}; // class SimpleBlockFunction	
	
} // namespace LaGO

#endif /*LAGOSIMPLEBLOCKFUNCTION_HPP_*/
