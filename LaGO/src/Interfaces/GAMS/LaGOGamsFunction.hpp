// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOGAMSFUNCTION_HPP_
#define LAGOGAMSFUNCTION_HPP_

#include "LaGObase.hpp"
#include "LaGOFunction.hpp"
#include "LaGOGamsReader.hpp"

namespace LaGO {

class GamsFunction : public Function {
private:
	SmartPtr<GamsReader::Data> data;
	
	int connr;
	
	vector<int> sparsity;
	
	void addRowName(string& message) const;

public:
	GamsFunction(int connr_, const SmartPtr<GamsReader::Data>& data_, const set<int>& sparsity_);

	double eval(const DenseVector& x) const;

	void gradient(DenseVector& grad, const DenseVector& x) const;

	void evalAndGradient(double& value, DenseVector& grad, const DenseVector& x) const;

	void hessianVectorProduct(DenseVector& product, const DenseVector& x, const DenseVector& factor) const;

	bool haveSparsity() const { return true; }

	const vector<int>& getSparsity() const { return sparsity; }
	
	void print(ostream& out) const;	
}; // class GamsFunction


	
} // namespace LaGO

#endif /*LAGOGAMSFUNCTION_HPP_*/
