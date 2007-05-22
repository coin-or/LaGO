// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id: base.hpp 94 2007-05-21 13:54:40Z stefan $

#ifndef LAGODENSEVECTOR_HPP_
#define LAGODENSEVECTOR_HPP_

#include "LaGObase.hpp"

#include "IpDenseVector.hpp"

namespace LaGO {
	
#define DenseVectorSpace Ipopt::DenseVectorSpace

class DenseVector : public Ipopt::DenseVector {
public:
	DenseVector(const DenseVectorSpace* owner_space)
	: Ipopt::DenseVector(owner_space)
	{ }

	DenseVector& operator*=(double factor) {
		Scal(factor);
		return *this;
	}
	
	double operator()(int index) const {
		assert(index>=0 && index<Dim());
		if (IsHomogeneous()) return Scalar();
		return Values()[index];
	}
};
	
	
} // namespace LaGO


#endif /*LAGODENSEVECTOR_HPP_*/
