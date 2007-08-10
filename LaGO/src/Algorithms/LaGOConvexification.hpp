// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id: LaGOQuadraticFunction.hpp 109 2007-06-03 19:31:14Z stefan $

#ifndef LAGOCONVEXIFICATION_HPP_
#define LAGOCONVEXIFICATION_HPP_

#include "LaGObase.hpp"
#include "LaGOMINLPData.hpp"

namespace LaGO {

class Convexification {
private:
	MINLPData& data;

public:
	Convexification(MINLPData& data_);
	
	/** Convexifies all nonconvex quadratic terms of a problem.
	 */
	void convexify();
	
	/** Convexifies all nonconvex quadratic terms of the objective or a constraint.
	 */
	void convexify(MINLPData::ObjCon& objcon, bool need_lower, bool need_upper);

	/** Convexifies all nonconvex quadratic terms of a block function.
	 */
	void convexify(BlockFunction& func, bool do_lower, bool do_upper);

	/** Computes convexification and concavification parameters for a given matrix.
	 * @param alpha_convexify If not NULL, then the convexification parameters are computed and stored in here.
	 * @param alpha_concavify If not NULL, then the concavification parameters are computed and stored in here.
	 */	
	void convexify(SymSparseMatrix& matrix, const DenseVector& diam, DenseVector* alpha_convexify, DenseVector* alpha_concavify);
	
}; // class Convexification	
	
} // namespace LaGO

#endif /*LAGOCONVEXIFICATION_HPP_*/
