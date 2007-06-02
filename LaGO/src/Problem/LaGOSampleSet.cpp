// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#include "LaGOSampleSet.hpp"

namespace LaGO {
	
double SamplePoint::tol=1e-4;

void SampleSet::addVector(const vector<DenseVector>& pointvector, const vector<int>& indices, bool set_startpoint_flag) {
	DenseVector newpoint;
	for (vector<DenseVector>::const_iterator it(pointvector.begin()); it!=pointvector.end(); ++it) {
		newpoint.setToBlock(*it, indices);
		pair<iterator, bool> ret=insert(newpoint);
		if (set_startpoint_flag) ret.first->is_startpoint=true;		
	}
}

void SampleSet::addVector(const vector<DenseVector>& pointvector, bool set_startpoint_flag) {
	for (vector<DenseVector>::const_iterator it(pointvector.begin()); it!=pointvector.end(); ++it) {
		pair<iterator, bool> ret=insert(*it);
		if (set_startpoint_flag) ret.first->is_startpoint=true;		
	}
}

} // namespace LaGO
