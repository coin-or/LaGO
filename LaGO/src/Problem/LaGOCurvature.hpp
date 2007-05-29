// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOCURVATURE_HPP_
#define LAGOCURVATURE_HPP_

#include "LaGObase.hpp"

namespace LaGO {

enum Curvature {
	UNKNOWN = 0,
	CONVEX = 1,
	CONCAVE = 2,
	CONVEXCONCAVE = CONVEX | CONCAVE, // =3
	INDEFINITE = 4	
};

inline Curvature addCurvatures(Curvature curv1, Curvature curv2) {
	if (curv1==UNKNOWN || curv2==UNKNOWN) return UNKNOWN;
	if (curv1==INDEFINITE || curv2==INDEFINITE) return INDEFINITE;
	switch (curv1) {
		case CONVEXCONCAVE: return curv2; break;
		case CONVEX: if (curv2==CONCAVE) return INDEFINITE; else return CONVEX; break;
		case CONCAVE: if (curv2==CONVEX) return INDEFINITE; else return CONCAVE; break;
		default: return UNKNOWN; // the program shouldn't get here 
	}	 
}

inline ostream& operator<<(ostream& out, const Curvature& curv) {
	switch(curv) {
		case UNKNOWN: out << "unknown"; break;
		case CONVEX: out << "convex"; break;
		case CONCAVE: out << "concave"; break;
		case CONVEXCONCAVE: out << "convex&concave"; break;
		case INDEFINITE: out << "indefinite"; break;		
	}
	return out;
}

} // namespace LaGO

#endif /*LAGOCURVATURE_HPP_*/
