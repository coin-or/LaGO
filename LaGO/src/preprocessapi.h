// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#ifndef PREPROCESSAPI_H
#define PREPROCESSAPI_H

#include "standard.h"
#include "tools.h"
#include "param.h"
#include "problem.h"

extern "C" struct dictRec;

/** Class for preprocessing a point before a LaGO local optimization is called.
*/
class LocOptPreprocessing {
	protected:
		Pointer<Param> param;

	public:
		LocOptPreprocessing(const Pointer<Param>& param_)
		: param(param_)
		{ }

		virtual ~LocOptPreprocessing() { };

		virtual int run(UserVector<double>& x)=0;
};

/** Creates a LocOptPreprocessing object.
    @param param Parameters.
		@param LaGO_prob The problem in LaGO.
*/
extern "C" typedef LocOptPreprocessing* (preprocess_create_t)(const Pointer<Param>& param, struct dictRec* LaGO_dict, const MinlpProblem& LaGO_prob);

/** Destroys a LocOptPreprocessing object.
*/
extern "C" typedef void (preprocess_destroy_t)(LocOptPreprocessing*);

#endif
