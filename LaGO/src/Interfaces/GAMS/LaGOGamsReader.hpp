// Copyright (C) Stefan Vigerske 2007
// All Rights Reserved.
// This code is published under the Common Public License.

// $Id$

#ifndef LAGOGAMSREADER_HPP_
#define LAGOGAMSREADER_HPP_

#include "LaGObase.hpp"

namespace LaGO {

class MINLPData;
class GamsSolver;

class GamsReader {
	friend class GamsFunction;
	friend class GamsSolver;
private:
	class Data : public ReferencedObject {
	public:
		/** Dimension in gams.
		*/
		int gamsdim;

		/** Dictionary for variable and constraint names.
		*/
		struct dictRec* dict;

		/** Stripped NL Instruction data. */
		unsigned int* instr;
		double* nlCons;
		int* startIdx;
		int* numInstr;

		/** Sum of Instruction lengths. */
		int lenins;
		/** Max instruction length. */
		int maxins;

		/** Used for function evaluation. */
		double* s;
		double* sbar;
		double* resstack;
		int resstacksize;

		/** Number of domain violations so far.
		*/
		int domain_violations;

		Data();

		~Data();
	};

	double obj_factor;
	int objcon, objvar;
	bool is_minimization;
	bool reformed;
	
	// used by GamsSolver
	int* con_types;
	double* con_rhs;

	SmartPtr<Data> data;

	/** get name of row i
	    @param dict
	    @param i row index, [0..nRows)
	    @param bufLen size of target buffer
			@param name target buffer for row name
      @return name on success, NULL on failure
	*/
	static char* getRowName (struct dictRec* dict, int i, char *name, int bufLen);

	/** get name of column j
	    @param dict
	    @param j column index, [0..nCols)
	    @param bufLen size of target buffer
			@param name target buffer for column name
      @return name on success, NULL on failure
	*/
	static char* getColName (struct dictRec* dict, int j, char *name, int bufLen);

#if defined(COIN_HAS_CPX) && defined(COIN_HAS_GAMSCPLEXLICE)
	void initCPLEXLicense(int connr, int varnr, int nnz, int nlnz, int ndisc) const;
#endif
public:
  GamsReader();

  ~GamsReader();

	SmartPtr<MINLPData> getProblem(char* cntr_file);
	
	string getParameterfilename();

	void writeSolutionFile(const DenseVector& sol_point, double obj_value, int model_status, int solver_status, int iter, double time);
	
}; // class GamsReader

} // namespace LaGO

#endif /*LAGOGAMSREADER_HPP_*/
