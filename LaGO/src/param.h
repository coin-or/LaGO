// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

/** Classes to read parameters.
*/
#ifndef PARAM_H
#define PARAM_H

#include "standard.h"

/** A binary tree to represent and quick find read parameters.
*/
class ParamTree {
  public:
    /** Left and right tree.
    */
    ParamTree *left, *right;

    /** Name of the parameter in this node.
    */
    const char *name;

    /** Value of the parameter in this node.
    */
    const char *value;

    /** Constructor for a new node.
        @param name_ Name of the parameter, not NULL !.
        @param value_ Value of the parameter, default is NULL.
        @see ParamTree(ParamTree&)
    */
    ParamTree(const char* name_, const char* value_=(char*)NULL)
    : left(NULL), right(NULL), name(name_), value(value_)
    { if (! name_) {
        out_err << "Can't handle NULL-pointer as argument for ParamTree::ParamTree(char*, char*=0)." << endl;
        exit(-1);
      }
    }

    /** Copy-Constructor.
        Copys name, value and the left and right tree.
        @param PT The ParamTree to copy.
        @see ParamTree()
    */
    ParamTree(ParamTree &PT)
    : name(strdup(PT.name)), value(PT.value ? strdup(PT.value) : NULL),
      left(PT.left ? new ParamTree(*PT.left) : NULL),
      right(PT.right ? new ParamTree(*PT.right) : NULL)
    { }

    /** Destructor.
        Deletes the left and right tree, name and value.
    */
		~ParamTree();

    /** Add's a new node in the tree.
        If a node for this name still exist's, the value will be overwritten.
        @param name_ The name of the parameter to add to the tree.
        @param value_ The value of the parameter to add to the tree, default is NULL.
    */
    void add(const char* name_, const char* value_=(const char*)NULL);

    /** Get's the value of a parameter from the tree.
        @param name_ The name of the parameter to look for.
        @return The value of the parameter. If it doesn't exist, returns NULL.
    */
    const char* get(const char* name_);

    /** Print's the whole tree.
        Print's the left tree, the root and the right tree.
        @param out The ostream to print to, default is out_log.
    */
    void dump(ostream& out);

};

/** Class for reading parameter-files and storing them in a tree (ParamTree).
*/
class Param {

  /** Print's the read parameters.
      @see ParamTree::dump(ostream&)
  */
  friend ostream& operator<<(ostream &out, Param& a) {
    if (a.head) a.head->dump(out);
    else out << "Empty parameter-tree";
    return out;
  }

private:
  /** Counter for the number of the files, which were read by now.
      -1, when no file was read.
  */
  int filenr;

  /** File to read.
  */
  ifstream *file;

  /** A list of parameter files to read.
  */
  vector<char*> paramfiles;

  /** The head of the parameter-tree.
  */
  ParamTree *head;

public:
  /** Base-directory for the parameter files.
  */
  char *basedir;

  /** Constructor for a list of parameter files.
      Copys the strings form the given argument.
      @param paramfiles_ The list of parameter files as vector of char*.
      @param basedir_ The base-directory, where the files can be found.
      @see Param(char*, char*)
      @see Param()
      @see Param(Param&)
  */
  Param(vector<char*> &paramfiles_, const char* basedir_="resource/")
  : paramfiles(paramfiles_.size()), filenr(-1), head(NULL), basedir(basedir_ ? strdup(basedir_) : new char(0))
  { for (int i=0; i<paramfiles_.size(); i++)
      paramfiles[i]=(paramfiles_[i] ? strdup(paramfiles_[i]) : NULL);
  }

  /** Constructor for the name of one parameter file.
      Copys the pointer to the filename.
      @param paramfile_ The name of the parameter file.
      @param basedir_ The base-directory, where the files can be found.
      @see Param(vector<char*>&, char*)
      @see Param()
      @see Param(Param&)
  */
  Param(const char* paramfile_, const char* basedir_="resource/")
  : paramfiles(paramfile_ ? 1 : 0), filenr(-1), head(NULL), basedir(basedir_ ? strdup(basedir_) : new char(0))
  {
    if (paramfile_) paramfiles[0]=strdup(paramfile_);
  }

  /** Standard-Constructor for no parameteter files.
      @see Param(vector<char*>&, char*)
      @see Param(char*, char*)
      @see Param(Param&)
  */
  Param() : filenr(-1), head(NULL), basedir(strdup("resource/")) { }

  /** Copy-Constructor.
      Call's Assign-Operator.
      @param P The Param to copy.
      @see operator=(const Param&)
      @see Param(vector<char*>&, char*)
      @see Param(char*, char*)
      @see Param()
  */
  Param(Param& P) {
    *this=P;
  }

  /** Destructor.
      Deletes the filename-list and frees the parameter-tree.
      @see ParamTree::~ParamTree
  */
  ~Param();
  
  /** Add's a file to the list of the parameter files.
      @param filename The name of the file.
  */
  void add_file(char* filename) {
    if (filename) paramfiles.push_back(strdup(filename));
  }

  /** Read's all parameter-files, starting with file filenr+1 and stores them in the parameter-tree.
      @return 0, if everything went right; -1, if there was an error when opening a file; or the number of parse-errors.
      @see head
  */
  int read();

  /** Get's the value for a parameter-name.
      @param name The name of the parameter to look for.
      @param def A default value for the parameter, if it wasn't set, default is NULL.
      @return A copy of the value of the parameter or NULL, if no parameter of the specific name was read.
  */
	Pointer<char> get(const char *name, const char* def=(const char*)NULL) const {
    const char* parvalue=(head ? head->get(name) : NULL);
    return (parvalue ? strdup(parvalue) : (def ? strdup(def) : NULL));
  }

  /** Get's the value of a parameter from the tree and converts it to a double.
      @param name_ The name of the parameter to look for.
      @param def The default value for the parameter, if it wasn't set, default is 0.
      @return The value of the parameter, converted to a double.
  */
  double get_d(const char* name_, double def=0) const {
    Pointer<char> parvalue=get(name_);
    double val=(parvalue ? (double)atof(parvalue) : def);

//    if (parvalue) free(parvalue);
    return val;
  };

  /** Get's the value of a parameter from the tree and converts it to an int.
      @param name_ The name of the parameter to look for.
      @param def The default value for the parameter, if it wasn't set, default is 0.
      @return The value of the parameter, converted to an integer.
  */
  int get_i(const char* name_, int def=0) const {
    Pointer<char> parvalue=get(name_);
    int val=(parvalue ? atoi(parvalue) : def);

//    if (parvalue) free(parvalue);
    return val;
  }

  /** Adds a (name, value) pair to the parameter tree.
      The char*'s are copied, using strdup, if not NULL.
      @param name_ The name of the parameter to look for.
      @param value_ The value of the parameter.
  */
  void add(const char* name_, const char* value_) {
    if (!name_) return;
    if (head) head->add(strdup(name_), value_ ? strdup(value_) : NULL);
    else head=new ParamTree(strdup(name_), value_ ? strdup(value_) : NULL);
  }

  /** Assing-Operator.
      Copy's the basedir, file list, the parameter-tree and filenr.
      @param p The SnoptParam to copy.
  */
  Param& operator = (const Param &p)
  {
    if (&p != this) {  // Check against self assignment
      basedir = (p.basedir ? strdup(p.basedir) : new char(0));

      for (int i=0; i<p.paramfiles.size(); i++)
        paramfiles.push_back(p.paramfiles[i] ? strdup(p.paramfiles[i]) : NULL);

      filenr=p.filenr;

      head=(p.head ? new ParamTree(*p.head) : NULL);
    }
    return *this;
  }
};

#endif // PARAM_H
