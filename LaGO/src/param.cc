// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Stefan Vigerske

#include "param.h"

ParamTree::~ParamTree() {
	if (left) delete left;
	if (right) delete right;
	free((char*)name);
	if (value) free((char*)value);
}

void ParamTree::add(const char* name_, const char* value_) {
  if (! name_) return;
  int cmp=strcmp(name_, name);
  if (cmp==0) { value=value_; return; }
  if (cmp<0) {
    if (left) left->add(name_, value_);
    else left=new ParamTree(name_, value_);
    return;
  }
  if (right) right->add(name_, value_);
  else right=new ParamTree(name_, value_);
}

const char* ParamTree::get(const char* name_) {
  if (! name_) return NULL;

  int cmp=strcmp(name_, name);
  if (cmp==0) return value;
  if (cmp<0) return (left ? left->get(name_) : NULL);
  return (right ? right->get(name_) : NULL);
}

void ParamTree::dump(ostream& out) {
  if (left) left->dump(out);
	out << name << ": " << value << endl;
  if (right) right->dump(out);
}

// --------------------------------------------------------------------

Param::~Param() {
	for (int i=0; i<paramfiles.size(); i++)
		if (paramfiles[i]) free(paramfiles[i]);
	if (head) delete head;
	if (basedir) free(basedir);
}

int Param::read() {
  bool readnext=true;
  char *line, *name, *value;
  int linenr=0;
  int ret=0;

  while (filenr+1<paramfiles.size()) {  // go through all files
    if (paramfiles[++filenr]) { // if filename != NULL
      if (basedir) {
        char *fullname=new char[strlen(basedir)+strlen(paramfiles[filenr])+1];
        strcpy(fullname, basedir);
        strcpy(fullname+strlen(basedir), paramfiles[filenr]);
        file = new ifstream(fullname, ios::in);
        delete[] fullname;
      }
      else file=new ifstream(paramfiles[filenr], ios::in);
      if (! file->good()) {  // if file was not correct opened
        cerr << "Param::read(): Error opening file " << paramfiles[filenr] << ":";
      	if (file->eof()) cerr << "eof ";
      	if (file->fail()) cerr << "fail ";
      	if (file->bad()) cerr << "bad ";
      	cerr << endl;
        delete file;
        return -1;
      }
			line=new char[255];
      while (!file->eof()) {
        if (file->getline(line, 255)) { // read's one line.
          linenr++;
          if (*line=='#') continue; // comment line

          name=line;
          while (*name==' ') name++; // go over whitespaces
          if (! *name) continue; // empty line

          value=name;
          while (*++value) if (*value==':') { *value++='\0'; break; } // search for ":" and cut the name there
          while (*value==' ') value++; // go over whitespaces

          if (*value && *name) {
            for (int i=strlen(name)-1; i>=0 && name[i]==' '; name[i--]='\0'); // delete trailing whitespaces
            for (int i=strlen(value)-1; i>=0 && value[i]==' '; value[i--]='\0'); // delete trailing whitespaces
            if (head) head->add(strdup(name), strdup(value));
            else head=new ParamTree(strdup(name), strdup(value));
          }
          else {
            cerr << "Param::read(): Error parsing file " << paramfiles[filenr] << ": line " << linenr << endl;
            ret++;
          }
        }
      }
			delete[] line;
      delete file;  // close file
    }
  }
  return ret;
}

