// Copyright (C) 2006 Ivo Nowak and Stefan Vigerske
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Author: Ivo Nowak, Stefan Vigerske

// tools

#ifndef TOOLS_H
#define TOOLS_H

#include <sys/resource.h> // for Timer-class
//#include <time.h>
#include <cassert>

// IO
#include <iomanip>
#include <iostream>
#include <fstream>

// data structures
#include <vector>
#include <map>
#include <list>
#include <set>

using namespace std;

#define rtol 1e-10          // 1e-12 tolerance of real numbers

#ifndef MIN
#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a) >= (b) ? (a) : (b))
#endif

// simply used as a really large number
#ifdef INFINITY
#undef INFINITY
#endif
#define INFINITY 1e99

// similarly for integers
#define INF 65535

/** The (global) counters for the objects, we point to and want to delete.
*/
extern map<const void*, int> refcount;

/** Prints all registered pointers and the number of Pointers to each of them.
    @param out The ostream to print to.
*/
extern void print_all_Pointer(ostream& out);

/** A Smartpointer with reference counter.
    When the last Pointer to an object is killed, the object itselfe will be killed.
*/
template <class Type> class Pointer {
  public:
    /** Gives the number of Pointer's to an object.
        @param obj_ The pointer to the object.
        @return The number of Pointer's to *obj_.
        @see count(const void&)
    */
    static int count(const Type* obj_) {
      map<const void*, int>::iterator p(refcount.find((const void*)obj_));
      if (p==refcount.end()) return 0;
      return p->second;
    }

  private:
    /** The object, this pointer points to.
    */
    Type* obj;

    /** Increase the counter for this object of adds it to the map.
        If obj doesn't exist in the map, it will be added.
        Else the counter for obj will be increased.
    */
    void inc_count() {
      pair<map<const void*, int>::iterator, bool> p(refcount.insert(pair<const void*, int>(obj, 1)));
      if (!p.second) p.first->second++;  // obj existed before
    }

    /** Decrease the counter for obj.
        @return The new counter.
    */
    int dec_count() {
      map<const void*, int>::iterator p(refcount.find(obj));
#ifndef NO_NULLPOINTER_CHECK
      assert(p!=refcount.end());
#endif
      int ret=--p->second;
      if (!ret) refcount.erase(p);
      return ret;
    }

    /** Indicates, whether the object should be deleted, when the last Pointer to it will be deleted.
        The default value is true.
    */
    bool delete_obj;

  public:


    /** Constructor for NULL-Pointer */
    Pointer()
    : obj(NULL), delete_obj(true)
    { }

    /** Constructor for an Type*.
        @param obj_ The pointer to the object, this pointer shoud point to.
        @param del_ Indicates, whether to object should be deleted, when the last Pointer to it will be deleted. Default is true.
    */
    Pointer(Type* obj_, bool del_=true)
    : obj(obj_), delete_obj(del_)
    { if (obj && delete_obj) inc_count();
    }

    /** Copy-Constructor for a Pointer<Type>.
        @param p The Pointer to copy.
    */
    Pointer(const Pointer<Type>& p)
    : obj(p.obj), delete_obj(p.delete_obj)
    { if (obj && delete_obj) inc_count();
    }

    /** Destructor.
        If this was the last Pointer to the object, the object is not NULL and delete_obj is true, the object will be deleted.
    */
    ~Pointer() {
      if (obj && delete_obj && (!dec_count())) delete obj;
    }

    /** Assign-Operator for another Pointer.
        @param p The Pointer to copy.
    */
    Pointer<Type>& operator=(const Pointer<Type>& p) {
      if (obj==p.obj) return *this;
      if (obj && delete_obj && (!dec_count())) delete obj;
      obj=p.obj; delete_obj=p.delete_obj;
      if (obj) inc_count();
      return *this;
    }

    /** Assign-Operator for a normal pointer.
        Sets delete_obj to true, if it's no self-assignment.
        @param obj_ The pointer to copy.
    */
    Pointer<Type>& operator=(Type* obj_) {
      if (obj==obj_) return *this;
      if (obj && delete_obj && (!dec_count())) delete obj;
      if (obj=obj_) inc_count();
      delete_obj=true;
      return *this;
    }

    /** Cast-operator to give the pointer to the object, this Pointer points to.
        @return obj.
    */
    operator Type*() const { return obj; }

    /** Cast-operator to give the pointer to the object, this Pointer points to as const.
        @return obj.
    */
//    operator const Type*() const { return (const Type*)obj; }

    /** Cast-operator to bool.
        @return True, if this is not a NULL-Pointer. False, else.
    */
//    operator bool() const { return (bool)obj; }

    /** Dereference-operator *.
        @return A reference to *obj (&*obj), if obj is not NULL.
    */
    Type& operator*() const {
#ifndef NO_NULLPOINTER_CHECK
      assert(obj);
#endif
      return *obj;
    }

    /** Dereference-operator ->.
        @return A reference to *obj (&*obj), if obj is not NULL.
    */
    Type* operator->() const {
#ifndef NO_NULLPOINTER_CHECK
      assert(obj);
#endif
      return obj;
    }

    /** Compare operator.
        @param obj_ The pointer to compare with.
        @return True, if obj_ points to the same object, obj does. False, else.
    */
/*    bool operator==(const Type* obj_) const {
      return (obj==obj_);
    }
*/
    /** Compare operator.
        @param obj_ The pointer to compare with.
        @return True, if obj_ doesn't point to the same object, obj does. False, else.
    */
    bool operator!=(const Type* obj_) const {
      return (obj!=obj_);
    }

    /** Gives the number of Pointers, pointing to *obj.
        @return The number of Pointers, pointing to *obj.
    */
    int count() const {
      return Pointer<Type>::count(obj);
    }

};

extern Pointer<ostream> out_out_p;  // defined in tools.cc
#define out_out if (out_out_p) (*out_out_p)

extern Pointer<ostream> out_log_p;  // defined in tools.cc
#define out_log if (out_log_p) (*out_log_p)

extern Pointer<ostream> out_err_p;  // defined in tools.cc
#define out_err if (out_err_p) (*out_err_p)

/** Gives a random integer.
    @param lb The lower bound.
    @param ub The upper bound.
    @return A random integer x \in {lb, ..., ub}.
*/
extern int random(int lb, int ub);

/** Gives a random double.
    @param lb The lower bound.
    @param ub The upper bound.
    @return A uniformly distributed random number from (lb, ub).
*/
extern double random(double lb, double ub);

/** Initialize the random-number-generator, using srand(time(NULL)).
*/
//extern int init_random();

/** A class to stop processor time.
    Uses getrusage.
*/
class Timer {
  /** Output-operator.
      Calls print.
      @param out The ostream to print to.
      @param w The Timer to print.
      @see print(ostream&)
  */
  friend ostream& operator<<(ostream& out, Timer& w) {
		print(out, (double)w);
//    w.print(out);
    return out;
  }

  private:
    /** Stores the starttime of this Timer.
    */
    struct timeval starttime;
    /** Stores the endtime of this Timer.
    */
    struct timeval endtime;

//		clock_t starttime;
//		clock_t endtime;

  public:
    /** Default-Constructor.
        Calls start().
        @see start()
    */
    Timer() {
      start();
    }

    /** Starts the timer.
        Calls getrusage to set starttime.
        @see stop()
    */
    void start();

    /** Stops the timer.
        Calls getrusage to set endtime.
        @return This Timer, so you can write "out << w.stop();" to print the time, which is used since this Timer was started.
        @see start()
    */
    Timer& stop();

    /** Gives the elapsed time as double.
        @return The time-difference between starttime and endtime in seconds with a precision of E-6.
    */
    operator double() const;

    /** Prints the elapsed time.
        Prints the time-difference endtime-starttime in the format hour:minute:sec, where sec is printed with a precision of E-3 (milliseconds).
        You should have set the endtime by a call to stop().
        @param out The ostream to print to.
    */
		static void print(ostream& out, double time);
//    void print(ostream& out) const;

};


unsigned int get_mem();

int start_process(const char* name, char*const* const args, int timelimit=0, const char* envvarname=NULL, const char* envvarvalue=NULL);

#endif // TOOLS_H
