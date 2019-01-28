/*  Side effect machine class of finite state automata
 *
 *  header file
 *
 *
 */#ifndef	_SEM_H
#define	_SEM_H

using namespace std;

//number of transitions out per state
#define TR 4

class aut {  //side effect machine class

 public:

  //constructors and destructors
  aut();            //create an empty automata
  aut(int nv);      //create an automata with nv states
  aut(const aut &); //copy constructor
  ~aut();           //deallocate the automata

  //management routines
  void blank();          //create an empty automata
  void clear();          //clear the automata
  void create(int nv);   //create an automata with nv states
  void recreate();       //randomize an allocated automata
  void copy(aut &other); //make a copy of another automata

  //use routines
  void reset();          //set state to zero and clear counters
  void run(int v);       //run on input value v
  void report(int *rp);  //report the state useage vectors

  //genetic operations
  void tpc(aut &other);  //perform two point crossover
  void mutate();         //perform a transition mutation

  //input and output
  void write(ostream &aus);
  void read(istream &inp);
  
 private:

  int n;          //number of states
  int cur;        //current state   
  int **states;   //state vector
  int *cntr;      //counter vector


};


#endif /* _SEM_H */
