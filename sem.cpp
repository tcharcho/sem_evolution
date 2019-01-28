/*  Side effect machine class of finite state automata
 *
 *  code file
 *
 *
 */
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <cmath>

#include "sem.h"


//constructors and destructors
aut::aut(){//create an empty automata

  blank();

}

aut::aut(int nv){//create an automata with nv states

  n=0;
  create(nv);

}

aut::aut(const aut &other){//copy constructor

int i,j;

  n=other.n;
  cur=other.cur;
  cntr=new int[n];
  states=new int*[n];
  for(i=0;i<n;i++){
    states[i]=new int[TR];
    for(j=0;j<TR;j++)states[i][j]=other.states[i][j];
    cntr[i]=other.cntr[i];
  }
}

aut::~aut(){//deallocate the automata

  clear();

}

//management routines
void aut::blank(){//create an empty automata

  n=0;

}

void aut::clear(){//clear the automata

  if(n){
    for(int i=0;i<n;i++)delete [] states[i];
    delete [] states;
    delete [] cntr;
  }
  n=0;

}

void aut::create(int nv){//create an automata with nv states

  if(n)clear();
  n=nv;
  states=new int*[n];
  for(int i=0;i<n;i++)states[i]=new int[TR];
  cntr=new int[n];
  recreate();

}

void aut::recreate(){//randomize an allocated automata

int i,j;

  if(n){//failsafe the recreation of empty automata
    for(i=0;i<n;i++){//loop over states
      for(j=0;j<TR;j++)states[i][j]=lrand48()%n;  //make transitions
    }  
  }
}

void aut::copy(aut &other){//make a copy of another automata

int i,j;

  if(n!=other.n)create(other.n);
  cur=other.cur;
  for(i=0;i<n;i++){
    for(j=0;j<TR;j++)states[i][j]=other.states[i][j];
    cntr[i]=other.cntr[i];
  }
}

//use routines
void aut::reset(){//set state to zero and clear counters

  cur=0;                         //set current state to zero
  for(int i=0;i<n;i++)cntr[i]=0; //zero the counters
  cntr[0]=1;                     //record the starting state

}

void aut::run(int v){//run on input value v

  cur=states[cur][v];  //make the state transition
  cntr[cur]++;         //record the entry into the state

}

void aut::report(int *rp){//report the state useage vectors

  //not that rp is presumed to be long enough
  for(int i=0;i<n;i++)rp[i]=cntr[i];

}

//genetic operations
void aut::tpc(aut &other){//perform two point crossover

int i,j,sw,cp1,cp2;

  cp1=lrand48()%n; //find first crossover point
  cp2=lrand48()%n; //find second crossover point
  if(cp1>cp2){sw=cp1;cp1=cp2;cp2=sw;} //order them
  for(i=cp1;i<=cp2;i++){
    for(j=0;j<TR;j++){//cross over the transitions
      sw=states[i][j];
      states[i][j]=other.states[i][j];
      other.states[i][j]=sw;
    }  
  }
}

void aut::mutate(){//perform a transition mutation

  states[lrand48()%n][lrand48()%TR]=lrand48()%n;

}

//input and output
void aut::write(ostream &aus){//output method

int i,j;

  aus << n << " states" << endl;
  for(i=0;i<n;i++){
    aus << i << ":";
    aus << states[i][0];   
    for(j=1;j<TR;j++)aus << "," << states[i][j];
    aus << endl;
  }

}

void aut::read(istream &inp){//input method

char buf[1000];
int i,j,k;
  inp.getline(buf,999);
  create(atoi(buf));
  for(i=0;i<n;i++){
    inp.getline(buf,999);
    j=0;
    while(buf[j]!=':')j++;
    for(k=0;k<TR;k++){
      j++;
      states[i][k]=atoi(buf+j);
      while(isdigit(buf[j]))j++;
    }
  }  

}
