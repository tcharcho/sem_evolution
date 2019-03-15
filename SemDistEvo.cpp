/*
 *  Distance based SEM evolver.  This code learns to distinguish
 *  categories of DNA data.  The input file has the format:
 *
 *  label SEQUENCE
 *  (repeat)
 *
 *  Where the label is a whole number.  Fitness is the total distance
 *  between members of different categories.  This fitness is to be
 *  maximized.
 *
 *
 */

//Include the relevant includes
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <sys/stat.h>

using namespace std;

#include "sem.h"
#include "stat.h"

//Random number seed
#define RNS 91207819

//Running statistical prints?  1=yes 0=no
#define verbose 1

//File with training examples
#define inputdata "DehydrinSet.dat"
//Input length - needs to be a bit longer than the longest sequence
#define IPL 5000
//Maximum data size -- an upper bound onthe number of sequences in the data set
#define MDS 1000

//Evolutionary algorithm details; population and tournament sizes,
//number of replicates to perform, number of mating events to use,
//reporting interval,and maximum number of mutations in a new
//structure.
// #define popsize 1000
#define tsize 7
#define runs 30
#define mevs 10000
#define RI 100
// #define MNM 3

//Population control details; number of STATES in each SEM
// #define STATES 6

void readdata();           //read in the data
void initalg();            //initialize the algorithm
double fitness(aut &A);    //report the fitness of an automata
void initpop();            //initialize a population
void matingevent();        //update the population
void report(ostream &aus); //report current statistics
void rbest(ostream &aus);  //report the best population member

//Here is the storage for the data objects
int NDI;                   //number of data items
int *data[MDS];            //stores the numerical versions of the data
int leng[MDS];             //lengths of the data items
int cate[MDS];             //data categories

//Population varaibles
aut pop[popsize];          //the population of automata
double fit[popsize];       //fitness values
int dx[popsize];           //sorting index

int main(int argc, char *argv[]) {
  fstream stat,crit;  //statistics and best SEM reporting channels
  char fn[60];        //file name construction buffer per run
  char fn2[60];       //file name construction buffer for all runs
  char dir[60];       // directory name buffer for storing best files
  int run,mev;        //run and mating event counters

  if (argc > 1) sprintf(dir, "%s", argv[1]);
  else sprintf(dir, "test_runs");

  // Creating a directory
  if (mkdir(dir, 0777) == -1) {
    // if the error isn't that the directory already exist, exit program
    if (errno != 17) {
      cerr << "Error :  " << strerror(errno) << endl;
      return 0;
    }
  }

  //construct file name for best file
  sprintf(fn2,"%s/best_%d_%d_%d.sem",dir, popsize, MNM, STATES);

  //initialize the algorithm -- including reading in the data
  initalg();

  //open the best structure reporting file
  crit.open(fn2,ios::out);

  //loop over requested replicates
  for(run=0;run<runs;run++) {
    // sprintf(fn,"run%02d.dat",run);  //construct file name for stat file
    // stat.open(fn,ios::out);        //open the stat file
    initpop();                     //initialize a new population
    // report(stat);                  //report initial statistics

    //loop over mating events
    for(mev=0;mev<mevs;mev++) {
      //update the population
      matingevent();

      //if its time for a report
      // if((mev+1)%RI==0) {
      //   if(verbose==1)cout << run << " " << (mev+1)/RI << "\n";  //eye candy?
      //   report(stat);   //report the statistics
      // }
    }
    rbest(crit);                   //report the best structure found
    // stat.close();                  //close the stat file
  }
  crit.close();  //close the best structure file
  cout << fn2 << " done\n";
  return(0);  //keep the system happy
}

//change DNA to numbers C=0 G=1 A=2 T=3
int numDNA(int c) {

  switch(c){//what letter is it?
    case 'c':
    case 'C': return(0);
    case 'g':
    case 'G': return(1);
    case 'a':
    case 'A': return(2);
    case 't':
    case 'T': return(3);
  default: return(-1);
  }
}

//read in the data
void readdata() {

  char buf[IPL];    //input buffer for characters
  int ibf[IPL];     //assembly buffer for bases
  fstream inp;      //input file
  int p;            //pointer into the data buffer
  int val;          //value of the current DNA base
  int i;            //loop index

  inp.open(inputdata,ios::in);  //open the input file
  NDI=0;                        //zero the data item counter
  inp.getline(buf,IPL-1);       //read in the first data item
  while(strlen(buf)>3){//while we still have data
    cate[NDI]=atoi(buf);  //capture the category
    p=0;                  //initialize the data pointer
    while(buf[p]!=' ')p++;//find the delimiter
    p++;                  //move past the delimiter
    leng[NDI]=0;          //zero the data length
    val=numDNA(buf[p]);   //get the value of the first base
    while(val!=-1){//while we have valid input
      ibf[leng[NDI]++]=val;     //capture the numericalized DNA base
      val=numDNA(buf[++p]);     //get the next base
    }
    data[NDI]=new int[leng[NDI]];  //get storage for this sequence
    for(i=0;i<leng[NDI];i++)data[NDI][i]=ibf[i]; //transfer the sequence
    NDI++;                         //register the new data item
    inp.getline(buf,IPL-1);        //get the next data item
  }
  inp.close();  //close the input file

}

//initialize the algorithm
void initalg() {
  int i;  //loop index

  srand48(RNS);  //seed the random number generator
  readdata();    //read in the data
  for(i=0;i<popsize;i++)pop[i].create(STATES);  //allocate the SEMs
}

//compute Euclidian distance
double EucDis(double *a,double *b) {
  double delta,ttl;   //coordinate distance and total
  int i;              //loop index

  ttl=0.0;  //zero the accumulator
  for(i=0;i<STATES;i++){//loop over the coordinates
    delta=a[i]-b[i];    //find the coordinate difference
    ttl+=(delta*delta); //sum in the squared different
  }

  return(sqrt(ttl));  //Square root to finish the distance and return it

}

//report the fitness of an automata
double fitness(aut &A) {

  static double inj[MDS][STATES];   //injected point buffer
  int cnt[STATES];                  //vector for reporting counts
  int i,j;                          //loop indices
  int ttl;                          //the total of the count vector
  double accu;                      //distance accumulator
  double accu2;                     //internal distance accumulator

  //compute the injection of the count vectors into Euclidian space
  //loop over the data items
  for(i=0;i<NDI;i++) {
    //reset the side effect machine
    A.reset();

    //traverse data item
    for(j=0;j<leng[i];j++)A.run(data[i][j]);

    //retrieve the count
    A.report(cnt);
    //zero the total
    ttl=0;

    for(j=0;j<STATES;j++)ttl+=cnt[j];  //total the counts
    for(j=0;j<STATES;j++)inj[i][j]=((double)cnt[j])/((double)ttl);  //save
  }

  // compute the between/withing fitness
  accu=0.0;    //zero the distance accumulator
  accu2=0.0;   //zero the second distance accumulator

  //loop over the pairs
  for(i=0;i<NDI-1;i++)for(j=i+1;j<NDI;j++) {
    //if the categories are different
    if(cate[i]!=cate[j]) {
      accu+=EucDis(inj[i],inj[j]);  //sum in the between distance
    } else accu2+=EucDis(inj[i],inj[j]);  //sum in the within distance
  }

  //this is the fitness
  return(accu/(accu2+1));
}

//initialize a population
void initpop() {

  int i;  //loop index

  for(i=0;i<popsize;i++){//loop over the population
    pop[i].recreate();      //put in new transitions
    fit[i]=fitness(pop[i]); //compute initial fitness
    //cout << fit[i] << " ";  //report the fitness
    dx[i]=i;                //refresh the sorting index
  }
  //cout << endl;  //terminate line

}

//update the population
void matingevent() {

  int i;    //loop index
  int nm;   //number of mutations

  tselect(fit,dx,tsize,popsize);       //select a tournament
  pop[dx[0]].copy(pop[dx[tsize-2]]);   //second best over worst
  pop[dx[1]].copy(pop[dx[tsize-1]]);   //best over second worst
  pop[dx[0]].tpc(pop[dx[1]]);          //two point crossover
  nm=lrand48()%MNM+1;                  //select number of mutations, child 1
  for(i=0;i<nm;i++)pop[dx[0]].mutate();//perform mutations
  nm=lrand48()%MNM+1;                  //select number of mutations, child 2
  for(i=0;i<nm;i++)pop[dx[1]].mutate();//perform mutations
  fit[dx[0]]=fitness(pop[dx[0]]);      //find fitness of child 1
  fit[dx[1]]=fitness(pop[dx[1]]);      //find fitness of child 2

}

//make a statistical report
void report(ostream &aus) {

  dset D;   //data set for statistical reporting

  D.add(fit,popsize);  //register the fitness data
  aus << D.Rmu() << " " << D.RCI95() << " "
      << D.Rsg() << " " << D.Rmax() << endl;  //report fitness data
  if(verbose==1){//if eye-candy requested
    cout << D.Rmu() << " " << D.RCI95() << " "
         << D.Rsg() << " " << D.Rmax() << endl;  //report fitness data
  }

}

//report the best population member
void rbest(ostream &aus) {

  int i,b;  //loop index and best pointer

  b=0;  //initialize best pointer
  for(i=1;i<popsize;i++)if(fit[i]>fit[b])b=i;  //find the best critter
  aus << fit[b] << " -fitness" << endl;        //report its fitness
  pop[b].write(aus);                           //write it out

}
