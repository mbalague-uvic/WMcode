// Mean Field Implementation
// Marco Loh, loh@in.tum.de

#ifndef _IF_BRAIN_H
#define _IF_BRAIN_H

#include "pool.h"

class brain {

friend class pool;

public:
    //Constructor
    brain();

int flagavr;
int flagtstat;
int trials;
double zeit;
double avr[maxpools][maxsteps];
double avrsumup[maxpools][maxsteps];


double pp0,pp1,pp2,pp3;
double numatrix[maxsteps];
double nu2matrix[maxsteps];
int nusteps;
    
    //Initialiser
    void initBrain(int pc, int c, int cext, double deltat);
    void initPool(int p, double poolsize, double startnu, double startavrv,
double extrate, poolinfo d);
    
    //Calculation Methods
    void MeanField(double eulersteps, double eulerdelta);
    void MeanFieldEffective(double , double , int , double, int, double);
    void Spiking(int, int , ostream* , ostream*, ostream*, int,int);
    void resetPools();
    void resetavr();

    //Adaptation Methods
    void setConnection(int Pool1,int Pool2, double ampa, double nmda, double 
gaba);
    void setStartPoolrate(int p, double nu);
    void setStartAvrV(int p, double v);
    void setExtrate(int p, double nu);
    
    double getNMDAerror();
    double getPoolrate(int p);
    double getConnection(int Pool1, int Pool2, char art);
    double getPoolsize(int p);

    double getFlow(int);
    
    void printPool(int p);


private:
    void setPoolrate(int p, double nu);
    void setAvrV(int p, double v);
        
    static double poolnu[maxpools];	
    static double poolf[maxpools];	

    static int C;  // connections to pyr cells
    static int Cext;// external connections
    static int PoolCount;

    static double dt;
    static double sumAMPArec[maxpools];
    static double sumNMDA[maxpools];
    static double sumGABA[maxpools];
    static int spikes[maxpools];
    static double sumup[maxpools];

    static double ggaba[maxpools];
    static double gampa[maxpools];
    static double gnmda[maxpools];

    pool pools[maxpools];   
    double nmdaerror;   
};
#endif
