// Mean Field Implementation
// Marco Loh, loh@in.tum.de

#ifndef _IF_POOL_H
#define _IF_POOL_H

#define VK		-80.   // (mV) reversal Ca2+
#define gAHP		7.5    // (nS , 0.015 mS/cm2) 7.5
#define alphaCa	0.1   // (0.2 muM) increment Ca2+  0.005
#define taoCa		50.   // (ms) decay Ca2+ 600
#define AHP_Ca    0   // flag

#define taoF      750 
#define taoD       50
#define U	  0.15


#define WINH      0.97

#define DILUTION  0.4

using namespace std;

const int maxpools = 11;
const int maxpoolneurons = 2600;
const int maxsteps = 20000;

const int NMAX = 7;
const double PI = 3.1415926535;

const double a1 = -1.26551223;
const double a2 = 1.00002368;
const double a3 = .37409196;
const double a4 = .09678418;
const double a5 = -.18628806;
const double a6 = .27886087;
const double a7 = -1.13520398;
const double a8 = 1.48851587;
const double a9 = -.82215223;
const double a10 = .17087277;

struct connweights {
    double ampa, nmda, gaba;
    };

struct poolinfo {
	double VL,Vthr,Vreset;	// Resting potential, firing threshold, reset potential
	double VI,VE;			//Constants

	double Cm;		// Membrane capacitance
	double gm;		// Membrane leak conductance
	double taurp;		// Refractory period
	double taum;		// time constant
	double calpha, cbeta, cgamma; //constants

//synaptic conductances
	double gAMPAext, gAMPArec, gNMDA, gGABA;
//gating variables
	double tauAMPA, tauGABA, tauNMDAdecay, tauNMDArise;
};

class pool {

friend class brain;

public:
	// configuration
    pool();
    void initPool(int number, poolinfo data);
    void setExternal(double nu);

	//Spiking Functions
	void resetSpikingPool();
	void calcSynapses();
	void calcPotentials(double,ostream*,int);

	// Meanfield Functions
    double Euler(double);
    double Flow();
    void newavrV();

private:

	// Surrouding World
    connweights w[maxpools];	// Incoming Connections
    int nr;						// Pool ID in Array of brain
	 int neurons;				// Number of neurons
    double nuext;				// External Input

	// Pool Specific Variables
   double VL,Vthr,Vreset;	// Resting potential, firing threshold, reset
	double VI,VE;			// constants
	double Cm;		// Membrane capacitance
	double gm;		// Membrane leak conductance
	double taurp;		// Refractory period
	double taum;		// time constant
	double calpha, cbeta, cgamma; //constants
	double gAMPAext, gAMPArec, gNMDA, gGABA;
	double tauAMPA, tauGABA, tauNMDAdecay, tauNMDArise;

	// Spiking Variables
	double sAMPAext[maxpoolneurons];
	double sAMPArec[maxpoolneurons];
	double sNMDA[maxpoolneurons];
	double xNMDA[maxpoolneurons];
	double sGABA[maxpoolneurons];
	double V[maxpoolneurons];
	double dV[maxpoolneurons];
	double Ca[maxpoolneurons];
	double u[maxpoolneurons];
	double x[maxpoolneurons];
        double sumu;

	double lastspiketime[maxpoolneurons];
	double expdttauAMPA,expdttauNMDArise,expdttauGABA;

	//Spiking functions
	int extspike();
	double drand();

	// Mean Field Functions
    void initVar();
    void generateintegraldata();

    double Phi();    //old
    double Phi2();  // new
	inline double nerf(double z);
    inline double falpha();
    inline double fbeta();
    inline double tauE();
    inline double muE();
    inline double sigmaE();
    inline double Sx();
    inline double rho1();
    inline double rho2();
    inline double nx();
    inline double Nx();
    inline double nix();
    inline double psi(double nu);    //old
    inline double psi2(double nu);   //new
	inline double J();
    inline double T(int n);

	// Mean Field Variables
	double startnu;
    double startavrV;
    double avrV;
	double fak2[NMAX+1];
    double bin[NMAX+1][NMAX+1];
    double alphatauNMDArisen[NMAX+1];
    double TEext, TEAMPA, TEI;
    double crho1,crho2;
    double tauNMDA, onenutauNMDA;
    double vnx, vNx, vnix;
    double vSx, vtauE, vmuE;
    double vsigmaE, csigmaE;

    double psidata[15000];
    double intdata[15000];
    double intdatabroad[300];
};

#endif
