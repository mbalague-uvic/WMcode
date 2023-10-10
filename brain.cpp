// Mean Field Implementation
// Marco Loh, loh@in.tum.de
#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "brain.h"
#include "pool.h"

int brain::PoolCount = 0;

int brain::C = 0;

int brain::Cext = 0;

double brain::poolnu[maxpools];

double brain::poolf[maxpools];

double brain::dt = 0;

double brain::sumAMPArec[maxpools];

double brain::sumNMDA[maxpools];

double brain::sumGABA[maxpools];

int brain::spikes[maxpools];

double brain::sumup[maxpools];

double brain::ggaba[maxpools];

double brain::gampa[maxpools];

double brain::gnmda[maxpools];

//constuctor

brain::brain(void)
{
}

void brain::resetavr(void)
{
   for (int i = 0; i < maxpools; i++)
      for (int j = 0; j < maxsteps; j++)
         avr[i][j] = 0;

   for (int i = 0; i < maxpools; i++)
   {
      ggaba[i] = 0;
      gampa[i] = 0;
      gnmda[i] = 0;
   }
}

void brain::initBrain(int pc, int c, int cext, double deltat)
{
   pp0 = 0;
   pp1 = 0;
   pp2 = 0;
   pp3 = 0;
   PoolCount = pc;
   C = c;
   Cext = cext;
   //deltat for spiking neurons-
   dt = deltat;
   int i, j;

   for (i = 0; i < maxpools; i++)
      for (j = 0; j < maxpools; j++)
      {
         pools[j].w[i].ampa = 0;
         pools[j].w[i].nmda = 0;
         pools[j].w[i].gaba = 0;
      }

   for (i = 0; i < maxpools; i++)
      for (j = 0; j < maxsteps; j++)
         avr[i][j] = 0;
}

void brain::initPool(int p, double poolsize, double startnu, double startavrv,
        double extrate, poolinfo d)
{
   poolf[p] = poolsize;
   pools[p].startnu = startnu;
   pools[p].startavrV = startavrv;
   pools[p].setExternal(extrate);
   pools[p].initPool(p, d);
}

void brain::resetPools()
{
   for (int j = 0; j < brain::PoolCount; j++)
   {
      setPoolrate(j, pools[j].startnu);
      setAvrV(j, pools[j].startavrV);
      pools[j].resetSpikingPool();
   }

}

void brain::setConnection(int Pool1, int Pool2, double ampa, double nmda,
        double gaba)
{
   pools[Pool2].w[Pool1].ampa = ampa;
   pools[Pool2].w[Pool1].nmda = nmda;
   pools[Pool2].w[Pool1].gaba = gaba;
}

double brain::getConnection(int Pool1, int Pool2, char art)
{
   switch (art)
   {
   case 'a':
      return (pools[Pool2].w[Pool1].ampa);
   case 'n':
      return (pools[Pool2].w[Pool1].nmda);
   case 'g':
      return (pools[Pool2].w[Pool1].gaba);
   default:
      return (99999);           
   }
}

void brain::setPoolrate(int p, double nu)
{
   poolnu[p] = nu;
}

void brain::setAvrV(int p, double v)
{
   pools[p].avrV = v;
}

void brain::setExtrate(int p, double nu)
{
   pools[p].setExternal(nu);
}


void brain::setStartPoolrate(int p, double nu)
{
   pools[p].startnu = nu;
}

void brain::setStartAvrV(int p, double v)
{
   pools[p].startavrV = v;
}

double brain::getPoolrate(int p)
{
   return poolnu[p];
}

double brain::getPoolsize(int p)
{
   return (poolf[p]);
}

double brain::getNMDAerror()
{
   return (nmdaerror);
}

void brain::Spiking(int offsettime, int timespan, ostream * rateoutput,
        ostream * spout, ostream * uout, int avtime, int flpr)
{
   int steps = int (timespan / dt + 0.5);

   for (int j = 0; j < steps; j++)
   {
      for (int i = 0; i < PoolCount; i++)
         pools[i].calcSynapses();
      for (int i = 0; i < PoolCount; i++)
         pools[i].calcPotentials((j + 1) * dt + offsettime, spout, flpr);

/////
      if (flagavr == 5)
      {
         if (((j + 1) % int ((avtime / dt) + 0.5)) == 0)
         {
            *rateoutput << (j + 1) * dt + offsettime << "  ";
            *uout << (j + 1) * dt + offsettime << "  ";
            for (int i = 0; i < PoolCount; i++)
            {
               *rateoutput << double (spikes[i] * (1000 / avtime)) /
                       pools[i].neurons << "  ";
               *uout << double (sumup[i] * (dt / avtime)) /
                       pools[i].neurons << "  ";
            }
            *rateoutput << endl;
            *uout << endl;

         }
      }
/////
      // Output
      if (flagavr == 1)
      {
         if (((j + 1) % int ((avtime / dt) + 0.5)) == 0)
         {
            *rateoutput << (j + 1) * dt + offsettime << "  ";
            *uout << (j + 1) * dt + offsettime << "  ";
            for (int i = 0; i < PoolCount; i++)
            {
               *rateoutput << (avr[i][(int) ((j + 1) * dt + offsettime)]
                       +
                       double (spikes[i] * (1000 / avtime)) /
                       pools[i].neurons) /trials << "  ";
               *uout << (avrsumup[i][(int) ((j + 1) * dt + offsettime)] +
                       double (sumup[i] * (dt / avtime)) /
                       pools[i].neurons) /trials << "  ";
               if (i == 1)
               {
                  numatrix[nusteps] =
                          (double (spikes[i] * (1000 / avtime)) /
                          pools[i].neurons);
               }
               if (i == 2)
               {
                  nu2matrix[nusteps] =
                          (double (spikes[i] * (1000 / avtime)) /
                          pools[i].neurons);
                  nusteps++;
               }
               spikes[i] = 0;
               sumup[i] = 0;
            }
            *rateoutput << endl;
            *uout << endl;
         }
      }
      else
      {
         if (((j + 1) % int ((avtime / dt) + 0.5)) == 0)
         {
            for (int i = 0; i < PoolCount; i++)
            {
               avr[i][(int) ((j + 1) * dt + offsettime)] =
                       avr[i][(int) ((j + 1) * dt + offsettime)] +
                       double (spikes[i] * (1000 / avtime)) /
                       pools[i].neurons;
               if (i == 1)
               {
                  numatrix[nusteps] =
                          (double (spikes[i] * (1000 / avtime)) /
                          pools[i].neurons);
               }
               if (i == 2)
               {
                  nu2matrix[nusteps] =
                          (double (spikes[i] * (1000 / avtime)) /
                          pools[i].neurons);
                  nusteps++;
               }
               spikes[i] = 0;
               avrsumup[i][(int) ((j + 1) * dt + offsettime)] =
                       avrsumup[i][(int) ((j + 1) * dt + offsettime)]
                       + double (sumup[i] * (dt / avtime)) / pools[i].neurons;
               sumup[i] = 0;
            }
         }
      }
   }
}

void brain::MeanField(double eulersteps, double eulerdelta)
{
   resetPools();
   char d;

   double tempnu[maxpools];

   int ok = 0;

   int eulercount = 0;

   nmdaerror = 0;

   ofstream aus;

   char dateiname[20];

   sprintf(dateiname, "Mfrate");
   aus.open(dateiname, ios::out);

   while (ok == 0)
   {

      
      eulercount++;
      if (eulercount == eulersteps)
         ok = 1;

      // Calculation
      aus << eulercount * eulerdelta;
      for (int j = 0; j < PoolCount; j++)
      {
         tempnu[j] = pools[j].Euler(eulerdelta);
         aus << " " << tempnu[j];
      }
      aus << endl;
      // Setting rates for new iteration + output
      for (int j = 0; j < PoolCount; j++)
      {
         if (tempnu[j] - poolnu[j] > 0.0000001)
         {
            d = '+';
         }
         else if (tempnu[j] - poolnu[j] < -0.0000001)
         {
            d = '-';
         }
         else
            d = 'o';
        

         setPoolrate(j, tempnu[j]);
      }

      // Updating avrV
      for (int j = 0; j < PoolCount; j++)
      {
         pools[j].newavrV();
         //Warning
         if (pools[j].avrV < -56.)
            nmdaerror += pow(-56. - pools[j].avrV, 2);
         if (pools[j].avrV > -48.)
            nmdaerror += pow(-48. - pools[j].avrV, 2);
         
      }
      
   }
}

void brain::MeanFieldEffective(double eulersteps, double eulerdelta, int fix1,
        double u1, int fix2, double u2)
{
   resetPools();
   char d;

   double tempnu[maxpools];

   int ok = 0;

   int eulercount = 0;

   nmdaerror = 0;

   ofstream aus;

   char dateiname[20];

   sprintf(dateiname, "Mfrate");
   aus.open(dateiname, ios::out);

   setPoolrate(fix1, u1);
   setPoolrate(fix2, u2);

   while (ok == 0)
   {

      // loop termination (now like a "for" loop)
      eulercount++;
      if (eulercount == eulersteps)
         ok = 1;

      // Calculation
      aus << eulercount * eulerdelta;
      for (int j = 0; j < PoolCount; j++)
      {
         tempnu[j] = pools[j].Euler(eulerdelta);
         aus << " " << tempnu[j];
      }
      aus << endl;
      // Setting rates for new iteration + output
      for (int j = 0; j < PoolCount; j++)
      {
         if (tempnu[j] - poolnu[j] > 0.0000001)
         {
            d = '+';
         }
         else if (tempnu[j] - poolnu[j] < -0.0000001)
         {
            d = '-';
         }
         else
            d = 'o';
      

         if (j != fix1 && j != fix2)
            setPoolrate(j, tempnu[j]);
         setStartPoolrate(j, tempnu[j]);
      }

      // Updating avrV
      for (int j = 0; j < PoolCount; j++)
      {
         pools[j].newavrV();
         //Warning
         if (pools[j].avrV < -56.)
            nmdaerror += pow(-56. - pools[j].avrV, 2);
         if (pools[j].avrV > -48.)
            nmdaerror += pow(-48. - pools[j].avrV, 2);
         
      }
      
   }
}

double brain::getFlow(int p)
{
   return pools[p].Flow();
}

// bunch of output for testing purposes
void brain::printPool(int p)
{
   cout << "====== Global Vars =======" << endl;
   cout << "C: " << C << endl;  // connections to pyr cells
   cout << "Cext: " << Cext << endl;    // external connections
   cout << "PoolCount: " << PoolCount << endl;
   cout << "Frequencies: ";
   for (int i = 0; i < PoolCount; i++)
      cout << poolnu[i] * 1000 << " ";
   cout << endl;

   cout << "Poolsizes: ";
   for (int i = 0; i < PoolCount; i++)
      cout << poolf[i] << " ";
   cout << endl;

   cout << "Extrates: ";
   for (int i = 0; i < PoolCount; i++)
      cout << pools[i].nuext << " ";
   cout << endl;

   cout << "====== AMPA Connections =======" << endl;
   printf("  ");
   for (int i = 0; i < PoolCount; i++)
      printf("%5d", i);
   printf("\n");
   for (int i = 0; i < PoolCount; i++)
   {
      printf("%2d ", i);
      for (int j = 0; j < PoolCount; j++)
         printf("%5.2f", pools[j].w[i].ampa);
      printf("\n");
   }

   cout << "====== NMDA Connections =======" << endl;
   printf("  ");
   for (int i = 0; i < PoolCount; i++)
      printf("%5d", i);
   printf("\n");
   for (int i = 0; i < PoolCount; i++)
   {
      printf("%2d ", i);
      for (int j = 0; j < PoolCount; j++)
         printf("%5.2f", pools[j].w[i].nmda);
      printf("\n");
   }

   cout << "====== GABA Connections =======" << endl;
   printf("  ");
   for (int i = 0; i < PoolCount; i++)
      printf("%5d", i);
   printf("\n");
   for (int i = 0; i < PoolCount; i++)
   {
      printf("%2d ", i);
      for (int j = 0; j < PoolCount; j++)
         printf("%5.2f", pools[j].w[i].gaba);
      printf("\n");
   }

   cout << "====== Local Vars =======" << endl;
   cout << "Pool Nr: " << p << endl;
   cout << "VL: " << pools[p].VL << endl;
   cout << "Vthr: " << pools[p].Vthr << endl;
   cout << "Vreset: " << pools[p].Vreset << endl;
   cout << "VI: " << pools[p].VI << endl;
   cout << "VE: " << pools[p].VE << endl;       // constants
   cout << "avrV: " << pools[p].avrV << endl;;
   cout << "Cm: " << pools[p].Cm << endl;;      // Membrane capacitance
   cout << "gm: " << pools[p].gm << endl;;      // Membrane leak conductance
   cout << "taurp: " << pools[p].taurp << endl;;        // Refractory period
   cout << "taum: " << pools[p].taum << endl;;  // time constant
   cout << "calpha: " << pools[p].calpha << endl;
   cout << "cbeta: " << pools[p].cbeta << endl;
   cout << "cgamma: " << pools[p].cgamma << endl;

//synaptic conductances
   cout << "gAMPAexp: " << pools[p].gAMPAext << endl;
   cout << "gAMPArec: " << pools[p].gAMPArec << endl;
   cout << "gNMDA: " << pools[p].gNMDA << endl;
   cout << "gGABA: " << pools[p].gGABA << endl;

//gating variables
   cout << "tauAMPA: " << pools[p].tauAMPA << endl;
   cout << "tauGABA: " << pools[p].tauGABA << endl;
   cout << "tauNMDAdecay: " << pools[p].tauNMDAdecay << endl;
   cout << "tauNMDArise: " << pools[p].tauNMDArise << endl;
   cout << "Nr: " << pools[p].nr << endl;
   cout << "nuext: " << pools[p].nuext << endl; // overall external input (HZ)

}
