// Mean Field Implementation
// Marco Loh, loh@in.tum.de

#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "brain.h"
#include "pool.h"

pool::pool(void)
{
}

void pool::setExternal(double nu)
{
   nuext = nu;
}

void pool::initPool(int number, poolinfo pi)
{

   int     sd = 1001;

   nr = number;
   neurons = int (brain::C * brain::poolf[nr] + 0.5);

   //to make the program better readable
   VL = pi.VL;
   Vthr = pi.Vthr;
   Vreset = pi.Vreset;
   VE = pi.VE;
   VI = pi.VI;
   Cm = pi.Cm;
   gm = pi.gm;
   taurp = pi.taurp;
   taum = pi.taum;
   calpha = pi.calpha;
   cbeta = pi.cbeta;
   cgamma = pi.cgamma;

   gAMPAext = pi.gAMPAext;
   gAMPArec = pi.gAMPArec;
   gNMDA = pi.gNMDA;
   gGABA = pi.gGABA;
   tauAMPA = pi.tauAMPA;
   tauGABA = pi.tauGABA;
   tauNMDAdecay = pi.tauNMDAdecay;
   tauNMDArise = pi.tauNMDArise;

   //precalculations
   bin[0][0] = bin[1][0] = bin[1][1] = 1.0;
   fak2[0] = fak2[1] = 1.0;
   alphatauNMDArisen[0] = 1.0;
   alphatauNMDArisen[1] = -1.0 * calpha * tauNMDArise;

   for (int n = 2; n <= NMAX; n++)
   {
      fak2[n] = fak2[n - 1] * (double) (n);

      bin[n][0] = 1.0;
      for (int k = 1; k < n; k++)
         bin[n][k] = bin[n - 1][k - 1] + bin[n - 1][k];
      bin[n][n] = 1.0;

      alphatauNMDArisen[n] = -1.0 * calpha * tauNMDArise *
              alphatauNMDArisen[n - 1];
   }

   TEext = (gAMPAext * brain::Cext * tauAMPA) / gm;
   TEAMPA = (gAMPArec * brain::C * tauAMPA) / gm;
   TEI = (gGABA * brain::C * tauGABA) / gm;
   tauNMDA = calpha * tauNMDArise * tauNMDAdecay;

   crho1 = (gNMDA * brain::C) / gm;
   crho2 = cbeta * crho1;
   csigmaE =
           (pow(gAMPAext, 2) * brain::Cext * pow(tauAMPA, 2)) / pow(gm * taum,
           2);

   for (int i = 0; i < 15000; i++)
   {
      psidata[i] = psi(i * 0.02 / 1000);
   }
   generateintegraldata();

   //Spiking
   time_t t1 = time(NULL);
   sd = t1 % 32000; // take the modulus
   if (nr == 0) // show the seed just for pool 0: it remains the same for each pool
   {
      cout << "seed " << sd << endl;
   }
   srand(sd);



   expdttauAMPA = exp(-brain::dt / tauAMPA);
   expdttauNMDArise = exp(-brain::dt / tauNMDArise);
   expdttauGABA = exp(-brain::dt / tauGABA);

   resetSpikingPool();
}

void pool::resetSpikingPool()
{
   for (int i = 0; i < maxpoolneurons; i++)
   {
      sAMPAext[i] = 0;
      sAMPArec[i] = 0;
      sNMDA[i] = 0;
      xNMDA[i] = 0;
      sGABA[i] = 0;
      Ca[i] = 0;
      u[i] = U;
      x[i] = 1;
      V[i] = Vreset + drand() * (Vthr - Vreset);
      dV[i] = 0;
      lastspiketime[i] = taurp + 1;
   }
}

double pool::drand()
{
   //pseudo drand
   int i = rand();

   int j = RAND_MAX;

   double result = double (i) / double (j);


   return (result);
}

int pool::extspike()
{
   if (drand() < (nuext * brain::Cext * brain::dt))
      return 1;
   else
      return 0;
}


void pool::calcSynapses()
{
   brain::sumAMPArec[nr] = 0;
   brain::sumNMDA[nr] = 0;
   brain::sumGABA[nr] = 0;

   for (int i = 0; i < neurons; i++)
   {                            // Loop over all neurons


      sAMPAext[i] = sAMPAext[i] * expdttauAMPA;
      if (extspike() == 1)
      {
         sAMPAext[i] += exp(-brain::dt * drand() / tauAMPA);
      }

      double dsold =
              (-sNMDA[i] / tauNMDAdecay) + calpha * xNMDA[i] * (1 - sNMDA[i]);
      double sNMDAeuler = sNMDA[i] + brain::dt * dsold;

      sAMPArec[i] = sAMPArec[i] * expdttauAMPA;
      xNMDA[i] = xNMDA[i] * expdttauNMDArise;
      sGABA[i] = sGABA[i] * expdttauGABA;
      sNMDA[i] =
              sNMDA[i] + brain::dt * 0.5 * (dsold -
              (sNMDAeuler / tauNMDAdecay) + calpha * xNMDA[i] * (1 -
                      sNMDAeuler));

      brain::sumAMPArec[nr] += sAMPArec[i] * u[i] * x[i]; // x[i] = 1; no synaptic depression
      brain::sumNMDA[nr] += sNMDA[i] * u[i] * x[i];
      brain::sumGABA[nr] += sGABA[i];
   }
}

void pool::calcPotentials(double zeit, ostream * spout, int flprint)
{
   //Calculate Incoming Synaptic Potentials
   double inputAMPArec = 0, inputNMDA = 0, inputGABA = 0;

   for (int i = 0; i < brain::PoolCount; i++)
   {
      inputAMPArec += w[i].ampa * brain::sumAMPArec[i];
      inputNMDA += w[i].nmda * brain::sumNMDA[i];
      inputGABA += w[i].gaba * brain::sumGABA[i];
   }
   sumu = 0;
   for (int i = 0; i < neurons; i++)
   {     // Loop over all neurons
      Ca[i] = Ca[i] - brain::dt * Ca[i] / taoCa;
      u[i] = u[i] + brain::dt * (U - u[i]) / taoF;
      sumu = sumu + u[i] * x[i];
      if (lastspiketime[i] > taurp)
      {
         //Heun Algorithm for new potential

         double Veuler = V[i] + brain::dt * dV[i];

         V[i] = V[i] + brain::dt * 0.5 * (dV[i] - (gm / Cm) * (Veuler - VL)
                 - (gAMPAext * (Veuler - VE) * sAMPAext[i]
                         + gAMPArec * (Veuler - VE) * inputAMPArec
                         + gNMDA * (Veuler - VE) * inputNMDA / (1. +
                                 cgamma * exp(-cbeta * Veuler)) +
                         gGABA * (Veuler - VI) * inputGABA +
                         AHP_Ca * gAHP * Ca[i] * (Veuler - VK)) / Cm);

         dV[i] = -(gm / Cm) * (V[i] - VL)
                 - (gAMPAext * (V[i] - VE) * sAMPAext[i]
                 + gAMPArec * (V[i] - VE) * inputAMPArec
                 + gNMDA * (V[i] - VE) * inputNMDA / (1. +
                         cgamma * exp(-cbeta * V[i])) + gGABA * (V[i] -
                         VI) * inputGABA + gAHP * (V[i] -
                         VK) * AHP_Ca * Ca[i]) / Cm;
         // Spike Generation
         if (V[i] > Vthr)
         {
            if (flprint == 1)
               if (i < 10 )
                  *spout << zeit << " " <<nr <<" " << i << endl;

            brain::spikes[nr]++;
            lastspiketime[i] = 0;
            sAMPArec[i] += 1.;
            xNMDA[i] += 1.;
            sGABA[i] += 1.;
            V[i] = Vreset;
            Ca[i] += alphaCa;
            u[i] += U * (1 - u[i]);

         }
      }
      lastspiketime[i] += brain::dt;
   }
   brain::sumup[nr] = brain::sumup[nr] + sumu;
}

//=======================
//---- meanfield section
//=======================

void pool::initVar()
{


   vnx = nx();
   vNx = Nx();
   vnix = nix();
   vSx = Sx();
   vtauE = tauE();
   vmuE = muE();
   vsigmaE = sigmaE();
}

double pool::Euler(double eulerdelta)
{


   initVar();

   return (brain::poolnu[nr] + eulerdelta * (-1 * brain::poolnu[nr] + Phi2())
           / vtauE);
}

double pool::Flow()
{


   initVar();

   return (Phi2());
}

double pool::Phi2()
{

   double b = fbeta();

   double a = falpha();

   if ((b < -20) || (a > 10))
   {
      cout << "Integral Range exceeds precalculated values" << endl;
      exit(1);
   }

   int xi = (int) ((b + 20) / 0.002);

   int yi = (int) ((a + 20) / 0.002);

   double xf = ((b + 20) - (xi * 0.002)) / 0.002;

   double yf = ((a + 20) - (yi * 0.002)) / 0.002;

   double sum = 0;

   int start, end;

   if (xf > 0.5)
   {
      sum += intdata[xi + 1] * (1.5 - xf);
      start = xi + 2;
   }
   else
   {
      sum += intdata[xi] * (0.5 - xf);
      start = xi + 1;
   }

   if (yf > 0.5)
   {
      sum += intdata[yi + 1] * (yf - 0.5);
      end = yi;
   }
   else
   {
      sum += intdata[yi] * (yf + 0.5);
      end = yi - 1;
   }


   int inter1 = start + (50 - (start % 50));

   int inter2 = end - (end % 50);

   int broad1 = (int) (start / 50) + 1;

   int broad2 = (int) (end / 50);

   for (int i = start; i < inter1; i++)
      sum += intdata[i];

   for (int i = broad1; i < broad2; i++)
      sum += intdatabroad[i];

   for (int i = inter2; i <= end; i++)
      sum += intdata[i];

   sum *= sqrt(PI) * 0.002;
   return (1 / (taurp + vtauE * sum));
}


double pool::Phi()
{
   double a = falpha();

   double b = fbeta();

   double z;

   int N = 1000;

   double sum = 0;

   for (int i = 0; i <= N; i++)
   {
      z = b + (a - b) * i / N;
      if (i == 0 | i == N)
         sum += 0.5 * nerf(z);
      else
         sum += nerf(z);
   }

   sum *= sqrt(PI) * fabs(a - b) / N;

   return (1 / (taurp + vtauE * sum));
}

//=======================
double pool::falpha()
{
   return (Vthr - vmuE) / vsigmaE * (1.0 + 0.5 * tauAMPA / vtauE) + 1.03 *
           sqrt(tauAMPA / vtauE) - 0.5 * tauAMPA / vtauE;
}

double pool::fbeta()
{
   return (Vreset - vmuE) / vsigmaE;
}

double pool::tauE()
{
   return Cm / (gm * vSx);
}

double pool::nerf(double z)
{
   double t, ef, at;

   double w;

   w = fabs(z);
   t = 1. / (1. + 0.5 * w);
   at = a1 + t * (a2 + t * (a3 + t * (a4 + t * (a5 + t * (a6 + t * (a7 +
                                                   t * (a8 + t * (a9 +
                                                                   t *
                                                                   a10))))))));
   ef = t * exp(at);
   if (z > 0.)
      ef = 2. * exp(w * w) - ef;
   return (ef);
}

//=======================
double pool::muE()
{
   double res;

   if (VE == 0)
      res = (rho2() * vNx * avrV + TEI * vnix * VI + VL) / vSx;
   else
      res = ((TEext * nuext + TEAMPA * vnx + rho1() * vNx) * VE +
              rho2() * vNx * avrV + TEI * vnix * VI + VL) / vSx;
   return res;
}

double pool::sigmaE()
{

   return sqrt(pow(avrV - VE, 2) * vtauE * csigmaE * nuext);
}

double pool::Sx()
{
   return 1 + TEext * nuext + TEAMPA * vnx + (rho1() + rho2()) * vNx +
           TEI * vnix;
}

//=======================

double pool::rho1()
{

   return crho1 / J();
}

double pool::rho2()
{

   return crho2 * (avrV - VE) * (J() - 1) / pow(J(), 2);
}

void pool::newavrV()
{

   avrV = vmuE - (Vthr - Vreset) * brain::poolnu[nr] * vtauE;

}

double pool::nx()
{
   double res = 0;

   for (int i = 0; i < brain::PoolCount; i++)
   {
      res += brain::poolf[i] * w[i].ampa * brain::poolnu[i];

   }
   return res;
}

double pool::Nx()
{
   double res = 0;

   for (int i = 0; i < brain::PoolCount; i++)
   {
      res += brain::poolf[i] * w[i].nmda * psi2(brain::poolnu[i]);
   }
   return res;
}

double pool::nix()
{
   double res = 0;

   for (int i = 0; i < brain::PoolCount; i++)
   {
      if (w[i].gaba != 0)
         res += brain::poolf[i] * w[i].gaba * brain::poolnu[i];

   }
   return res;
}

//=======================
double pool::J()
{
   return 1 + cgamma * exp(-cbeta * avrV);
}

double pool::psi2(double nu)
{
   int a = (int) (nu * 50000);

   if (nu < 0.3)
      return (psidata[a] + (psidata[a + 1] - psidata[a]) * (nu -
                      (a * 0.00002)) * 50000);
   else
      return psi(nu);
}

double pool::psi(double nu)
{
   double sum1 = 0;

   onenutauNMDA = 1 + nu * tauNMDA;
   for (int n = 1; n < 7; n++)
   {

      sum1 += alphatauNMDArisen[n] * T(n) / fak2[n + 1];
   }
   return (nu * tauNMDA) / onenutauNMDA * (1 + (1 / onenutauNMDA) * sum1);
}

//=======================

double pool::T(int n)
{

   double sum = 0;

   for (int k = 0; k <= n; k++)
   {
      if (k % 2 == 0)
         sum += bin[n][k] * (tauNMDArise * onenutauNMDA) / (tauNMDArise *
                 onenutauNMDA + k * tauNMDAdecay);
      else
         sum -= bin[n][k] * (tauNMDArise * onenutauNMDA) / (tauNMDArise *
                 onenutauNMDA + k * tauNMDAdecay);
   }
   return sum;
}

void pool::generateintegraldata()
{
   int c = 0;

   int d = 0;

   for (int i = 0; i < 15000; i++)
   {
      intdata[i] = nerf(-20. + i * 0.002);
      if (c < 50)
         intdatabroad[d] += intdata[i];
      else
      {
         c = 0;
         d++;
         intdatabroad[d] = intdata[i];
      }
      c++;
   }
}
