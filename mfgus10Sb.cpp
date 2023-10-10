#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "brain.h"
#include "pool.h"
#include <sstream>
#include <string>

#define N 1 // default 1; multiple for number of thousands of neurons in a simulation

brain   pfc;

int     poolcount;

struct poolinfo poolinh, poolexc;

void setpooldata()
{
   poolinh.VL = -70;
   poolinh.Vthr = -50;
   poolinh.Vreset = -55;
   poolinh.VI = -70;
   poolinh.VE = 0;

   poolinh.Cm = 200;
   poolinh.gm = 20;
   poolinh.taurp = 1;
   poolinh.taum = 10;
   poolinh.calpha = 0.5;
   poolinh.cbeta = 0.062;
   poolinh.cgamma = 0.2801120448;

//synaptic conductances

   poolinh.gAMPAext = 1.62;
   poolinh.gAMPArec = 0.081 / N;
   poolinh.gNMDA = 0.258 / N;
   poolinh.gGABA = 0.973 / N;

//gating variables
   poolinh.tauAMPA = 2;
   poolinh.tauGABA = 10;
   poolinh.tauNMDAdecay = 100;
   poolinh.tauNMDArise = 2;

   poolexc.VL = -70;
   poolexc.Vthr = -50;
   poolexc.Vreset = -55;
   poolexc.VI = -70;
   poolexc.VE = 0;

   poolexc.Cm = 500;
   poolexc.gm = 25;
   poolexc.taurp = 2;
   poolexc.taum = 20;
   poolexc.calpha = 0.5;
   poolexc.cbeta = 0.062;
   poolexc.cgamma = 0.2801120448;

//synaptic conductances
   poolexc.gAMPAext = 2.08;
   poolexc.gAMPArec = 0.104 / N;
   poolexc.gNMDA = 0.327 / N;
   poolexc.gGABA = 1.25 / N;

//gating variables
   poolexc.tauAMPA = 2;
   poolexc.tauGABA = 10;
   poolexc.tauNMDAdecay = 100;
   poolexc.tauNMDArise = 2;
}

// Set connections with normalisation
// Pool 0 is Inhib; all other pools as Specific
void setconnections(brain & pfc, double wplus, double inh)
{
   double  wminus;

   wminus = 1. - 0.1 * (wplus - 1.) / (1. - 0.1);

// Set everything to 1
   for (int i = 0; i < poolcount; i++)
   {
      pfc.setConnection(0, i, 0, 0, 1);
      for (int j = 1; j < poolcount; j++)
      {
         pfc.setConnection(j, i, 1.0, 1.0, 0);
      }
   }

   pfc.setConnection(0, 1, 0, 0, inh); 
   pfc.setConnection(0, 2, 0, 0, inh);
   pfc.setConnection(0, 3, 0, 0, inh);
   pfc.setConnection(0, 4, 0, 0, inh);
   pfc.setConnection(0, 5, 0, 0, inh);
   pfc.setConnection(0, 6, 0, 0, inh);
   pfc.setConnection(0, 7, 0, 0, inh);
   pfc.setConnection(0, 8, 0, 0, inh);
   pfc.setConnection(0, 9, 0, 0, inh);
   pfc.setConnection(0, 10, 0, 0, inh); 

//Selective pools
//AMPA
   for (int i = 1; i < poolcount; i++)
   {

      for (int j = 1; j < (poolcount); j++) 
      {
         pfc.setConnection(i, j, wminus, wminus, 0);
      }
      if (i < (poolcount)) 
         pfc.setConnection(i, i, wplus, wplus, 0);
   }

}



inline std::string stringify_extension(double i)
{
	
	ostringstream temp;
	temp.width(4);
	temp.fill('0');
	temp << i;
	return temp.str();
}


//////////////////////////
// 
//   Parameter space 
//
//////////////////////////



void run_spikes(brain & pfc, double aa)
{
   int     flag, i;

   ofstream dateiausgabe;
   string rate_file = "rate";
   string rate_file_ext = rate_file + "_" +
		    stringify_extension(aa);


   dateiausgabe.open(rate_file_ext.c_str(), ios::out);


   ofstream spout;
   string spikes_file = "spike";
   string spikes_file_ext = spikes_file + "_" +
		    stringify_extension(aa);

   spout.open(spikes_file_ext.c_str(), ios::out);

   ofstream fsumu;
   string sumu_file = "sumu";
   string sumu_file_ext = sumu_file + "_" +
		    stringify_extension(aa);
   
   fsumu.open(sumu_file_ext.c_str(), ios::out);


   flag = 1;



   pfc.flagavr = 5;
   pfc.resetavr();
   pfc.resetPools();

   for (i = 0; i < poolcount; i++)
   {

      pfc.setExtrate(i, 0.003 + 0.04 / 800);
   }
   pfc.Spiking(0, 500, &dateiausgabe, &spout, &fsumu, 100, flag); // spontaneous period

   for (int k = 1; k <2 ; k++) 
   {
      pfc.setExtrate(k, 0.003 + 0.25 / 800); 
   }
   pfc.Spiking(500, 1000, &dateiausgabe, &spout, &fsumu, 100, flag); // cue application period


   for (int k = 0; k < poolcount; k++) 
   {
      pfc.setExtrate(k, 0.003 + 0.04 / 800.); 
   }
   pfc.Spiking(1500, 1000, &dateiausgabe, &spout, &fsumu, 100, flag); 
   
   
   for (int k = 2; k <3 ; k++) 
   {
      pfc.setExtrate(k, 0.003 + 0.25 / 800); 
   }
   pfc.Spiking(2500, 1000, &dateiausgabe, &spout, &fsumu, 100, flag);
   
   
   for (int k = 0; k < poolcount; k++) 
   {
      pfc.setExtrate(k, 0.003 + 0.04 / 800.); 
   }
   pfc.Spiking(3500, 1000, &dateiausgabe, &spout, &fsumu, 100, flag); 

for (int k = 3; k <4 ; k++) 
   {
      pfc.setExtrate(k, 0.003 + 0.25 / 800); 
   }
   pfc.Spiking(4500, 1000, &dateiausgabe, &spout, &fsumu, 100, flag); 
   
   
   for (int k = 0; k < poolcount; k++) 
   {
      pfc.setExtrate(k, 0.003 + 0.04 / 800.); 
   }
   pfc.Spiking(5500, 1000, &dateiausgabe, &spout, &fsumu, 100, flag); 
   
   for (int k = 4; k <5 ; k++) 
   {
      pfc.setExtrate(k, 0.003 + 0.25 / 800); 
   }
   pfc.Spiking(6500, 1000, &dateiausgabe, &spout, &fsumu, 100, flag);
   
   for (int k = 0; k < poolcount; k++) 
   {
      pfc.setExtrate(k, 0.003 + 0.04 / 800.); 
   }
   pfc.Spiking(7500, 1000, &dateiausgabe, &spout, &fsumu, 100, flag); 
   
   for (int k = 5; k <6 ; k++) 
   {
      pfc.setExtrate(k, 0.003 + 0.25 / 800); 
   }
   pfc.Spiking(8500, 1000, &dateiausgabe, &spout, &fsumu, 100, flag);
   
   for (int k = 0; k < poolcount; k++) 
   {
      pfc.setExtrate(k, 0.003 + 0.04 / 800.); 
   }
   pfc.Spiking(9500, 1000, &dateiausgabe, &spout, &fsumu, 100, flag); 
   
   for (int k = 6; k <7 ; k++) 
   {
      pfc.setExtrate(k, 0.003 + 0.25 / 800); 
   }
   pfc.Spiking(10500, 1000, &dateiausgabe, &spout, &fsumu, 100, flag);
   
   for (int k = 0; k < poolcount; k++) 
   {
      pfc.setExtrate(k, 0.003 + 0.04 / 800.); 
   }
   pfc.Spiking(11500, 1000, &dateiausgabe, &spout, &fsumu, 100, flag); 
   
   for (int k = 7; k <8 ; k++) 
   {
      pfc.setExtrate(k, 0.003 + 0.25 / 800); 
   }
   pfc.Spiking(12500, 1000, &dateiausgabe, &spout, &fsumu, 100, flag);
   
   for (int k = 0; k < poolcount; k++) 
   {
      pfc.setExtrate(k, 0.003 + 0.04 / 800.); 
   }
   pfc.Spiking(13500, 1000, &dateiausgabe, &spout, &fsumu, 100, flag); 
   
   for (int k = 8; k <9 ; k++)
   {
      pfc.setExtrate(k, 0.003 + 0.25 / 800); 
   }
   pfc.Spiking(14500, 1000, &dateiausgabe, &spout, &fsumu, 100, flag);
   
   for (int k = 0; k < poolcount; k++)
   {
      pfc.setExtrate(k, 0.003 + 0.04 / 800.); 
   }
   pfc.Spiking(15500, 1000, &dateiausgabe, &spout, &fsumu, 100, flag); 
   
   for (int k = 9; k <10 ; k++) 
   {
      pfc.setExtrate(k, 0.003 + 0.25 / 800); 
   }
   pfc.Spiking(16500, 1000, &dateiausgabe, &spout, &fsumu, 100, flag);
   
   for (int k = 0; k < poolcount; k++) 
   {
      pfc.setExtrate(k, 0.003 + 0.04 / 800.); 
   }
   pfc.Spiking(17500, 3000, &dateiausgabe, &spout, &fsumu, 100, flag); 
   
   dateiausgabe.close();
   spout.close();
   fsumu.close();
}

int main()
{
   //Initialize pfc
   poolcount = 11;
   setpooldata();
   pfc.initBrain(poolcount, (int) (1000 * N), 800, 0.05);
   pfc.initPool(0, 0.2, 0.009, -52.5, 0.003, poolinh);
   pfc.initPool(1, 0.08, 0.003, -52.5, 0.003, poolexc);
   pfc.initPool(2, 0.08, 0.003, -52.5, 0.003, poolexc);
   pfc.initPool(3, 0.08, 0.003, -52.5, 0.003, poolexc);
   pfc.initPool(4, 0.08, 0.003, -52.5, 0.003, poolexc);
   pfc.initPool(5, 0.08, 0.003, -52.5, 0.003, poolexc);
   pfc.initPool(6, 0.08, 0.003, -52.5, 0.003, poolexc);
   pfc.initPool(7, 0.08, 0.003, -52.5, 0.003, poolexc);
   pfc.initPool(8, 0.08, 0.003, -52.5, 0.003, poolexc);
   pfc.initPool(9, 0.08, 0.003, -52.5, 0.003, poolexc);
   pfc.initPool(10, 0.08, 0.003, -52.5, 0.003, poolexc);

   
     pfc.printPool(6);
   
  setconnections(pfc, 2.3, 0.97); 
  run_spikes(pfc,0.97);  

}                              
