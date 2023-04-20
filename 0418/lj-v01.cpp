//
// A Simple MD Code for Lennard-Jones System
//     mass=1; epsilon=1; sigma=1;
// Version 01: Prototype
//
#include <stdio.h>
#include <math.h>

#define NUM_LATTICE   5
#define NUM_ATOM    (NUM_LATTICE*NUM_LATTICE*NUM_LATTICE)

#define TOTAL_STEP  20000
#define SAVE_STEP   10

double DEL_T= 0.001;
double CELL_X= 10.0;
double CELL_Y= 10.0;
double CELL_Z= 10.0;

void initial( );            // Function Prototypes
void force( );
void move( );
void statistics(int );
void initplot( );
void finalplot( );

double posx[NUM_ATOM], posy[NUM_ATOM], posz[NUM_ATOM];
double momx[NUM_ATOM], momy[NUM_ATOM], momz[NUM_ATOM];
double frcx[NUM_ATOM], frcy[NUM_ATOM], frcz[NUM_ATOM];
double eng_kin, eng_pot;

FILE *fgnuplot, *fsave, *fout;

int main( )
{
   int step;

      initial( );
      initplot( );
      fout=fopen("lj.dat","w");

   for (step=0; step<=TOTAL_STEP; step++) {
      force( );
      move( );
      if (step%SAVE_STEP==0) {
         statistics(step);
      }
   }

      finalplot( );
      fclose(fout);
   return 0;
}

//------------------------------------------------
//   Initialize Gnuplot Command File
//
void initplot( )
{
   fgnuplot=fopen("lj.plt","w");

   fprintf(fgnuplot,"set size square\n");
   fprintf(fgnuplot,"set key box font \",24\"\n");
   fprintf(fgnuplot,"set xran [0:%f]\n",CELL_X);
   fprintf(fgnuplot,"set yran [0:%f]\n",CELL_Y);
}
//------------------------------------------------
//   Finalize Gnuplot Command File
//
void finalplot( )
{
   fclose(fgnuplot);
}
//------------------------------------------------
//   Make Simple Cubic Lattice with Random Momenta
//
void initial( )
{
   int i=0;
   int ix,iy,iz;
   double dx,dy,dz;
      dx=CELL_X/NUM_LATTICE;
      dy=CELL_Y/NUM_LATTICE;
      dz=CELL_Z/NUM_LATTICE;

   for (ix=0;ix<NUM_LATTICE;ix++) {
   for (iy=0;iy<NUM_LATTICE;iy++) {
   for (iz=0;iz<NUM_LATTICE;iz++) {
      posx[i] = dx*ix;
      posy[i] = dy*iy;
      posz[i] = dz*iz;
      momx[i] = (rand( )/(double)RAND_MAX -0.5);
      momy[i] = (rand( )/(double)RAND_MAX -0.5);
      momz[i] = (rand( )/(double)RAND_MAX -0.5);
      i++;
   } } }
}
//------------------------------------------------
//   Calculate Forces
//
void force( )
{
   int i, j;
   double dx,dy,dz,r2,ri06,ri12;
   double eng, frc;

      eng_pot=0.0;
   for (i=0; i<NUM_ATOM; i++) {
      frcx[i]=0.0;
      frcy[i]=0.0;
      frcz[i]=0.0;
   }

   for (i=0; i<NUM_ATOM; i++) {
   for (j=i+1; j<NUM_ATOM; j++) {
      dx=posx[i]-posx[j];
      dy=posy[i]-posy[j];
      dz=posz[i]-posz[j];
      r2=dx*dx+dy*dy+dz*dz;
      ri06=1.0/(r2*r2*r2);
      ri12=ri06*ri06;

      eng=4*(ri12-ri06);
      frc=4*(12*ri12-6*ri06)/r2;

      eng_pot += eng;
      frcx[i] += frc*dx;
      frcy[i] += frc*dy;
      frcz[i] += frc*dz;
      frcx[j] -= frc*dx;
      frcy[j] -= frc*dy;
      frcz[j] -= frc*dz;
   } }
}
//------------------------------------------------
//   Update Particle Positon & Momentum: Repulsive Boundaries
//
void move( )
{
   int i;

      eng_kin=0.0;

   for (i=0; i<NUM_ATOM; i++) {
      momx[i] += frcx[i]*DEL_T;
      momy[i] += frcy[i]*DEL_T;
      momz[i] += frcz[i]*DEL_T;
         eng_kin += (momx[i]*momx[i]+momy[i]*momy[i]+momz[i]*momz[i]);

      posx[i] += momx[i]*DEL_T;
         if (posx[i]<0.0) {            // Reflective Boundary
            posx[i]=-posx[i];
            momx[i]=-momx[i];}
         else if (posx[i]>CELL_X) {
            posx[i]=2*CELL_X-posx[i];
            momx[i]=-momx[i];}
      posy[i] += momy[i]*DEL_T;
         if (posy[i]<0.0) {
            posy[i]=-posy[i];
            momy[i]=-momy[i];}
         else if (posy[i]>CELL_Y) {
            posy[i]=2*CELL_Y-posy[i];
            momy[i]=-momy[i];}
      posz[i] += momz[i]*DEL_T;
         if (posz[i]<0.0) {
            posz[i]=-posz[i];
            momz[i]=-momz[i];}
         else if (posz[i]>CELL_Z) {
            posz[i]=2*CELL_Z-posz[i];
            momz[i]=-momz[i];}
   }
      eng_kin /= 2;
}
//------------------------------------------------
//   Output Statistical Data & Configuration
//
void statistics(int step)
{
   int i;
   char fname[100];

   printf("%8d %10.5f %10.5f %10.5f\n",
      step, eng_kin/NUM_ATOM, eng_pot/NUM_ATOM, (eng_kin+eng_pot)/NUM_ATOM);
   fprintf(fout,"%8d %10.5f %10.5f %10.5f\n",
      step, eng_kin/NUM_ATOM, eng_pot/NUM_ATOM, (eng_kin+eng_pot)/NUM_ATOM);

      sprintf(fname,"lj%8.8d.dat",step);
      fsave=fopen(fname,"w");
   for (i=0;i<NUM_ATOM;i++) {
      fprintf(fsave,"%8.3f %8.3f %8.3f\n",posx[i],posy[i],posz[i]);
   }
      fclose(fsave);

      fprintf(fgnuplot,"plot \"%s\" using 1:2 t \"%6.6d\" w po pt 6 ps 3 lw 2\n",fname,step);

      fprintf(fgnuplot,"pause 0.01\n");


}