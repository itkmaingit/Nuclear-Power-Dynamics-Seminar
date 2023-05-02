// vaf.cpp : calculate velocity autocorrelation function
//
#include <stdio.h>
#include <math.h>
#include <boost/timer/progress_display.hpp>

#define NUM_ATOM 512
#define TOTAL_STEP 20000
#define SAVE_STEP 10
#define DEL_T 0.001
#define NUM_HIST 100
#define NUM_DATA (TOTAL_STEP / 2 / SAVE_STEP + 1)

void getdata(int), vaf(); // Function Prototype

double momx[NUM_ATOM][NUM_DATA], momy[NUM_ATOM][NUM_DATA], momz[NUM_ATOM][NUM_DATA];
int ndata;

double CELL_X;
double CELL_Y;
double CELL_Z;
double TemperatureArray[4] = {0.7, 1.0, 1.3, 2.0};
double CELL_SIZE[2] = {4.0, 8.0};

int main() //===========================================================
{
  for (int j = 0; j < 2; j++)
  {
    CELL_X = CELL_SIZE[j];
    CELL_Y = CELL_X;
    CELL_Z = CELL_X;
    for (int i = 0; i < 4; i++)
    {
      t_target = TemperatureArray[i];
      int step;

      ndata = 0;
      for (step = TOTAL_STEP / 2; step <= TOTAL_STEP && ndata < NUM_DATA; step += SAVE_STEP)
      {
        getdata(step);
        ndata++;
      }
      vaf();
    }
  }
  return 0;
}
void getdata(int step) //================================================
{
  int i;
  double posx, posy, posz;
  double cx, cy, cz;
  char fname[100];
  FILE *fsave;

  sprintf(fname, "lj%8.8d.dat", step);
  printf("%s\n", fname);
  fsave = fopen(fname, "r");
  if (fsave == NULL)
  {
    printf("Input File NOT Exist\n");
    return;
  }

  fscanf(fsave, "%lf %lf %lf", &cx, &cy, &cz);
  for (i = 0; i < NUM_ATOM; i++)
  {
    fscanf(fsave, "%lf %lf %lf", &posx, &posy, &posz);
    fscanf(fsave, "%lf %lf %lf", &momx[i][ndata], &momy[i][ndata], &momz[i][ndata]);
  }
  fclose(fsave);
  return;
}
void vaf() //====================================================
{
  int i, t1, t2, tmax, h, count[NUM_HIST];
  double vaf[NUM_HIST];
  FILE *fout;

  for (h = 0; h < NUM_HIST; h++)
  {
    vaf[h] = 0;
    count[h] = 0;
  }

  for (i = 0; i < NUM_ATOM; i++)
  {
    for (t1 = 0; t1 < ndata; t1++)
    {
      tmax = t1 + NUM_HIST;
      if (tmax > ndata)
        tmax = ndata;
      for (t2 = t1; t2 < tmax; t2++)
      {
        vaf[t2 - t1] += momx[i][t2] * momx[i][t1] + momy[i][t2] * momy[i][t1] + momz[i][t2] * momz[i][t1];
        count[t2 - t1]++;
      }
    }
  }
  sprintf(datfname, "temperature_%lf_cellsize_%lf.dat", t_target, CELL_X);
  fout = fopen(datfname, "w");
  for (h = 0; h < NUM_HIST; h++)
  {
    fprintf(fout, "%8.3f %10.4f\n", h * DEL_T * SAVE_STEP, vaf[h] / count[h]);
  }
  fclose(fout);
  return;
}