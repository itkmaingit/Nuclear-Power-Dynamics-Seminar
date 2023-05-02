// rdf.cpp : calculate radial distribution function using save files
//
#include <stdio.h>
#include <math.h>

#define NUM_ATOM 512
#define R_MAX 5.0
#define TOTAL_STEP 20000
#define SAVE_STEP 10
#define NUM_HIST 100

void getdata(int); // Function prototype
void ini_rdf(), calc_rdf(), out_rdf();

double posx[NUM_ATOM], posy[NUM_ATOM], posz[NUM_ATOM];
double momx[NUM_ATOM], momy[NUM_ATOM], momz[NUM_ATOM];
int nsum, hist[NUM_HIST];
double delta_r, cellx, celly, cellz, chx, chy, chz;

int main() //=================================
{
  int step;

  ini_rdf();
  for (step = TOTAL_STEP / 2; step <= TOTAL_STEP; step += SAVE_STEP)
  {
    getdata(step);
    calc_rdf();
  }
  out_rdf();

  return 0;
}

void getdata(int step) //===============================
{
  int i;
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

  fscanf(fsave, "%lf %lf %lf", &cellx, &celly, &cellz);
  chx = cellx / 2;
  chy = celly / 2;
  chz = cellz / 2;
  for (i = 0; i < NUM_ATOM; i++)
  {
    fscanf(fsave, "%lf %lf %lf", &posx[i], &posy[i], &posz[i]);
    fscanf(fsave, "%lf %lf %lf", &momx[i], &momy[i], &momz[i]);
  }
  fclose(fsave);
  return;
}
void ini_rdf() //=============== Initialize Counter
{
  int h;

  nsum = 0;
  for (h = 0; h < NUM_HIST; h++)
    hist[h] = 0;
  delta_r = R_MAX / NUM_HIST;

  return;
}
void calc_rdf() //============== Accumulate Histogram
{
  int i, j, h;
  double dx, dy, dz;

  for (i = 1; i < NUM_ATOM; i++)
  {
    for (j = 0; j < i; j++)
    {
      dx = posx[j] - posx[i]; // Periodic Boundary taken into account
      if (dx < -chx)
        dx += cellx;
      else if (dx > chx)
        dx -= cellx;
      dy = posy[j] - posy[i];
      if (dy < -chy)
        dy += celly;
      else if (dy > chy)
        dy -= celly;
      dz = posz[j] - posz[i];
      if (dz < -chz)
        dz += cellz;
      else if (dz > chz)
        dz -= cellz;
      h = (int)(sqrt(dx * dx + dy * dy + dz * dz) / delta_r);
      if (h < NUM_HIST)
        hist[h]++;
    }
  }
  nsum++;

  return;
}
void out_rdf() //============== Output RDF
{
  int h;
  double density, volume, r, rdf;
  FILE *fout;

  fout = fopen("rdf.dat", "w");

  density = NUM_ATOM / (cellx * celly * cellz);
  for (h = 0; h < NUM_HIST; h++)
  {
    volume = 4.0 * 3.141592654 / 3.0 * (3 * h * h + 3 * h + 1) * delta_r * delta_r * delta_r;
    rdf = hist[h] / (double)nsum / (0.5 * NUM_ATOM) / (density * volume);
    r = (h + 0.5) * delta_r;
    fprintf(fout, "%8.4f %8.4f %10d\n", r, rdf, hist[h]);
  }
  fclose(fout);

  return;
}