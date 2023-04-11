#include <stdio.h>
#include <math.h>

#define MAX_TIME 20.0
#define MASS 1.0   // Particle Mass
#define SPRING 1.0 // Spring Constant
#define INIT_X 1.0 // Initial Position
#define INIT_P 0.0 // Initial Momentum

int main()
{
  double x = INIT_X;
  double p = INIT_P;
  double t = 0.0;
  double dt, force, energy, omega;
  double exact; // exact solution x=cos(omega*t)
  FILE *fout;

  fout = fopen("leapflog.dat", "w");
  printf("Input dt: ");
  scanf("%lf", &dt); // Input from Keyboard

  omega = sqrt(SPRING / MASS);

  while (t < MAX_TIME)
  {
    t += dt;
    x += dt * p / MASS;
    force = -SPRING * x;
    p += dt * force;
    exact = cos(omega * t);
    energy = 0.5 * p * p / MASS + 0.5 * SPRING * x * x;
    fprintf(fout, "%8.5f %8.5f %8.5f %8.5f\n", t, x, exact, energy);
  }

  fclose(fout);
  return 0;
}