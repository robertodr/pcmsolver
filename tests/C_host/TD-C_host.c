#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcmsolver.h"
#include "PCMInput.h"

#include "C_host-functions.h"

#define NR_NUCLEI 1

FILE * output;

void host_writer(const char * message, size_t UNUSED(message_length))
{
  fprintf(output, "%s\n", message);
}

int main()
{

  output = fopen("TD-C_host.log", "w+");
  if (!pcmsolver_is_compatible_library())
  {
    fprintf(stderr, "%s\n", "PCMSolver library not compatible");
    exit(EXIT_FAILURE);
  }

  fprintf(output, "%s\n", "Starting a PCMSolver calculation");
  // Use C2H4 in D2h symmetry
  double charges[NR_NUCLEI] = {1.0};
  double coordinates[3 * NR_NUCLEI] = { 0.0, 0.0, 0.0 };
  // This means the molecular point group has three generators:
  // the Oxy, Oxz and Oyz planes
  int symmetry_info[4] = {0, 0, 0, 0};
  struct PCMInput host_input = pcmsolver_input();

  pcmsolver_context_t * pcm_context = pcmsolver_new(PCMSOLVER_READER_OWN,
      NR_NUCLEI, charges, coordinates, symmetry_info, &host_input);

  pcmsolver_print(pcm_context);

  size_t grid_size = pcmsolver_get_cavity_size(pcm_context);
  double * grid = (double *) calloc(3*grid_size, sizeof(double));
  pcmsolver_get_centers(pcm_context, grid);

  double * mep = (double *) calloc(grid_size, sizeof(double));
  for (size_t j = 0; j < grid_size; j++) {
    // Column-major ordering. Offsets: col_idx * nr_rows + row_idx
    double dist = pow((coordinates[0] - grid[j*3]), 2)
      + pow((coordinates[1] - grid[j*3 + 1]), 2)
      + pow((coordinates[2] - grid[j*3 + 2]), 2);
    dist = sqrt(dist);
    mep[j] += grid[j*3 + 2] / pow(dist, 3);
  }
  // Calculate potential from a point dipole along the z direction
  const char * mep_lbl = {"MEP"};
  pcmsolver_set_surface_function(pcm_context, grid_size, mep, mep_lbl);
  // This is the totally symmetric irreducible representation
  int irrep = 0;
  const char * mept_lbl = {"MEP_t"};
  const char * meptdt_lbl = {"MEP_tdt"};
  const char * asc_lbl = {"ASC"};
  const char * asct_lbl = {"ASC_t"};
  const char * asctdt_lbl = {"ASC_tdt"};
  pcmsolver_initialize_propagation(pcm_context, mep_lbl, asc_lbl, mept_lbl, asct_lbl, meptdt_lbl, asctdt_lbl, irrep);
  const double secondsToAU = 2.418884326509e-17;
  const double convertBohrToAngstrom = 0.52917721092;
  double dt = 0.2; // 4.838 as
  double total_time = 100e-15 / secondsToAU; // Total simulation time: 100 fs
  int steps = (int)total_time/dt + 1;
  // Data to calculate reference values according to Onsager model
  double radius = 1.181 * 1.10 / convertBohrToAngstrom;
  double e_0 = 35.69;
  double e_d = 1.807;
  double tau = 2000; // 48.38 fs
  // Propagate
  double t_0 = 0.0, t = 0.0;
  double energy = pcmsolver_compute_polarization_energy(pcm_context, mept_lbl, asct_lbl);
  double rf_energy = 0.0;
  for (int i = 0; i < steps; i++) {
    t = t_0 + i * dt;
    rf_energy = -2.0 * energy;
    fprintf(output, "%20.12f    %20.12f      %20.12f\n", (t * secondsToAU / 1.0e-15), reactionField(radius, e_0, e_d, tau, t), rf_energy);
    energy = pcmsolver_propagate_asc(pcm_context, mept_lbl, asct_lbl, meptdt_lbl, asctdt_lbl, dt, irrep);
  }

  pcmsolver_write_timings(pcm_context);

  pcmsolver_delete(pcm_context);

  free(grid);
  free(mep);

  fclose(output);

  return EXIT_SUCCESS;
}
