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

  output = fopen("D2h-VPCM-C_host.log", "w+");
  if (!pcmsolver_is_compatible_library())
  {
    fprintf(stderr, "%s\n", "PCMSolver library not compatible");
    exit(EXIT_FAILURE);
  }

  fprintf(output, "%s\n", "Starting a PCMSolver calculation");
  // Use C2H4_D2h in D2h symmetry
  double charges[NR_NUCLEI] = {6.0};
  double coordinates[3 * NR_NUCLEI] = { 0.0, 0.0, 0.0};
  // This means the molecular point group has three generators:
  // the Oxy, Oxz and Oyz planes
  int symmetry_info[4] = {3, 4, 2, 1};
  struct PCMInput host_input = pcmsolver_input();

  pcmsolver_context_t * pcm_context = pcmsolver_new(PCMSOLVER_READER_OWN,
      NR_NUCLEI, charges, coordinates, symmetry_info, &host_input);

  pcmsolver_print(pcm_context);

  size_t grid_size = pcmsolver_get_cavity_size(pcm_context);
  size_t irr_grid_size = pcmsolver_get_irreducible_cavity_size(pcm_context);
  double * grid = (double *) calloc(3*grid_size, sizeof(double));
  pcmsolver_get_centers(pcm_context, grid);

  double * mep = nuclear_mep(NR_NUCLEI, charges, coordinates, grid_size, grid);
  const char * mep_lbl = {"NucMEP"};
  pcmsolver_set_surface_function(pcm_context, grid_size, mep, mep_lbl);
  pcmsolver_print_surface_function(pcm_context, mep_lbl);
  const char * asc_lbl = {"NucASC"};
  // This is the Ag irreducible representation (totally symmetric)
  int irrep = 0;
  // Get initial guess
  double nuc_chg = 6.0;
  pcmsolver_compute_initial_guess_asc(pcm_context, mep_lbl, asc_lbl, nuc_chg, irrep);
  pcmsolver_print_surface_function(pcm_context, asc_lbl);

  double * asc = (double *) calloc(grid_size, sizeof(double));
  pcmsolver_get_surface_function(pcm_context, grid_size, asc, asc_lbl);
  double energy = pcmsolver_compute_polarization_energy(pcm_context, mep_lbl, asc_lbl);

  // Compute error (aka the residual changed of sign)
  const char * err_lbl = {"ErrASC"};
  pcmsolver_compute_error_asc(pcm_context, mep_lbl, asc_lbl, err_lbl, irrep);
  pcmsolver_print_surface_function(pcm_context, err_lbl);
  double * error = (double *) calloc(grid_size, sizeof(double));
  pcmsolver_get_surface_function(pcm_context, grid_size, error, err_lbl);
  double residual_norm = norm(grid_size, error);

  int iteration = 0;
  int max_it = 5;
  fprintf(output, "  Iteration          Residual norm          Polarization energy\n");
  while(residual_norm >= 1.0e-05)
  {
    fprintf(output, "    %2i       %20.10f          %20.12f\n", iteration, residual_norm, energy);
    if (iteration >= max_it) break;
    // Compute update
    pcmsolver_compute_update_asc(pcm_context, asc_lbl, err_lbl, irrep);
    pcmsolver_print_surface_function(pcm_context, asc_lbl);
    // Compute error
    pcmsolver_compute_error_asc(pcm_context, mep_lbl, asc_lbl, err_lbl, irrep);
    pcmsolver_print_surface_function(pcm_context, err_lbl);
    // Compute residual norm
    pcmsolver_get_surface_function(pcm_context, grid_size, error, err_lbl);
    residual_norm = norm(grid_size, error);
    energy = pcmsolver_compute_polarization_energy(pcm_context, mep_lbl, asc_lbl);
    iteration += 1;
  }

  // Check that everything calculated is OK
  // Cavity size
  const size_t ref_size = 576;
  if (grid_size != ref_size) {
    fprintf(stderr, "%s\n", "Error in the cavity size, please file an issue on: https://github.com/PCMSolver/pcmsolver");
    exit(EXIT_FAILURE);
  } else {
    fprintf(output, "%s\n", "Test on cavity size: PASSED");
  }
  // Irreducible cavity size
  const size_t ref_irr_size = 72;
  if (irr_grid_size != ref_irr_size) {
    fprintf(stderr, "%s\n", "Error in the irreducible cavity size, please file an issue on: https://github.com/PCMSolver/pcmsolver");
    exit(EXIT_FAILURE);
  } else {
    fprintf(output, "%s\n", "Test on irreducible cavity size: PASSED");
  }
  // Polarization energy
  const double ref_energy = -0.437960027982;
  if (!check_unsigned_error(energy, ref_energy, 1.0e-7)) {
    fprintf(stderr, "%s\n", "Error in the polarization energy, please file an issue on: https://github.com/PCMSolver/pcmsolver");
    exit(EXIT_FAILURE);
  } else {
    fprintf(output, "%s\n", "Test on polarization energy: PASSED");
  }
  // Surface functions
  //test_surface_functions(output, grid_size, mep, asc_Ag, asc_B3g, asc_neq_B3g);

  pcmsolver_write_timings(pcm_context);

  pcmsolver_delete(pcm_context);

  free(grid);
  free(mep);
  free(asc);
  free(error);

  fclose(output);

  return EXIT_SUCCESS;
}
