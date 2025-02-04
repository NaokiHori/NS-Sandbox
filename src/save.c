#include <stdio.h> // snprintf
#include <stdbool.h> // false
#include <string.h> // strlen
#include <errno.h> // errno, EEXIST
#include <sys/stat.h> // mode_t, S_IRWXU, S_IRWXG, S_IRWXO
#include "memory.h"
#include "logger.h"
#include "domain.h"
#include "flow_field.h"
#include "./save.h"
#include "./save/snpyio.h"

#define ROOT_DIRECTORY "output/save/"

#define NDIMS 2

static int concat_dir_name(
    const size_t id,
    char ** const dir_name
){
  const char prefix[] = ROOT_DIRECTORY;
  const int ndigits = 10;
  const int nchars = strlen(prefix) + ndigits + 1;
  *dir_name = memory_alloc(nchars, sizeof(char));
  (*dir_name)[nchars - 1] = '\0';
  if (nchars - 1 != snprintf(*dir_name, nchars, "%s%0*zu", prefix, ndigits, id)) {
    LOGGER_FAILURE("snprintf returns unexpected result");
    goto abort;
  }
  return 0;
abort:
  return 1;
}

static int create_directory(
    const char dir_name[]
){
  errno = 0;
  if (0 != mkdir(dir_name, S_IRWXU | S_IRWXG | S_IRWXO)) {
    perror(dir_name);
    // treat EEXIST as expected to override present files
    if (EEXIST != errno) {
      LOGGER_FAILURE("failed to create a directory (other than 'already exist' error)");
      goto abort;
    }
  }
  return 0;
abort:
  return 1;
}

static int write_npy_file(
    const char dir_name[],
    const char dset_name[],
    const size_t ndims,
    const size_t * shape,
    const char dtype[],
    const size_t size,
    const void * data
){
  int error_code = 0;
  char * file_name = NULL;
  FILE * fp = NULL;
  size_t header_size = 0;
  // assign file_name
  {
    const char slash[] = {"/"};
    const char suffix[] = {".npy"};
    const int nchars =
      + strlen( dir_name)
      + strlen(    slash)
      + strlen(dset_name)
      + strlen(   suffix)
      + 1;
    file_name = memory_alloc(nchars, sizeof(char));
    file_name[nchars - 1] = '\0';
    if (nchars - 1 != snprintf(file_name, nchars, "%s%s%s%s", dir_name, slash, dset_name, suffix)) {
      error_code = 1;
      LOGGER_FAILURE("snprintf returns unexpected result");
      goto abort;
    }
  }
  // write npy header
  {
    errno = 0;
    fp = fopen(file_name, "w");
    if (NULL == fp) {
      perror(file_name);
      error_code = 1;
      LOGGER_FAILURE("failed to open file (attempted to write NPY header)");
      goto abort;
    }
    size_t header_size = 0;
    if (0 != snpyio_w_header(ndims, shape, dtype, false, fp, &header_size)) {
      error_code = 1;
      LOGGER_FAILURE("failed to write NPY header");
      goto abort;
    }
    fclose(fp);
  }
  // append main data
  {
    errno = 0;
    fp = fopen(file_name, "a");
    if (NULL == fp) {
      perror(file_name);
      error_code = 1;
      LOGGER_FAILURE("failed to open file (attempted to write data)");
      goto abort;
    }
    if (0 != fseek(fp, (long)header_size, SEEK_SET)) {
      error_code = 1;
      LOGGER_FAILURE("failed to move file pointer after NPY header");
      goto abort;
    }
    size_t nitems = 1;
    for (size_t dim = 0; dim < ndims; dim++) {
      nitems *= shape[dim];
    }
    if (nitems != fwrite(data, size, nitems, fp)) {
      error_code = 1;
      LOGGER_FAILURE("failed to write data");
      goto abort;
    }
  }
abort:
  memory_free(file_name);
  if (NULL != fp) {
    fclose(fp);
  }
  return error_code;
}

int save(
    const size_t id,
    const size_t step,
    const double time,
    const domain_t * const domain,
    const flow_field_t * const flow_field
) {
  int error_code = 0;
  char * dir_name = NULL;
  if (0 != concat_dir_name(id, &dir_name)) {
    error_code = 1;
    LOGGER_FAILURE("failed to concatenate directory name");
    goto abort;
  }
  if (0 != create_directory(dir_name)) {
    error_code = 1;
    LOGGER_FAILURE("failed to create a directory");
    goto abort;
  }
  const size_t nx = domain->nx;
  const size_t ny = domain->ny;
  write_npy_file(dir_name, "step", 0, NULL, "'<u8'", sizeof(size_t), &step);
  write_npy_file(dir_name, "time", 0, NULL, "'<f8'", sizeof(double), &time);
  write_npy_file(dir_name, "ux", NDIMS, (size_t [NDIMS]){ny + 2, nx + 2}, "'<f8'", sizeof(double), &flow_field->ux[0][0]);
  write_npy_file(dir_name, "uy", NDIMS, (size_t [NDIMS]){ny + 2, nx + 2}, "'<f8'", sizeof(double), &flow_field->uy[0][0]);
  write_npy_file(dir_name,  "p", NDIMS, (size_t [NDIMS]){ny + 2, nx + 2}, "'<f8'", sizeof(double), &flow_field-> p[0][0]);
abort:
  memory_free(dir_name);
  return error_code;
}

