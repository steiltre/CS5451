#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * * @brief Write an array of integers to a file.
 * *
 * * @param filename The name of the file to write to.
 * * @param numbers The array of numbers.
 * * @param nnumbers How many numbers to write.
 * */
static void print_numbers(
            char const * const filename,
                uint32_t const * const numbers,
                    uint32_t const nnumbers)
{
      FILE * fout;

        /* open file */
        if((fout = fopen(filename, "w")) == NULL) {
                fprintf(stderr, "error opening '%s'\n", filename);
                    abort();
                      }

          /* write the header */
          fprintf(fout, "%d\n", nnumbers);

            /* write numbers to fout */
            for(uint32_t i = 0; i < nnumbers; ++i) {
                    fprintf(fout, "%d\n", numbers[i]);
                      }

              fclose(fout);
}

/* Gives us high-resolution timers. */
#define _POSIX_C_SOURCE 200809L
#include <time.h>

/* OSX timer includes */
#ifdef __MACH__
  #include <mach/mach.h>
  #include <mach/mach_time.h>
#endif

/**
 * * @brief Return the number of seconds since an unspecified time (e.g., Unix
 * *        epoch). This is accomplished with a high-resolution monotonic timer,
 * *        suitable for performance timing.
 * *
 * * @return The number of seconds.
 * */
static inline double monotonic_seconds()
{
#ifdef __MACH__
      /* OSX */
      static mach_timebase_info_data_t info;
        static double seconds_per_unit;
          if(seconds_per_unit == 0) {
                  mach_timebase_info(&info);
                      seconds_per_unit = (info.numer / info.denom) / 1e9;
                        }
            return seconds_per_unit * mach_absolute_time();
#else
              /* Linux systems */
              struct timespec ts;
                clock_gettime(CLOCK_MONOTONIC, &ts);
                  return ts.tv_sec + ts.tv_nsec * 1e-9;
#endif
}

/**
 * * @brief Output the seconds elapsed while sorting. This excludes input and
 * *        output time. This should be wallclock time, not CPU time.
 * *
 * * @param seconds Seconds spent sorting.
 * */
static void print_time(
            double const seconds)
{
      printf("Sort Time: %0.04fs\n", seconds);
}
