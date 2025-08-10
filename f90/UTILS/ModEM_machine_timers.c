#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

/* In order to use these C routines in ModEM, place the following manually in
 your makefile:

 1. Define -DUSE_C_TIMERS
 2. Set CC env variable to your desired C compiler: e.g.: `export CC=mpicc`
 3. Manually add rules/prerequisit to your makefile by doing the following:

```
$(OBJDIR)/ModEM_machine_timers.o: UTILS/ModEM_machine_timers.c
	 $(CC) -c UTILS/ModEM_machine_timers.c -o $(OBJDIR)/ModEM_machine_timers.o
```

You will then need add `$(OBJDIR)/ModEM_machine_timers.o` to both the `OBJ`
variable (append to the end) *and* to the prerequisit of the
`$(OBJDIR)/ModEM_timing.o` rule in the makefile:

```
$(OBJDIR)/ModEM_timing.o:UTILS/ModEM_timing.f90 $(OBJDIR)/ModEM_machine_timers.o
	 $(F90) -c $(MODULE) $(FFLAGS) $(MPIFLAGS) UTILS/ModEM_timing.f90 -o $(OBJDIR)/ModEM_timing.o
```

*/

/* This file contatins the following timing functions:

    rdtscp(i)
    : Read Time-Step Counter and Processor ID

    timer_start(timer_id)
    timer_end(timer_id, sec, nsec)
    : Uses the RTC of the system to measure time

    cpseed(polls)
    : Uses rdtscp and timer_start and timer_end to caclulate the
    : speed of the system's cpu in hertz.

*/

/* In Fortran, use the following as an interface for rdtscp:
    use iso_c_binding, only : c_long
 
    interface
        subroutine rdtscp(i) bind(C)
           use iso_c_binding, only : c_long
           integer (c_long), intent(out) :: i
        end subroutine rdtscp
    end interface
 
    integer (c_long) :: tsc_start, tsc_end 
*/
    
/* In Fortran, use the following as an interfaces for rtc timer:
    use iso_c_binding, only : c_int
 
    interface
        subroutine timer_start(timer_id) bind(C)
           use iso_c_binding, only : c_int
           integer (c_int), intent(in), value :: timer_id
        end subroutine timer_start

        subroutine timer_stop(timer_id, sec, nsec) bind(C)
           use iso_c_binding, only : c_int
           integer (c_int), intent(in), value :: timer_id
           integer (c_int), intent(out) :: sec, nsec
        end subroutine timer_stop
    end interface
 
    integer (c_int) :: timer_id, sec, nsec 
*/

/* In Fortran, use the following as an interface for cspeed:

   use iso_c_binding, only : c_float, c_int
    
    interface
        real(c_float) function cspeed(polls) BIND(C)
            use iso_c_binding, only : c_float, c_int
            integer(c_int) :: polls
        end function cspeed
    end interface

    integer (c_int) :: polls
    real (c_float) :: clockSpeed
*/

#ifdef __LINUX__
void rdtscp( long *i )
{
	unsigned rax, rdx;
	asm volatile ("RDTSCP\n\t"
			"mov %%edx, %0\n\t"
			"mov %%eax, %1\n\t": "=r" (rdx), "=r" (rax));
	*i = ((unsigned long)rdx << 32) + rax;
}
#endif

#define MAX_TIMERS 10

#ifdef GETTIMEOFDAY
#include <sys/time.h>
#endif

#ifdef __MACH__
#include <mach/mach.h>
#include <mach/mach_time.h>
#include <unistd.h>
#endif

#ifdef __linux__
#include <time.h>
#endif

#ifdef GETTIMEOFDAY
struct timeval start_time[MAX_TIMERS];
struct timeval end_time[MAX_TIMERS];
#endif

#ifdef __MACH__
uint64_t start_time[MAX_TIMERS];
uint64_t end_time[MAX_TIMERS];
#endif

#ifdef AIX
timebasestruct_t start_time[MAX_TIMERS];
timebasestruct_t end_time[MAX_TIMERS];
#endif

#ifdef __linux__
struct timespec start_time[MAX_TIMERS];
struct timespec end_time[MAX_TIMERS];
#endif

void timer_start(int n)
{
#ifdef GETTIMEOFDAY
   gettimeofday(&start_time[n], NULL);
#endif

#ifdef __MACH__
   start_time[n] = mach_absolute_time();
#endif

#ifdef AIX
   read_real_time(&start_time[n], TIMEBASE_SZ);
#endif

#ifdef __linux__
   clock_gettime(CLOCK_MONOTONIC_RAW, &start_time[n]);
#endif
}

void timer_stop(int n, int *secs, int *n_secs)
{
#ifdef GETTIMEOFDAY
   gettimeofday(&end_time[n], NULL);
  
   *secs   = (int)(end_time[n].tv_sec - start_time[n].tv_sec);
   *n_secs = (int)(end_time[n].tv_usec - start_time[n].tv_usec) * 1000;

   if (*n_secs < 0)  {
      *secs   -= 1;
      *n_secs += 1000000000;
   }
#endif

#ifdef __MACH__
   uint64_t elapsed, elapsedNano;
   static mach_timebase_info_data_t sTimebaseInfo;

   end_time[n] = mach_absolute_time();

   elapsed = end_time[n] - start_time[n];

    if ( sTimebaseInfo.denom == 0 ) {
        (void) mach_timebase_info(&sTimebaseInfo);
    }

    // Do the maths. We hope that the multiplication doesn't 
    // overflow; the price you pay for working in fixed point.

    elapsedNano = elapsed * sTimebaseInfo.numer / sTimebaseInfo.denom;


   *secs   = (int)(elapsedNano / 1000000000);
   *n_secs = (int)(elapsedNano % 1000000000);
#endif

#ifdef AIX
   read_real_time(&end_time[n], TIMEBASE_SZ);
   time_base_to_time(&start_time[n], TIMEBASE_SZ);
   time_base_to_time(&end_time[n], TIMEBASE_SZ);

   *secs = end_time[n].tb_high - start_time[n].tb_high;
   *n_secs = end_time[n].tb_low - start_time[n].tb_low;

   if (*n_secs < 0)  {
      *secs   -= 1;
      *n_secs += 1000000000;
   }
#endif

#ifdef __linux__
   clock_gettime(CLOCK_MONOTONIC_RAW, &end_time[n]);

   *secs = (int)(end_time[n].tv_sec - start_time[n].tv_sec);
   *n_secs = (int)(end_time[n].tv_nsec - start_time[n].tv_nsec);

   if (*n_secs < 0)  {
      *secs   -= 1;
      *n_secs += 1000000000;
   }
#endif
}
