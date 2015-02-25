#ifndef _TIME_BENCH_
#define _TIME_BENCH_

#include <ctime>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>

long unsigned int diff_timeval(struct timeval start, struct timeval end){
  /* struct timeval の差をmicro秒の単位で返す */
  time_t diff_sec = end.tv_sec - start.tv_sec;
  suseconds_t diff_usec = end.tv_usec - start.tv_usec;
  return ((long unsigned int)diff_sec*1000000 + diff_usec);
}

#if 0 /* usage */
{
  struct timeval t_start, t_end;
  gettimeofday(&t_start, NULL);
    {
      /* 計測したい処理 */
    }
  gettimeofday(&t_end, NULL);
  diff_timeval(t_start, t_end);
}
#endif

#endif
