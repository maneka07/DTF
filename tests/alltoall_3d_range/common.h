#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#define TEST_DEBUG
#define WARMUP 5
#define NITER 100

#define TIME_DIFF(t1, t2) (((t2).tv_sec - (t1).tv_sec)*1e9 + ((t2).tv_nsec - (t1).tv_nsec))

#endif
