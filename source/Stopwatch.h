#ifndef STOPWATCH_H
#define	STOPWATCH_H

#include <time.h>

class Stopwatch {
public:
   double elapsedTime;

   Stopwatch() {
      reset();
   }

   ~Stopwatch() {
   }

   void reset() {
      elapsedTime = 0;
   }

   void restart() {
      reset();
      start();
   }

   void start() {
      startTime = clock();
   }

   void stop() {
      endTime = clock();
      elapsedTime += ((double) (endTime - startTime)) / CLOCKS_PER_SEC;
   }

private:
   clock_t startTime, endTime;

};

#endif	/* STOPWATCH_H */
