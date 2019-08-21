#ifndef logger_h
#define logger_h
#include <sys/time.h>
#include "types.h"

//! Timer and memory logger
class Logger {
private:
  Event         tic;                                            //!< Timer (map from string to double)

//! Timer function
  double get_time() {
    struct timeval tv;                                          // Time value
    gettimeofday(&tv, NULL);                                    // Get time of day in seconds and microseconds
    return double(tv.tv_sec+tv.tv_usec*1e-6);                   // Combine seconds and microseconds and return
  }

public:
  bool  printNow;                                               //!< Switch to print timings
  Event timer;                                                  //!< Stores timings for all events

//! Constructor
  Logger() {
    printNow = false;                                           // Don't print by default
  }
//! Destructor
  ~Logger() {}

//! Start timer for given event
  void startTimer(std::string event) {
    tic[event] = get_time();                                    // Get time of day and store in tic
  }

//! Stop timer for given event
  double stopTimer(std::string event, bool print=false) {
    double toc = get_time();                                    // Get time of day and store in toc
    timer[event] += toc - tic[event];                           // Accumulate event time to timer
    if(print) std::cout << event << " : " << timer[event] << std::endl;// Print event and timer to screen
    return toc - tic[event];                                    // Return the event time
  }

//! Erase entry in timer
  void eraseTimer(std::string event) {
    timer.erase(event);                                         // Erase event from timer
  }
};

#endif
