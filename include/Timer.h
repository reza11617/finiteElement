#ifndef TIMER_H
#define TIMER_H

#include <iostream>
#include <chrono>

class Timer {
public:
  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<float> duration;
  Timer();
  void endHere();
};

#endif
