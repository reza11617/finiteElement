#ifndef TIMER_H
#define TIMER_H

#include <iostream>
#include <chrono>
#include <string>

class Timer {
public:
  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::duration<float> duration;
  std::string message;
  Timer();
  Timer(const std::string &);
  ~Timer();
  //void endHere();
};

#endif
