#include "Timer.h"

Timer::Timer()
{
  start = std::chrono::high_resolution_clock::now();
  message = "Time spend in this function is: ";
};

Timer::Timer(const std::string &mes)
  : message(mes)
{
  start = std::chrono::high_resolution_clock::now();
};

Timer::~Timer()
{
  end = std::chrono::high_resolution_clock::now();
  duration = end - start;
  double ms = duration.count() * 1000.0f;
  std::cout << message << ms << "ms " << std::endl;
};
