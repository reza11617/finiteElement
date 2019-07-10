#ifndef LOG_H
#define LOG_H

#include <iostream>
#include <string>
#include <vector>


class Log {
public:
  static const unsigned int LevelError = 0;
  static const unsigned int LevelWarning = 1;
  static const unsigned int LevelInfo = 2;
private:
  static unsigned int m_LogLevel;

public:
  static Log& Logger()
  {
    static Log instanceOfLog;
    return instanceOfLog;
  };
  Log();
  void setLevel(unsigned int);
  void Error(const std::string&);
  void Warning(const std::string&);
  void Info(const std::string&);
  void Info(const double);
  void Info(const unsigned int *, unsigned int);
  void Info(const float *, unsigned int);
  void Info(const double *, unsigned int);
};
  
#endif
