#ifndef LOG_H
#define LOG_H

#include <iostream>
#include <string>
#include <vector>


class Log {
public:
  enum Level
    {
      LevelError = 0, LevelWarning, LevelInfo
    };
private:
  Level m_LogLevel;

public:
  Log();
  Log(Level);
  void setLevel(Level);
  void Error(const std::string&);
  void Warning(const std::string&);
  void Info(const auto& message)
  {
    if (m_LogLevel >= LevelInfo)
      std::cout << "\033[1;37m[Info]: \033[0m" << message << std::endl;
  };

  void PrintArray(auto arrayToPrint[], auto size)
  {
    if (m_LogLevel >= LevelInfo)
      {
	auto s = size > 20 ? 20 : size;
	std::cout << "\033[1;34m[Array]: \033[0m" << std::endl;
	for (unsigned int i = 0; i < s; i++)
	  std::cout << "a[" << i << "]= " << arrayToPrint[i] << std::endl;
      }
  };
  
  void PrintVector(std::vector<auto> &vectorToPrint)
  {
    unsigned int printSize = 20;
    if (m_LogLevel >= LevelInfo)
    {
      int size = vectorToPrint.size()>printSize ? printSize: vectorToPrint.size();
      std::cout << "\033[1;34m[Vector]: \033[0m " << vectorToPrint.size() << std::endl;
      for (unsigned int i = 0; i < size; i++)
        std::cout << "v[" << i << "]= " << vectorToPrint.at(i) << std::endl;
    }
  };

  void PrintVector(std::vector<auto> &vectorToPrint, auto printSize)
  {
    if (m_LogLevel >= LevelInfo)
    {
      int size = vectorToPrint.size()>printSize ? printSize: vectorToPrint.size();
      std::cout << "\033[1;34m[Vector]: \033[0m " << vectorToPrint.size() << std::endl;
      for (unsigned int i = 0; i < size; i++)
        std::cout << "v[" << i << "]= " << vectorToPrint.at(i) << std::endl;
    }
  };

  void PrintVector2(std::vector<std::vector<auto>> &vectorToPrint)
  {
    if (m_LogLevel >= LevelInfo)
    {
      unsigned int printSize = 4;
      int size = vectorToPrint.size()>printSize ? printSize: vectorToPrint.size();
      std::cout << "\033[1;34m[Vector2]: \033[0m" << vectorToPrint.size() << std::endl;
      for (unsigned int i = 0; i < size; i++)
        PrintVector(vectorToPrint[i]);
    }
  };

  void PrintVector2(std::vector<std::vector<auto>> &vectorToPrint,  auto printSize)
  {
    if (m_LogLevel >= LevelInfo)
    {
      int size = vectorToPrint.size()>printSize ? printSize: vectorToPrint.size();
      std::cout << "\033[1;34m[Vector2]: \033[0m" << vectorToPrint.size() << std::endl;
      for (unsigned int i = 0; i < size; i++)
        PrintVector(vectorToPrint[i]);
    }
  };

};

#endif
