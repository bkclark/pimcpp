#include "Common.h"

bool fileExists(const std::string& fileName)
{
  std::fstream fin;
  fin.open(fileName.c_str(),std::ios::in);
  if( fin.is_open() )
    {
      fin.close();
      return true;
    }
  fin.close();
  return false;
}
