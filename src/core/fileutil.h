#include <string>
#include <string.h>

#include "cpbrt.h"

bool IsAbsolutePath(const std::string &filename);

std::string AbsolutePath(const std::string &filename);

std::string ResolveFilename(const std::string &filename);