#ifndef CPBRT_CORE_FILEUTIL_H
#define CPBRT_CORE_FILEUTIL_H

#include <string>
#include <string.h>

#include "cpbrt.h"

bool IsAbsolutePath(const std::string &filename);

std::string AbsolutePath(const std::string &filename);

std::string ResolveFilename(const std::string &filename);

#endif // CPBRT_CORE_FILEUTIL_H