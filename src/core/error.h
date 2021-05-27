#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef CPBRT_CORE_ERROR_H
#define CPBRT_CORE_ERROR_H

#include "cpbrt.h"

// Error Reporting Declarations.

#ifdef __GNUG__
#define PRINTF_FUNC __attribute__((__format__(__printf__, 1, 2)))
#else
#define PRINTF_FUNC
#endif  // __GNUG__
void Warning(const char *, ...) PRINTF_FUNC;
void Error(const char *, ...) PRINTF_FUNC;

#endif // CPBRT_CORE_ERROR_H