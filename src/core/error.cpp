#include <mutex>
#include <stdarg.h>

#include "error.h"
#include "parser.h"
#include "stringprint.h"

template <typename... Args> static std::string StringVaprintf(const std::string &fmt, va_list args) {

    va_list argsCopy;
    va_copy(argsCopy, args);
    size_t size = vsnprintf(nullptr, 0, fmt.c_str(), args) + 1;
    std::string str;
    str.resize(size);
    vsnprintf(&str[0], size, fmt.c_str(), argsCopy);
    str.pop_back();
    return str;
}

static void processError(
    Loc *loc,
    const char *format,
    va_list args,
    const char *errorType
) {
    std::string errorString;

    if (loc) {
        errorString = StringPrintf("%s:%d:%d: ", loc->filename.c_str(),
                                   loc->line, loc->column);
    }

    errorString += errorType;
    errorString += ": ";
    errorString += StringVaprintf(format, args);

    static std::string lastError;
    static std::mutex mutex;
    std::lock_guard<std::mutex> lock(mutex);
    if (errorString != lastError) {
        // LOG(INFO) << errorString;
        fprintf(stderr, "%s\n", errorString.c_str());
        lastError = errorString;
    }
}

void Warning(const char *format, ...) {
    if (CpbrtOptions.quiet) return;
    va_list args;
    va_start(args, format);
    processError(parserLoc, format, args, "Warning");
    va_end(args);
}

void Error(const char *format, ...) {
    va_list args;
    va_start(args, format);
    processError(parserLoc, format, args, "Error");
    va_end(args);
}