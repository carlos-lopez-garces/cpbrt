#include "fileutil.h"

static std::string searchDirectory;

bool IsAbsolutePath(const std::string &filename) {
    if (filename.empty()) return false;
    return (filename[0] == '\\' || filename[0] == '/' ||
            filename.find(':') != std::string::npos);
}

std::string AbsolutePath(const std::string &filename) {
    char full[_MAX_PATH];
    if (_fullpath(full, filename.c_str(), _MAX_PATH)) {
        return std::string(full);
    } else {
        return filename;
    }
}

std::string ResolveFilename(const std::string &filename) {
    if (searchDirectory.empty() || filename.empty()) {
        return filename;
    } else if (IsAbsolutePath(filename)) {
        return filename;
    }

    char searchDirectoryEnd = searchDirectory[searchDirectory.size() - 1];
    if (searchDirectoryEnd == '\\' || searchDirectoryEnd == '/') {
        return searchDirectory + filename;
    } else {
        return searchDirectory + "\\" + filename;
    }
}