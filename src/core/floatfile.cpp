#include "floatfile.h"
#include <ctype.h>
#include <stdlib.h>

bool ReadFloatFile(const char *filename, std::vector<Float> *values) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        Error("Unable to open file \"%s\"", filename);
        return false;
    }

    int c;
    bool inNumber = false;
    char curNumber[32];
    int curNumberPos = 0;
    int lineNumber = 1;
    while ((c = getc(f)) != EOF) {
        if (c == '\n') ++lineNumber;
        if (inNumber) {
            // TODO: define.
            // CHECK_LT(curNumberPos, (int)sizeof(curNumber))
            //     << "Overflowed buffer for parsing number in file: " << filename
            //     << ", at line " << lineNumber;
            if (isdigit(c) || c == '.' || c == 'e' || c == '-' || c == '+')
                curNumber[curNumberPos++] = c;
            else {
                curNumber[curNumberPos++] = '\0';
                values->push_back(atof(curNumber));
                inNumber = false;
                curNumberPos = 0;
            }
        } else {
            if (isdigit(c) || c == '.' || c == '-' || c == '+') {
                inNumber = true;
                curNumber[curNumberPos++] = c;
            } else if (c == '#') {
                while ((c = getc(f)) != '\n' && c != EOF)
                    ;
                ++lineNumber;
            } else if (!isspace(c)) {
                Warning("Unexpected text found at line %d of float file \"%s\"",
                        lineNumber, filename);
            }
        }
    }
    fclose(f);
    return true;
}