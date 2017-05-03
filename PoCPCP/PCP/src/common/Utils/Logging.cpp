#include "Logging.hpp"
#include <string>

#include <set>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <time.h>

namespace PCP_Project {

    Logger::Logger(const std::string& msg, const std::string& fileName) {
#ifdef _DEBUG
        static std::set<std::string> fileSet;
        if (fileSet.find(fileName) == fileSet.end()) { // File has not been accessed yet.
            fileSet.insert(fileName);
            std::remove(fileName.c_str());
        }
        std::ofstream logFile(fileName, std::ios_base::out | std::ios_base::app);
        if (logFile.fail()) {
		std::cout << "Warning: unable to open \"" << fileName << "\". Following log message not written:" << std::endl;
		std::cout << "  \"" << msg << "\"" << std::endl;
            return;
        }
        logFile << msg << std::endl;
#endif
    }

    void Logger::logTime(const std::string& fileName) {
#ifdef _DEBUG
        time_t time = std::time(0); // get time in seconds (since epoch)
        struct tm* now = localtime(&time);
        std::stringstream s;
        s   << "Time: " 
            << std::setfill('0') << std::setw(2) << now->tm_hour << ":" 
            << std::setfill('0') << std::setw(2) << now->tm_min << ":" 
            << std::setfill('0') << std::setw(2) << now->tm_sec << ", " 
            << std::setfill('0') << std::setw(2) << now->tm_mday << "/" 
            << std::setfill('0') << std::setw(2) << now->tm_mon + 1 << "/" 
            << now->tm_year + 1900 ; 
        Logger(s.str(), fileName);
#endif
    }

} // namespace PCP_Project 
