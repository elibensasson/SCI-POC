#include <string>

#ifndef  __Logging_HPP
#define  __Logging_HPP

namespace PCP_Project {

/********************************************************/
/*******************     Logging    *********************/
/********************************************************/
/**
 * The Logging class implements the functionality of writing to a logfile
 */

//TODO: Very strange to put the code in the constructor. Discuss. --Eran
class Logger {
	public:
		Logger(const std::string& msg, const std::string& fileName = "_COMMON_log.txt");
		static void logTime(const std::string& fileName = "_COMMON_log.txt");
}; // class Logger

} // namespace PCP_Project 

#endif   // __Logging_HPP
