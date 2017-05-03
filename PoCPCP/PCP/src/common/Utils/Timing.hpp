#ifndef  __Timing_HPP
#define  __Timing_HPP

#include <ctime>

namespace PCP_Project {

/********************************************************/
/*******************     Timing    **********************/
/********************************************************/

/**
 * A lightweight class for elapsed-time measurements
 */
class Timer {
	private:
#ifdef __linux__
		struct timespec startTimeSpec;
#else
		clock_t startTime;
#endif
	public:
		Timer();

		/** Elapsed time in seconds */
		double getElapsed();
};

} // namespace PCP_Project 

#endif   // __Timing_HPP
