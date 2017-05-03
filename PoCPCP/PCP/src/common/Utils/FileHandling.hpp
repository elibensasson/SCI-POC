#include "common/Infrastructure/Infrastructure.hpp"

#include <fstream>
#include <sstream>

#ifndef  __FileHandling_HPP
#define  __FileHandling_HPP

namespace PCP_Project {

/********************************************************/
/******************* File Handling **********************/
/********************************************************/

/**
 * Descendant of ifstream that verifies that the file was opened correctly, and aborts otherwise.
 */
class checkedIfstream : public std::ifstream {
	private:
		void assertOpen(std::string filename) {
			if (!*this || !this->is_open()) {
				std::stringstream msg;
				msg << "Cannot open input file " << filename;
				_COMMON_FATAL(msg.str());
			}
		}
	public:
		checkedIfstream(const char* filename, std::ios_base::openmode mode = std::ios_base::in) : std::ifstream(filename, mode) {
			assertOpen(filename);
		}

		checkedIfstream(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in) : std::ifstream(filename.c_str(), mode) {
			assertOpen(filename);
		}


		void open(const char* filename, std::ios_base::openmode mode = std::ios_base::in) {
			std::ifstream::open(filename, mode);
			assertOpen(filename);
		}
		void open(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in) {
			std::ifstream::open(filename.c_str(), mode);
			assertOpen(filename);
		}
};

/**
 * Descendant of ofstream that verifies that the file was opened correctly, and aborts otherwise.
 */
class checkedOfstream : public std::ofstream {
	private:
		void assertOpen(std::string filename) {
			if (!*this || !this->is_open()) {
				std::stringstream msg;
				msg << "Cannot open output file " << filename;
				_COMMON_FATAL(msg.str());
			}
		}
	public:
		checkedOfstream(const char* filename, std::ios_base::openmode mode =  std::ios_base::out | std::ios_base::trunc) : std::ofstream(filename, mode) {
			assertOpen(filename);
		}

		checkedOfstream(const std::string& filename, std::ios_base::openmode mode =  std::ios_base::out | std::ios_base::trunc) : std::ofstream(filename.c_str(), mode) {
			assertOpen(filename);
		}


		void open(const char* filename, std::ios_base::openmode mode =  std::ios_base::out | std::ios_base::trunc) {
			std::ofstream::open(filename, mode);
			assertOpen(filename);
		}
		void open(const std::string& filename, std::ios_base::openmode mode =  std::ios_base::out | std::ios_base::trunc) {
			std::ofstream::open(filename.c_str(), mode);
			assertOpen(filename);
		}
};
} // namespace PCP_Project 

#endif   // __FileHandling_HPP
