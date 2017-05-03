#include "common/Infrastructure/Infrastructure.hpp"

#include <boost/lexical_cast.hpp>
#include <string>

#ifndef  __EnvParams_HPP
#define  __EnvParams_HPP

namespace PCP_Project {
/**
 * Read a parameter from an environment variable. If the environment variable isn't given, use defaultVal.
 * Aborts (fatal error) if the environment variable exists but cannot be parsed.
 * The supported parameter types are those supported by boost::lexical_cast. In particular, this includes strings and ints.
 */
template<class T>
	T getEnvParam(std::string envName, const T& defaultVal) {
		const char* env = getenv(envName.c_str());
		if (env) {
			try {
				return boost::lexical_cast<T>(env);
			} catch ( boost::bad_lexical_cast ) {
				_COMMON_FATAL("Environment variable cannot be parsed: "+envName+"="+env);
			}
		} else {
			return defaultVal;
		}
	}

/**
 * Read a parameter from an environment variable.
 * The supported parameter types are those supported by boost::lexical_cast. In particular, this includes std::strings and ints.
 * Aborts (fatal error) if the environment does not exist, or exists but cannot be parsed.
 */template<class T>
T getMandatoryEnvParam(std::string envName) {
	const char* env = getenv(envName.c_str());
	if (env) {
		try {
			return boost::lexical_cast<T>(env);
		} catch ( boost::bad_lexical_cast ) {
			_COMMON_FATAL("Environment variable cannot be parsed: "+envName+"="+env);
		}
	} else {
		_COMMON_FATAL("Mandatory environment variable not set: "+envName);
	}
}

} // namespace PCP_Project 

#endif   // __EnvParams_HPP
