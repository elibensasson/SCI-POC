#ifndef __CONFIGURATION_HPP__
#define __CONFIGURATION_HPP__

#include <stdio.h>
#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <limits.h>
#include <boost/tokenizer.hpp>

const int UNSPECIFIED_TEST_NUM				= INT_MIN;

const std::string PRIMARY_TAPE_DEFAULT		= "DefaultPrimaryInput.txt";
const std::string AUXILIARY_TAPE_DEFAULT	= "DefaultAuxInput.txt";
const unsigned int REGISTER_NUM_DEFAULT		= (16);
const unsigned int REGISTER_SIZE_DEFAULT	= (16);
const bool INSN_DEC_DBG						= false;

class Configuration {
public:
	///////////////////////////////////////////////////////////
	// Methods
	///////////////////////////////////////////////////////////
	static Configuration& getInstance(void);
	void initEnv(void);
	void initSwitches(int argc, char *argv[]);
	void usage(void);
	std::string			getPrimaryTape()			{ return m_primaryTape; }
	std::string			getAuxiliaryTape()			{ return m_auxiliaryTape; }
	unsigned int		getRegistersNum()			{ return m_registersNum; }
	unsigned int		getRegistersSize()			{ return m_registersSize; }
	std::string			getGtests()					{ return m_gtests; }
	bool				hasGtests()					{ return m_hasGtests; }
	unsigned int		getOmpThreadsNum()			{ return m_ompThreadsNum; }
	std::vector<std::string> getRandomArgs()		{ return m_randArgs; }

		///////////////////////////////////////////////////////////
		// Attributes
		///////////////////////////////////////////////////////////
		
	private:
		///////////////////////////////////////////////////////////
		// Methods
		///////////////////////////////////////////////////////////
		Configuration(void);
		Configuration(const Configuration &that) {}
		bool isNumber(const std::string& s);

		///////////////////////////////////////////////////////////
		// Attributes
		///////////////////////////////////////////////////////////
		bool				m_hasGtests;
		std::string			m_gtests;
		std::string			m_primaryTape;
		std::string			m_auxiliaryTape;
		unsigned int		m_registersNum;
		unsigned int		m_registersSize;
		unsigned int		m_ompThreadsNum;
		std::vector<std::string> m_randArgs;

	protected:
		///////////////////////////////////////////////////////////
		// Methods
		///////////////////////////////////////////////////////////

		///////////////////////////////////////////////////////////
		// Attributes
		///////////////////////////////////////////////////////////

};

#endif	// #ifndef __TRANSCRIPT_HPP__
