#include <boost/program_options.hpp>

#include "Configuration.hpp"
#include "common/Utils/commonUtils.hpp"
#include <omp.h>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

Configuration::Configuration(void)
{
	m_primaryTape				= PRIMARY_TAPE_DEFAULT;
	m_auxiliaryTape				= AUXILIARY_TAPE_DEFAULT;
	m_registersNum				= REGISTER_NUM_DEFAULT;
	m_registersSize				= REGISTER_SIZE_DEFAULT;
	m_ompThreadsNum				= omp_get_max_threads();
	m_ompThreadsNum				= m_ompThreadsNum%2 ? (m_ompThreadsNum+1)/2 : m_ompThreadsNum/2;
}

Configuration& Configuration::getInstance(void)
{
	static Configuration instance;
	return instance;
}

void Configuration::initEnv(void)
{
	m_primaryTape			= PCP_Project::getEnvParam("TR_INPUT1",			PRIMARY_TAPE_DEFAULT);
	m_auxiliaryTape			= PCP_Project::getEnvParam("TR_INPUT2",			AUXILIARY_TAPE_DEFAULT);
	m_registersNum			= PCP_Project::getEnvParam("TR_NUM_REGISTERS",	REGISTER_NUM_DEFAULT);
	m_registersSize			= PCP_Project::getEnvParam("TR_REGISTERS_SIZE",	REGISTER_SIZE_DEFAULT);
}

void Configuration::initSwitches(int argc, char *argv[])
{
	std::stringstream msg;
	boost::program_options::variables_map vm;
	boost::program_options::options_description description(					
		"Usage: _COMMON "
		"--tr-input1 <file name> --tr-input2 <file name> "
		"[--gtest <arg>] "
		"[--tr-num-registers <num>] "
		"[--tr-registers-size <num>] "
		"[--omp-max-threads <num>] "
		"[--add-extra-debug-instructions] \n"
		);
	description.add_options()
		("help,h",				"Display this help message")
		("version,v",			"Display simulator version")
		("gtest",
			boost::program_options::value<std::string>(),
			"Run gtests tests (to run all tests use --gtest \"\")")
		("tr-input1",
			boost::program_options::value<std::string>(),
			"Set the name of the primary tape file (override the TR_INPUT1 environment variable]")
		("tr-input2",
			boost::program_options::value<std::string>(), 
			"Set the name of the auxilury tape file (override the TR_INPUT2 environment variable]")
		("tr-num-registers",
			boost::program_options::value<unsigned int>(),
			"Set the number of registers to use (override the TR_NUM_REGISTERS environment variable]")
		("tr-registers-size",
			boost::program_options::value<unsigned int>(),
			"Set the size of registers in bits (override the TR_REGISTERS_SIZE environment variable]")
		("omp-max-threads",
			boost::program_options::value<unsigned int>(),
			"Run FFT lib with <omp-max-threads> threads")
		("add-extra-debug-instructions",
			"Add debug print instructions support")
		("extra-args",
			boost::program_options::value<std::string>(),
			"Pass extra arguments to gtests")
		;
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(description).run(), vm);
		boost::program_options::notify(vm);
	} catch (std::exception &e) {
		msg << std::endl << e.what() << std::endl;
		msg << description << std::ends;
		PCP_Project::ErrorHandling::fatalError(msg);
	}
	
	if (vm.count("help")){
		msg << description << std::ends;
		PCP_Project::ErrorHandling::fatalError(msg);
	}
	m_hasGtests = vm.count("gtest");
	if (vm.count("gtest")) { m_gtests = vm["gtest"].as<std::string>(); }
	if (vm.count("tr-input1")) { m_primaryTape = vm["tr-input1"].as<std::string>(); }
	if (vm.count("tr-input2")) { m_auxiliaryTape = vm["tr-input2"].as<std::string>(); }
	if (vm.count("tr-num-registers")) { m_registersNum = vm["tr-num-registers"].as<unsigned int>(); }
	if (vm.count("tr-registers-size")) { m_registersSize = vm["tr-registers-size"].as<unsigned int>(); }
	if (vm.count("omp-max-threads")) { m_ompThreadsNum = vm["omp-max-threads"].as<unsigned int>(); }
	if (vm.count("extra-args")) {
		std::string gtestRand = vm["extra-args"].as<std::string>();
		boost::split(m_randArgs, gtestRand, boost::is_any_of(" \t"));
	}
}
