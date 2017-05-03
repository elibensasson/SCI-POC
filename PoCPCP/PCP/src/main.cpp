/**
 * Executable's main file, containing main() which performs some initialization and invokes the desired tests.
 */
#include "PCP/SoundnessParams.hpp"
#include "languages/TinyRAM/TinyRAMDefinitions.hpp"
#include "boost/lexical_cast.hpp"
#include "common/Configuration.hpp"
#include "gtest/gtest.h"
#include "common/Utils/TaskReporting.hpp"
#include "common/EnvParams.hpp"
#include <gtest/gtest.h>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <FFT.h>
#include <omp.h>

using namespace std;
using namespace NTL;
using namespace PCP_Project;

#define MAIN_IRR_DEGREE 64

// switch to disable transcript output for
extern "C" {
	bool swPrintTranscript = true;
	bool swSimTypeHv = true;
};

int main(int argc, char *argv[]) {

    //int tmp = _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG);
    //tmp |= _CRTDBG_DELAY_FREE_MEM_DF | _CRTDBG_CHECK_ALWAYS_DF | _CRTDBG_CHECK_CRT_DF | _CRTDBG_LEAK_CHECK_DF;
    //tmp = (tmp & 0x0000FFFF) | _CRTDBG_CHECK_EVERY_16_DF;
    //_CrtSetDbgFlag(tmp);
    //_CrtSetReportMode( _CRT_ERROR, _CRTDBG_MODE_FILE );
    //_CrtSetReportFile( _CRT_ERROR, _CRTDBG_FILE_STDERR );

	//seed randomness
	srand(time(0)); //TODO: Consider changing to a better seed for better randomness
	
	//execute program

	TASK("Main program");
	
	::Algebra::details::initGF2E(MAIN_IRR_DEGREE);
	SoundnessParameters::setupSoundnessParams();

	
	///Verifying the correctness of the soundness parameters (Mu, Eta, k0 etc.)

	if (!areSoundnessParametersValid()){
		cout << "invalid soundness params";
		_COMMON_FATAL("Soundness parameters are incorrect");
	}

	::Configuration &conf = ::Configuration::getInstance();
	conf.initEnv();
	conf.initSwitches(argc, argv);

	FFF::set_omp_state(conf.getOmpThreadsNum());
	int retval;
	if (conf.hasGtests()) {
		int i;
		std::string gtests = conf.getGtests();
		std::vector<std::string> tokens;
		boost::split(tokens, gtests, boost::is_any_of(" \t"));
		int argc_ = argc + tokens.size();
		char **argv_ = (char**)alloca(argc_*sizeof(char*));
		for (i = 0; i < tokens.size(); i++) { argv_[i] = &tokens.at(i)[0]; }
		for (; i < argc; i++) { argv_[i] = argv[i]; }
		argv_[i] = '\0';
		::testing::InitGoogleTest(&argc_, argv_);
		retval = RUN_ALL_TESTS();
	}

	//release NTL global info
	delete ::NTL::GF2EInfo;
	return retval;
}
