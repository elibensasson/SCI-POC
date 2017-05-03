#include "TinyRAMDefinitions.hpp"
#include "common/EnvParams.hpp"
#include "common/Configuration.hpp"

namespace PCP_Project {

int trNumRegisters = -1;
int trRegisterLen = -1;
int trNumRegistersInclPC = -1;
int trRegisterOrPCIndexLen = -1;
int trRegisterNoPCIndexLen = -1;
int trRegisterIndexLen = -1;
int trPCLength = -1;
int trInstructionLength = -1;
int trConfigurationLength = -1;
int trFirstReadInstr = -1;
int trInputReadSteps = -1;
int trInputRegisterNum = -1;

static bool trParametersInited = false;

void clearTinyRAMParams(){
    trParametersInited = false;
    trNumRegisters = -1;
    trRegisterLen = -1;
    trNumRegistersInclPC = -1;
    trRegisterOrPCIndexLen = -1;
    trRegisterNoPCIndexLen = -1;
    trRegisterIndexLen = -1;
    trPCLength = -1;
    trInstructionLength = -1;
    trConfigurationLength = -1;
    trFirstReadInstr = -1;
    trInputReadSteps = -1;
    trInputRegisterNum = -1;
    cout<< "TinyRAM parameters are cleared"<<endl;
}

void initTinyRAMParams(int numRegisters, int registerLen) {
    trParametersInited = true;

    _COMMON_ASSERT(numRegisters>=1, "Invalid number of registers");
    _COMMON_ASSERT(registerLen>=1, "Invalid register length");
    cout << "TinyRAM parameters: K=" << numRegisters << " registers, W=" << registerLen << " bits" << endl;

    trNumRegisters = numRegisters;
    trRegisterLen = registerLen;
    trNumRegistersInclPC = trNumRegisters + 1;
    trRegisterOrPCIndexLen = (int)ceil( ::Infrastructure::Log2((double)(trNumRegistersInclPC)) );
    trRegisterNoPCIndexLen = (int)ceil( ::Infrastructure::Log2((double)(trNumRegisters)) );
    trRegisterIndexLen = trRegisterOrPCIndexLen;
    trPCLength = trRegisterLen;
    trInstructionLength = trOpcodeLen + 1 + 2*trRegisterIndexLen + max(trRegisterIndexLen,trRegisterLen);
    trConfigurationLength = trPCLength + trNumRegisters * trRegisterLen + 1;

    //Input Reading:
    trFirstReadInstr = 1;
    trInputReadSteps = 3;
    trInputRegisterNum = 0;
}

void initTinyRAMParamsFromEnvVariables() {
    initTinyRAMParams(getEnvParam("TR_NUM_REGISTERS", DEFAULT_NUM_REGISTERS),
                      getEnvParam("TR_REGISTER_LEN", DEFAULT_REGISTER_LENGTH));
}

/*******************************************************/
/**************** class PreConfiguration ***************/
/*******************************************************/

PreConfiguration::PreConfiguration() {
    PCValue = 0;
    for (int i = 0; i < trNumRegisters; i++)
        registersValues.push_back(0);
    flag = 0;
}

unsigned int PreConfiguration::getRegisterValue(int index) const {
    assert(index>=0);
    assert(index<trNumRegisters);
    return registersValues.at(index);
}

void PreConfiguration::setRegisterValue(unsigned int index, unsigned int newValue) {
    _COMMON_ASSERT(newValue>=0 && newValue< ::Infrastructure::POW2(trRegisterLen),
                 "register value out of bounds");
    registersValues.at(index) = newValue;
}

int PreConfiguration::getPC() {
    return PCValue;
}

void PreConfiguration::setPC(int newPCValue) {
    assert(newPCValue>=0);
    PCValue = newPCValue;
}

bool PreConfiguration::getFlag() const {
    return flag;
}

void PreConfiguration::setFlag(bool newVal) {
    flag = newVal;
}

void PreConfiguration::write(ostream& file) {
    file << PCValue << " ";
    for (int i = 0; i < trNumRegisters; i++)
        file << registersValues.at(i) << " ";
    file << flag << " " << endl;
}


/*******************************************************/
/**************** class Configuration ******************/
/*******************************************************/

/**
* Initialize all registers to zero. PC is initialized to "-1" so in after the first command he will point to line zero.
*/
Configuration::Configuration () {
    assert(trNumRegisters >= 0);
    registerValues.resize(trNumRegisters);
    for (int i=0; i<trNumRegisters;i++)
        for (int j=0;j<trRegisterLen;j++)
            registerValues[i][j] = 0;
    for (int j=0;j<trPCLength;j++)
        PCValue[j] = 0;
    flag = 0;
};

/**
* Copy Constructor.
*/
Configuration::Configuration(const Configuration& other) {
    this->setPCValue(other.getPCValue());

    registerValues.resize(trNumRegisters);
    for (int i=0;i<trNumRegisters;i++)
        this->setRegisterValue(i,other.getRegisterValue(i));
    this->flag = other.getFlag();
}

/**
 * Write a configuration to stream.
 * @param outputStream is the stream to which we write.
 */
void Configuration::print(ostream& outputStream)
{
    for (int i=0; i<trPCLength;i++)
        outputStream << (PCValue[i]);
    outputStream<<" ";

    for (int i=0; i<trNumRegisters;i++){
        for (int j=0; j<trRegisterLen;j++)
            outputStream << registerValues[i][j];
        outputStream<<" ";
    }

    outputStream << flag << endl;
}

/**
* Constructs the configuration from a string.
* The expected format is a binary list in the following order: PC -> R0..Rn -> flag.
*/
void Configuration::constructFromLine (string line){
    for (int i=0; i<trPCLength;i++) this->PCValue[i] = ((atoi(line.substr(i,1).c_str()))==1);
    for (int i=0; i<trNumRegisters;i++){
        for (int j=0; j<trRegisterLen;j++){
            this->registerValues[i][j] = ((atoi(line.substr(trPCLength+1+(trRegisterLen+1)*i+j,1).c_str()))==1);
        }
    }
    this->flag = ((atoi(line.substr(trPCLength+1+(trRegisterLen+1)*trNumRegisters,1).c_str()))==1);
}

/**
* Configuration Equality operator
*/
bool operator==(const Configuration& first, const Configuration& second) {
    if (first.getPCValue() != second.getPCValue()) return false;
    if (first.getFlag() != second.getFlag()) return false;
    if (first.getNumRegisters() != second.getNumRegisters()) return false;
    for (size_t i = 0; i < first.getNumRegisters(); ++i) {
        if (first.getRegisterValue(i) != second.getRegisterValue(i)) {
            return false;
        }
    }
    return true;
}


/*******************************************************/
/**** Yucky functions that produce ugly filenames ******/
/*******************************************************/

//TODO: getAssemblyFiles(): Remove all hard-coded paths. Fix the strange and inconsistent filenames. --Eran

string getAsmPath() {
	::Configuration &conf = ::Configuration::getInstance();
	const char* defaultAsmPath = "../../AssemblyFiles";
    const char* asmPathEnv = "_COMMON_ASM_PATH";

    // Yucky hardcoded path can be orrriden by setting the _COMMON_ASM_PATH environment variable to the path to asm files
	string path = defaultAsmPath;

    assert(path.length()>0);
    path += "/";
    //string path = "C:\\Users\\sohadb\\Documents\\Visual Studio 2010\\Projects\\_COMMON_NEW\\code\\AssemblyFiles\\";
    return path;
}

RAMExecutionFilenames getAssemblyFile(string fileName) {
    RAMExecutionFilenames toReturn;
	::Configuration &conf = ::Configuration::getInstance();
	const string path = getAsmPath();
    toReturn.codeFile = path + fileName;
	toReturn.firstInput = path + conf.getPrimaryTape();
	toReturn.secondInput = path + conf.getAuxiliaryTape();
    return toReturn;
}

/**
 * Create input instances conveniently for the RAMParams constructor.
 * Assumes a hard-coded directory where the file resides.
 * Use this version if you only wish to create GCP instance with some t, with a default name.
 */
RAMExecutionFilenames getAssemblyFile(const int t) {
    stringstream ss;
    ss << "assemblyT" << t << ".S";
    return getAssemblyFile(ss.str());
}



} // of namepsace
