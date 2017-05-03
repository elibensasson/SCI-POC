#ifndef TINY_RAM_DEFINITIONS_HPP
#define TINY_RAM_DEFINITIONS_HPP
/************************************************************************TinyRAMDefinitions.hpp****************************************/
/** @file
 * Basic definitions of parameters and constants for TinyRAM, together with data types needed in the reductions from TinyRAM.
 * It also contains functions for creating input instances for the RAMParams constructor.
 **/
//TODO: Keep only basic TinyRAM defs here, and move reduction-related stuff (deBruin etc.) to a separate file.


#include <cstdint>
#include <iostream>
#include <bitset>
#include <fstream>
#include <sstream>
#include <boost/filesystem.hpp>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include "common/Algebra/details/AlgebraUtils.hpp"
#include "common/Infrastructure/Infrastructure.hpp"

#define DEFAULT_NUM_REGISTERS 15
#define DEFAULT_REGISTER_LENGTH 16
#define BIN_EOF_STRING "EOF"

using namespace std;

namespace PCP_Project {

// TinyRAM flexible parameters (set via parameters to initTinyRAMParams)

/** Number of registers that the TinyRAM has (excluding PC) */
extern int trNumRegisters;

/** The bit-length of TinyRAM registers */
extern int trRegisterLen;

// TinyRAM derived/inflexible parameters (mostly set inside initTinyRAMParams))

/** Numbers of registers including the program counter (used for code clarity)
 * Set to trNumRegisters + 1. */
extern int trNumRegistersInclPC;

/** Length of an index into the register file, including PC.
 * Set to ceil( ::Infrastructure::Log2(trNumRegistersInclPC)).
 */
extern int trRegisterOrPCIndexLen;

/** Length of an index into the register file, excluding PC.
 * Set to ceil( ::Infrastructure::Log2(trNumRegistersPC) ).
 */
extern int trRegisterNoPCIndexLen;

/** Length of an index into the register file, including PC.
 * @note Our legacy TinyRAM binary format uses indices that allow PC, and in particular uses this constant unnecessarily.
 * @deprecated Unless interfacing with the legacy TinyRAMcode, avoid this constant and use the explicit trRegisterOrPCIndexLen or trRegisterNoPCIndexLen.
 */
//TODO: The TinyRAM instruction format wastes a register index on the PC (contrary to spec!). Should be fixed. Any occurance of trRegisterInclPCIndexLen in that code indicates a problem spot. --Eran
extern int trRegisterIndexLen;


/** Length of PC register.
 * Currently same as the general-purpose registers so used just for code readibility.
 * (There are probably some places that say trRegisterLen instead of trPCLength, and breaking the equivalence would require changes in the K-Circuit structure) */
extern int trPCLength;

/** The bit size of an op's binary representation.
 * Set to: trOpcodeLen + 1 + 2*trRegisterIndexLen + max(trRegisterIndexLen,trRegisterLen)
 * */
extern int trInstructionLength;

/** The bit size of a TinyRAM configuration in binary representation.
 * Set to: trPCLength + trNumRegisters * trRegisterLen + 1
 */
extern int trConfigurationLength;

/** 
 * index of the 1st reada instruction. 
 * We assume the first instruction is "store r0 r0", and immediately after we start reading the input. 
 * Thus the value we use for this parameter is 1 (meaning the 2nd instruction)
 */ 
extern int trFirstReadInstr; 

/** 
 * Each TinyRAM program starts with "store r0 r0", followed by "reada R0" instructions that appear every k'th step, until the 
 * input is fully read. The parameter trInputReadSteps defines k.
 */
extern int trInputReadSteps;

/**
 * The register number to which we read the input.
 * Currently set to R0.
 */
extern int trInputRegisterNum;

/**
 * The enum TinyRAMOpcode_depracated gives the binary opcode values for TinyRAM ops, as produced by the
 * simulator and verified by the circuit.
 */
//TODO: This enum is should be depracated, still here for Ohad's legacy code for testing purposes.
//      Use TinyRAM_opcode enum class instead
enum TinyRAMOpcode_depracated {
    AdditionOpcode		= 0,
    NopOpcode			= 1,
    ploadOpcode			= 2,
    storeOpcode			= 4,
    loadOpcode			= 5,
    cJmpOpcode			= 6,
    JmpOpcode			= 7,
    AndOpcode			= 8,
    notOpcode			= 9,
    xorOpcode			= 10,
    orOpcode			= 11,
    mulOpcode			= 12,
    shiftLeftOpcode 	= 13,
    shiftRightOpcode	= 14,
    cmpeOpcode			= 15,
    cmpgOpcode			= 16,
    cmpgeOpcode			= 17,
    readAOpcode			= 19,
    readBOpcode			= 20,
    movOpcode			= 21,
    cmovOpcode			= 22,
    subOpcode			= 23,
    doesNothingOpcode	= 24,
    cmpsgOpcode			= 25,
    cmpsgeOpcode		= 26,
    cnJmpOpcode			= 27,
    acceptOpcode		= 31,
    umullOpcode			= 32,
    umulhOpcode			= 33,
    smullOpcode			= 34,
    smulhOpcode			= 35,
    udivOpcode			= 36,
    umodOpcode			= 37,
    storeByteOpcode		= 38,
    loadByteOpcode		= 39,
    storeWordOpcode		= 40,
    loadWordOpcode		= 41,
};

/**
    The TinyRAMOpcode::Opcode enum class gives the binary opcode encodings.
    Last update: 16/01/13 TinyRAM Architecture Specification v1.995
    It is namespced in order to get scoping and readable syntax. Strongly typed enums where not
    used in order to allow implicit conversion to int.
*/
namespace TinyRAMOpcode {
enum Opcode {
    AND     = 0,
    OR      = 1,
    XOR     = 2,
    NOT     = 3,
    ADD     = 4,
    SUB     = 5,
    MULL    = 6,
    UMULH   = 7,
    SMULH   = 8,
    UDIV    = 9,
    UMOD    = 10,
    SHL     = 11,
	SHR		= 12,
	SHAR	= 13,

    CMPE    = 14,
    CMPA    = 15,
    CMPAE   = 16,
    CMPG    = 17,
    CMPGE   = 18,

    MOV     = 19,
    CMOV    = 20,
    JMP     = 21,

    CJMP    = 22,
    CNJMP   = 23,

    RESERVED_OPCODE_24 = 24,
    RESERVED_OPCODE_25 = 25,

    STOREB  = 26,
    LOADB   = 27,
    STOREW  = 28,
    LOADW   = 29,
    READ    = 30,
    ANSWER  = 31,
    NUM_OPCODES = 32
}; // enum Opcode

const Opcode controlFlowInstructions[] = {
    JMP,
    CJMP,
    CNJMP
};

const Opcode stallInstructions[] = {
    RESERVED_OPCODE_24,
    RESERVED_OPCODE_25,
    ANSWER
};

const Opcode registerInstructions[] = {
    AND,
    OR,
    XOR,
    NOT,
    ADD,
    SUB,
    MULL,
    UMULH,
    SMULH,
    UDIV,
    UMOD,
    SHL,
	SHR,
	SHAR,

    CMPE,
    CMPA,
    CMPAE,
    CMPG,
    CMPGE,

    MOV,
    CMOV,

    LOADB,
    LOADW,
    READ
};



} // namespace TinyRAMOpcode


/// Holds the number of possible opcodes
const int trNumOpcodes = TinyRAMOpcode::NUM_OPCODES;

/** 
 * Holds the length of an opcode representation. Should be Log2ceil(trNumOpcodes). Must currently
 * be const and not a function call for legacy code compatibility (trOpcodeLength is used as a 
 * bitset<> template argument). When MSVC starts supporting 'constexpr' second format can be used.
 */
const int trOpcodeLen = 5;
// constexpr const int trOpcodeLen = ::Infrastructure::Log2ceil(trOpcodeLen);

/**
 * Holds the filenames needed by the TinyRAM simulator (the code file and two input files).
 */
struct RAMExecutionFilenames {
    string codeFile;
    string firstInput;
    string secondInput;
};


/************ Bit vectors for various TinyRAM values **************/

// Dynamic allocation. Overhead for keeping track of length is empirically insigificant in the grand scheme of things.
#define DEFINE_TINYRAM_BIT_VECTOR(_type, _len) \
        class _type : public boost::dynamic_bitset<> { \
        public: \
        _type() : boost::dynamic_bitset<>(_len) {} \
        };

DEFINE_TINYRAM_BIT_VECTOR(InstructionBits, trInstructionLength);
DEFINE_TINYRAM_BIT_VECTOR(RegisterBits, trRegisterLen);
DEFINE_TINYRAM_BIT_VECTOR(MemoryAddressBits, trRegisterLen);
DEFINE_TINYRAM_BIT_VECTOR(CodeAddressBits, trPCLength);

/**
 * Convert register value to memory address.
 * Explicit, in case we change TinyRAM memory addressing (e.g., from word-addressible to byte-addressible).
 * @param regAddr - address word from register (or immediate)
 * @return Address bits as index into memory words
 */
inline MemoryAddressBits convert_RegisterBits_to_MemoryAddressBits(const RegisterBits & regAddr) {
    MemoryAddressBits addr;
    for (int i=0; i<trRegisterLen; ++i)
        addr[i] = regAddr[i];
    return addr;
}

/**
 * Convert register value to code address.
 * Explicit, in case we change TinyRAM code addressing.
 * @param regAddr - address word from register (or immediate)
 * @return Address bits as index into instruction list
 */
inline CodeAddressBits convert_RegisterBits_to_CodeAddressBits(const RegisterBits & regAddr) {
    CodeAddressBits addr;
    for (int i=0; i<<trRegisterLen; ++i)
        addr[i] = regAddr[i];
    return addr;
}


/**
 * Represent a TinyRAM configuration, with values stored as integers.
 * Includes the PC, registers and flag.
 * Used in the TinyRAM interpreter.
 */
class PreConfiguration{
private:
    int PCValue;
    vector<int> registersValues;
    bool flag;
public:
	PreConfiguration();
	/**
	 *Write a pre-configuration to file.
	 *Used by the transcript producer, that calculate all the pre-configurations and outputs them to a file.
	 **/
	void write(ostream& file);
	int getPC();
	void setPC(int newPCValue);
	unsigned int getRegisterValue(int index) const;
	void setRegisterValue(unsigned int index, unsigned int newValue);
	bool getFlag() const;
	void setFlag(bool newVal);
	;
};

/**
 * Represent a TinyRAM configuration, with values stored as bit vectors.
 * Includes the PC, registers and flag.
 */
class Configuration{
protected:
    //an array of registers and pc.
    CodeAddressBits PCValue;
    vector<RegisterBits> registerValues;
    bool flag;
public:

    /**
     * Initialize all registers to zero. PC is initialized to "-1" so in after the first command he will point to line zero.
     */
    Configuration();

    /**
     * Copy Constructor.
     */
    Configuration(const Configuration& other);

    //getters and setters in the bit level and in the register level.
    bool getFlag() const {return flag;}
    void setFlag(bool newVal) {flag = newVal;}
    bool getIthBitOfRegisterJ(int i,int j) const {return registerValues[j][i];}
    void setIthBitOfRegisterJ(int i, int j,bool value) {registerValues[j][i] = value;}
    size_t getNumRegisters() const {return registerValues.size();}
    RegisterBits getRegisterValue(int registerNumber) const {return registerValues[registerNumber];}
    void setRegisterValue(int registerNumber, const RegisterBits & value) {registerValues[registerNumber] = value;}
    void setPCValue(CodeAddressBits value) {PCValue = value;}
    void setIthBitOfPC(int i,bool value) {PCValue[i] = value;}
    bool getIthBitOfPC(int i) const {return PCValue[i];}
    CodeAddressBits getPCValue() const {return PCValue;}
    
    /**
     * Write a configuration to stream.
     * Used by the toBinary module, that calculate all the configurations and output them to a file.
     **/
    void print(ostream& outputStream);
    void print() { print(cout); }		//Print configration to stdout.

    /**
     * Constructs the configuration from a string.
     * The expected format is a binary list in the following order: PC -> R0..Rn -> flag.
     */
    void constructFromLine (string line);

};

bool operator==(const Configuration& first, const Configuration& second);

/**
 * The color of each vertex at the de Bruijn graph. Each color consists of a configuration and an index.
 */
struct NodeColor{
    int Tau; // Timestamp
    Configuration configuration;
    void print(ofstream & outputFile){
        outputFile<<"Tau is: ";
        outputFile<<Tau<< " .\n";
        outputFile<<" Configuration is: .\n ";
        configuration.print(outputFile);
    }
    void print (){
        cout<<"Tau is: ";
        cout<<Tau<< " .\n";
        cout<<" Configuration is: .\n ";
        configuration.print();
    }
};

// A different name to differentiate between RAM context and Graph context.typedef NodeColor ComputationStep;
typedef NodeColor ComputationStep;
// Representation of a computation transcript as a vector of computation steps
typedef vector<ComputationStep*> BinaryRamTranscript;
// Binary representation of an assembler instruction
  // TODO (shaul,ohad) Code Spaghetti -> this is a parent class for timyRAM but is using TinyRAM's trInstructionLength.
// Binary representation of a vector of instructions
typedef vector<InstructionBits> BinaryRamCode;


/**
 * Abstract class for object that access the memory.
 * Provides a comparion operator< between memory related objects, based on their addresses (to be provided by concrete derived classes).
 * The address is given as a bit vector, most-significant bit first.
 * Used for creating the memory permutation.
 */
class MemoryRelatedObject{
protected:
    MemoryAddressBits addressAccessed;
public:
    MemoryRelatedObject(MemoryAddressBits ourAddress) {addressAccessed = ourAddress;};
    virtual MemoryAddressBits getAddress () {return addressAccessed;};
    bool operator< (const MemoryRelatedObject & other) const {
        for (int i=0; i<trRegisterLen; i++){
            if (this->addressAccessed[i] < other.addressAccessed[i]) return true;
            if (this->addressAccessed[i] > other.addressAccessed[i]) return false;
        }
        return false;
    };
};


/**
 * A TinyRAM configuration that "access the memory", i.e. executes a STORE or PLOAD op.
 * Inherits from both Configuration and MemoryRelatedObject.
 * Used for creating the memory permutation.
 */
class MemoryConfiguration : public Configuration, public MemoryRelatedObject{
private:
    int timeStamp;

public:
    /** Class Constructor */
    MemoryConfiguration (const Configuration & x, int time , const MemoryAddressBits & address) : 
      Configuration(x), MemoryRelatedObject(address), timeStamp(time) {}

    /** Timestamp getter */
    const int GetTime () const {return timeStamp;};
};


/**
 * Set the various TinyRAM parameters.
 * Must be called before using any TinyRAM-related variables or classes.
 * Must not be called more than once.
 */
void initTinyRAMParams(int numRegisters, int registerLen);
void initTinyRAMParamsFromEnvVariables();

/**
 * Clears the TinyRAM params
 * sets the state as before it was initialized
 */
void clearTinyRAMParams();



/************************ Filenames *************************/

// Filename extensions
const string asmFilenameExt = ".S";
const string canonicalAsmFilenameExt = ".code.S";          //TODO: Change to ".canS"        --Eran
const string binaryCodeExt = ".code.bin";            //TODO: Change to ".code.bin"    --Eran

const string textualTranscriptExt = ".conf-trans.txt";            //TODO: Change to ".trans.txt"   --Eran
const string binaryTranscriptExt = ".conf-trans.bin";  //TODO: Change to ".trans.bin"   --Eran
const string textualIntructionTranscriptExt = ".inst-trans.txt";
const string binaryIntructionTranscriptExt = ".inst-trans.bin";
const string memoryPermutationExt = ".mem-perm.txt";	
const string inputTimestampsExt = ".input-trans.txt";	//Contains a list of timestamps in which reada was done.

//TODO: "first input" and "second input" should be added here and consistently renamed to "input" and "witness". --Eran


/**
 * Ugly hardcoded function which gets path to assembly files
 */
::std::string getAsmPath();
     
/**
 * Create input instances conveniently for the RAMParams constructor.
 * Assumes a hard-coded directory where the file resides.
 * Use this version if you have a specific filename to run.
 */
RAMExecutionFilenames getAssemblyFile(string fileName);

/**
 * Create input instances conveniently for the RAMParams constructor.
 * Assumes a hard-coded directory where the file resides.
 * Use this version if you only wish to create GCP instance with some t, with a default name.
 */
RAMExecutionFilenames getAssemblyFile(const int t);



typedef unsigned int ordinal_t; // ??? --Eran


} // namespace



#endif

