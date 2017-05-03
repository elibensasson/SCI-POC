/************************************** TinyRAMInstance.hpp **************************************/
/**
 * @file
 *
 * @brief The file TinyRAMInstance.hpp contains the interface for a TinyRAM instance.
 * 
 * Main classes defined in the file:
 *     TinyRAMPartialInstance
 *     TinyRAMFullInstance
 */
  /***********************************************************************************************/

#ifndef __TINYRAM_INSTANCE_HPP
#define __TINYRAM_INSTANCE_HPP

#include <vector>
#include <string>
#include <cstdint>
#include "languages/TinyRAM/TinyRAMDefinitions.hpp"
#include <NTL/GF2E.h>
#include <algorithm>
#include <math.h>
#include "../../common/Infrastructure/Infrastructure.hpp"
#include <algebraLib/FieldElement.hpp>
#include <algebraLib/CircuitPolynomial.hpp>
#include <algebraLib/variable.hpp>
#include <gadgetlib/infrastructure.hpp>
#include <gadgetlib/common_use.hpp>
#include <gadgetlib/protoboard.hpp>
using gadgetlib::Opcode;
using gadgetlib::Log2ceil;

namespace PCP_Project {
	//enum class Opcode : int {
	//	MEMORY = -2,
	//	NONE = -1,
	//	AND = 0,
	//	OR = 1,
	//	XOR = 2,
	//	NOT = 3,
	//	ADD = 4,
	//	SUB = 5,
	//	MULL = 6,
	//	UMULH = 7,
	//	SMULH = 8,
	//	UDIV = 9,
	//	UMOD = 10,
	//	SHL = 11,
	//	SHR = 12,
	//	SHAR = 13,

	//	CMPE = 14,
	//	CMPA = 15,
	//	CMPAE = 16,
	//	CMPG = 17,
	//	CMPGE = 18,

	//	MOV = 19,
	//	CMOV = 20,
	//	JMP = 21,

	//	CJMP = 22,
	//	CNJMP = 23,

	//	RESERVED_OPCODE_24 = 24,
	//	RESERVED_OPCODE_25 = 25,

	//	STOREB = 26,
	//	LOADB = 27,
	//	STOREW = 28,
	//	LOADW = 29,
	//	READ = 30,
	//	ANSWER = 31,
	//	NUM_OPCODES = 32
	//}; // enum Opcode



	/*************************************************************************************************/
	/*************************************************************************************************/
	/****************************                                         ****************************/
	/****************************         class TinyRAMArchParams         ****************************/
	/****************************                                         ****************************/
	/*************************************************************************************************/
	/*************************************************************************************************/

	struct TinyRAMArchParams {
		size_t numRegisters;
		size_t registerLength;

		bool operator==(const TinyRAMArchParams& rhs) const;
	};

	/*************************************************************************************************/
	/*************************************************************************************************/
	/****************************                                         ****************************/
	/****************************         class TinyRAMProgram            ****************************/
	/****************************                                         ****************************/
	/*************************************************************************************************/
	/*************************************************************************************************/
	/// A data object which holds a TinyRAM code and auxilary information.
	struct MachineInstruction {
		Opcode opcode_ = Opcode::ANSWER;
		bool arg2isImmediate_ = true;
		size_t destIdx_ = 0;
		size_t arg1Idx_ = 0;
		size_t arg2IdxOrImmediate_ = 1;

		MachineInstruction(
			const Opcode& opcode,
			const bool arg2isImmediate,
			const size_t destIdx,
			const size_t arg1Idx,
			const size_t arg2IdxOrImmediate);
	};

	class TinyRAMProgram {
	public:
		typedef ::std::vector<MachineInstruction> TinyRAMMachineCode;
	private:
		::std::string name_;
		TinyRAMArchParams archParams_;
		TinyRAMMachineCode code_;
	public:
		TinyRAMProgram(const ::std::string& name,
			const TinyRAMArchParams& archParams,
			const TinyRAMMachineCode& code) :
			name_(name), archParams_(archParams), code_(code) {}

		TinyRAMProgram(const ::std::string& name,
			size_t numRegisters,
			size_t wordSize) :
			name_(name), archParams_(TinyRAMArchParams{ numRegisters, wordSize }) {
		}


		::std::string name() const { return name_; }
		const TinyRAMMachineCode& code() const { return code_; }
		const size_t size() const { return code_.size(); }
		const TinyRAMArchParams& archParams() const { return archParams_; }
		const MachineInstruction& getInstructionAtPc(const size_t pc) const { return code_[pc]; }
		void addInstruction(const MachineInstruction& instruction) {
			code_.emplace_back(instruction);
		}
		void loadProgram();
		unsigned int pcLength() const {
			int codeSize = code_.size();
			if (codeSize == 0){ _COMMON_FATAL("TinyRAMProgram : The code is not initialized"); };
			if (codeSize == 1) { return 1; }
			return  gadgetlib::Log2ceiled(codeSize);
		}


	};

	::std::ostream& operator<<(::std::ostream& os, const TinyRAMProgram& prog);

	


	///*************************************************************************************************/
	///*************************************************************************************************/
	///*******************                                                            ******************/
	///*******************					Memory Info				 				******************/
	///*******************                                                            ******************/
	///*************************************************************************************************/
	///*************************************************************************************************/

	//class MemoryInfo{
	//private:
	//	int serialNumber_;
	//	bool isMemOp_;
	//	bool isLoadOp_;
	//	Algebra::FElem timestamp_;
	//	int timestampDegree_;
	//	Algebra::FElem address_;
	//	Algebra::FElem value_;
	//public:
	//	MemoryInfo();
	//	MemoryInfo(int serialNumber, bool isMemOp, bool isLoadOp,
	//		const Algebra::FElem& timestamp, int timestampDegree,
	//		const Algebra::FElem& address, const Algebra::FElem& value);
	//	int getSerialNumber() const{ return serialNumber_; }
	//	bool getIsMemOp() const{ return isMemOp_; }
	//	bool getIsLoadOp() const { return isLoadOp_; }
	//	Algebra::FElem getTimestamp() const{ return timestamp_; }
	//	int getTimestampDegree() const{ return timestampDegree_; }
	//	Algebra::FElem getAddress() const { return address_; }
	//	Algebra::FElem getValue() const{ return value_; }
	//	void updateIsMemOp(bool isMemOp){ isMemOp_ = isMemOp; }
	//	void updateIsLoadOp(bool isLoadOp){ isLoadOp_ = isLoadOp; }
	//	void updateAddress(const Algebra::FElem& address){ address_ = address; }
	//	void updateSerialNumber(int serialNumber){ serialNumber_ = serialNumber; }
	//	void updateValue(const Algebra::FElem& value){ value_ = value; }
	//	void updateTimestamp(Algebra::FElem timestamp, int timestampDegree);
	//};

	//inline bool operator==(const MemoryInfo& a, const MemoryInfo& b){
	//	return (a.getSerialNumber() == b.getSerialNumber() && a.getIsMemOp() == b.getIsMemOp() && a.getIsLoadOp() == b.getIsLoadOp() &&
	//		a.getTimestamp() == b.getTimestamp() && a.getAddress() == b.getAddress() && a.getValue() == b.getValue());
	//}

	//bool sortMemoryInfo(MemoryInfo a, MemoryInfo b);



	/*************************************************************************************************/
	/*************************************************************************************************/
	/****************************                                         ****************************/
	/****************************         class TinyRAMProtoboardParams   ****************************/
	/****************************                                         ****************************/
	/*************************************************************************************************/
	/*************************************************************************************************/

	class TinyRAMProtoboardParams : public gadgetlib::ProtoboardParams {
	private:
		TinyRAMArchParams archParams_;
		size_t opcodeWidth_;
		size_t timeBound_;
		size_t pcIncrement_;
	public:
		TinyRAMProtoboardParams(unsigned int numRegisters, unsigned int registerLength,
			size_t opcodeWidth, size_t timeBound, size_t pcIncrement)
			: archParams_(TinyRAMArchParams{ numRegisters, registerLength }),
			opcodeWidth_(opcodeWidth),
			timeBound_(timeBound), pcIncrement_(pcIncrement) {}
		TinyRAMProtoboardParams()
			: archParams_(TinyRAMArchParams{ 0, 0 }), opcodeWidth_(0), timeBound_(0), pcIncrement_(0) {}
		TinyRAMArchParams archParams() const { return archParams_; }
		size_t opcodeWidth() const { return opcodeWidth_; }
		size_t numRegisters() const { return archParams_.numRegisters; }
		size_t registerLength() const { return archParams_.registerLength; }
		size_t registerIndexLength() const { return Log2ceil(numRegisters()); }
		size_t arg2length() const { return std::max({ registerIndexLength(), registerLength() }); }
		size_t numOpcodes() const { return 1u << (opcodeWidth()); }
		size_t timeBound() const { return timeBound_; }
		size_t pcIncrement() const { return pcIncrement_; }
	}; // class TinyRAMProtoboardParams



} // namespace PCP_Project
#endif // __TINYRAM_INSTANCE_HPP
