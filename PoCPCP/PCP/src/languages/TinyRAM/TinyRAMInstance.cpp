#include "TinyRAMInstance.hpp"
#include "gadgetlib/common_use.hpp"
namespace PCP_Project{
	MachineInstruction::MachineInstruction(const Opcode& opcode, const bool arg2isImmediate,
		const size_t destIdx, const size_t arg1Idx, const size_t arg2IdxOrImmediate) :
		opcode_(opcode), arg2isImmediate_(arg2isImmediate), destIdx_(destIdx), arg1Idx_(arg1Idx),
		arg2IdxOrImmediate_(arg2IdxOrImmediate){}
}