#include <gadgetlib/common_use.hpp>
#include <gadgetlib/infrastructure.hpp>
#include <algebraLib/variable_operators.hpp>
#include <NTL/GF2XFactoring.h>
#include <NTL/GF2.h>
#include <stdint.h>

namespace gadgetlib{

	//useful for constructing SelectorSum objects
	std::vector<Algebra::Variable> getPCVars(const Algebra::UnpackedWord& pc){
		std::vector<Algebra::Variable> retVal;
		for (auto i = 0; i < pc.size(); i++)
			retVal.push_back(pc[i]);
		return retVal;
	}


	//MachineInstruction::MachineInstruction(const Opcode& opcode, const bool arg2isImmediate,
	//	const size_t destIdx, const size_t arg1Idx, const size_t arg2IdxOrImmediate) :
	//	opcode_(opcode), arg2isImmediate_(arg2isImmediate), destIdx_(destIdx), arg1Idx_(arg1Idx),
	//	arg2IdxOrImmediate_(arg2IdxOrImmediate){}

::NTL::GF2X getPrimitivePoly(){
	::NTL::GF2X res;
	uint64_t polyRepresentation;
	//We assume the polynomial x^64 + x^4 + x^3 + x + 1
	for (int digitIndex = 0; digitIndex <= 64; digitIndex++) {
		if (digitIndex == 0 || digitIndex == 1 || digitIndex == 3 || digitIndex == 4 || digitIndex == 64){
			SetCoeff(res, digitIndex, 1);
			continue;
		}
		SetCoeff(res, digitIndex, 0);
	}
	//return res;
	::NTL::GF2X Irr;
	::NTL::BuildIrred(Irr, 64);
	return Irr;
}

void initGF2EwithPrimitivePoly() {
	const NTL::GF2X primitive = getPrimitivePoly();
	NTL::GF2E::init(primitive);
}
NTL::GF2E getGF2E_X() {
	return NTL::to_GF2E(NTL::GF2X(1, 1));
}

Algebra::FieldElement getBit(Algebra::FieldElement elem, int i){
	std::cout << elem << std::endl;
	const ::NTL::GF2X& elem2 = rep(NTL::GF2E(elem));
	const ::NTL::GF2& coeff = ::NTL::coeff(elem2, i);
	std::cout << coeff << std::endl;
	return ::NTL::to_GF2E(coeff) == NTL::GF2E(Algebra::zero()) ? Algebra::zero() : Algebra::one();
}

/*::NTL::GF2E& getGenerator2toTheNthPower(const unsigned int n){
	static std::vector<NTL::GF2E> g_2_i = { getGF2E_X() };
	const unsigned int degree = ::NTL::GF2EInfo->p.n;
	static unsigned int prevDegree = degree;

	if (degree != prevDegree){
		g_2_i = { getGF2E_X() };
		prevDegree = degree;
	}

	if (g_2_i.size() <= n){
		for (int i = g_2_i.size(); i <= n; ++i){
			g_2_i.push_back(NTL::GF2EInfo->power(g_2_i[i - 1], 2));
		}
	}
	return g_2_i[n];
}*/
const Algebra::FieldElement& getGenerator2toTheNthPower(const unsigned int n){
	static std::vector<Algebra::FieldElement> g_2_i = { Algebra::FieldElement(getGF2E_X()) };
	for (int i = g_2_i.size(); i <= n; ++i)
			g_2_i.push_back(Algebra::power(g_2_i[i - 1], 2));
	return g_2_i[n];
}

::NTL::GF2E getGeneratorPower(const unsigned int n){
	const int degree = ::NTL::GF2EInfo->p.n;
	const ::NTL::GF2E generator = getGF2E_X();
	return ::NTL::GF2EInfo->power(generator, n);
}



#define tempSelector
Algebra::CircuitPolynomial getSelector(int programLength, int instructionLine, Algebra::UnpackedWord unpakedPC){
	GADGETLIB_ASSERT(unpakedPC.size() >= Log2ceiled(programLength), "Number of unpacked bits used should be at least log length of the program");
	Algebra::FElem value = (instructionLine & 1U) ? Algebra::one() : Algebra::zero();
#ifndef tempSelector

	Algebra::CircuitPolynomial selectorPoly(Algebra::one() + value + unpakedPC[0]);

	instructionLine >>= 1U;
	int i = 1;
	while (instructionLine || i < unpakedPC.size()){
		Algebra::FElem value = (instructionLine & 1U) ? Algebra::one() : Algebra::zero();
		selectorPoly = selectorPoly * (Algebra::one() + value + unpakedPC[i]);
		i++;
		instructionLine >>= 1U;
	}
	return selectorPoly;
#else
	std::vector<Algebra::LinearCombination> lcVec;
	Algebra::LinearCombination lc(Algebra::one() + value + unpakedPC[0]);
	lcVec.push_back(lc);
	instructionLine >>= 1U;
	int i = 1;
	while (instructionLine || i < unpakedPC.size()){
		Algebra::FElem value = (instructionLine & 1U) ? Algebra::one() : Algebra::zero();
		lcVec.push_back(Algebra::one() + value + unpakedPC[i]);
		i++;
		instructionLine >>= 1U;
	}
	return Algebra::CircuitPolynomial(lcVec);

#endif
}

/*************************************************************************************************/
/*************************************************************************************************/
/*******************                                                            ******************/
/*******************					Memory Info				 				******************/
/*******************                                                            ******************/
/*************************************************************************************************/
/*************************************************************************************************/
MemoryInfo::MemoryInfo(){
	serialNumber_ = 0;
	isMemOp_ = 0;
	isLoadOp_ = 0;
	timestamp_ = Algebra::zero();
	timestampDegree_ = -1;
	address_ = Algebra::zero();
	value_ = Algebra::zero();
}

MemoryInfo::MemoryInfo(int serialNumber, bool isMemOp, bool isLoadOp,
	const Algebra::FElem& timestamp, int timestampDegree,
	const Algebra::FElem& address, const Algebra::FElem& value){
	serialNumber_ = serialNumber;
	isMemOp_ = isMemOp;
	isLoadOp_ = isLoadOp;
	timestamp_ = timestamp;
	timestampDegree_ = timestampDegree;
	address_ = address;
	value_ = value;
	GADGETLIB_ASSERT(timestamp == Algebra::FElem(getGeneratorPower(timestampDegree_)), "g^timestampDegree != timestamp");
}

void MemoryInfo::updateTimestamp(Algebra::FElem timestamp, int timestampDegree){
	GADGETLIB_ASSERT(timestamp == Algebra::FElem(getGeneratorPower(timestampDegree)), "g^timestampDegree != timestamp");
	timestamp_ = timestamp;
	timestampDegree_ = timestampDegree;
}

bool sortMemoryInfo(MemoryInfo a, MemoryInfo b){
	if (!a.getIsMemOp() && !b.getIsMemOp()){
		return a.getSerialNumber() < b.getSerialNumber();
	}
	else if (a.getIsMemOp() && !b.getIsMemOp()){
		return true;
	}
	else if (!a.getIsMemOp() && b.getIsMemOp()){
		return false;
	}
	else{ // MemOp = True -> Sort By Address
		int aAddr = mapFieldElementToInteger(0, 16, a.getAddress());
		int bAddr = mapFieldElementToInteger(0, 16, b.getAddress());
		if (aAddr < bAddr){
			return true;
		}
		else if (bAddr < aAddr){
			return false;
		}
		else{ // MemOp = True and Same Address -> Sort By TimeStamp
			return a.getTimestampDegree() < b.getTimestampDegree();
		}
	}
}



} // namespace
