#include "cs2BREX.hpp"
#include "../RamToContraintSystem/MemoryConsraints.hpp"


namespace PCP_Project{

/*************************************************************************************************/
/*************************************************************************************************/
/*******************                                                            ******************/
/*******************                      cs2BREX				                ******************/
/*******************                                                            ******************/
/*************************************************************************************************/
/*************************************************************************************************/

cs2BREX::cs2BREX(ProtoboardPtr pb,const TinyRAMProgram& program, const int transcriptLen, bool constructWitness) : pb_(pb),
	program_(program),transcript_len(transcriptLen),
	followingTraceVariable_(program.pcLength(), std::dynamic_pointer_cast<const TinyRAMProtoboardParams>(pb_->params())->numRegisters()),
	memoryFollowingTraceVariables_(followingTraceVariable_.first_.timeStamp_, followingTraceVariable_.second_.timeStamp_),
	doesProgramUsesMemory_(false){
		this->init();
		this->generateConstraints();
		this->boundaryConstraints();
        if(constructWitness){
            generateWitness();
            generateMemoryWitness();
        }
}

void cs2BREX::init() {
	// Init the First Values of timestamp
	pb_->addDegreeTranslation(Algebra::FieldElement(getGF2E_X()), 1);
	pb_->addDegreeTranslation(Algebra::one(), 0);
	// GadgetPtr
	transitionFunction_ = TransitionFunction::create(pb_, followingTraceVariable_,memoryFollowingTraceVariables_ , program_);
	checkMemoryUse();
	if (doesProgramUsesMemory_){
		memoryConstraints_ = MemoryConstraints::create(pb_, memoryFollowingTraceVariables_);
	}
}

void cs2BREX::boundaryConstraints() const{
	::std::shared_ptr<const TinyRAMProtoboardParams> params = std::dynamic_pointer_cast<const TinyRAMProtoboardParams>(pb_->params());
	pb_->addBoundaryConstraint(followingTraceVariable_.first_.flag_, 0, Algebra::zero());
	pb_->addBoundaryConstraint(followingTraceVariable_.first_.timeStamp_, 0, Algebra::one());
	for (int i = 0; i < program_.pcLength(); i++){
		pb_->addBoundaryConstraint(followingTraceVariable_.first_.pc_[i], 0, Algebra::zero()); 
	}
	//for (int i = 0; i < params->numRegisters(); i++){
	//	pb_->addBoundaryConstraint(followingTraceVariable_.first_.registers_[i], 0, Algebra::zero());
	//}
}

void cs2BREX::initInitialVars(){
	pb_->val(followingTraceVariable_.first_.flag_) = Algebra::zero();
	pb_->val(followingTraceVariable_.first_.timeStamp_) = Algebra::one();
	int pcLength = program_.pcLength();
	int numRegisters = std::dynamic_pointer_cast<const TinyRAMProtoboardParams>(pb_->params())->numRegisters();
	for (int i = 0; i < pcLength; i++){
		pb_->val(followingTraceVariable_.first_.pc_[i]) = Algebra::zero();
	}
	for (int i = 0; i < numRegisters; i++){
		pb_->val(followingTraceVariable_.first_.registers_[i]) = Algebra::zero();
	}
}

void cs2BREX::checkMemoryUse(){
	for (int i = 0; i < program_.code().size(); ++i){
		Opcode opcode = program_.code()[i].opcode_;
		if (opcode == Opcode::STOREB || opcode == Opcode::STOREW ||
			opcode == Opcode::LOADB || opcode == Opcode::LOADW){
			doesProgramUsesMemory_ = true;
			break;
		}
	}
}

std::vector<Variable> cs2BREX::variablesToVector(TraceVariables traceVariables){
	std::vector<Variable> v;
	v.emplace_back(traceVariables.flag_);
	v.emplace_back(traceVariables.timeStamp_);
	::std::shared_ptr<const TinyRAMProtoboardParams> params = std::dynamic_pointer_cast<const TinyRAMProtoboardParams>(pb_->params());
	for (int i = 0; i < traceVariables.pc_.size(); i++){
		v.emplace_back(traceVariables.pc_[i]);
	}
	for (int i = 0; i < params->numRegisters(); i++){
		v.emplace_back(traceVariables.registers_[i]);
	}
	return v;
}

std::vector<Variable> cs2BREX::memoryVariablesToVector(MemoryTraceVariables traceVariables){
	std::vector<Variable> v;
	v.emplace_back(traceVariables.isMemOp_);
	v.emplace_back(traceVariables.isLoad_);
	v.emplace_back(traceVariables.timeStamp_);
	v.emplace_back(traceVariables.address_);
	v.emplace_back(traceVariables.value_);
	return v;
}

void cs2BREX::copyTraceOutputValuesToTraceInput(){
	pb_->val(followingTraceVariable_.first_.timeStamp_) = pb_->val(followingTraceVariable_.second_.timeStamp_);
	pb_->val(followingTraceVariable_.first_.flag_) = pb_->val(followingTraceVariable_.second_.flag_);
	::std::shared_ptr<const TinyRAMProtoboardParams> params = std::dynamic_pointer_cast<const TinyRAMProtoboardParams>(pb_->params());
	int pcLength = program_.pcLength();
	for (int i = 0; i < pcLength; i++){
		pb_->val(followingTraceVariable_.first_.pc_[i]) = pb_->val(followingTraceVariable_.second_.pc_[i]);
	}
	for (int i = 0; i < params->numRegisters(); i++){
		pb_->val(followingTraceVariable_.first_.registers_[i]) = pb_->val(followingTraceVariable_.second_.registers_[i]);
	}
}


void cs2BREX::createTranslationVector(){
	Algebra::Variable::set traceFirstVariables;
	Algebra::Variable::set traceSecondVariables;
	
	std::vector<Variable> trace1 = variablesToVector(followingTraceVariable_.first_);
	std::vector<Variable> trace2 = variablesToVector(followingTraceVariable_.second_);
	std::vector<Variable> memoryTrace1 = memoryVariablesToVector(memoryFollowingTraceVariables_.first_);
	std::vector<Variable> memoryTrace2 = memoryVariablesToVector(memoryFollowingTraceVariables_.second_);
	
	traceFirstVariables.insert(trace1.begin(), trace1.end());
	traceFirstVariables.insert(memoryTrace1.begin(), memoryTrace1.end());
	traceSecondVariables.insert(trace2.begin(), trace2.end());
	traceSecondVariables.insert(memoryTrace2.begin(), memoryTrace2.end());

	Algebra::Variable::set auxVars = pb_->constraintSystem(Opcode::NONE).getUsedVariables();
	Algebra::Variable::set memoryauxVars = pb_->constraintSystem(Opcode::MEMORY).getUsedVariables();
	auxVars.insert(memoryauxVars.begin(), memoryauxVars.end());
	for (const auto& var : trace1) {
		auxVars.erase(var);

	}
	for (const auto& var : trace2) {
		auxVars.erase(var);
	}
	for (const auto& var : memoryTrace1) {
		auxVars.erase(var);
	}
	for (const auto& var : memoryTrace2) {
		auxVars.erase(var);
	}
 	translation_.insert(translation_.end(), traceFirstVariables.cbegin(), traceFirstVariables.cend());
	translation_.insert(translation_.end(), auxVars.cbegin(), auxVars.cend());
	translation_.insert(translation_.end(), traceSecondVariables.cbegin(), traceSecondVariables.cend());
	Algebra::Variable unUsed("Unused aux Variable");
	for (size_t i = 0; i < auxVars.size(); i++){
		translation_.emplace_back(unUsed);
	}
	// Jenya - pretty ugly - we leave it like this for now - we'll change it as soon as I have time.
	pb_->setNewIndicesForTranslation(translation_);
	translation_ = pb_->getTranslationVector();
}

void cs2BREX::generateConstraints(){
	transitionFunction_->generateConstraints();
	if (doesProgramUsesMemory_){
		memoryConstraints_->generateConstraints();
	}
	createTranslationVector();
}

//VariableAssignment <-> vector<FieldElement> mappings
std::vector<Algebra::FieldElement> cs2BREX::assignmentToVec(const VariableAssignment& ass)const{
    std::vector<Algebra::FieldElement> res(translation_.size()/2);
    for(size_t i=0; i< res.size(); i++){
        res[i] = ass.at(translation_[i]);
    }
    return res;
}

VariableAssignment cs2BREX::vectorToAssignment(const std::vector<Algebra::FieldElement>& vec)const{
    VariableAssignment res;
    for(size_t i=0; i< translation_.size()/2; i++){
        res[translation_[i]] = vec[i];
    }
    return res;
}

//#define falseWitness//checking how good PCP is at finding small errors
void cs2BREX::generateWitness(){
	//const unsigned int transcript_len = POW2(TRANSCIPT_LEN_LOG) - 1;
	// First Assignment should be zero
    initInitialVars();
    for (int i = 0; i < transcript_len; ++i){
		::std::dynamic_pointer_cast<TransitionFunction>(transitionFunction_)->generateWitness(i);
        Algebra::VariableAssignment assignment = pb_->assignment();
        traceAssignmentTable_.push_back(assignmentToVec(assignment));
#ifdef falseWitness
		if (i == 3)
			traceAssignmentTable_[i][4] = Algebra::xFE();
#endif
		copyTraceOutputValuesToTraceInput();
	}
}

ConstraintSystem cs2BREX::getMemoryConstratins() const{
	if (doesProgramUsesMemory_){
		return pb_->constraintSystem(Opcode::MEMORY); 
	}
	return ConstraintSystem();
}

void cs2BREX::generateMemoryWitness(){

    Variable::set usedMemoryVariables = pb_->constraintSystem(Opcode::MEMORY).getUsedVariables();
	usedMemoryVariables.erase(memoryFollowingTraceVariables_.first_.timeStamp_);
	usedMemoryVariables.erase(memoryFollowingTraceVariables_.first_.isMemOp_); 
	usedMemoryVariables.erase(memoryFollowingTraceVariables_.first_.isLoad_);
	usedMemoryVariables.erase(memoryFollowingTraceVariables_.first_.value_);
	usedMemoryVariables.erase(memoryFollowingTraceVariables_.first_.address_);
	usedMemoryVariables.erase(memoryFollowingTraceVariables_.second_.timeStamp_);
	usedMemoryVariables.erase(memoryFollowingTraceVariables_.second_.isMemOp_);
	usedMemoryVariables.erase(memoryFollowingTraceVariables_.second_.isLoad_);
	usedMemoryVariables.erase(memoryFollowingTraceVariables_.second_.value_);
	usedMemoryVariables.erase(memoryFollowingTraceVariables_.second_.address_);

	std::vector<MemoryInfo> memoryTrace = pb_->getMemoryTrace();
	std::sort(memoryTrace.begin(), memoryTrace.end(), sortMemoryInfo);
	GADGETLIB_ASSERT(memoryTrace.size() == traceAssignmentTable_.size(), "memoryInfo size should be the same as the coloring");
	for (int i = 0; i < memoryTrace.size(); i++){
		size_t serialNumber1 = memoryTrace[i].getSerialNumber();
		size_t serialNumber2 = memoryTrace[(i + 1) % memoryTrace.size()].getSerialNumber();
		memoryPermutation_[serialNumber1] = serialNumber2;
		if (doesProgramUsesMemory_){
			VariableAssignment currAssignment = vectorToAssignment(traceAssignmentTable_[serialNumber1]);
			VariableAssignment nextAssignment = vectorToAssignment(traceAssignmentTable_[serialNumber2]);
			pb_->clearAssignment();
			pb_->val(memoryFollowingTraceVariables_.first_.timeStamp_) = currAssignment[memoryFollowingTraceVariables_.first_.timeStamp_];
			pb_->val(memoryFollowingTraceVariables_.first_.isMemOp_) = currAssignment[memoryFollowingTraceVariables_.first_.isMemOp_];
			pb_->val(memoryFollowingTraceVariables_.first_.isLoad_) = currAssignment[memoryFollowingTraceVariables_.first_.isLoad_];
			pb_->val(memoryFollowingTraceVariables_.first_.value_) = currAssignment[memoryFollowingTraceVariables_.first_.value_];
			pb_->val(memoryFollowingTraceVariables_.first_.address_) = memoryTrace[i].getAddress();

			pb_->val(memoryFollowingTraceVariables_.second_.timeStamp_) = nextAssignment[memoryFollowingTraceVariables_.first_.timeStamp_];
			pb_->val(memoryFollowingTraceVariables_.second_.isMemOp_) = nextAssignment[memoryFollowingTraceVariables_.first_.isMemOp_];
			pb_->val(memoryFollowingTraceVariables_.second_.isLoad_) = nextAssignment[memoryFollowingTraceVariables_.first_.isLoad_];
			pb_->val(memoryFollowingTraceVariables_.second_.value_) = nextAssignment[memoryFollowingTraceVariables_.first_.value_];
			pb_->val(memoryFollowingTraceVariables_.second_.address_) = memoryTrace[(i+1) % memoryTrace.size()].getAddress();
			if (Algebra::one() == pb_->val(memoryFollowingTraceVariables_.first_.isMemOp_)){

				memoryConstraints_->generateWitness();

				if (i != memoryTrace.size() - 1){
					GADGETLIB_ASSERT(pb_->isSatisfied(Opcode::MEMORY), "MemoryConstraints are not satisfied");
					for (const Variable& v : usedMemoryVariables){
						currAssignment[v] = pb_->val(v);
					}
				}
				else{
					for (const Variable& v : usedMemoryVariables){
						currAssignment[v] = Algebra::zero();
					}
				}
				traceAssignmentTable_[serialNumber1] = assignmentToVec(currAssignment);
			}
		}
	}
}


size_t cs2BREX::varsAmount() const{
	std::vector<Algebra::Variable> translation = pb_->getTranslationVector();
	GADGETLIB_ASSERT(translation.size() % 2 == 0, "translation vector size is expected to be even.");
	/*
	return translation_.size() / 2;
	*/
	// Michael said that he wants the size of all the variables.
	return translation.size();
}


BREXInstance::boundaryConstraints_t cs2BREX::getBoundaryConstraints() const{
	BREXInstance::boundaryConstraints_t boundaryConstraints;

	BoundaryVariables boundaryVariables = pb_->boundaryVariables();
	BoundaryTimestamps boundaryTimestamps = pb_->boundaryTimestamps();
	BoundaryAssignments boundaryAssignment = pb_->boundaryAssignments();

	_COMMON_ASSERT(boundaryVariables.size() == boundaryTimestamps.size(),
		"Number of variables should be equal to the number of timestamps in boundary constaints");
	_COMMON_ASSERT(boundaryAssignment.size() == boundaryVariables.size(),
		"Number of variabled should be equal to the number os assignments in boundary constraints");

	for (int i = 0; i < boundaryVariables.size(); ++i){
		for (int j = 0; j < translation_.size(); ++j){
			if (translation_[j] == boundaryVariables[i]){
				BREXInstance::point_t location(boundaryTimestamps[i], j);
				boundaryConstraints[location] = boundaryAssignment[i];
			}
		}
	}
	return boundaryConstraints;
}

Algebra::Variable::set cs2BREX::getStateVariables() const {
	Algebra::Variable::set retSet;
	retSet.insert(followingTraceVariable_.first_.flag_);
	retSet.insert(followingTraceVariable_.second_.flag_);
	retSet.insert(followingTraceVariable_.first_.timeStamp_);
	retSet.insert(followingTraceVariable_.second_.timeStamp_);
	GADGETLIB_ASSERT(followingTraceVariable_.first_.pc_.size() == followingTraceVariable_.second_.pc_.size(), 
							"CS2BREX: unpacked pc should have the exact same size in both of the states");
	for (int i = 0; i < followingTraceVariable_.first_.pc_.size(); i++) {
		retSet.insert(followingTraceVariable_.first_.pc_[i]);
		retSet.insert(followingTraceVariable_.second_.pc_[i]);
	}
	GADGETLIB_ASSERT(followingTraceVariable_.first_.registers_.size() == followingTraceVariable_.second_.registers_.size(),
									"CS2BREX: number of registers should be the same in both of the states");
	for (int i = 0; i < followingTraceVariable_.first_.registers_.size(); i++) {
		retSet.insert(followingTraceVariable_.first_.registers_[i]);
		retSet.insert(followingTraceVariable_.second_.registers_[i]);
	}
	return retSet;
}

void cs2BREX::printReductionData() const{
	std::cout << std::endl;
	std::cout << "================== RAM DATA ==================" << std::endl;
	std::cout << "PROGRAM LENGTH: " << program_.size() << endl;
	GADGETLIB_ASSERT(followingTraceVariable_.first_.pc_.size() == followingTraceVariable_.second_.pc_.size(),
									"CS2Brex: unpacked PC size should be equal in both of the states");
	std::cout << "UNPACKED PC LENGTH: " << followingTraceVariable_.first_.pc_.size() << endl;
	GADGETLIB_ASSERT(followingTraceVariable_.first_.registers_.size() == followingTraceVariable_.second_.registers_.size(),
									"CS2BREX: number of registers should be the same in both of the states");
	std::cout << "#REGISTERS (without PC): " << followingTraceVariable_.second_.registers_.size() << endl;
	GADGETLIB_ASSERT(followingTraceVariable_.first_.size() == followingTraceVariable_.second_.size(),
									"CS2BREX: number of registers should be the same in both of the states");
	std::cout << "#STATE VARIABLES: " << followingTraceVariable_.first_.size() << endl;
	ConstraintSystem regularConstraints = pb_->constraintSystem(Opcode::NONE);
	int numConstraints = regularConstraints.getConstraints().size();
	std::cout << "#REGULAR constraints: " << numConstraints << std::endl;
	int numMemoryConstraints = doesProgramUsesMemory_ ? pb_->constraintSystem(Opcode::MEMORY).getConstraints().size() : 0;
	std::cout << "#MEMORY constraints: " << numMemoryConstraints << std::endl;
	Algebra::Variable::set regularUsedVariables = regularConstraints.getUsedVariables();
	std::cout << "#VARIABLES REGULAR constraints: " << regularUsedVariables.size() << std::endl;
	Algebra::Variable::set memoryUsedVariables = doesProgramUsesMemory_ ? pb_->constraintSystem(Opcode::MEMORY).getUsedVariables() : Algebra::Variable::set(); 
	std::cout << "#VARIABLES MEMORY constraints: " << memoryUsedVariables.size() << std::endl;
	Algebra::Variable::set allUsedVariables = Algebra::Variable::set();
	allUsedVariables.insert(regularUsedVariables.begin(), regularUsedVariables.end());
	allUsedVariables.insert(memoryUsedVariables.begin(), memoryUsedVariables.end());
	std::cout << "#VARIABLES: " << allUsedVariables.size() << std::endl;
	int registerLength = std::dynamic_pointer_cast<const TinyRAMProtoboardParams>(pb_->params())->registerLength();
	std::cout << "REGISTER LENGTH (except PC): " << registerLength << endl;
	std::cout << "================ END RAM DATA ================" << std::endl;
	std::cout << std::endl;
}

/*************************************************************************************************/
/*************************************************************************************************/
/*******************                                                            ******************/
/*******************                      cs2BREXConstraints	                ******************/
/*******************                                                            ******************/
/*************************************************************************************************/
/*************************************************************************************************/


cs2BREXConstaraints::cs2BREXConstaraints(const cs2BREX& cs2brex) : cs2brex_(cs2brex){
	ConstraintSystem cs = cs2brex_.getConstraints();
	Variable::set usedVars = cs.getUsedVariables();
	constraints_ = cs.getConstraintPolynomials();

};
const Algebra::FiniteField cs2BREXConstaraints::contextField() const{
	return  Algebra::FiniteField(NTL::GF2EInfo->p);
}

size_t cs2BREXConstaraints::varsAmount() const{
	return cs2brex_.varsAmount();
}
const ConstraintSys::polySet_t& cs2BREXConstaraints::constraints() const{
	return constraints_;
}

/*************************************************************************************************/
/*************************************************************************************************/
/*******************                                                            ******************/
/*******************                      cs2BREXColoring		                ******************/
/*******************                                                            ******************/
/*************************************************************************************************/
/*************************************************************************************************/

cs2BREXColoring::cs2BREXColoring(const cs2BREX& cs2brex) : traceAssignments_(cs2brex.getTraceAssignmet()){
	translation_ = cs2brex.getTranslationVector();
};

BREXWitness::color_t cs2BREXColoring::getElementByIndex(index_t index) const{
	GADGETLIB_ASSERT(index < traceAssignments_.size(), "Attempted to get an illegal element");
	return traceAssignments_[index];
}

/*************************************************************************************************/
/*************************************************************************************************/
/*******************                                                            ******************/
/*******************                      cs2BREXMemory			                ******************/
/*******************                                                            ******************/
/*************************************************************************************************/
/*************************************************************************************************/

cs2BREXMemory::cs2BREXMemory(const cs2BREX& cs2brex) : memoryPermutation_(cs2brex.getMemoryPermutation()){}

size_t cs2BREXMemory::getElementByIndex(index_t index) const{
	return memoryPermutation_.at(index);
}

/*************************************************************************************************/
/*************************************************************************************************/
/*******************                                                            ******************/
/*******************                      cs2BREXCS				                ******************/
/*******************                                                            ******************/
/*************************************************************************************************/
/*************************************************************************************************/

cs2BrexMemoryCS::cs2BrexMemoryCS(const cs2BREX& cs2brex) : cs2brex_(cs2brex){
	ConstraintSystem cs = cs2brex_.getMemoryConstratins();
	translation_ = cs2brex.getTranslationVector();
	constraints_ = cs.getConstraintPolynomials();
};

const Algebra::FiniteField cs2BrexMemoryCS::contextField() const {
	// check NTL global degree and return field
	const size_t extensionDegree = ::NTL::deg(::NTL::GF2EInfo->p);
	return Algebra::FiniteField(extensionDegree);
}

size_t cs2BrexMemoryCS::varsAmount() const{
	return cs2brex_.varsAmount();
}

const ConstraintSys::polySet_t& cs2BrexMemoryCS::constraints() const {
	return constraints_;
}



} //namespace
