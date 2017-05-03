#ifndef _COMMON_CONSTRAINTSLIB2_CONSTRAINTSYSTEMTOBREX_cs2BREX_HPP_
#define _COMMON_CONSTRAINTSLIB2_CONSTRAINTSYSTEMTOBREX_cs2BREX_HPP_

#include <algebraLib/FiniteField.hpp>

#include "gadgetlib/common_use.hpp"
#include "reductions/TinyRAMtoBREX/RamToContraintSystem/transitionFunction.hpp"
#include "reductions/TinyRAMtoBREX/RamToContraintSystem/transitionFunction.hpp"
#include "languages/BREX/ConstraintsSys.hpp"
#include "languages/BREX/BrexWitness.hpp"
#include "languages/BREX/BrexInstance.hpp"



namespace PCP_Project{

/*************************************************************************************************/
/*************************************************************************************************/
/*******************                                                            ******************/
/*******************                      cs2BREX				                ******************/
/*******************                                                            ******************/
/*************************************************************************************************/
/*************************************************************************************************/

class cs2BREX{
private:
	// Input Variables
	ProtoboardPtr pb_;
	TinyRAMProgram program_;

	// Internal Variables
	int transcript_len;
	// Is there any Memory Use in the program
	bool doesProgramUsesMemory_;
	// Memory and Regular constraints
	GadgetPtr transitionFunction_;
	GadgetPtr memoryConstraints_;
	// Trace Variables for regular and memory constraints.
	FollowingTraceVariables followingTraceVariable_;
	MemoryFollowingTraceVariables memoryFollowingTraceVariables_;
	// Translation vector between CS variables to BREX
	std::vector<Algebra::Variable> translation_;
    
	std::vector<std::vector<Algebra::FieldElement>> traceAssignmentTable_;
    std::map<size_t,size_t> memoryPermutation_;

	// Help Functions
	void initInitialVars();
	void checkMemoryUse();
	void copyTraceOutputValuesToTraceInput();
	void createTranslationVector();
	std::vector<Variable> variablesToVector(TraceVariables traceVariables);
	std::vector<Variable> memoryVariablesToVector(MemoryTraceVariables traceVariables);
	void generateMemoryConstraints() const;
	void boundaryConstraints() const;
	Algebra::Variable::set getStateVariables() const;

	//Functions
	void init();
	void generateConstraints();
	void generateWitness();
	void generateMemoryWitness();

    //VariableAssignment <-> vector<FieldElement> mappings
    std::vector<Algebra::FieldElement> assignmentToVec(const VariableAssignment& ass)const;
    VariableAssignment vectorToAssignment(const std::vector<Algebra::FieldElement>& vec)const;

public:
	cs2BREX(ProtoboardPtr pb, const TinyRAMProgram& program, int transcriptLen, bool constructWitness);

	ConstraintSystem getConstraints() const { return pb_->constraintSystem(Opcode::NONE); }
	ConstraintSystem getMemoryConstratins() const;

	const std::vector<std::vector<Algebra::FieldElement>>& getTraceAssignmet()const { return traceAssignmentTable_; }
	std::vector<Variable> getTranslationVector() const { return pb_->getTranslationVector(); }
	const std::map<size_t, size_t>& getMemoryPermutation()const { return memoryPermutation_; }
	size_t varsAmount() const; // Number of Variables in one traceLine
	BREXInstance::boundaryConstraints_t getBoundaryConstraints() const; 
	void printReductionData() const;
}; // cs2BREX


/*************************************************************************************************/
/*************************************************************************************************/
/*******************                                                            ******************/
/*******************                      cs2BREXConstraints	                ******************/
/*******************                                                            ******************/
/*************************************************************************************************/
/*************************************************************************************************/

class cs2BREXConstaraints :public ConstraintSys {
private:
	const cs2BREX& cs2brex_;
	polySet_t constraints_;
public:
	cs2BREXConstaraints(const cs2BREX& cs2brex);
	size_t varsAmount() const; // amount of variables in 1 traceline (including aux variables)
	const Algebra::FiniteField contextField() const;
	const polySet_t& constraints() const;
};

/*************************************************************************************************/
/*************************************************************************************************/
/*******************                                                            ******************/
/*******************                      cs2BREXColoring		                ******************/
/*******************                                                            ******************/
/*************************************************************************************************/
/*************************************************************************************************/

class cs2BREXColoring : public BREXWitness::assignment_t{
private:
	const std::vector<std::vector<FieldElement>>& traceAssignments_;
	std::vector<Variable> translation_;
public:
	cs2BREXColoring(const cs2BREX& cs2brex);
	BREXWitness::color_t getElementByIndex(index_t index) const;
};

/*************************************************************************************************/
/*************************************************************************************************/
/*******************                                                            ******************/
/*******************                      cs2BREXMemory			                ******************/
/*******************                                                            ******************/
/*************************************************************************************************/
/*************************************************************************************************/

// For now - stub. We don't use memory permutation
class cs2BREXMemory : public Sequence<size_t>{
private:
	const std::map<size_t, size_t>& memoryPermutation_;
public:
	cs2BREXMemory(const cs2BREX& cs2brex);
	size_t getElementByIndex(index_t index) const;

};

/*************************************************************************************************/
/*************************************************************************************************/
/*******************                                                            ******************/
/*******************                      cs2BREXCS				                ******************/
/*******************                                                            ******************/
/*************************************************************************************************/
/*************************************************************************************************/

//ForNow we use empty constraint system
class cs2BrexMemoryCS : public ConstraintSys {
private:
	const cs2BREX& cs2brex_;
	std::vector<Variable> translation_;
	polySet_t constraints_;
	
public:
	cs2BrexMemoryCS(const cs2BREX& cs2brex);

	const Algebra::FiniteField contextField() const;

	size_t varsAmount() const;

	const polySet_t& constraints() const;

};


} // namespace
#endif // _COMMON_CONSTRAINTSLIB2_CONSTRAINTSYSTEMTOBREX_cs2BREX_HPP_
