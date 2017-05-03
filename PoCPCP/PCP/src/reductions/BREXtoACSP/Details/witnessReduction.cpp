#include "witnessReduction.hpp"
#include "../Routing/LongSymmetricDeBruijnNetwork.hpp"
#include "common/Utils/TaskReporting.hpp"

#include <vector>

namespace PCP_Project{
namespace BREXtoACSP{

using std::vector;
using std::unique_ptr;
using Algebra::UnivariatePolynomialGeneral;
using Algebra::FieldElement;
using Algebra::zero;
using Algebra::one;

typedef LongSymmetricDeBruijnNetwork permNetwork_t;

unique_ptr<ACSPWitness> witnessReduction::reduceWitness( const BREXInstance& instance, const BREXWitness& witness){
    
    TASK("Reducing BREX witness to ACSP");
    
    //get common information
    common commonDef(instance);
    witnessMappings witnessMapping(commonDef);

    //get the embedding
    evaluation_t mapping = getEmbeddingMapping(instance, witness, commonDef,witnessMapping);

    //ordered basis for witness space
    const auto& basis = witnessMapping.getImageSpaceOrderedBasis();
//#define ZERO_WITNESS  //sometimes for debugging it's convenient to use an identically zero witness
#ifdef ZERO_WITNESS 
	for (auto& p : mapping)
		p = zero();
#endif
    //construct witness
    //and return it
    unique_ptr<const ACSPWitness::polynomial> Assignment_ptr(new UnivariatePolynomialGeneral(mapping,basis,zero()));
	
//#ifdef IDENTITY_WITNESS
//	unique_ptr<UnivariatePolynomialGeneral> idPoly(new UnivariatePolynomialGeneral(std::vector<FieldElement>({ Algebra::zero(), Algebra::one() })));
//	return unique_ptr<ACSPWitness>(new ACSPWitness(move(idPoly)));
//#endif //IDENTITY_WITNESS

	unique_ptr<ACSPWitness> witness_ptr(new ACSPWitness(move(Assignment_ptr)));
    return move(witness_ptr);
}

class exponentsPermutation : public Sequence<size_t>{
public:
    exponentsPermutation(const BREXInstance& instance):instance_(instance),commonDef_(instance){};

	size_t getElementByIndex(index_t index)const{
        const size_t singeltonIndex = instance_.domainSize();
        const auto numBits = commonDef_.heightSpaceDimension();
           
        //singleton case 
        if(index == singeltonIndex) return 0;

        //general case
        return expModulu(index);
    }

private:
    const BREXInstance& instance_;
    const common commonDef_;

    size_t expModulu(const size_t exp)const{
        using NTL::GF2X;
        using NTL::SetCoeff;

        const GF2X primitivePoly = commonDef_.rowsModulus();
        const size_t numBits = commonDef_.heightSpaceDimension();
        GF2X x_power;
        SetCoeff(x_power,exp);
        const GF2X resPoly = x_power % primitivePoly;
        size_t res=0;
        for(int i=numBits-1; i>=0; i--){
            short bitVal = NTL::IsOne(NTL::coeff(resPoly,i));
            res*=2;
            res+=bitVal;
        }
        return res;
    }
};

class inversePermutation : public Sequence<size_t>{
    public :
        inversePermutation(const Sequence<size_t>& src, const size_t numElements): seq_(numElements){
            for(size_t val=0; val< numElements; val++){
                const size_t index = src.getElementByIndex(val);
                seq_[index] = val;
            }
        }
		size_t getElementByIndex(index_t index)const{
            if (index < seq_.size()) return seq_[index];
            else {
                _COMMON_FATAL("Access to such index is unexpected");
            }
        }
    private:
        vector<size_t> seq_;
};

//represents the permutation \f$ g \circ h \circ g^{-1} \f$
class conjugatePermutation : public Sequence<size_t>{
    public:
        conjugatePermutation(const Sequence<size_t>& g, const Sequence<size_t>& h, const size_t numElements):
            g_(g), h_(h), g_inv_(g,numElements){};
        
		size_t getElementByIndex(index_t index)const{
            const size_t v1 = g_inv_.getElementByIndex(index);
            const size_t v2 = h_.getElementByIndex(v1);
            const size_t v3 = g_.getElementByIndex(v2);

            return v3;
        } 

    private:
        const Sequence<size_t>& g_;
        const Sequence<size_t>& h_;
        const inversePermutation g_inv_;
};

class addSingeltonToPermutation : public Sequence<size_t>{
    public:
        addSingeltonToPermutation(const Sequence<size_t>& orig, const size_t singletoneIndex): origPerm_(orig) ,singletoneIndex_(singletoneIndex){};
        size_t getElementByIndex(index_t index)const{
            if (index == singletoneIndex_) return singletoneIndex_;
            return origPerm_.getElementByIndex(index);
        }
    
    private:
        const Sequence<size_t>& origPerm_;
        const size_t singletoneIndex_;
};

//Mapping
witnessReduction::evaluation_t witnessReduction::getEmbeddingMapping( const BREXInstance& instance, const BREXWitness& witness, const common& commonDef, const witnessMappings& witnessMapping){
 
    //set the context field
    commonDef.contextField().setContext(); 

    //define the mapping on which the assignment polynomial will be interpolated on
    //This mapping is the arithmetization of the witness
    const size_t mappingSize = Infrastructure::POW2(witnessMapping.getImageSpaceOrderedBasis().size());
    evaluation_t mapping(mappingSize,Algebra::zero());
    
    //Map the coloring (assignment)
    mapChi(instance,witness,mapping,commonDef,witnessMapping);

    //Map the routing network (Pi) if needed
    if(commonDef.hasRoutingNetwork()){
        mapNetwork(instance,witness,mapping,commonDef,witnessMapping);
    }
    
    //return the mapping
    return mapping;
}

void witnessReduction::mapChi(const BREXInstance& instance, const BREXWitness& witness, evaluation_t& mapping, const common& commonDef, const witnessMappings& witnessMapping){
    const size_t cyclicDomainSize = instance.domainSize();
    
    const vector<size_t> unroutedVars = commonDef.variablesNonPerm();

    //ordered basis for witness space
    const auto& basis = witnessMapping.getImageSpaceOrderedBasis();

    //Map the coloring of the circle
    {
    auto currRow_spaceIndex = witnessMapping.getFirstRow_spaceIndex();
    for( size_t vecId =0; vecId < cyclicDomainSize; vecId++){
        const auto& assignment = witness.get_color(vecId);
        for (const size_t& varId :  unroutedVars){
            const size_t varIndex = commonDef.getVarLocation(varId).index;
            const size_t indicator = witnessMapping.mapIndexOfNonPermutationVariable_spaceIndex(currRow_spaceIndex,varIndex);
            mapping[indicator] = assignment[varId];
        }
        currRow_spaceIndex = witnessMapping.getNextRow_spaceIndex(currRow_spaceIndex);
    }

    }

}

void witnessReduction::mapNetwork(const BREXInstance& instance, const BREXWitness& witness, evaluation_t& mapping, const common& commonDef, const witnessMappings& witnessMapping){
    
    /// We want the "log" permutation, the one that maps:
    /// \f$ g^i \mapsto i \f$, and for the special case: \f$ 0 \mapsto singletoneIndex \f$
    /// The log permutation defines the mapping from the label in the first column
    /// of the DeBruijn to the vectorId it should hold,
    /// thus it defines the mapping from the values in the network to vectors ids
    const exponentsPermutation expPerm(instance);
    const inversePermutation logPerm(expPerm, instance.domainSize()+1);

    //Add another element to permutation as a singleton, so we can rout it using DeBruijn
    const addSingeltonToPermutation expendedPerm(witness.permutation(),instance.domainSize());
    
    /// Because we mess with the indexes, we have to change the permutation as well
    /// If the good permutation was \f$ i \mapsto \pi(i) \f$
    /// after changing the numbering, using some permutation \f$ \sigma \f$
    /// if the routing network will (originally) rout \f$ \pi \f$ we
    /// would get the permutation \f$ \sigma(i) \mapsto \sigma \circ \pi(i) \f$
    /// But what we want is \f$ \sigma(i) \mapsto \pi \circ \sigma(i) \f$,
    /// So instead we rout the conjugate permutation \f$ \sigma^{-1} \circ \pi \circ \sigma \f$
    /// This way we get: 
    /// \f$ \sigma(i) \mapsto \sigma \circ \sigma^{-1} \circ \pi \circ \sigma(i) = \pi \circ \sigma(i) \f$
    const conjugatePermutation permToRout(expPerm,expendedPerm, instance.domainSize()+1);

    //map routing network
    permNetwork_t net(instance.domainSizeIndicator());
    net.rout(permToRout);
   
    const vector<size_t> routedIndexes = commonDef.variablesPerm();

    //define the index for the additional vector
    const size_t addedVecId = instance.domainSize();
    
    //ordered basis for witness space
    const auto& basis = witnessMapping.getImageSpaceOrderedBasis();

    //map the main routing network
    {
    
    for (size_t rowId=0; rowId < net.height(); rowId++){
        for(size_t columnId = 0; columnId < net.getWingWidth() ; columnId++){
            for(short netPartId = 0; netPartId < 2; netPartId++){//loop over the two halves

                const RoutingNetwork::dataID_t dataId = net.getDataID(netPartId, columnId, rowId);
                const size_t vecId = logPerm.getElementByIndex(dataId);
                const auto& coloring =  (vecId == addedVecId? instance.paddingPi() : witness.get_color(vecId));
                
                //map routing network
                for (size_t packetId = 0; packetId < routedIndexes.size(); packetId++){
                    const auto currPacketIndex = routedIndexes[packetId];
                    const FieldElement val = coloring[currPacketIndex];
                    const auto indicator_index = witnessMapping.mapNetworkElement_spaceIndex(rowId,columnId,2*packetId + netPartId);
                    mapping[indicator_index] = val;
                }

                //map routing bits
                if( columnId < net.getWingWidth()-1){
                    const short bitVal = net.routingBit(netPartId,columnId,rowId);
                    const auto indicator_index = witnessMapping.mapNetworkRoutingBit_spaceIndex(rowId,columnId,netPartId);
                    
                    switch(bitVal){
                        case 0: mapping[indicator_index] = zero(); break;
                        case 1: mapping[indicator_index] = one(); break;
                        default : _COMMON_FATAL("Bad value of bit");
                    }
                }
            }
        }
    }
    
    }
}
    
} //namespace BREXtoACSP
} //namespace PCP_Project
