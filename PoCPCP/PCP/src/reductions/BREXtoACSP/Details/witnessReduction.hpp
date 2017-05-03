#ifndef BREX_TO_ACSP_WITNESS_REDUCTION_HPP__
#define BREX_TO_ACSP_WITNESS_REDUCTION_HPP__

#include "common.hpp"
#include "commonMappings.hpp"
#include "languages/ACSP/ACSPWitness.hpp"
#include "languages/BREX/BrexWitness.hpp"
#include "languages/BREX/BrexInstance.hpp"
#include "witnessMappings.hpp""

#include <algebraLib/UnivariatePolynomialGeneral.hpp>
#include <memory>

namespace PCP_Project{
namespace BREXtoACSP{

class witnessReduction{
public:
    static std::unique_ptr<ACSPWitness> reduceWitness(const BREXInstance& instance, const BREXWitness& witness);
    
protected:
    typedef std::vector<Algebra::FieldElement> evaluation_t;
    
    static evaluation_t getEmbeddingMapping( const BREXInstance& instance, const BREXWitness& witness, const common& commonDef, const witnessMappings& witnessMapping);
    static void mapChi(const BREXInstance& instance, const BREXWitness& witness, evaluation_t& mapping, const common& commonDef, const witnessMappings& witnessMapping);
    static void mapNetwork(const BREXInstance& instance, const BREXWitness& witness, evaluation_t& mapping, const common& commonDef, const witnessMappings& witnessMapping);
};
    
} //namespace BREXtoACSP
} //namespace PCP_Project


#endif // BREX_TO_ACSP_WITNESS_REDUCTION_HPP__
