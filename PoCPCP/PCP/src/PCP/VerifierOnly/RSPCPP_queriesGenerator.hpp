#ifndef RSPCPP_QUERIES_GEN_HPP__
#define RSPCPP_QUERIES_GEN_HPP__

#include "languages/ACSP/ACSPInstance.hpp"
#include <algebraLib/FieldElement.hpp>

#include <vector>

namespace PCP_Project{
namespace Verifier{

//
// Common
//

class ResultLocation{
public:
    void addAnswerPtr(Algebra::FieldElement* ptr);
    void answer(const Algebra::FieldElement& res)const;

private:
    std::vector<Algebra::FieldElement*> answerLocations_;
};

//
// Reed-Solomon PCPP related
//

class RS_PCPP_result{
public:
    std::vector<Algebra::FieldElement> results;
    
    void init(const std::vector<Algebra::FieldElement>& basis, const short degBound);
    bool verify()const;
private:
    std::vector<Algebra::FieldElement> basis_;
    short degBound_;
};

struct RS_PCPP_queriesTree{
    std::map<size_t,ResultLocation> localQueries;
    std::map<size_t,RS_PCPP_queriesTree> rowsToExtract;
    std::map<size_t,RS_PCPP_queriesTree> columnsToExtract;
};

std::map<size_t,Algebra::FieldElement*> addRandomQuery(const short degBound_logCeil, const std::vector<Algebra::FieldElement>& evaluationBasis, const Algebra::FieldElement& affineShift, RS_PCPP_queriesTree& queries, RS_PCPP_result& results, const size_t depth, const bool isRoot = true);

//
// ACSP consistency related
//

class ACSP_CONSISTENCY_result{
public:
    Algebra::FieldElement compositionPoly_res;
    std::vector<Algebra::FieldElement> boundaryPoly_res;

    void init(const ACSPInstance& instance, const Algebra::FieldElement& consistencyPoint);
    bool verify()const;

private:
    const ACSPInstance* instance_;
    Algebra::FieldElement consistencyPoint_;
    std::vector<Algebra::FieldElement> boundaryVanishingVals_;
    std::vector<Algebra::FieldElement> boundaryPolyVals_;
    Algebra::FieldElement vanishingSpacePolyVal_;
};

class ACSP_CONSISTENCY_query{
public:
    ACSP_CONSISTENCY_query(const Algebra::FieldElement& consistencyPoint_, ACSP_CONSISTENCY_result& resLocation);
    const Algebra::FieldElement consistencyPoint;
    ResultLocation compositionPoly_res;
    std::vector<ResultLocation> boundaryPoly_res;
};

} //namespace Verifier
} // namespace PCP_Project

#endif // RSPCPP_QUERIES_GEN_HPP__
