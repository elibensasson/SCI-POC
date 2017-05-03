#include "RSPCPP_queriesGenerator.hpp"
#include "PCP/PCP_common.hpp"
#include "common/Algebra/ShiftedSubspacePolynomial.hpp"
#include "common/Algebra/details/Polynomials.hpp"
#include <algebraLib/UnivariatePolynomialGeneral.hpp>

namespace PCP_Project{
namespace Verifier{

using std::pair;
using std::max;
using std::map;
using Algebra::FieldElement;
using Algebra::zero;
using Infrastructure::POW2;

void RS_PCPP_result::init(const vector<FieldElement>& basis, const short degBound){
    results.resize(POW2(basis.size()));
    basis_ = basis;
    degBound_ = degBound;
}

bool RS_PCPP_result::verify()const{
    const Algebra::UnivariatePolynomialGeneral P_Poly(results,basis_,Algebra::zero());
    bool res = (Algebra::PolynomialDegree::integral_t(P_Poly.getDegree()) < degBound_);
    return res;
}

void ResultLocation::addAnswerPtr(FieldElement* ptr){
    answerLocations_.push_back(ptr);
}

void ResultLocation::answer(const FieldElement& res)const{
    for (const auto& e : answerLocations_){
        *e = res;
    }
}

namespace{

size_t bivariatePointToIndex(const size_t rowSize, const size_t rowIndex, const size_t columnIndex){
    return rowSize*rowIndex + columnIndex;
}

size_t queryUnivariatePolynomial_index(const size_t rowSize, const short cosetSizeLog, const size_t queryPointIndex){
    const size_t numCommonColumns = POW2(cosetSizeLog + Mu);
    const size_t rowIndex = queryPointIndex>>cosetSizeLog;
    const size_t columnIndex = (queryPointIndex < numCommonColumns? queryPointIndex : (numCommonColumns + queryPointIndex%POW2(cosetSizeLog)));
	
    return bivariatePointToIndex(rowSize, rowIndex, columnIndex);
}

pair<bool,size_t> getParentIndexOfPoint(const vector<FieldElement>& evaluationBasis, const size_t rowIdx, const size_t columnIdx){
    
    const short L0BasisSize = PCP_common::getL0Basis(evaluationBasis).size();
    const short L0PrimeBasisSize = PCP_common::getL0PrimeBasis(evaluationBasis).size();
    const short cosetIndex = columnIdx>>L0BasisSize;
    const short cosetsInL0Prime = POW2(L0PrimeBasisSize - L0BasisSize);
    
    //check if the side case of the first few segments that intersect the rectangle
    //(and those right after them)
    //In this case the original index is the same as the columnIdx
    if ((rowIdx < cosetsInL0Prime) && (cosetIndex == rowIdx)){
        return pair<bool,size_t>(true,columnIdx);
    }

    //if the column is in the rectangle, and we are still here,
    //this point is not a parent point
    if((cosetIndex < cosetsInL0Prime) || (rowIdx < cosetsInL0Prime)){
        return pair<bool,size_t>(false,0);
    }

    //Check for the cosets outside the rectangle
    if(cosetIndex == cosetsInL0Prime){
        return pair<bool,size_t>(true, rowIdx*POW2(L0BasisSize) + (columnIdx%POW2(L0BasisSize)));
    }
    
    //It is not a parent point
    return pair<bool,size_t>(false,0);
}

}

std::map<size_t,Algebra::FieldElement*> addRandomQuery(const short degBound_logCeil, const vector<FieldElement>& evaluationBasis, const FieldElement& affineShift, RS_PCPP_queriesTree& queries, RS_PCPP_result& results, const size_t depth, const bool isRoot){

    //
    //Global constants
    //
    const size_t MIN_DEPTH = 1;

    //
    //some constants
    //
    const vector<FieldElement> columnsBasis(PCP_common::getL1Basis(evaluationBasis));
    const vector<FieldElement> L0Basis = PCP_common::getL0Basis(evaluationBasis);
    const size_t numRows = POW2(columnsBasis.size());
    const size_t numColumns = POW2(PCP_common::getL0PrimeBasis(evaluationBasis).size());
    const short sizeOfRowBasis = PCP_common::dimOfLBeta(evaluationBasis.size());
    const size_t sizeOfRow = POW2(sizeOfRowBasis);

    vector<FieldElement> basisForColumnsProof(columnsBasis);
    FieldElement columnsAffineShift;
    {
        const Algebra::elementsSet_t rowsBasis(L0Basis.begin(), L0Basis.end());

        const Algebra::ShiftedSubspacePolynomial q(rowsBasis,affineShift);
        const FieldElement q_on_ZERO = q.eval(zero());
        columnsAffineShift = q.eval(affineShift);
        for(short i=0; i<columnsBasis.size(); i++){
            basisForColumnsProof[i] = q.eval(columnsBasis[i] + q_on_ZERO);
        }
    }

    //
    //Decide if to query row or column - if BIASED_RECURSION flag is off, simply choose uniformly
    //
 //   typedef enum{ROW,COLUMN} queryType_t; -now defined in SoundnessParameters class
#ifdef BIASED_RECURSION
	const queryType_t queryType = SoundnessParameters::chooseRowOrCol(depth);
#else
	const queryType_t queryType = static_cast<queryType_t>(rand() % 2);
#endif
    //
    //Decide index of ROW/COLUMN to query
    //
    size_t queryIdx = rand();
    if(queryType == ROW){
        queryIdx %= numRows;
    }
    else{
        queryIdx %= numColumns;
    }

    const bool passToParent = !isRoot;
    map<size_t,FieldElement*> queriesToParent;
    map<size_t,FieldElement*> queriesFromChild;

    if(queryType == ROW){
        const vector<FieldElement> currRowBasis(PCP_common::getLBeta(evaluationBasis, queryIdx));
        const short next_degBound_log = L0Basis.size();
        if(depth == MIN_DEPTH){
            results.init(currRowBasis,POW2(next_degBound_log));
            const size_t queriesIndexOffset = sizeOfRow*queryIdx;
            
            for(size_t i=0; i < sizeOfRow; i++){
                const auto isParentPoint = getParentIndexOfPoint(evaluationBasis,queryIdx,i);
                if(passToParent && isParentPoint.first){
                    queriesToParent[isParentPoint.second] = &(results.results[i]);
                }
                else{
                    const size_t currIndex = queriesIndexOffset+i;
                    auto query_iter = queries.localQueries.insert(pair<size_t,ResultLocation>(currIndex,ResultLocation())).first;
                    query_iter->second.addAnswerPtr(&(results.results[i]));
                }
            }
        }
        else{
            const vector<FieldElement> currRowBasis_shuffled(PCP_common::changeRowBasisBeforeRecursion(currRowBasis));
            auto rowQuery = queries.rowsToExtract.find(queryIdx);
            if(rowQuery == queries.rowsToExtract.end()){
                rowQuery = queries.rowsToExtract.insert(pair<size_t,RS_PCPP_queriesTree>(queryIdx,RS_PCPP_queriesTree())).first;
            }
            queriesFromChild = addRandomQuery(next_degBound_log,currRowBasis_shuffled,affineShift,rowQuery->second,results,depth-1,false);
        }
    }
    else{ // hendle column query
        const short curr_degree_log=max(0,short(basisForColumnsProof.size())-short(evaluationBasis.size()-degBound_logCeil));
        
        if(depth == MIN_DEPTH){
            results.init(basisForColumnsProof,POW2(curr_degree_log));
            for(size_t i = 0; i < numRows; i++){
                const auto isParentPoint = getParentIndexOfPoint(evaluationBasis,i,queryIdx);
                if( passToParent && isParentPoint.first){
                    queriesToParent[isParentPoint.second] = &(results.results[i]);
                }
                else{
                    const size_t currIndex = i*sizeOfRow + queryIdx;
                    auto query_iter = queries.localQueries.insert(pair<size_t,ResultLocation>(currIndex,ResultLocation())).first;
                    query_iter->second.addAnswerPtr(&(results.results[i]));
                }
            }
        }
        else{
            auto columnQuery = queries.columnsToExtract.find(queryIdx);
            if(columnQuery == queries.columnsToExtract.end()){
                columnQuery = queries.columnsToExtract.insert(pair<size_t,RS_PCPP_queriesTree>(queryIdx,RS_PCPP_queriesTree())).first;
            }
            queriesFromChild = addRandomQuery(curr_degree_log,basisForColumnsProof,columnsAffineShift,columnQuery->second,results,depth-1,false);
        }
    }

    //
    //Handle queries from child
    //
    for(const auto& qChild : queriesFromChild){
        pair<bool,size_t> isParentPoint;
        if(queryType == ROW){
            //shuffle back to original row index
            const size_t origRowIndex = PCP_common::shuffledRowIndexToOriginal(qChild.first,sizeOfRow);
            isParentPoint = getParentIndexOfPoint(evaluationBasis,queryIdx,origRowIndex);
            if(passToParent && isParentPoint.first){
                queriesToParent[isParentPoint.second] = qChild.second;
            }
            else{
                const size_t currIndex = queryIdx*sizeOfRow + origRowIndex;
                auto query_iter = queries.localQueries.insert(pair<size_t,ResultLocation>(currIndex,ResultLocation())).first;
                query_iter->second.addAnswerPtr(qChild.second);
            }
        }
        else{ //queryType == COLUMN
            isParentPoint = getParentIndexOfPoint(evaluationBasis,qChild.first,queryIdx);
            if(passToParent && isParentPoint.first){
                queriesToParent[isParentPoint.second] = qChild.second;
            }
            else{
                const size_t currIndex = queryIdx + sizeOfRow*qChild.first;
                auto query_iter = queries.localQueries.insert(pair<size_t,ResultLocation>(currIndex,ResultLocation())).first;
                query_iter->second.addAnswerPtr(qChild.second);
            }
        }
    }

    //
    //Pass queries to parent
    //
    return queriesToParent;

}

void ACSP_CONSISTENCY_result::init(const ACSPInstance& instance, const FieldElement& consistencyPoint){
    instance_ = &instance;
    boundaryVanishingVals_.resize(instance.neighborPolys().size());
    boundaryPolyVals_.resize(instance.neighborPolys().size());
    boundaryPoly_res.resize(instance.neighborPolys().size());
    consistencyPoint_ = consistencyPoint;
    
    //get the neighbors values for the consistency test
    {
        Algebra::details::UnivariatePolynomial Ax;
        Algebra::details::UnivariatePolynomial Z_X;
        {
            std::vector<FieldElement> S_x;
            Algebra::details::UniPolynomialEvaluation f_eval;
            for (const auto& p : instance_->boundaryConstraints()) {
                f_eval.addPoint(p.first, p.second);
                S_x.push_back(p.first);
            }
            Ax.interpolation(f_eval, S_x);
            Z_X.findNonSparseVanishingPoly(S_x);
        }

        const auto& neighbors = instance_->neighborPolys();
        for (int i=0; i < neighbors.size(); i++) {
            const FieldElement currNeighborVal = neighbors[i]->eval(consistencyPoint_);
            boundaryVanishingVals_[i] = Z_X.queryAtPoint(currNeighborVal);
            boundaryPolyVals_[i] = Ax.queryAtPoint(currNeighborVal);
        }
    }
    
    //get the vanishing space polynomial value for consistency test (aka Z_H)
    {
        vanishingSpacePolyVal_ = instance_->vanishingSet().vanishingPoly()->eval(consistencyPoint_);
    }
}

bool ACSP_CONSISTENCY_result::verify()const{
    vector<FieldElement> inputValuesVector;
    inputValuesVector.push_back(consistencyPoint_);

    const auto& neighbors = instance_->neighborPolys();
    for (int i=0; i<neighbors.size(); i++) {
        const FieldElement inputValue(boundaryPoly_res[i]);
        inputValuesVector.push_back(inputValue*boundaryVanishingVals_[i] + boundaryPolyVals_[i]);
    }
    const auto& acspPoly = instance_->constraintPoly();
    const FieldElement res1 = compositionPoly_res * vanishingSpacePolyVal_;
    const FieldElement res2 = acspPoly.eval(inputValuesVector);

    return res1 == res2;
}

ACSP_CONSISTENCY_query::ACSP_CONSISTENCY_query(const Algebra::FieldElement& consistencyPoint_, ACSP_CONSISTENCY_result& resLocation):
    consistencyPoint(consistencyPoint_){
        compositionPoly_res.addAnswerPtr(&resLocation.compositionPoly_res);
        
        boundaryPoly_res.resize(resLocation.boundaryPoly_res.size());
        for(size_t i=0; i< resLocation.boundaryPoly_res.size(); i++){
            boundaryPoly_res[i].addAnswerPtr(&resLocation.boundaryPoly_res[i]);
        }
    }

} //namespace Verifier
} // namespace PCP_Project
