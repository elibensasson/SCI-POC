#include "QueriesResFetcher.hpp"
#include "PCP/PCP_common.hpp"
#include "common/Algebra/ShiftedSubspacePolynomial.hpp"

namespace PCP_Project{
namespace QueriesResFetcher{

using Infrastructure::POW2;
using Algebra::FieldElement;
using Algebra::zero;
using Algebra::elementsSet_t;
using Algebra::ShiftedSubspacePolynomial;
using std::vector;
using std::pair;

namespace{
omp_lock_t path_manipulation_lock;
}

//Michael: can you explain this method?
FieldElement readValAndKeepPath(const BivariateExpansionProof& bivariate, const size_t index, paths_t& paths){
    const size_t blockIndex = bivariate.merkleTree().getBlockIndex(index);
    const short offsetInDualBlock = bivariate.merkleTree().getOffsetInDualBlock(index);
    FieldElement res;

omp_set_lock(&path_manipulation_lock);
{
    pathIndicator_t currIndicator(bivariate.merkleTree().getRoot() , blockIndex & ~1UL);
    const auto it = paths.find(currIndicator);
    if(it != paths.end()){
        res = it->second.data[offsetInDualBlock];
    }
    else{
        //if got here, need to actually do the query
        queryCommitmentPath& currPath = paths[currIndicator];
        currPath.blockIndex = blockIndex;
        currPath.root = bivariate.merkleTree().getRoot();
        currPath.path = bivariate.merkleTree().getPathToBlock(blockIndex);
        for(size_t i=0; i< bivariate.merkleTree().getDualBlockSize(); i++){
            currPath.data.push_back(bivariate.readValue(index - offsetInDualBlock + i));
        }

        //return result
        res =  currPath.data[offsetInDualBlock];
    }
}
omp_unset_lock(&path_manipulation_lock);
    return res;
}


//question to Michael: Please explain method
FieldElement readValAndKeepPath(const vector<FieldElement>& row,const size_t rowId, const size_t columnId, const MerkleTreeInterface& currRowTree, const MerkleTreeInterface& finalTree, paths_t& paths){
    const size_t globalIndex = (rowId * row.size()) + columnId;
    const size_t blockIndex = MerkleTree::getBlockIndex(globalIndex);
    const short offsetInDualBlock = MerkleTree::getOffsetInDualBlock(globalIndex);
    FieldElement res;

omp_set_lock(&path_manipulation_lock);
{
    pathIndicator_t currIndicator(currRowTree.getRoot() , blockIndex & ~1UL);
    const auto it = paths.find(currIndicator);
    if(it != paths.end()){
        res = it->second.data[offsetInDualBlock];
    }
    else{
        //if got here, need to actually do the query
        queryCommitmentPath& currPath = paths[currIndicator];
        currPath.blockIndex = blockIndex;
        currPath.root = finalTree.getRoot();

        currPath.path =  currRowTree.getPathToBlock(MerkleTree::getBlockIndex(columnId));
        const auto secondPathPart = finalTree.getPathToBlock(2*rowId);
        currPath.path.insert(currPath.path.end(), secondPathPart.begin(), secondPathPart.end());

        for(size_t i=0; i< MerkleTree::getDualBlockSize(); i++){
            currPath.data.push_back(row[columnId - offsetInDualBlock + i]);
        }

        //return result
        res = currPath.data[offsetInDualBlock];
    }
}
omp_unset_lock(&path_manipulation_lock);
    return res;
}


//takes an index in space L, on which we want to evaluate the univariate, and finds matching index (in one dim format) in the bivariate table
size_t queryUnivariatePolynomial_index(const BivariateExpansionProof& bivariate, const short cosetSizeLog, const size_t queryPointIndex){
    const size_t numCommonColumns = POW2(cosetSizeLog + Mu);
    const size_t rowIndex = queryPointIndex>>cosetSizeLog;
    const size_t columnIndex = (queryPointIndex < numCommonColumns? queryPointIndex : (numCommonColumns + queryPointIndex%POW2(cosetSizeLog)));
	
    return bivariate.getValIndex(rowIndex,columnIndex);
}

void queryPCPP_recursive(const BivariateExpansionProof& proof,const Prover::ProverAlgorithm& prover, const vector<FieldElement>& evaluationBasis, const FieldElement& affineShift,
    const Verifier::RS_PCPP_queriesTree& queries,
    paths_t& all_paths,std::function<void(const BivariateExpansionProof& proof)> collect_proof_spec){
    
    //query current bivariate
    for (const auto& q : queries.localQueries) {
        const FieldElement result = readValAndKeepPath(proof, q.first, all_paths);
        q.second.answer(result);
    }
    
    //query from proofs of rows and columns - should go in here only if recursion Depth is larger than 2
    {
        //
        //some constants
        //
        //const vector<FieldElement> columnsBasis(evaluationBasis.begin()+HALF_K(evaluationBasis.size())-Gamma+1,evaluationBasis.end());
		const vector<FieldElement> columnsBasis(PCP_common::getL1Basis(evaluationBasis));
        vector<FieldElement> basisForColumnsProof(columnsBasis);
        FieldElement columnsAffineShift;
#pragma omp critical
        {
            //const elementsSet_t rowsBasis(evaluationBasis.begin(),evaluationBasis.begin()+HALF_K(evaluationBasis.size())-Gamma+1);
			vector<FieldElement> L0Basis = PCP_common::getL0Basis(evaluationBasis);
			const elementsSet_t rowsBasis(L0Basis.begin(), L0Basis.end());

            const ShiftedSubspacePolynomial q(rowsBasis,affineShift);
            const FieldElement q_on_ZERO = q.eval(zero());
            columnsAffineShift = q.eval(affineShift);
            for(short i=0; i<columnsBasis.size(); i++){
                basisForColumnsProof[i] = q.eval(columnsBasis[i] + q_on_ZERO);
            }
        }
        const size_t numRows = POW2(columnsBasis.size());
        //const size_t sizeOfRowBasis = 2 + HALF_K(evaluationBasis.size()) - Gamma + Mu;
		const size_t sizeOfRowBasis = PCP_common::dimOfLBeta(evaluationBasis.size());

        
        for(const auto& q : queries.rowsToExtract){
            const size_t rowIndex = q.first;

            //Here we are creating a subproof for a row. Basis needs to be shuffled to guarantee each row in subproof has intersection with original row
            vector<FieldElement> currRowBasis(PCP_common::changeRowBasisBeforeRecursion(PCP_common::getLBeta(evaluationBasis, rowIndex)));
            const auto currSubProof = prover.proofRow_lowDegree(proof,rowIndex,currRowBasis,affineShift);
            collect_proof_spec(*currSubProof);
            queryPCPP_recursive(*currSubProof, prover, currRowBasis, affineShift, q.second, all_paths, collect_proof_spec);
        }

        for(const auto& q : queries.columnsToExtract){
            const size_t columnIndex = q.first;
            const auto currSubProof = prover.proofColumn_lowDegree(proof,columnIndex,basisForColumnsProof,columnsAffineShift);
            collect_proof_spec(*currSubProof);
            queryPCPP_recursive(*currSubProof, prover,basisForColumnsProof,columnsAffineShift,q.second,all_paths,collect_proof_spec);
        }
    }
}

typedef struct columnQueryInfo_new{
    size_t columnId;
    const Verifier::RS_PCPP_queriesTree& subQueries;
    vector<FieldElement> columnVals;
    columnQueryInfo_new(const size_t colId ,const Verifier::RS_PCPP_queriesTree& subQueries_) : columnId(colId) , subQueries(subQueries_){};
};

typedef struct queryToFill{
    size_t rowId;
    size_t colId;
    std::function<void(const FieldElement&)> setVal;
    queryToFill(size_t rowId_, size_t colId_, std::function<void(const FieldElement&)> setVal_) : rowId(rowId_), colId(colId_), setVal(setVal_){};
};

bool cmp_1(const queryToFill &a,const queryToFill &b){
    if (a.rowId < b.rowId) return true;
    if (a.rowId > b.rowId) return false;

    //now we are sure a.x == a.b
    return a.colId < b.colId;
}

//this method generates the first recursion layer, which always has to be fully generated and comitted to.
void queryPCPP_recursive_firstLayer(const BivariateExpansionProof& proof,const Prover::ProverAlgorithm& prover, const vector<FieldElement>& evaluationBasis, const FieldElement& affineShift,
    const Verifier::RS_PCPP_queriesTree& queries,
    paths_t& all_paths,vector<vector<queryToFill>>& dataQueries,std::function<void(const BivariateExpansionProof& proof)> collect_proof_spec){
    
    //
    //some constants
    //
//    const vector<FieldElement> columnsBasis(evaluationBasis.begin()+HALF_K(evaluationBasis.size())-Gamma+1,evaluationBasis.end());
	const vector<FieldElement> columnsBasis(PCP_common::getL1Basis(evaluationBasis));
	vector<FieldElement> basisForColumnsProof(columnsBasis);
    FieldElement columnsAffineShift;
    {
        //const elementsSet_t rowsBasis(evaluationBasis.begin(),evaluationBasis.begin()+HALF_K(evaluationBasis.size())-Gamma+1);
		vector<FieldElement> L0Basis = PCP_common::getL0Basis(evaluationBasis);
		const elementsSet_t rowsBasis(L0Basis.begin(),L0Basis.end());
		const ShiftedSubspacePolynomial q(rowsBasis,affineShift);
        const FieldElement q_on_ZERO = q.eval(zero());
        columnsAffineShift = q.eval(affineShift);
        for(short i=0; i<columnsBasis.size(); i++){
            basisForColumnsProof[i] = q.eval(columnsBasis[i] + q_on_ZERO);
        }
    }
    const size_t numRows = POW2(columnsBasis.size());
	const size_t sizeOfRowBasis = PCP_common::dimOfLBeta(evaluationBasis.size());
	
	//
    //convert queries to comfortable structs
    //
    vector<columnQueryInfo_new> columnsToQuery;

    {
        for (const auto& localQ : queries.localQueries) {
            auto setValFunc = [=](const FieldElement& res){ localQ.second.answer(res); };
            const size_t rowId = localQ.first >> sizeOfRowBasis;
            const size_t columnId = localQ.first & (POW2(sizeOfRowBasis)-1);
            dataQueries[rowId].push_back(queryToFill(rowId,columnId,setValFunc));
        }
        
        for(const auto& col : queries.columnsToExtract){
            columnsToQuery.push_back(columnQueryInfo_new(col.first,col.second));
        }

        for(auto& c : columnsToQuery){
            c.columnVals.resize(numRows);
        }
    }

    //
    //get the final Merkle tree
    //
    const MerkleTreeInterface& finalTree = proof.merkleTree();
    
    //
    //Answer queries row by row - 
    //
    {
#pragma omp parallel for num_threads(SCIPR_NUM_THREADS)
        for (long long rowIdx=0; rowIdx< numRows; rowIdx++){
            
            //get current row
            const auto currRow = proof.getRow(rowIdx);

            //read relevant data queries
            if(!dataQueries[rowIdx].empty()){
                    const MerkleTree currRowTree(&currRow[0],sizeOfRowBasis);
                    
                    for(auto& currentQueryInfo : dataQueries[rowIdx]){
                        currentQueryInfo.setVal(readValAndKeepPath(currRow, rowIdx, currentQueryInfo.colId, currRowTree, finalTree, all_paths));//THIS is the line where actual query results are stored, though hard to track where with all these lambda functions
                    }
            }

            //prove row if needed. If running the PCPP at recursion depth 1, will never go in here
            {
                const auto currRowQuery = queries.rowsToExtract.find(rowIdx);
                if(currRowQuery != queries.rowsToExtract.end()){
					//Here we are creating a subproof for a row. Basis needs to be shuffled to guarantee each row in subproof has intersection with original poly
					vector<FieldElement> currRowBasis(PCP_common::changeRowBasisBeforeRecursion(PCP_common::getLBeta(evaluationBasis, rowIdx)));
					vector<FieldElement> shuffledRow(PCP_common::changeRowValuesBeforeRecursion(currRow));
                    const auto currSubProof = prover.expandUnivariateToBivariate(&shuffledRow[0],currRowBasis,affineShift);
                    collect_proof_spec(*currSubProof);
					queryPCPP_recursive(*currSubProof,prover, currRowBasis, affineShift, currRowQuery->second , all_paths, collect_proof_spec);

                }
            }

            //collect relevant columns data In *this* *first* recursion level
            for(auto& col : columnsToQuery){
                col.columnVals[rowIdx] = currRow[col.columnId];
            }
        }
    }

    //
    //Answer column queries
    //
    {
#pragma omp parallel for num_threads(SCIPR_NUM_THREADS)
        for(long long i=0; i< columnsToQuery.size(); i++){
            const auto& col = columnsToQuery[i];
            const auto currSubProof = prover.expandUnivariateToBivariate(&col.columnVals[0],basisForColumnsProof,columnsAffineShift);
            collect_proof_spec(*currSubProof);
            queryPCPP_recursive(*currSubProof,prover,basisForColumnsProof,columnsAffineShift,col.subQueries,all_paths,collect_proof_spec);
        }
    }

}

paths_t resultsFillingAlgorithm::fillResults(const Prover::proof_t& proof,const Prover::ProverAlgorithm& prover,
                        const ACSPInstance& instance, const queries_t& queries,
                        std::function<void(const BivariateExpansionProof& proof)> collect_proof_spec)const{
   
    //lock intialization
    omp_init_lock(&path_manipulation_lock);
    
    //get some constants
    const Algebra::details::Basis basisPCPP(PCP_common::basisForPCPP(instance));
    const size_t numRows = POW2(PCP_common::getL1Basis(basisPCPP.asVector()).size());
    
    //queries lists
    vector<vector<queryToFill>> dataQueries_composition(numRows);
    vector<vector<queryToFill>> dataQueries_boundary(numRows);
    
    //Merkle paths
    paths_t all_paths;

    //
    //answer consistency check
    //
    {
		/* Iterating over the consistency queries of all the verifiers - which are *univariate* inputs, and preparing the appropriate *bivariate* inputs
		 * according to index calcluations in queryUnivariatePolynomialIndex, storing queries in dataQueries_composition and dataQueriesBoundary
		 * afterwards, using queryPCPP method, we fill in answers in the results variable.
		 */
		const short cosetSizeLog = PCP_common::L0Dim(basisPCPP.asVector().size());

		for(const auto& q : queries.ACSP_consistency){
            const FieldElement alpha = q.consistencyPoint;
            const size_t alpha_index = Algebra::mapFieldElementToInteger(0,64,alpha); //getSpaceIndexOfElement(basisPCPP.asVector(), zero(),alpha);
            const size_t res_index = queryUnivariatePolynomial_index(*proof.compositionRS_proof,cosetSizeLog,alpha_index);

            auto setValFunc = [=,&q](const FieldElement& res){ 
                q.compositionPoly_res.answer(res); 
            };

            {
                const size_t rowId = proof.compositionRS_proof->getIndexRow(res_index);
                dataQueries_composition[rowId].push_back(queryToFill(rowId,proof.compositionRS_proof->getIndexColumn(res_index),setValFunc));
            }


            const auto& neighbours = instance.neighborPolys();
            size_t nsize = neighbours.size();
            for (int affine_num = 0; affine_num < nsize; affine_num++){
                const FieldElement currNeighborVal = neighbours[affine_num]->eval(alpha);
                const size_t neighborValIndex = Algebra::mapFieldElementToInteger(0,64,currNeighborVal); //getSpaceIndexOfElement(basisPCPP.asVector(),zero(),currNeighborVal);
                const size_t proofRes_index = queryUnivariatePolynomial_index(*proof.boundaryRS_proof,cosetSizeLog,neighborValIndex);

                auto setValFunc = [=,&q](const FieldElement& res){ 
                    q.boundaryPoly_res[affine_num].answer(res); 
                };
                
                {
                    const size_t rowId = proof.boundaryRS_proof->getIndexRow(proofRes_index);
                    dataQueries_boundary[rowId].push_back(queryToFill(rowId,proof.boundaryRS_proof->getIndexColumn(proofRes_index),setValFunc));
                }
            }
        }
    }

    //
    //Fetched composition RS queries
    //
    queryPCPP_recursive_firstLayer(*proof.compositionRS_proof,prover,basisPCPP.asVector(),zero(), queries.RS_composition ,all_paths, dataQueries_composition,collect_proof_spec);
	
    //
    //Fetch boundary RS queries
    //
    queryPCPP_recursive_firstLayer(*proof.boundaryRS_proof,prover,basisPCPP.asVector(),zero(), queries.RS_boundary,all_paths, dataQueries_boundary,collect_proof_spec);
    
    
    omp_destroy_lock(&path_manipulation_lock);

    //return paths
    return all_paths;
}

paths_t resultsFillingAlgorithm_dummy::fillResults(const Prover::proof_t& proof,const Prover::ProverAlgorithm& prover,
                        const ACSPInstance& instance, const queries_t& queries,
                        std::function<void(const BivariateExpansionProof& proof)> collect_proof_spec)const{
   
    //lock intialization
    omp_init_lock(&path_manipulation_lock);
    
    //get some constants
    const Algebra::details::Basis basisPCPP(PCP_common::basisForPCPP(instance));
    paths_t all_paths;

    //
    //answer consistency check
    //
    {
		/* Iterating over the consistency queries of all the verifiers - which are *univariate* inputs, and preparing the appropriate *bivariate* inputs
		 * according to index calcluations in queryUnivariatePolynomialIndex, storing queries in dataQueries_composition and dataQueriesBoundary
		 * afterwards, using queryPCPP method, we fill in answers in the results variable.
		 */
		const short cosetSizeLog = PCP_common::L0Dim(basisPCPP.asVector().size());

		for(const auto& q : queries.ACSP_consistency){
            const FieldElement alpha = q.consistencyPoint;
            const size_t alpha_index = Algebra::mapFieldElementToInteger(0,64,alpha); //getSpaceIndexOfElement(basisPCPP.asVector(), zero(),alpha);
            const size_t res_index = queryUnivariatePolynomial_index(*proof.compositionRS_proof,cosetSizeLog,alpha_index);
            const FieldElement res = proof.compositionRS_proof->readValue(res_index);
            q.compositionPoly_res.answer(res); 
            
            const auto& neighbours = instance.neighborPolys();
            size_t nsize = neighbours.size();
            for (int affine_num = 0; affine_num < nsize; affine_num++){
                const FieldElement currNeighborVal = neighbours[affine_num]->eval(alpha);
                const size_t neighborValIndex = Algebra::mapFieldElementToInteger(0,64,currNeighborVal); //getSpaceIndexOfElement(basisPCPP.asVector(),zero(),currNeighborVal);
                const size_t proofRes_index = queryUnivariatePolynomial_index(*proof.boundaryRS_proof,cosetSizeLog,neighborValIndex);
                const FieldElement res = proof.boundaryRS_proof->readValue(proofRes_index);
                q.boundaryPoly_res[affine_num].answer(res); 
            }
        }
    }

    //
    //Fetched composition RS queries
    //
    queryPCPP_recursive(*proof.compositionRS_proof,prover,basisPCPP.asVector(),zero(), queries.RS_composition ,all_paths, collect_proof_spec);
	
    //
    //Fetch boundary RS queries
    //
    queryPCPP_recursive(*proof.boundaryRS_proof,prover,basisPCPP.asVector(),zero(), queries.RS_boundary,all_paths, collect_proof_spec);
    
    
    omp_destroy_lock(&path_manipulation_lock);

    //return paths
    return all_paths;
}

} //namespace QueriesResFetcher
} //namespace PCP_Project
