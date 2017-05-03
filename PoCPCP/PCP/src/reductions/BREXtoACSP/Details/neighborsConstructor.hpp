#ifndef NEIGHBORS_CONSTRUCTOR_HPP__
#define NEIGHBORS_CONSTRUCTOR_HPP__

#include "common.hpp"
#include "instanceMappings.hpp"
#include "constraintSystemsTestLocations.hpp"
#include <algebraLib/LinearPolynomial.hpp>

#include <array>
#include <vector>
#include <map>
#include <memory>

namespace PCP_Project{
namespace BREXtoACSP{

class AcspNeighbors{
public:
    typedef size_t polynomialIndicator_t;

    AcspNeighbors(const BREXInstance& brexInstance,const common& commonDef, const instanceMappings& instanceMapping, const CS_testLocations& testLocations);
    std::vector<std::unique_ptr<const Algebra::UnivariatePolynomialInterface> > getNeighborPolynomials()const;
    size_t polynomialsNum()const;

    size_t locationOfId()const;
    size_t locationOfTwinLayer()const;
    size_t locationOfDeBruijn(const short dbNeighborId, const short affineCosetId ,const size_t layerId)const;
    size_t locationOfRoutingBit(const size_t layerId)const;
    size_t locationOfPermCS(polynomialIndicator_t poly, const size_t varId)const;
    bool   existsPermCS(polynomialIndicator_t poly, const size_t varId)const;
    size_t locationOfAssignmentCS(polynomialIndicator_t poly, const size_t varId, const size_t neighborVersion)const;
    bool   existsAssignmentCS(polynomialIndicator_t poly, const size_t varId)const;

    
private:

    struct NeighborLocation_stc{
        bool exists;
        size_t index;

        NeighborLocation_stc(): exists(false), index(0){};
        NeighborLocation_stc(size_t index_): exists(true), index(index_){};
    };

    //reduction common
    const common& commonDef_;
    const instanceMappings instanceMapping_;
    const CS_testLocations testLocations_;
    
    //usage map of constraint systems polynomials
    const std::map< polynomialIndicator_t, std::vector <bool> > permutationCS_usageVector_;
    const std::map< polynomialIndicator_t, std::vector <bool> > assignmentCS_usageVector_;

    //
    //neighbor polynomials
    //

    std::vector<Algebra::LinearPolynomial> neighbors_;


    /******************************************
    * Neighbors locations
    ******************************************/

    //
    //locations for common neighbors
    //
    
    // ID polynomial (aka x)
    NeighborLocation_stc locationOfId_;

    // Permutation constraint system neighbors
    std::map<polynomialIndicator_t , std::vector<NeighborLocation_stc>> locationOfPermCS_;
    
    // Assignment constraint system neighbors
    std::map<polynomialIndicator_t , std::vector<std::array<NeighborLocation_stc,2>>> locationOfAssignmentCS_;
    
    //
    //locations for routing network neighbors
    //

    NeighborLocation_stc locationOfTwinLayer_;

    //accessed using the index [layer ID][DeBruijn neighbor][coset ID]
    std::vector<NeighborLocation_stc> locationOfDeBruijn_[2][4];
    
    //accessed using the layer ID
    std::vector<NeighborLocation_stc> locationOfRoutingBit_;

    /********************************************
    * private methods
    *********************************************/

    struct NeighborLocation_stc addNeighbor(const Algebra::LinearPolynomial& neighbor);
    static size_t retLocation(const struct NeighborLocation_stc& loc);
    void initRoutingNetworkNeighbors();
    Algebra::FieldElement getGenerator()const;
    static std::map<polynomialIndicator_t , std::vector<bool>> getCsUsageVector(const std::vector<std::unique_ptr<Algebra::PolynomialInterface>>& cs, size_t varsAmount);

    //neighbors constructors
    static Algebra::LinearPolynomial constructIdNeighbor();
    Algebra::LinearPolynomial constructTwinLayerNeighbor()const;
    void constructDeBruijn();
    void constructRoutingBitAccess();
    void constructPermCS();
    void constructAssignmentCS();
    
    //neighbors constructor helper functions
    Algebra::LinearPolynomial applyShiftedLinearOperation(const size_t layer_id, const Algebra::LinearPolynomial operation)const;
    Algebra::LinearPolynomial constructLinearDeBruijn_N0()const;
    Algebra::LinearPolynomial constructLinearDeBruijn_N1()const;
    Algebra::FieldElement DeBruijn_fixRowIdCarry()const;
    Algebra::FieldElement DeBruijn_fixColumnIdCarry()const;
    Algebra::FieldElement DeBruijn_getFixByCoset(const short cosetId)const;
    Algebra::LinearPolynomial constructPermCS(const size_t nonPermElemId, const size_t varId)const;
    Algebra::LinearPolynomial constructAssignmentCS(const size_t nonPermElemId, const size_t varId, const bool withCarry)const;
    Algebra::LinearPolynomial moveFromPointToVarId(const Algebra::FieldElement src, const size_t varId)const;
};


} //namespace BREXtoACSP
} //namespace PCP_Project


#endif //#ifdef NEIGHBORS_CONSTRUCTOR_HPP__
