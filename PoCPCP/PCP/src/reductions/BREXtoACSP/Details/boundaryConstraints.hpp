#ifndef BREXTOACSP_INPUTINDICATOR_HPP__
#define BREXTOACSP_INPUTINDICATOR_HPP__

namespace PCP_Project{
namespace BREXtoACSP{

class boundaryMapping : private instanceMappings{
public:
    boundaryMapping(const common& commonInfo):
        instanceMappings(commonInfo){};

    Algebra::FieldElement mapPoint(const size_t rowId, const size_t varId)const{
        return mapRowId(rowId) + mapVariable(varId);
    }

private:
    //Methods
    Algebra::FieldElement mapRowId(const size_t rowId)const{
        return map_x_power_modulu_poly(rowId, commonInfo_.rowsModulus());
    }
};

ACSPInstance::boundaryConstraints_t reduceBoundaryConstraints(const BREXInstance::boundaryConstraints_t& brexBoundary, const common& commonDef){
    boundaryMapping mapping(commonDef);
    ACSPInstance::boundaryConstraints_t boundaryConstraints;
    for(const auto& p : brexBoundary){
        const auto currPoint = p.first;
        const size_t rowId = currPoint.first;
        const size_t varId = currPoint.second;
        const auto x = mapping.mapPoint(rowId,varId);
        boundaryConstraints[x] = p.second;
    }

    return boundaryConstraints;
}

} //namespace BREXtoACSP
} //namespace PCP_Project

#endif //#ifndef BREXTOACSP_INPUTINDICATOR_HPP__
