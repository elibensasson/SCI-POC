#ifndef BREXTOACSP_CONSTRAINTS_HPP__
#define BREXTOACSP_CONSTRAINTS_HPP__

#include "neighborsConstructor.hpp"
#include "spaces.hpp"
#include "common.hpp"
#include "instanceMappings.hpp"

namespace PCP_Project{
namespace BREXtoACSP{

class constraints{
public:
    static std::unique_ptr<Algebra::PolynomialInterface> getConstraintsPoly(
const BREXInstance& instance,
const AcspNeighbors& neighbors,
const common& commonDef,
const instanceMappings& instanceMapping,
const spaces& spacesGenerator,
const CS_testLocations& testLocations);
};

} //namespace BREXtoACSP
} //namespace PCP_Project

#endif //#ifndef BREXTOACSP_CONSTRAINTS_HPP__
