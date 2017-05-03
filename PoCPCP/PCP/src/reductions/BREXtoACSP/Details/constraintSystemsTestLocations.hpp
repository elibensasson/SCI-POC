#ifndef CONSTRAINTSYSTEMS_TEST_LOCATIONS_HPP__
#define CONSTRAINTSYSTEMS_TEST_LOCATIONS_HPP__

#include "common.hpp"
#include <map>

namespace PCP_Project{
namespace BREXtoACSP{

class CS_testLocations{
public:
    typedef size_t polynomialIndicator_t;
    CS_testLocations(const common& commonDef);
    size_t indexOfConstraint_Assignment(polynomialIndicator_t poly)const;
    size_t indexOfConstraint_Permuation(polynomialIndicator_t poly)const;

private:
    std::map<polynomialIndicator_t, size_t> indexesAssignment_;
    std::map<polynomialIndicator_t, size_t> indexesPermutation_;
};

} //namespace BREXtoACSP
} //namespace PCP_Project


#endif // #ifndef CONSTRAINTSYSTEMS_TEST_LOCATIONS_HPP__
