#include "algebraLib/FiniteField.hpp"
#include "algebraLib/ErrorHandling.hpp"
#include <NTL/GF2XFactoring.h>
#include <NTL/GF2E.h>



namespace Algebra {

FiniteField::FiniteField():modulus_(::NTL::BuildIrred_GF2X(1)){};
FiniteField::FiniteField(const size_t extensionDegree)
:modulus_(::NTL::BuildIrred_GF2X(extensionDegree)) {}
	
FiniteField::FiniteField(const FiniteField& other):modulus_(other.modulus_){}

size_t FiniteField::degree()const{
	return ::NTL::deg(modulus_);
}

elementsSet_t FiniteField::getBasis()const {
	return getStandartBasis(degree());
}
	
void FiniteField::setContext()const{
	//verify it is the constant irreducible Matans library (aka FFF) works with
    ALGEBRALIB_ASSERT( NTL::deg(modulus_) == 64 ,"The context field irreducible does not match FFF ireducible");
    
    for(int i=0; i<65 ; i++){
        switch(i){
            case 0:
            case 1:
            case 3:
            case 4:
            case 64:
                ALGEBRALIB_ASSERT( NTL::coeff(modulus_,i) == 1 ,"The context field irreducible does not match FFF ireducible");
                break;
            default:
                ALGEBRALIB_ASSERT( NTL::coeff(modulus_,i) == 0 ,"The context field irreducible does not match FFF ireducible");
                break;
        }
    }

    NTL::GF2E::init(modulus_);
}

void FiniteField::releaseContextField(){
    delete NTL::GF2EInfo;
}

} // namespace Algebra 
