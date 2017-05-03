/**
 *       @file  ACSPInstance.hpp
 *      @brief  Header file for ACSP Partial Instance
 *
 *	Contains ACSP Partiial instance definitions
 *
 *     @author  Michael Riabzev, RiabzevMichael@gmail.com
 * =====================================================================================
 */
#pragma once // This is not needed, here just to get rid of an annoying VS2010 warning.
             // Thanks again, Microsoft.

#ifndef  __ACSP_Instance_HPP
#define  __ACSP_Instance_HPP

#include "common/langCommon/Sequence.hpp"
#include <algebraLib/FiniteField.hpp>
#include <algebraLib/FieldElement.hpp>
#include <algebraLib/UnivariatePolynomialGeneral.hpp>
#include "common/Algebra/FiniteSetInterface.hpp"
#include "common/Algebra/details/FiniteFields.hpp" //ARIEL
#include "common/CXX11_macros.hpp"
#include <vector>
#include <memory>
#include <cassert>
#include <functional>

namespace PCP_Project {

/**
 * @class ACSPInstance
 * @brief class for ACSP Instance
 *
 * An ACSP partial instance is a tuple \f$(\mathbb{F},H,\vec{N},P,WitnessDegreeBound,B)\f$ such that: \n
 * \f$\mathbb{F}\f$ is a field
 *
 * \f$H \subset \mathbb{F}\f$ is a subset of \f$ \mathbb{F} \f$
 *
 * \f$\vec{N} \in \mathbb{F}[x]^k\f$ is a vector of univariate polynomials over \f$ \mathbb{F} \f$
 *
 * \f$P \in \mathbb{F}[x_0,x_1,x_2,\dots,x_k]\f$ is a multivariate polynomial over \f$ \mathbb{F} \f$
 *
 * \f$WitnessDegreeBound \in \mathbb{N}\f$ is an apper bound for satisfying witness
 *
 * \f$B \subset \mathbb{F} \times \amthbb{F}\f$ is the boundry constraints set
 *
 *
 * For a polynomial \f$Q \in \mathbb{F}[x]\f$ we define \f$P \circ (x \vert Q \circ \vec{N}) \in \mathbb{F}[x]\f$ by
 * \f$(P \circ (x \vert Q \circ \vec{N}))(x) := P(x,Q(N_1(x)),Q(N_2(x)),\dots,Q(N_k(x)))\f$
 *
 * An ACSP partial instance \f$(\mathbb{F},H,\vec{N},P,\lambda,I)\f$ is satisfiable partial instance if and only if 
 * there exists \f$A \in \mathbb{F}[x]\f$ and \f$ \vec{x} \in \mathbb{F}^n \f$ (for some \f$ n \in \mathbb{N} \f$) such that:
 * \f{eqnarray*}{
 * \forall z \in H: (P \circ (x \vert A \circ \vec{N}))(z) = 0 \\
 * \deg A <= WitnessDegreeBound \\
 * \forall (x,y) \in B : A(x)=y \\
 * \f}
 *
 * In the code we give more descriptive names:
 *  
 * we name \f$ \mathbb{F} \f$ as 'contextField' \n
 * we name \f$H\f$ as 'vanishingSet' \n 
 * we name \f$N\f$ as 'neighborPolys' \n
 * we name \f$P\f$ as 'constraintPoly' \n
 * we name \f$B\f$ as 'boundaryConstraints' \n
 *
 *
 * Methods:\n
 * Instance class contains only getters, 
 * constructor and a destructor.
 */

class ACSPInstance {
public:
	typedef Algebra::FiniteField field;
	typedef Algebra::FiniteSetInterface set;
	typedef Algebra::UnivariatePolynomialInterface uniPoly;
	typedef std::vector<std::unique_ptr<const uniPoly>> polynomialsVec;
	typedef Algebra::PolynomialInterface polynomial;
	typedef Algebra::FieldElement fieldElement;
    typedef std::map<Algebra::FieldElement, Algebra::FieldElement, Algebra::classCompElements> boundaryConstraints_t;
    typedef std::function< std::unique_ptr<uniPoly>(const ACSPInstance&, const uniPoly& witness)> compAlgorithm_t;

	ACSPInstance(
		const field& contextField,
		std::unique_ptr<const set>&& vanishingSet,
		polynomialsVec&& neighborPolys,
		std::unique_ptr<const polynomial>&& constraintPoly,
        const Algebra::PolynomialDegree& witnessDegreeBound, 
		const boundaryConstraints_t& boundaryConstraints,
        const compAlgorithm_t& compositionAlgorithm = naiveComposition_N_division_Alg
		)
		:
			contextField_(contextField),
			vanishingSet_(std::move(vanishingSet)),
			neighborPolys_(std::move(neighborPolys)),
			constraintPoly_(std::move(constraintPoly)),
            witnessDegreeBound_(witnessDegreeBound),
			boundaryConstraints_(boundaryConstraints),
            compositionAlgorithm_(compositionAlgorithm)
			{};
	
	/**
	 * @brief   Move constructor
	 * should be default, but not supported
	 * by MS VS2010, so just implements 'default'
	 */

	ACSPInstance(ACSPInstance&& src): // = default;
					contextField_(std::move(src.contextField_)),
					vanishingSet_(std::move(src.vanishingSet_)),
					neighborPolys_(std::move(src.neighborPolys_)),
                    witnessDegreeBound_(std::move(src.witnessDegreeBound_)),
					constraintPoly_(std::move(src.constraintPoly_)),
                    boundaryConstraints_(std::move(src.boundaryConstraints_)),
                    compositionAlgorithm_(std::move(src.compositionAlgorithm_))
	{};

	ACSPInstance() _COMMON_CXX11_DELETED; // Workaround. In MSVC we need the default
                                               // constructor to exist in order for the
                                               // ACSPFullInstance copy constructor to compile.
                                               // This constructor can be deleted if MSVC ever
                                               // supports C++11 =delete
    
	inline const field& contextField()const {
		return contextField_;
	}
	
	inline const set& vanishingSet()const {
		return *(vanishingSet_.get());
	}

	inline const polynomialsVec& neighborPolys()const {
		return neighborPolys_;
	}

	inline const polynomial& constraintPoly()const {
		return *(constraintPoly_.get());
	}

    inline const boundaryConstraints_t& boundaryConstraints()const{
        return boundaryConstraints_;
    }
    
    inline const Algebra::PolynomialDegree witnessDegreeBound()const{
        return witnessDegreeBound_;
    }

    inline std::unique_ptr<uniPoly> composeWithWitness_and_divideByVanishingSpacePoly(const uniPoly& witness)const{
        return compositionAlgorithm_(*this,witness);
    }
	
    private:

	field contextField_; 
	std::unique_ptr<const set> vanishingSet_;
	polynomialsVec neighborPolys_;
	std::unique_ptr<const polynomial> constraintPoly_;
	Algebra::PolynomialDegree witnessDegreeBound_;
    boundaryConstraints_t boundaryConstraints_;

    //A hint for fast construction of the composition polynomial
    compAlgorithm_t compositionAlgorithm_; 

	/**
	 * @brief  Copy constructor
	 * should be deleted, but not supported
	 * by MS VS2010, so it is private and throws exception 
	 */
	ACSPInstance(const ACSPInstance& src) _COMMON_CXX11_DELETED;

    static std::unique_ptr<uniPoly> naiveComposition_N_division_Alg(const ACSPInstance& instance, const uniPoly& witness){
        using std::vector;
        using std::min;
        using std::unique_ptr;
        using Algebra::FieldElement;
        using Algebra::PolynomialDegree;
        using Algebra::zero;
        using Algebra::UnivariatePolynomialGeneral;
       
        //
        //get the composition degree bound
        //
        vector<PolynomialDegree> constraintsInputDegrees;

        // first input is "x" which has degree 1
        constraintsInputDegrees.push_back(PolynomialDegree(1));

        // rest are composition of neighbor with witness
        const auto witnessDegree = witness.getDegree();
        for (const auto& n : instance.neighborPolys()){
            constraintsInputDegrees.push_back(n->getDegreeBound(witnessDegree));
        }

        // get the composition degree bound
        const PolynomialDegree degBound = instance.constraintPoly().getDegreeBound(constraintsInputDegrees);

        //
        //define interpolation space
        //
        const size_t degForPoly = degBound.isInteger()? ceil(log2(1+PolynomialDegree::integral_t(degBound))) : 0;
        const auto interpolationBasis = Algebra::details::buildStandardBasis(min(instance.contextField().degree(), degForPoly));

        //construct evaluation
        vector<FieldElement> evaluation(Infrastructure::POW2(interpolationBasis.asVector().size()));
        for(size_t i=0; i< evaluation.size(); i++){
            const FieldElement x = getSpaceElementByIndex(interpolationBasis.asVector(),zero(),i);
            
            //construct the assignment for the constraints poly
            vector<FieldElement> assignment;
            assignment.push_back(x);
            for(const auto& n : instance.neighborPolys()){
                assignment.push_back(witness.eval(n->eval(x)));
            }

            //evaluate and return
            evaluation[i] = instance.constraintPoly().eval(assignment);
        }

        UnivariatePolynomialGeneral compositionPoly(evaluation,interpolationBasis.asVector(),zero());

        //build denominator and divide
        return instance.vanishingSet().vanishingPoly()->divideByMe(compositionPoly);
    }

};

}// namespace PCP_Project

#endif   // __ACSP_Instance_HPP
