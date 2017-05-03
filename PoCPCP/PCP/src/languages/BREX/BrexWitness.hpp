#ifndef BREX_WITNESS_HPP__
#define BREX_WITNESS_HPP__

#include <algebraLib/FieldElement.hpp>
#include "common/langCommon/Sequence.hpp"
#include "common/CXX11_macros.hpp"

#include <memory>
#include <vector>

namespace PCP_Project{

/**
 * @class BREXWitness
 * @brief class for BREX Witness
 *
 * A BREX witness for BREX instance \f$(\mathbb{F},d,\mathcal{V},\mathcal{C}_\mathcal{A}, \mathcal{C}_\pi, B)\f$
 * is a pair \f$ (\mathcal{A},\pi) \f$ such that:
 *
 * \f$ \mathcal{A} : \mathcal{D} \to (\mathcal{V} \to \mathbb{F}) \f$ which is called the assignment
 *
 * \f$ \pi : \mathcal{D} \to \mathcal {D} \f$ which is called the permutation
 *
 * For two assignments \f$ \alpha , \beta : \mathcal{V} \to \mathbb{F} \f$ we define \f$ (\alpha,\beta):\mathcal{V} \times\{0,1\} \to \mathbb{F} \f$ by:
 * \f[ 
 * [(\alpha,\beta)](v,i) =
 * \begin{cases}
 * \alpha\left(v\right) & i=0\\
 * \beta\left(v\right) & i=1
 * \end{cases}
 * \f]
 *
 * It said that a BREX instance \f$(\mathbb{F},d,\mathcal{V},\mathcal{C}_\mathcal{A}, \mathcal{C}_\pi, B)\f$
 * is satisfied by a BREX witness \f$ (\mathcal{A},\pi) \f$ the following constraints hold:
 *
 * Assignment constraints satisfaction:
 * \f[
 * \forall n \in \mathcal{D} \setminus \{ 2^d - 2 \}, \forall p \in \mathcal{C}_\mathcal{A} : p(\mathcal{A}(n), \mathcal{A}(n+1)) = 0
 * \f]
 *
 * Permutation constraints satisfaction:
 * \f[
 * \begin{align*}
 * \left[\mathcal{V}\times\left\{ 1\right\} \right]\left(\mathcal{C}_\pi\right)\neq\emptyset\Rightarrow & \pi \in S_\mathcal{D}\text{ (meaning it is a permutation)}\\
 * \forall n \in \mathcal{D} \setminus \{2^d -2 \},\forall p \in \mathcal{C}_\pi : & p\left(\mathcal{A} \left(n\right),\mathcal{A} \left( \pi \left(n\right)\right)\right)=0
 * \end{align*}
 * \f]
 * 
 * Boundary constraints:
 * \f[
 * \forall ((x,y),\alpha) \in B : [\mathcal{A}(x)](y) = \alpha
 * \f]
 *
 * In the code we give more descriptive names:
 * 
 * we name \f$ \mathcal{A} \f$ as 'get_assignment' \n
 * we name \f$ \pi \f$ as 'permutation' \n
 *
 *
 * Methods:\n
 * Instance class contains only getters, 
 * constructor and a destructor.
 */
class BREXWitness {
public:
	typedef std::vector<Algebra::FieldElement> color_t;
	typedef Sequence<color_t> assignment_t;
	typedef std::unique_ptr<assignment_t> assignment_ptr;

	typedef Sequence<size_t> permutation_t;
	typedef std::unique_ptr<permutation_t> permutation_ptr;
	
	BREXWitness(
		assignment_ptr&& assignment,
		permutation_ptr&& permutation)
	:
		assignment_(std::move(assignment)),
		permutation_(std::move(permutation)){};
	
	/**
	 * @brief   Move constructor
	 * should be default, but not supported
	 * by MS VS2010, so just implements 'default'
	 */
	BREXWitness(BREXWitness&& src) : // = default;
		assignment_(std::move(src.assignment_)),
		permutation_(std::move(src.permutation_)){};

	inline color_t get_color(size_t vecIndex)const {
		return assignment_->getElementByIndex(vecIndex);
	}
	
	inline const permutation_t& permutation()const {
		return *permutation_;
	}
	
private:
	assignment_ptr assignment_;
	permutation_ptr permutation_;
	
	/**
	 * @brief  Copy constructor
	 * should be deleted, but not supported
	 * by MS VS2010, so it is private and throws exception 
	 */
	BREXWitness(const BREXWitness& src) _COMMON_CXX11_DELETED;
};

} // namespace PCP_Project

#endif // BREX_WITNESS_HPP__
