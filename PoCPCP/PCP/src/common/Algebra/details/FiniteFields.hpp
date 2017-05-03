/********************************************** FiniteFields.hpp ************************************************/
/**
 * @file.
 *
 * The file FiniteFields.hpp contains the implementation of finite fields functionality, including basis, linear
 * spaces, randomness etc.
 * 
 * Main classes defined in the file:
 * Basis - Implements the  functionality of a a field's basis, represented as a vector of field elements.
 * LinearSubset - Implements the functionality of fields & subsets, without actually keeping the field's elements in memory.
 * ExplicitLinearSubset - Implements the functionality of fields & subsets, by computing all the elements and saving them into memory.
 * 
 * Main date types defined in the file:
 * RandomnessString - A random binary string is stored as a vector of boolean variables.
 * basisShiftPair - A struct that contains a basis and an affine shift.
 *
 * Main functions defined in the file:
 * initGF2E - Finds an irreducible polynomial of the given degree and initialize GF2E elements.
 * Basis Functions - Generate a basis, Add or remove element, build a standard basis etc.
 * LinearSubset Functions - Construct a subset (without computing the elements), get an element by index etc.
 * ExplicitLinearSubset Functions - Construct a subset (compute and save the elements), get an element by index in O(1) time etc.
 */
  /************************************************************************************************************/
#ifndef FINITE_FIELDS_HPP_
#define FINITE_FIELDS_HPP_

#include <algebraLib/FieldElement.hpp>

#include <vector>
#include <memory>
#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>

#include <NTL/mat_GF2.h>
#include <NTL/mat_GF2E.h>
#include <NTL/vec_vec_GF2.h>

namespace Algebra {
namespace details {

/* Forward Declarations */
class UnivariatePolynomial;
class UniPolynomialEvaluation;
class BivariatePolynomialEvaluationMap;
class BivariatePolynomialEvaluationNaive;

/** 
 * Find an irreducible polynomial of the given degree and initialize GF2E elements. 
 * @param degree is the required degree of the irrducible polynomial (i.e. the field dimension).
 */
void initGF2E(int degree);	
typedef long BasisIndex;		///An index of a basis element.
typedef uint64_t FieldIndex;	///An index of a field element.

/** Univariate polynomial over GF(2) type */
typedef NTL::GF2X GF2Poly;


/********************************************************/
/*************** Finite Fields Behavior *****************/
/********************************************************/

/****************************/
/********** Basis ***********/
/****************************/

/**
 * The class Basis contains basic functionality of a a field's basis.
 * The basis is actually a vector of vectors over GF(2), namely a vector
 * of field points.
 */
class Basis {
private:

	/* Basis Serialization: */
	friend class boost::serialization::access;
	template <class Archive>
	void save (Archive & ar,const unsigned int version)const{
		ar & basisElements;
		ar & basisSize;
	}
	template <class Archive>
	void load (Archive & ar,const unsigned int version){
		ar & basisElements;
		ar & basisSize;
	}
	BOOST_SERIALIZATION_SPLIT_MEMBER()


protected:

	/* Class members: */
	std::vector<FieldElement> basisElements;		///The representation is a vector of field points.
	/*vec_GF2E basisElements;		///The representation is a vector of field points.*/
	BasisIndex basisSize;		///The number of vectors in the basis.

public:

    const std::vector<FieldElement>& asVector()const;
    
	/**
	 * The class's default C'tor
	 */
	Basis();

	/**
	 * Class C'tor - Builds a Basis object according to specific elements.
	 */
	Basis(const std::vector<FieldElement> b);
	//Basis(const vec_GF2E b);

	/**
	 * Class C'tor - Builds a Basis object of a field represented by the input irreducible polynomial.
	 */
	Basis(const GF2Poly& Irr);

	/**
	 *  The class's Copy C'tor
	 */
	Basis(const Basis& copy_from);

	/**
	* Returns an element of the basis according to a given index.
	* param idx is the index of the element.
	*/
	FieldElement getBasisElementByIndex(const BasisIndex idx) const;

	/**
	* Returns an element of the basis according to a given index.
	* param idx is the index of the element.
	*/
	FieldElement& getBasisElementByIndexByRef(const BasisIndex idx) const;

	/**
	 * Returns a Basis object which contains part of the elements in "this" object.
	 * @params from and to specify the range of this vector from which the result Basis
	 * object is built.
	 */
	Basis getPartialBasisFromTo(const BasisIndex from, const BasisIndex to) const;

	/**
	 * Add an element (of type GF2E) to the vector of elements.
	 * @param elem - element to be added.
	 */
	void addElement(const FieldElement& elem);

	/**
	* Add an element (of type GF2E) to the vector of elements, at the front, and shift all elements on place forward.
	* @param elem - element to be added.
	*/
	void addElementInFront(const FieldElement& elem);


	/**
	* returns this basis in reversed order.
	*/
	Basis reverseBasis();

	/**
	* Removes the last element inserted to the basis.
	* The function is used by the RS= prover, which builds the basis of L_beta by starting with L0', and going
	* over all beta such that for each beta it adds the element to the basis by calling addElement, and later removes
	* it by calling removeLastElement, before moving on the next beta.
	*/
	void removeLastElement();

	/**
	* Returns the index of the given element in the field spanned by the basis. 
	* @param elem is the field element whose index we're interested in.
	* @Implementation - The function finds the coefficients vector using NTL's mat_GF2, thus
	* computing the index of the element.
	*/
	FieldIndex getIndexOfElement(const FieldElement& elem, const bool naive = false) const;

	
	/**
	 * The function returns True iff the element is a part of the basis.
	 * @param elem is the element we're after.
	 */
	virtual bool containsElement(const FieldElement& elem) const;

	/**
	 * Returns the size of the basis (number of vectors).
	 */
	BasisIndex getSizeOfBasis() const;

	/**
	* Returns a matrix of GF2 elements.
	* @param resultVector - the vector on the last column that represent the result that we are interested in.
	* The returned matrix size is (max,basisSize+1), when max is the length of the longest base word.
	* On the first basisSize columns of these matrix we put the basis elements.
	* On the remaining columns we put the result vector.
	* We use this function from the "findIndexByElement" functions, in order to find the index (=the right combination of the basis elements
	* that should give us the resultVector
	*/
	NTL::mat_GF2 turnBasisToMatrix(NTL::vec_GF2& resultVector) const;
	
	/**
	 * Testing purposes only
	 */
	void print() const;

    

	/**
	 * The class destructor
	 */
	virtual ~Basis();
};


/** Returns the standard basis of the form 1,x,x^2,...,x^size */
extern Basis buildStandardBasis(int size);

/** A struct that contains a basis and an affine shift */
struct basisShiftPair {
	Basis basis;
	FieldElement affineShift;
    
    basisShiftPair():affineShift(Algebra::zero()){};
    basisShiftPair(const Basis& b, const FieldElement& s):basis(b),affineShift(s){};
};

/****************************/
/******* LinearSubset *******/
/****************************/
/**
 * The class LinearSubset Implements the functionality of fields & subsets, without actually keeping 
 * the field's elements in memory.
 */
class LinearSubset : public Basis {

protected:

	FieldIndex fieldSize;		///Number of elements in the subset. Added for efficiency 

public:

	/**
	 * The class's default C'tor
	 */
	LinearSubset();

	/**
	 * Class C'tor - Builds a LinearSubset object whose basis is param b.
	 */
	LinearSubset(const Basis b);

	/**
	 *  The class's Copy C'tor
	 */
	LinearSubset(const LinearSubset& copy_from);

	/**
	 * Returns an element of the set according to a given index.
	 * param idx is the index of the element.
	 */
	virtual FieldElement getFieldElementByIndex(const FieldIndex idx) const;

	/**
	 * Returns the size of the set \ field (number of elements).
	 * This method really should have been called getSizeOfSpace, but messy to change now.
	 */
	FieldIndex getSizeOfField() const;

};

/****************************/
/*** ExplicitLinearSubset ***/
/****************************/
/**
 * The class ExplicitLinearSubset handles behavior involved with fields & subsets.
 * The elements are calculated and saved when the object is constructed, and are later
 * available with an O(1) access.
 */
class ExplicitLinearSubset : public LinearSubset {

private:

	/* The elements are saved as a vector of field points. */
	std::vector<FieldElement> elements;	

	/* The affine shift of the space (if exists) */
	FieldElement affineShift;

	/**
	 * A non-recursive span algorithm. Spans the basis and fills in the vector of elements.
	 */
	void span(const Basis& basis);

public:

	/**
	 * The class's default C'tor
	 */
	ExplicitLinearSubset();

	/**
	 * Class C'tor - Builds a LinearSubset object whose basis is param b.
	 */
	ExplicitLinearSubset(const Basis b, const FieldElement& affineShift = Algebra::zero());

	/**
	 *  The class's Copy C'tor
	 */
	ExplicitLinearSubset(const ExplicitLinearSubset& copy_from);

	/**
	 * Returns an element of the set according to a given index.
	 * The function overrides LinearSubset's function for a better
	 * complexity of O(1).
	 * param idx is the index of the element.
	 */
	virtual FieldElement getFieldElementByIndex(const FieldIndex idx) const;

	/** Sets the affine shift of the space */
	void setAffineShift(const FieldElement& newAffineShift);

	/**
	 * The function returns True iff the element is a part of the field.
	 * @param elem is the element we're after.
	 */
	virtual bool containsElement(const FieldElement& elem) const;

	/**
	 * The function returns a union between the current span and the span of the basis element.
	 * For example: L.addBasisAndSpan(a) = L U (L+a).
	 * The way it works is by shifting all the elements by element.
	 * @param element is the basis element to be added.
	 */
	void addBasisAndSpan(const FieldElement& element);

	/**
	* The function removes the last element of the basis and spans the field accordingly.
	* The function is used by the RS= prover, which builds L_beta for every beta by starting 
	* with L0', and calling addBasisAndSpan for a specific beta, using it, and removing it back
	* to L0' by calling removeLastBasisAndSpan. 
	*/
	void removeLastBasisAndSpan();

	/**
	* The function returns a new ExplicitLinearSubset object that contains the 2nd half of the original's ExplicitLinearSubset elements.
	* For example, if this object contains element a_1, a_2...a_2n, the return object will include a_n+1,...,a_2n
	*/
	//ExplicitLinearSubset getShiftedSubset() const;

	/**
	 * Testing purposes only.
	 */
	void printElements() const;

	/**
	 * The class destructor
	 */
	~ExplicitLinearSubset();
};

} // of namespace details
} // of namespace Algebra

#endif /* FINITE_FIELDS_HPP_ */
