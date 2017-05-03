/********************************************** FiniteFields.cpp ************************************************/
/**
 * @file.
 *
 * The file FiniteFields.cpp contains the implementation of finite fields functionality, including basis, linear
 * spaces, randomness etc.
 * 
 * For more information - Read the documentation of the header file FiniteFields.hpp.
 */
  /************************************************************************************************************/
#include <ctime>
#include "FiniteFields.hpp"
#include "Polynomials.hpp"
#include <NTL/GF2XFactoring.h>
#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>

namespace Algebra {
namespace details {

/**Find an irreducible polynomial of the given degree and initialize the extension field. */
void initGF2E(int degree) {
	::NTL::GF2X Irr;
	::NTL::BuildIrred(Irr, degree);
	UnivariatePolynomial::IrreduciblePolynomial = Irr;
	::NTL::GF2E::init(Irr);
}

/********************************************************/
/**************** Basis Implementation ******************/
/********************************************************/

/** Default C'tor */
Basis::Basis() {
	basisElements.resize(0);
	basisSize = 0;
}

/** Builds a Basis object according to specific elements. */
Basis::Basis(const vector<FieldElement> b) {
	this->basisElements = b;				///assignment operator of vector<FieldElement> is O(n)
	basisSize = ::Infrastructure::safeConvert(basisElements.size());
}

//Basis::Basis(const vec_GF2E b) {
//	this->basisElements = b;				///assignment operator of vector<FieldElement> is O(n)
//	basisSize = basisElements.length();
//}

/** The class's copy C'tor */
Basis::Basis(const Basis& copy_from) {
	this->basisElements = copy_from.basisElements;	//assignment operator of vector<FieldElement> is O(n)
	basisSize = copy_from.basisSize;
}

/** Constructs an object that is the standard basis with deg(Irr) elements. */
Basis::Basis(const GF2Poly& Irr) {
	int degree = deg(Irr);
	::NTL::GF2X temp;

	basisElements.resize(degree);
	basisSize = degree;
	for (int i=0; i<degree; i++) {
		clear(temp);
		SetCoeff(temp, i);
		basisElements[i] = FieldElement(::NTL::to_GF2E(temp));
	}
}

/** Returns an element of the basis according to a given index. */
FieldElement Basis::getBasisElementByIndex(const BasisIndex idx) const {
	return getBasisElementByIndexByRef(idx);
}

/** Returns an element of the basis according to a given index. */
FieldElement& Basis::getBasisElementByIndexByRef(const BasisIndex idx) const {
	_COMMON_ASSERT(idx >= 0 && idx<basisSize ,"Bad indices in getBasisElementByIndex function");
	return (FieldElement&)basisElements[idx];			//O(1) access to vector<FieldElement>
}

/**
 * Returns a Basis object which contains part of the elements in "this" object.
 * Complexity is O(to-from).
 */
Basis Basis::getPartialBasisFromTo(const BasisIndex from, const BasisIndex to) const {

	if ((from > to) || (to > basisSize-1)){
		_COMMON_FATAL("Bad indices in getPartialBasisFromTo function");
	}
	vector<FieldElement> newBasis;
	newBasis.resize(to-from+1);
	BasisIndex i,j;
	for (j=0, i=from; i<=to; j++, i++) {
		newBasis[j] = this->basisElements[i];
	}
	return Basis(newBasis);
}

/** Add a field element to the vector of elements. */
void Basis::addElement(const FieldElement& elem) {
#ifdef _DEBUG //when not running in debug mode assuming no one is readding element, and so don't waste time checking for it
	for (BasisIndex i=0; i<basisSize; i++) {
		_COMMON_ASSERT((basisElements[i] != elem), "Element already exist in Basis::addElement function");
		if (basisElements[i] == elem) {
			std::cout << elem << std::endl;
			_COMMON_FATAL("Element already exist in Basis::addElement function");
		}
	}
#endif
	if (basisElements.size() < basisSize+1) 
		basisElements.resize(basisSize+1);

	basisElements[basisSize] = elem;
	basisSize++;
}


/**
* Add an element (of type GF2E) to the vector of elements, at the front, and shift all elements on place forward.
* @param elem - element to be added.
*/
void Basis::addElementInFront(const FieldElement& elem){
#ifdef _DEBUG //when not running in debug mode assuming no one is readding element, and so don't waste time checking for it
for (BasisIndex i = 0; i<basisSize; i++) {
	_COMMON_ASSERT((basisElements[i] != elem), "Element already exist in Basis::addElement function");
	if (basisElements[i] == elem) {
		std::cout << elem << std::endl;
		_COMMON_FATAL("Element already exist in Basis::addElement function");
	}
}
#endif
//
//if (basisElements.size() < basisSize + 1)
//	basisElements.resize(basisSize + 1);

basisElements.insert(basisElements.begin(), elem);
basisSize++;
}


/**
* reverses the order of basis elements
*/
Basis Basis::reverseBasis(){
	Basis res;
	for (BasisIndex i = basisSize - 1; i >= 0; i--)
		res.addElement(this->getBasisElementByIndex(i));
		
	return res;
}


/* 
* Removes the basis element that was added last.
* @Implementation - Since we use vectors, it is enough to decrease the vector size by 1. 
*/
void Basis::removeLastElement() {	
	basisSize--;
}

/**
 * The function returns True iff the element is a part of the basis.
 */
bool Basis::containsElement(const FieldElement& elem) const {
	for (int i=0; i<basisSize; i++) {
		if (basisElements[i] == elem)
			return true;
	}
	return false;
}

/** 
 * This function is an auxiliary function used by the function getIndexOfElement().
 * The function solves a linear equations system and outputs the vector of bits that form the index.
 * @param indexVector is the result vector of the index bits.
 * @param base_matrix contains the coefficients of the relevant field element.
 */
void solveEq(::NTL::vec_GF2 & indexVector, ::NTL::mat_GF2 base_matrix){
	long n = base_matrix.NumRows();
	long m = base_matrix.NumCols()-1;
	
	//ranking the matrix
	gauss(base_matrix);

	//If everything ok, after the ranking we have m non-zero rows, corresponding to the m variants that we have.
	//So we create a square matrix with the first m rows and solve the equation with the solve method 
	::NTL::mat_GF2 solve_mat;
	solve_mat.SetDims(m,m);
	for (long i=0;i<m;i++){
		for (long j=0;j<m;j++){
			//its [j][i] and not [i][j] because NTL solves index_Vector*matrix = result and not matrix*index_Vector = result
			solve_mat[i][j] = base_matrix[j][i];
		}
	}
	::NTL::vec_GF2 result_vector;
	result_vector.SetLength(m);
	for (long i=0;i<m;i++){
		result_vector[i] = base_matrix[i][m];
	}
	::NTL::GF2 gfone = ::NTL::to_GF2(1);
	solve (gfone,indexVector,solve_mat,result_vector);
	return;
}

/**
 * Returns the index of the given element in the field spanned by the basis. 
 * @param element is the field element whose index we seek.
 * @param naive is flag that tells whether the computation should be done in a naive manner or not.
 */
FieldIndex Basis::getIndexOfElement(const FieldElement& element, const bool naive) const {

	///Naive implementation - span the linear space and look for the element.
	if (naive) {
		ExplicitLinearSubset set(*this);
		for (FieldIndex i=0; i<set.getSizeOfField(); i++)
			if (set.getFieldElementByIndex(i) == element)
				return i;

		std::cout << "Index translation: Could not find element in field" << std::endl;
		return 0;
	}
	
	///Non-naive and efficient way - solve linear equations:
	::NTL::vec_GF2 resultVector = ::NTL::to_vec_GF2(::NTL::GF2E(element)._GF2E__rep);
	
	//create an appropriate matrix with the equation and results, for ranking
	::NTL::mat_GF2 basisMatrix = this->turnBasisToMatrix(resultVector);
	::NTL::vec_GF2 indexVector;
	indexVector.SetLength(basisMatrix.NumCols());
	
	//solve the equations, the result will be in indexVector
	solveEq(indexVector,basisMatrix);
	//calculate the index
	long Index = 0;
	for (int i=0; i<indexVector.length(); i++)	
		if ((((::NTL::GF2)(indexVector[i]))._GF2__rep)) Index += ::Infrastructure::safeConvert(::Infrastructure::POW2(i));

	return Index;	
}


/** Returns the size of the basis (number of vectors). */
BasisIndex Basis::getSizeOfBasis() const {
	return basisSize;
}

const std::vector<FieldElement>& Basis::asVector()const{
    return basisElements;
}


/**
 * Converts a vector of GF2 elements to a m X m matrix of GF2 elements.
 */
::NTL::mat_GF2 Basis::turnBasisToMatrix(::NTL::vec_GF2 &resultVector) const{

	::NTL::mat_GF2 result;
	long max = ::NTL::GF2E::degree();

	result.SetDims(max,this->basisSize+1);
	resultVector.SetLength(max);
	if (max<this->basisSize){
		std::cout<<"WRONG BASIS"<<std::endl;
	} 
	int k = this->basisSize;
	
	for (int j=0;j<k;j++){		
		long length = ::NTL::GF2E(basisElements[j])._GF2E__rep.xrep.length();
		for (long i=0; i<max; i++){
			long wi = i/NTL_BITS_PER_LONG;
			if (wi >= length) result[i][j] = ::NTL::to_GF2(0);
			long bi = i - wi*NTL_BITS_PER_LONG;
			result[i][j] = ::NTL::to_GF2((::NTL::GF2E(this->basisElements[j])._GF2E__rep.xrep[wi] & (1UL << bi)) != 0);
		}

	}
	for (long i=0;i<max;i++){
		result[i][k] = resultVector[i];
	}
	return result;
}

/** Testing purposes only */
void Basis::print() const {
	std::cout << "Basis elements are: ";
	for (BasisIndex i=0; i<basisSize; i++) {
		std::cout << basisElements[i];
	}
	std::cout << std::endl;
}

/** Class Destructor */
Basis::~Basis() {
}

/** Returns the standard basis of the form 1,x,x^2,...,x^{size-1} */
Basis buildStandardBasis(int size) {
	Basis result;
	for (int i=0; i<size; i++) {
		::NTL::GF2X elem_rep(i,1);
		result.addElement(FieldElement(::NTL::to_GF2E(elem_rep)));
	}
	return result;
}

/********************************************************/
/************* LinearSubset Implementation **************/
/********************************************************/

/** Default C'tor */
LinearSubset::LinearSubset() : Basis() {
	fieldSize = 0;
}

/**
 * Builds a LinearSubset object whose basis is param b.
 * @param basis is the basis elements that span the linear subset.
 * Complexity - Linear in the number of vectors in b.
 */
LinearSubset::LinearSubset(const Basis b) : Basis(b) {

	///basisElements are initialized in the superclass object
	fieldSize = ::Infrastructure::POW2( this->basisSize);
}

/**
 * The class's Copy C'tor.
 * Complexity - Linear in the number of vectors in copy_from's basis.
 */
LinearSubset::LinearSubset(const LinearSubset& copy_from) : Basis(copy_from) {

	//basisElements are initialized in the superclass object
	fieldSize = copy_from.fieldSize;
}

/**
 * Returns an element of the set according to a given index.
 * The function maps a number to a field element by looking at the
 * number as a string of bits (number of bits is the number of vectors in
 * field's basis), and summing the vectors whose bits are on.
 * Complexity - Linear in the size of the basis.
 */
FieldElement LinearSubset::getFieldElementByIndex(const FieldIndex idx) const {

	if ((idx < 0) || (idx > fieldSize))
		_COMMON_FATAL("Bad indices in getFieldElementByIndex function");

	FieldElement result = zero();
	BasisIndex i;
	FieldIndex index = idx;
	for (i=0; i<basisSize; i++) {
		if (index&1) {
			result += basisElements[i];
		}
		index >>= 1;
	}

	return result;
}


/**
 * Returns the size of the set
 */
FieldIndex LinearSubset::getSizeOfField() const {
	return fieldSize;
}

/********************************************************/
/******** ExplicitLinearSubset Implementation ***********/
/********************************************************/

/**
 * A non-recursive span algorithm. Spans the basis and fills in the vector of elements.
 * @param basis is the basis of the spanned field.
 */
void ExplicitLinearSubset::span(const Basis& basis) {
	long len;

	elements.resize(::Infrastructure::safeConvert(fieldSize));
	elements[0] = zero();

	for(BasisIndex i=0; i<basis.getSizeOfBasis(); i++) {
		len = ::Infrastructure::safeConvert(::Infrastructure::POW2(i));	///Can improve by just doubling i

		for(long j=0; j<len; j++)
			elements[j+len] = elements[j] + basis.getBasisElementByIndex(i);
	}
}

/** Default Constructor */
ExplicitLinearSubset::ExplicitLinearSubset() : LinearSubset() {
	elements.resize(0);
	affineShift = zero();
}

/**
 * Builds a LinearSubset object whose basis is param b.
 * @param basis is the basis elements that span the linear subset.
 * @param affineShift is the affine shift in case the field is shifted.
 * Complexity - Exponential in the size of the vector b.
 */
ExplicitLinearSubset::ExplicitLinearSubset(const Basis b, const FieldElement& affineShift) : LinearSubset(b) {
	///LinearSubset object is initialized in the initialization part.

	span(b);
	this->affineShift = affineShift;

	if(affineShift != zero())
		for (FieldIndex i=0; i<fieldSize; i++) 
			elements[(long)i] = elements[(long)i] + affineShift;
}

/** The class's copy C'tor */
ExplicitLinearSubset::ExplicitLinearSubset(const ExplicitLinearSubset& copy_from) {

	basisElements = copy_from.basisElements;
	basisSize = copy_from.basisSize;
	fieldSize = copy_from.fieldSize;
	elements = copy_from.elements;
	affineShift = copy_from.affineShift;
}

/**
 * Returns an element of the set according to a given index.
 * @param idx is the index of the required element.
 * Overrides LinearSubset's function.
 * Complexity - O(1).
 */
FieldElement ExplicitLinearSubset::getFieldElementByIndex(const FieldIndex idx) const {

	assert(((idx >= 0) && (idx <= fieldSize))
			&& "MSG: ExplicitLinearSubset::getFieldElementByIndex - Bad argument for index");	
	assert(idx < INT_MAX);	///Make sure the conversion from 64 bit to 32 is safe

	return elements[(long)idx];
}

/** Sets the affine shift of the space */
void ExplicitLinearSubset::setAffineShift(const FieldElement& newAffineShift) {
	
	FieldElement oldShift = this->affineShift;
	FieldElement sumShifts = oldShift + newAffineShift;

	for (FieldIndex i=0; i<fieldSize; i++) 
		elements[(long)i] += sumShifts;		//Add old shift to cancel it + new shift.
	
	this->affineShift = newAffineShift;
}

/**
 * The function returns True iff the element is a part of the field.
 */
bool ExplicitLinearSubset::containsElement(const FieldElement& elem) const {

	for (FieldIndex i=0; i<fieldSize; i++) {
		if (getFieldElementByIndex(i) == elem)
			return true;
	}
	return false;
}

/**
 * The function returns a union between the current span and the span of the basis element.
 * For example: L.addBasisAndSpan(a) = L U (L+a).
 * The way it works is by shifting all the elements by element.
 */
void ExplicitLinearSubset::addBasisAndSpan(const FieldElement& element) {
	FieldIndex oldFieldSize = fieldSize;
	if (oldFieldSize==0) {
		vector<FieldElement> newBasis;
		newBasis.resize(1);
		newBasis[0] = element;
		*this = ExplicitLinearSubset(Basis(newBasis), affineShift);
		return;
	}

	fieldSize = oldFieldSize*2;
	basisSize += 1;
	basisElements.resize(basisSize);
	basisElements[basisSize - 1] = element;

	elements.resize(::Infrastructure::safeConvert(fieldSize));
	FieldIndex i;
	for (i=0; i<oldFieldSize; i++) {
		elements[::Infrastructure::safeConvert(oldFieldSize + i)] = elements[::Infrastructure::safeConvert(i)] + element;
	}
}

/** 
 * Removes the last element of the basis and removes the 2nd half of the subspace accordingly.
 */
void ExplicitLinearSubset::removeLastBasisAndSpan() {
	
	assert(basisSize > 0);
	this->removeLastElement();
	fieldSize /= 2;
}

/** Testing purposes only */
void ExplicitLinearSubset::printElements() const {
	FieldIndex i;
	for (i=0; i<fieldSize; i++) {
		std::cout << "Element " << i << " is " << getFieldElementByIndex(i) << std::endl;
	}
}

/** Class Destructor */
ExplicitLinearSubset::~ExplicitLinearSubset() {
}

} // of namespace details
} // of namespace Algebra
