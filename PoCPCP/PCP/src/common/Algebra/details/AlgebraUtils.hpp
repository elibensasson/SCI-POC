/********************************************** commonUtils.hpp ************************************************/
/** @file.
 *
 * The file commonUtils.hpp contains various utility functions and date type definitions used by the different 
 * classes in the project.
 * 
 * Main classes defined in the file:
 * ElementsVector - Represents a vector of field elements, which can store a large number of elements (larger than 2^32).
 * 
 * Main date types defined in the file:
 * BasisIndex - An index of a basis element.
 * FieldIndex - An index of a field element.
 * FieldPoint - An element in an extension field of GF(2).
 * FieldPointPair - A two-dimensional point (x,y) in the field.
 * GF2Poly - Univariate polynomial over GF(2) type.
 *
 * Main functions defined in the file:
 * Basic Math Functions - Power macro, log with basis 2, extract a bit from an integer etc.
 * GF2E Functionality - Convert an integer to a vector of GF2 bits and vice versa, compare two GF2 extension field elements etc.
 */
/************************************************************************************************************/
#pragma once

#ifndef ALGEBRA_UTILS_HPP_
#define ALGEBRA_UTILS_HPP_

#include <vector>
#include <cstdint>

#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>

namespace Algebra {
namespace details {

	/********************************************************/
	/*************** Data Types Definitions *****************/
	/********************************************************/
	typedef long BasisIndex;		///An index of a basis element.
	typedef uint64_t FieldIndex;	///An index of a field element.

	typedef NTL::GF2E FieldPoint;		///An element in an extension field of GF(2)				
	typedef struct {				///A 2 dimensional point (x,y).
		FieldPoint x;
		FieldPoint y;
	} FieldPointPair;

	/** Univariate polynomial over GF(2) type */
	typedef NTL::GF2X GF2Poly;

	/********************************************************/
	/***************** GF2E Functionality *******************/
	/********************************************************/

	/** Returns a constant finite field element initialized by the given integer. */
	FieldPoint constFieldElement(const int i);

	/** 
	 * Conversions between integer and a list of field elements. 
	 * The result vector of elements is ordered LSB...MSB, i.e. index 0 contains the LSB.
	 */
	void convertIntToGF2s(const int bitsNum, FieldIndex index, std::vector<NTL::GF2>& indexBits);
	void convertIntToGF2sMSBFirst(const int bitsNum, FieldIndex index, std::vector<NTL::GF2>& indexBits);
	int convertGF2sToInt(const std::vector<NTL::GF2>& indexBits);
	int convertG2sToIntMSBFirst (const std::vector<NTL::GF2>& indexBits);
	void convertIntToGF2s(const int bitsNum, FieldIndex index, std::vector<FieldPoint>& indexBits);
	int convertGF2sToInt(const std::vector<FieldPoint>& indexBits);

	/** Converts a GF2E vector to a GF2 vector */
	//vector<GF2> convertNTLVectors(vector<GF2E> inp);

	/**
	 * Comparison function for NTL's field element type GF2E.
	 * @param elem1 and elem2 are the two field elements to compare.
	 */
	inline bool operator< (const NTL::GF2E& elem1, const NTL::GF2E& elem2) {

		long wlen1 = elem1._GF2E__rep.xrep.length();
		long wlen2 = elem2._GF2E__rep.xrep.length();

		if (wlen1 != wlen2)
			return (wlen1 < wlen2);

		const unsigned long* elem1_rep = elem1._GF2E__rep.xrep.elts();
		const unsigned long* elem2_rep = elem2._GF2E__rep.xrep.elts();

		for (long i = (wlen1)-1; i>=0; i--) {
			if (elem1_rep[i] != elem2_rep[i])
				return elem1_rep[i] < elem2_rep[i];
		}

		return false;

		/*
		 * GF2E is represented in the following manner:
		 * elem1, elem2	 - GF2E
		 * _GF2E__rep	 - GF2X
		 * x_rep		 - WordVector
		 * rep			 - _ntl_ulong* (unsigned long. every bit is a GF2 element)
		 */
	}

	/**
	 * GF2E_Compare Compares 2 elements in the field. Used for stl map data structure that represents
	 * a univariate evaluation.
	 * @param elem1 and elem2 are the two field elements to compare.
	 */
	struct GF2E_Compare {
		bool operator() (const NTL::GF2E& elem1, const NTL::GF2E& elem2) const {
			return elem1 < elem2;
		}
	};

	/**
	 * GF2E_Pair_Compare Compares two 2-dimensional elements in the field. Used for stl map data structure that
	 * represents a bivariate evaluation.
	 * @param elem1 and elem2 are the two 2D field elements to compare.
	 */
	struct GF2E_Pair_Compare {
		bool operator() (const FieldPointPair& elem1, const FieldPointPair& elem2) const {
			return ((elem1.x < elem2.x) || ((elem1.x == elem2.x) && (elem1.y < elem2.y)));
		}
	};

	/********************************************************/
	/**************** ElementsVector Class ******************/
	/********************************************************/
	/**
	 * The following class implements a vector of field elements.
	 * The class is actually a wrapper for NTL's vec_GF2E, which can allocate and access elements very efficiently.
	 * HOWEVER, NTL's vector has a limited size of 2^24 elements, which is not enough. Therefore, the class is currently
	 * used for testing purposes only.
	 */
	//class ElementsVector_NTL {
	//
	//private:
	//
	//	vec_GF2E vecRep;
	//
	//public:
	//	
	//	ElementsVector_NTL();
	//	ElementsVector_NTL(const int64_t len);
	//
	//	int64_t length() const;
	//	void SetLength(const int64_t len);
	//
	//	FieldPoint& operator[](const int64_t& idx);
	//	const FieldPoint& operator[](const int64_t& idx) const;
	//	
	//	void kill();
	//};

#define BUCKET_SIZE Infrastructure::POW2(24)
	/**
	 * The class implements a vector of field elements, which can store a large number of elements (larger than 2^32).
	 * The implementation is done by using "buckets" of NTL's vec_GF2E. Each bucket has its size, and the number of buckets
	 * is set according to the requested length. 
	 */
	class ElementsVector_Buckets {

		private:

			std::vector<NTL::vec_GF2E> elementsRep;	///A vector of NTL vectors, each contains BUCKET_SIZE elements.
			int numBuckets;					///The number of allocated buckets.
			int64_t size;				///Keeps the length (total number of elements), saved for performace.

		public:

			/** Class Constructors */
			ElementsVector_Buckets();
			ElementsVector_Buckets(const int64_t len);

			/** Setters and getters of the vector's length. */
			int64_t length() const;
			void SetLength(const int64_t len);

			/** O(1) access to the vector elements. */
			FieldPoint& operator[](const int64_t& idx);
			const FieldPoint& operator[](const int64_t& idx) const;

			/** Destroys the buckets, leavng a single empty bucket. */
			void kill();
	};

	typedef ElementsVector_Buckets ElementsVector;

} // of namespace details
} // of namespace Algebra

#endif	//ALGEBRA_UTILS_HPP_
