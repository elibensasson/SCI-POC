/********************************************** commonUtils.cpp ************************************************/
/**
 * @file.
 *
 * The file commonUtils.cpp contains various utility functions and date type definitions used by the different 
 * classes in the project.
 * 
 * For more information - Read the documentation of the header file commonUtils.hpp.
 */
  /************************************************************************************************************/
#include "AlgebraUtils.hpp"
#include "../../Infrastructure/Infrastructure.hpp"

#include <NTL/GF2XFactoring.h>
#include <NTL/GF2EXFactoring.h>

#include <cassert>

using namespace std;
using namespace NTL;

namespace Algebra {
namespace details {


/*****************************************************************************/
/*************************** GF2E Functionality ******************************/
/*****************************************************************************/

/** Returns a constant finite field element initialized by the given integer. */
FieldPoint constFieldElement(const int i) {
	return to_GF2E(i);
}

/** 
* Converts a given integer to a vector of GF2 elements containing its bits. 
* GF2 indexed 0 is the LSB, GF2 indexed bitsNum-1 is the MSB.
*/
void convertIntToGF2s(const int bitsNum, FieldIndex index, vector<GF2>& indexBits) {

	GF2 zero = to_GF2(0);
	GF2 one = to_GF2(1);
	if (indexBits.size() < bitsNum) indexBits.resize(bitsNum);

	for (int j=0; j<bitsNum; j++) {
		indexBits[j] = (index&1) ? one : zero;
		index >>= 1;
	}
}
void convertIntToGF2sMSBFirst(const int bitsNum, FieldIndex index, vector<GF2>& indexBits) {

	GF2 zero = to_GF2(0);
	GF2 one = to_GF2(1);
	if (indexBits.size() < bitsNum) indexBits.resize(bitsNum);

	for (int j=0; j<bitsNum; j++) {
		indexBits[bitsNum-1-j] = (index&1) ? one : zero;
		index >>= 1;
	}
}


void convertIntToGF2s(const int bitsNum, FieldIndex index, vector<FieldPoint>& indexBits) {

	FieldPoint zero = constFieldElement(0);
	FieldPoint one = constFieldElement(1);
	if (indexBits.size() < bitsNum) indexBits.resize(bitsNum);

	for (int j=0; j<bitsNum; j++) {
		indexBits[j] = (index&1) ? one : zero;
		index >>= 1;
	}
}

/** 
* Converts a given vector of GF2 elements to a single integer.  
* Note - Order of bits in the vector is MSB...LSB.
*/
int convertGF2sToInt(const vector<GF2>& indexBits) {

	GF2 zero = to_GF2(0);
	GF2 one = to_GF2(1);
	int bitsNum = (int)indexBits.size();
	long res = 0, expo = 1;

	for (int j=0; j<bitsNum; j++) {
		if (indexBits[j] == zero)
			expo *= 2;
		else if (indexBits[j] == one) {
			res += expo;
			expo *= 2;
		}
		else
			assert(false && "cannot convert a field element that is not 0,1 to integer");
	}
	return res;
}
int convertG2sToIntMSBFirst (const vector<GF2>& indexBits){
	GF2 zero = to_GF2(0);
	GF2 one = to_GF2(1);
	int bitsNum = (int)indexBits.size();
	long res = 0, expo = 1;

	for (int j=bitsNum-1; j>=0; j--) {
		if (indexBits[j] == zero)
			expo *= 2;
		else if (indexBits[j] == one) {
			res += expo;
			expo *= 2;
		}
		else
			assert(false && "cannot convert a field element that is not 0,1 to integer");
	}
	return res;

}
int convertGF2sToInt(const vector<FieldPoint>& indexBits) {

	FieldPoint zero = constFieldElement(0);
	FieldPoint one = constFieldElement(1);
	int bitsNum = (int)indexBits.size();
	long res = 0, expo = 1;

	for (int j=0; j<bitsNum; j++) {
		if (indexBits[j] == to_GF2(0))
			expo *= 2;
		else if (indexBits[j] == to_GF2(1)) {
			res += expo;
			expo *= 2;
		}
		else
			assert(false && "cannot convert a field element that is not 0,1 to integer");
	}
	return res;
}

/*****************************************************************************/
/******************* ElementsVector_NTL Implementation ***********************/
/*****************************************************************************/

//ElementsVector_NTL::ElementsVector_NTL() {}
//
//ElementsVector_NTL::ElementsVector_NTL(const int64_t len) {
//	vecRep.SetLength(::Infrastructure::safeConvert(len));
//}
//
//int64_t ElementsVector_NTL::length() const {
//	return vecRep.length();
//}
//
//void ElementsVector_NTL::SetLength(const int64_t len) {
//	vecRep.SetLength(::Infrastructure::safeConvert(len));
//}
//
//FieldPoint& ElementsVector_NTL::operator[](const int64_t& idx) {
//	return vecRep[::Infrastructure::safeConvert(idx)];
//}
//
//const FieldPoint& ElementsVector_NTL::operator[](const int64_t& idx) const {
//	return vecRep[::Infrastructure::safeConvert(idx)];
//}
//
//void ElementsVector_NTL::kill() {
//	vecRep.kill();
//}

/*****************************************************************************/
/****************** ElementsVector_Buckets Implementation ********************/
/*****************************************************************************/

/** Class default constructor. Creates a single bucket of length 0 */
ElementsVector_Buckets::ElementsVector_Buckets() {
	vec_GF2E vec;
	elementsRep.push_back(vec);
	numBuckets = 1;
	size = 0;
}

/** Class constructor. Creates a vector of len elements using as many buckets as needed. */
ElementsVector_Buckets::ElementsVector_Buckets(const int64_t len) {
	
	numBuckets = (int)ceil((double)len / (double)BUCKET_SIZE);
	size = len;

	/** One bucket is enough */
	if (numBuckets <= 1) {
		vec_GF2E vec;
		elementsRep.push_back(vec);
		elementsRep[0].SetLength(::Infrastructure::safeConvert(len));
		return;
	}

	/** Requested length is larger than a bucket size. Allocating several buckets. */	
	for (int i=0; i<numBuckets-1; i++) {	///First n-1 buckets are full
		vec_GF2E currVec;
		elementsRep.push_back(currVec);
		elementsRep[i].SetLength(::Infrastructure::safeConvert(BUCKET_SIZE));
	}

	//Last bucket might not be full.
	long lastBucketSize = ::Infrastructure::safeConvert(len - (numBuckets-1)*BUCKET_SIZE);
	assert(lastBucketSize>0 && lastBucketSize<=BUCKET_SIZE);
	vec_GF2E lastVec;
	elementsRep.push_back(lastVec);
	elementsRep[numBuckets-1].SetLength(lastBucketSize);
}

/** Returns the total number of elements in the vector. */
int64_t ElementsVector_Buckets::length() const {
	return size;
}

/** Sets the vector length. */
void ElementsVector_Buckets::SetLength(const int64_t len) {
	
	///Case 1 - No need for changes.
	if (len == size) 
		return;
	
	///Case 2 - New length is 0. Need one empty bucket.
	if (len==0) {
		for (int i=0; i<numBuckets; i++) elementsRep[i].kill();		///Killing the buckets that won't be needed anymore.
		elementsRep.resize(1);
		elementsRep[0].SetLength(0);
		size = len;
		numBuckets = 1;
		return;
	}

	int newNumBuckets = (int)ceil((double)len / (double)BUCKET_SIZE);

	///Case 3 - Number of buckets stays the same. Use NTL's SetLength on the last bucket.
	if (newNumBuckets == numBuckets) {
		long lastBucketSize = ::Infrastructure::safeConvert(len - (newNumBuckets-1)*BUCKET_SIZE);
		elementsRep[numBuckets-1].SetLength(lastBucketSize);
		size = len;
		return;
	}

	///Case 4 - We're currently using more buckets than needed. Need to deallocate the unncessary buckets.
	if (len < size) {
		assert(newNumBuckets<numBuckets && numBuckets>1);
		for (int i=newNumBuckets; i<=numBuckets-1; i++)		///Killing the buckets that won't be needed anymore.
			elementsRep[i].kill();

		long lastBucketSize = ::Infrastructure::safeConvert(len - (newNumBuckets-1)*BUCKET_SIZE);
		elementsRep[newNumBuckets-1].SetLength(lastBucketSize);
		numBuckets = newNumBuckets;
		size = len;
		return;
	}

	//Case 5 - New length is larger and requires allocation of more buckets.
	elementsRep.resize(newNumBuckets);

	///The newly buckets added (plus the bucket that is currently not full) should be of full.
	for (int i=numBuckets-1; i<newNumBuckets-1; i++)
		elementsRep[i].SetLength(::Infrastructure::safeConvert(BUCKET_SIZE));	

	long lastBucketSize = ::Infrastructure::safeConvert(len - (newNumBuckets-1)*BUCKET_SIZE);
	elementsRep[newNumBuckets-1].SetLength(lastBucketSize);	
	numBuckets = newNumBuckets;
	size = len;
}

/** O(1) access to the vector's elements. */
FieldPoint& ElementsVector_Buckets::operator[](const int64_t& idx) {
	_COMMON_ASSERT((idx >= 0 && idx<size), "vector index out of bounds");
	if (numBuckets == 1)
		return elementsRep[0][::Infrastructure::safeConvert(idx)];

	int relevantBucket = ::Infrastructure::safeConvert(idx / BUCKET_SIZE);
	int relevantIndex = ::Infrastructure::safeConvert(idx % BUCKET_SIZE);
	return elementsRep[relevantBucket][relevantIndex];
}

/** O(1) access to the vector's elements. */
const FieldPoint& ElementsVector_Buckets::operator[](const int64_t& idx) const {
	_COMMON_ASSERT((idx>=0 && idx<size),"vector index out of bounds");
	if (numBuckets == 1)
		return elementsRep[0][::Infrastructure::safeConvert(idx)];

	int relevantBucket = ::Infrastructure::safeConvert(idx / BUCKET_SIZE);
	int relevantIndex = ::Infrastructure::safeConvert(idx % BUCKET_SIZE);
	return elementsRep[relevantBucket][relevantIndex];
}

/** Destroys the buckets, leavng a single empty bucket. */
void ElementsVector_Buckets::kill() {
	for (int i=0; i<numBuckets; i++)
		elementsRep[i].kill();
	
	numBuckets = 1;
	size = 0;
}
} // of namespace details
} // of namespace Algebra
