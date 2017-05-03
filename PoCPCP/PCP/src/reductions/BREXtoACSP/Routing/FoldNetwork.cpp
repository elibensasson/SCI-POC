#include "FoldNetwork.hpp"

namespace PCP_Project{
namespace BREXtoACSP{

typedef RoutingNetwork::layerID_t layerID_t;
typedef RoutingNetwork::labelID_t labelID_t;
typedef RoutingNetwork::dataID_t dataID_t;

dataID_t foldLeftNetwork::getDataID(const layerID_t& l, const labelID_t& w) const{
	dataID_t mask = (1<<l) - 1; //the mask in binary is just 'l' least segnificant bit set
	mask = mask<<(k_-l); //the mask is now al the bits from $k-l-1$ to $k-1$ set

	//calculate the bits that were falded
	//in this layer (relativly to first layer)
	labelID_t foldedData = (w & mask)>>(k_-l);

	//calculate the reverse of the folded data
	labelID_t revData = reverseBits(foldedData)>>(8*sizeof(labelID_t) - l);

	//calculate mask for middle bit:
	//if $k$ is odd, there is a bit in the middle that
	//should not change, without a special
	//tritment it will be allwayes xored with itself, so
	//it will be alwayes cleared.
	//We calculate a mask that sets it if needed
	const dataID_t middleMask = (k_%2)<<(k_/2);
	return (w ^ revData) | (w & middleMask) ;
}

dataID_t foldRightNetwork::getDataID(const layerID_t& l, const labelID_t& w) const{
	dataID_t mask = (1<<l) - 1; //the mask in binary is just 'l' least segnificant bit set
	mask = mask<<(k_/2 - l); //the mask is 'l' set bits, just lower that the center of a 'k' length word

	//calculate the bits that were falded
	//in this layer (relativly to first layer)
	labelID_t foldedData = (w & mask)>>(k_/2 - l);

	//calculate the reverse of the folded data
	labelID_t revData = reverseBits(foldedData)>>(8*sizeof(labelID_t) - l - (k_/2 + k_%2));

	//calculate mask for middle bit:
	//if $k$ is odd, there is a bit in the middle that
	//should not change, without a special
	//tritment it will be allwayes xored with itself, so
	//it will be alwayes cleared.
	//We calculate a mask that sets it if needed
	const dataID_t middleMask = (k_%2)<<(k_/2);
	return (w ^ revData) | (w & middleMask) ;
}


} //namespace BREXtoACSP
} //namespace PCP_Project
