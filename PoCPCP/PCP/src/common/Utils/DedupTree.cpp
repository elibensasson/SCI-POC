/** @file
 * Implementation of deduplicated trees.
 */

#include "DedupTree.hpp"

using namespace std;

namespace PCP_Project {

/********************************************************/
/******************* Deduplicated tree ******************/
/********************************************************/


/**
 * Create a new tree node (public constructor).
 */
DedupTree::DedupTree(Payload* _payload) :
		isRooted(true),
		isDone(false),
		parent(NULL),
		prior(NULL),
		payload(_payload)
{
	activate();
}

/**
 * Create a new tree node, with the given parent and prior(private constructor, called by parent.)
 */
DedupTree::DedupTree(Payload* _payload, DedupTree* _parent, DedupTree* _prior) :
		isRooted(false),
		isDone(false),
		parent(_parent),
		prior(_prior),
		payload(_payload)
{
}

DedupTree::~DedupTree() {
	//assert(isDone);
}

/**
 * Called either internally or by parent (depending on how constructed) as soon as callbacks to parent (if any) are allowed
 */
void DedupTree::activate() {
	if (!prior || !payload->equivalent(prior->payload.get()))
		deviateFromPrior();
}

void DedupTree::becomeRooted() {
	if (isRooted)
		return;
	isRooted = true;
	payload->begin();
	for (size_t i=0; i<children.size(); ++i)
		children[i]->becomeRooted();
	if (isDone) {
		payload->end();
		//if (numRepeats)
		//	payload->repeated(numRepeats);
	}
}

/**
 * Handle the consequences of deviating from our current prior (if any).
 */
void DedupTree::deviateFromPrior() {
	//cout << "DEVIATED : " << name << endl;
	assert(!embryo); // If we recently had an embryo, it should have been born during the childDeviatedFromPrior that got us here.
	if (prior) {
		assert(!deeplyEquivalent(prior));
		prior = NULL;
	}

	if (parent) {
		// Tell ancestors we've deviated from our prior.
		// (If father considers us am embryo, he will now be convinced to promote us to its list of born children).
		// This may they deviated from theirs, and so on.
		// They may then call our setPrior() before returning.
		parent->childDeviatedFromPrior(this);
	} else {
		assert(isRooted);
		// an orphan won't get a becomeRooted() callback, so has to take care of itself
		payload->begin();
	}
}

/**
 * Compare payload-equivalence of trees, recurively.
 */
bool DedupTree::deeplyEquivalent(const DedupTree* other) const {
	assert(other);
	if (!payload->equivalent(other->payload.get()))
		return false;
	if (children.size() != other->children.size())
		return false;
	for (size_t i=0; i<children.size(); ++i)
		if (! children[i]->deeplyEquivalent(other->children[i].get()) )
			return false;
	return true;
}

/**
 * Callback from child/embryo, saying it's identical to its prior.
 * If embryo, discard it since it was found to be equivalent to its predecessor sibling.
 * May be called only as a callback by the child/embryo.
 */
void DedupTree::childIdenticalToPrior(DedupTree* _child) {
	assert(_child);
	if (_child==embryo.get()) {
		assert(!children.empty() && _child->deeplyEquivalent(children.back().get()));

		// Inform the prior that it's been duplicated:
		children.back()->payload->foundSiblingDuplicate(_child->payload.get());

		embryo = NULL; // deallocates the embryo and its children, due to reference counting
	} else {
		assert(_child==children.back().get());
		// Youngest child confirms it's identical to its prior. We've been speculatively assuming this, so great, just proceed.
	}
}

/**
 * Callback from child/embryo, saying it's different from its prior.
 * Happens as early as the child/embryo discovers this.
 * If it's am embryo, then promote it to a child.
 * In any case, handle consequences for our own equivalence to our prior.
 * May be called only as a callback by the child/embryo.
 */
void DedupTree::childDeviatedFromPrior(DedupTree* child) {
	assert(child);
	if (child == embryo.get()) {
		// Mazl tov! Promote embryo to a child.

		// The embryo isn't a dup of the the youngest already-born child, so the latter isn't going to be duplicated by siblings anymore.
		if (!children.empty())
			children.back()->payload->doneRepeating();

		// Promote embryo to a child.
		DedupTreePtr newborn = embryo;
		embryo = NULL;
		children.push_back(newborn);
		size_t newbornIdx = children.size()-1;

		// The newborn should stop comparing itself to its youngest sibling, and instead
		// compare itself to the analogous child in our prior.
		if (prior) {
			if (prior->children.size() < newbornIdx+1) {
				// Prior had fewer children, so *we* deviate:
				deviateFromPrior(); // will cause setPrior() calls to us before returning
			} else {
				// Tell newborn about its new prior:
				DedupTree* newbornPrior = prior->children[newbornIdx].get();
				newborn->setPrior(newbornPrior); // may immediate generate more childDeviatedFromPrior() callbacks. Viva la stack!
			}
		}
		if (isRooted)
			newborn->becomeRooted(); // (child may have already become rooted during a callback generated above)
	} else {
		// The (already-born) youngest child says it's different from its prior, so we deviate from our prior too.
		assert(!children.empty() && child==children.back().get());
		deviateFromPrior(); // may cause setPrior() calls to us before returning
	}
}

/**
 * Change the prior of this node (to a more distant cousin, or NULL if we ran out of cousins).
 * Any mismatch will (indirectly) cause a childDeviatedFromPrior callback before returning, and also cause a return of false.
 * @return true if we managed to recursively set the new prior all the way through; false if a mismatch found.
 */
bool DedupTree::setPrior(DedupTree* _prior) {
	assert(!embryo); // this should happen only by a call from a parent invoked (indirectly) by us calling deviateFromPrior
	if (prior == _prior)
		return true;
	prior = _prior;
	if (prior) {
		if (isDone) {
			// New prior for an already-done node, so just do a deep-equivalence test.
			if (!deeplyEquivalent(prior)) {
				deviateFromPrior();
				return false;
			}
		} else {
			// New prior for an under-construction node, so check if OK so far.
			if (!payload->equivalent(prior->payload.get())) { // Is our payload equivalent to the prior's?
				deviateFromPrior();
				return false;
			} else if (children.size() > prior->children.size()) { // Does the clone has enough children?
				deviateFromPrior();
				return false;
			} else {
				for (size_t i=0; i<children.size(); ++i) { // Inform the (known) children of their new prior
					DedupTree* childPrior = prior->children[i].get();
					bool deeplyGood = children[i]->setPrior(childPrior);
					if (!deeplyGood)
						return false; // they didn't like it, abort. Another setPrior() was generated by now with a new prior.
				}
			}
		}
	} else {
		// New prior is NULL. Bother telling only children that are still under construction.
		if (!children.empty() && !(children.back()->isDone))
			children.back()->setPrior(NULL);
	}
	return true;
}

/**
 * Create a new potential child, but don't yet accept it as a permanent child since it may turn out
 * to be identical to its predecessor sibling (if any).
 */
DedupTreePtr DedupTree::conceiveChild(Payload* childPayload) {
	assert(!embryo);
	DedupTree* embryoSibling = NULL;
	if (!children.empty())
			embryoSibling = children.back().get();
	DedupTreePtr newnode(new DedupTree(childPayload, this, embryoSibling));
	embryo = newnode;
	// If clone is null then the next line will cause a a birthEmbryo() callback to us, immediately, from the embryo.
	// In that case the object will move from the 'embryo' field to the list of children.
	// Otherwise we'll eventually get either a birthEmbyo() or an abotyEmryo() callback.
	newnode->activate();
	return newnode;
}

/**
 * Tell this tree that it will not have any more children added to it.
 */
void DedupTree::done() {
	assert(!embryo); // did you forget to call done() on some descendant node?
	assert(!isDone);
	isDone = true;

	// The youngest child (if any) is not going to be duplicated by siblings anymore
	if (!children.empty())
		children.back()->payload->doneRepeating();

	if (prior) {
		assert(parent);
		if (children.size()==prior->children.size()) {
			// Made it so far with non-NULL prior, so deeply equivalent to prior.
			assert(prior && deeplyEquivalent(prior));

			// Tell parent we're identical to our prior.
			// (If father considers us an embryo, he'll abort us now.)
			parent->childIdenticalToPrior(this);
			return; // return immediately, the parent may have just deleted us by reducing our ref count to 0
		} else {
			deviateFromPrior();
		}
	}
	if (isRooted) {
		payload->end();
	}
}




DedupTree::HelpfulPayload::HelpfulPayload() :
		visible(false), numRepeats(0) {
}

void DedupTree::HelpfulPayload::begin() {
	visible = true;
}

void DedupTree::HelpfulPayload::end() {
	if (numRepeats > 0)
		visiblyDoneRepeating();
}

void DedupTree::HelpfulPayload::foundSiblingDuplicate(Payload* dup) {
	++numRepeats;
	if (visible)
		visiblyRepeating();
}

void DedupTree::HelpfulPayload::doneRepeating() {
	if (visible)
		visiblyDoneRepeating();
}

} // of namespace PCP_Project
