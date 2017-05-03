/** @file
 * Class declarations of deduplicated trees.
 */

#ifndef DEDUPTREE_HPP_
#define DEDUPTREE_HPP_

#include <memory>
#include <vector>
#include <cassert>

namespace PCP_Project {

/********************************************************/
/******************* Deduplicated tree ******************/
/********************************************************/

class DedupTree;
typedef std::shared_ptr<DedupTree> DedupTreePtr;

/**
 * A deduplicated-tree data structure.
 * Intended to aid concise logging of nested, repeatitive operations.
 * Within each node, the children are ordered. Duplicate subtrees rooted in consecutive siblings are culled.
 * The nodes contain Payload objects, whose equivalent() member defines what it means to be equivalent
 * (it is applied recursively, in a deep comparison).
 * The tree must be filled in in depth-first order.
 * Each node's Payload object is informed when the node is "entered", i.e., ascertain to belong to the tree
 * rather than being culled. The Payload is also informed when the node "exists", i.e., when the tree-building
 * exits a committed node (this can happen after the node is both committed and done(), whichever is later, if any).
 * Nodes are also informed when they are duplicated due to culling. This happens only when duplicated by
 * siblings; if a whole subtree is duplicated, only the subtree's root is informed.
 */
class DedupTree {
public:
	/**
	 * Base class for payload object attached to DedupTree nodes.
	 * In the following, "visible" means the node has been ascertained to be non-culled.
	 */
	class Payload {
	public:
		/** Are two nodes' payloads equivalent? */
		virtual bool equivalent(const Payload* other) { return false; }

		/** The node has just become visible. */
		virtual void begin() {}

		/** The node has become visible and fully-constructed (whichever happens last, if any) */
		virtual void end() {}

		/**
		 * This node has just been duplicated by a sibling with the given payload.
		 * The node may be visible or not.
		 */
		virtual void foundSiblingDuplicate(Payload* dup) {}

		/**
		 * Called when it's known that the node will not be further duplicated by siblings.
		 * Called exactly once per node.
		 * The node may be visible or not.
		 */
		virtual void doneRepeating() {}

		virtual ~Payload() {}

	};

	/**
	 * Slightly more helpful payload, that keeps track of visibility and the number
	 * of duplicates, and gives events at the right times for progress-of-duplication events.
	 */
	class HelpfulPayload : public Payload {
	protected:
		bool visible;   /** Is this node currently visible? */
		int numRepeats; /** Number of times this node was duplicated (excluding initial) */

	public:
		HelpfulPayload();

		virtual void begin(); //override
		virtual void end(); //override
		virtual void foundSiblingDuplicate(Payload* dup); //override
		virtual void doneRepeating(); //override

	protected:
		/**
		 * Called whenever a duplicate sibling of an already-visible node is found.
		 */
		virtual void visiblyRepeating() {}

		/**
		 * Called as soon as the node is done repeating and is visible, whichever happens last (if any).
		 */
		virtual void visiblyDoneRepeating() {}
	};


private:
	std::vector<DedupTreePtr> children;
	DedupTreePtr embryo;

	bool isRooted;  /** Attached to the tree's root by a path of non-culled nodes */
	bool isDone;    /** Fully populated (i.e., done() was called) */

	DedupTree* parent;   /** Parent node to inform, by callbacks, of whether to cull this node */
	DedupTree* prior;    /** Prior node to which this node should be compared to. (Initially youngest sibling, later some cousins.) */

	std::unique_ptr<Payload> payload;   /** Payload of this node, informed of its status changes */

private: // to be called only by self
	void deviateFromPrior();

private: // to be called by self or parent
	bool deeplyEquivalent(const DedupTree* other) const;

private: // to be called only by parent
	DedupTree(Payload* payload, DedupTree* _parent, DedupTree* _prior);
	void activate();
	void becomeRooted();
	bool setPrior(DedupTree* _prior);

private: // callbacks to be called only by self's embryos and children
	void childDeviatedFromPrior(DedupTree* child);
	void childIdenticalToPrior(DedupTree* child);

public:
	DedupTree(Payload* payload);
	~DedupTree();
	DedupTreePtr conceiveChild(Payload* childPayload);
	void done();
};

} // of namespace PCP_Project

#endif /* DEDUPTREE_HPP_ */
