#include "TaskReporting.hpp"

#include <iostream>
#include <boost/format.hpp>

namespace PCP_Project {

/********************************************************/
/******************* Task reporting *********************/
/********************************************************/

class Task::TaskDedupPayload : public DedupTree::HelpfulPayload {
private:
	std::string name;
	std::string indent;
	double elapsedTime;
public:
	TaskDedupPayload(std::string _name, std::string _indent) : name(_name), indent(_indent), elapsedTime(0) {}

	virtual bool equivalent(const DedupTree::Payload* other) {
		return name == dynamic_cast<const Task::TaskDedupPayload*>(other)->name;
	}

	virtual void begin() {
		DedupTree::HelpfulPayload::begin();
		std::cout << indent << "/ " << name << std::endl;
	}

	virtual void end() {
		assert(elapsedTime>=0);
		std::cout << indent << boost::format("\\ [%.3f sec     ")%elapsedTime << name << "]" << std::endl;
		DedupTree::HelpfulPayload::end();
	}

	virtual void visiblyRepeating() {
		if (numRepeats==1)
			std::cout << indent << "* repeating";
		std::cout << "." << std::flush;
	}

	virtual void visiblyDoneRepeating() {
		if (numRepeats>0)
			std::cout << " " << (numRepeats+1) << " times " << std::endl;
	}

	void setElapsedTime(double elap) {
		elapsedTime = elap;
	}
};

int Task::currentNesting = 0;
std::list<DedupTreePtr> Task::taskStack;

//TODO: Use Logger? --Eran

Task::Task(std::string _name) :
	name(_name),
	isDone(false),
	nesting(currentNesting)
{
	if (taskStack.empty()) {
		DedupTree::Payload* rootPayload = new DedupTree::Payload();
		taskStack.push_back( DedupTreePtr(new DedupTree(rootPayload)) );
	}
	payload = new TaskDedupPayload(name, indent());
	taskTreeNode = taskStack.back()->conceiveChild( payload );
	taskStack.push_back( taskTreeNode );
	currentNesting++;
}

Task::~Task() {
	if (!isDone)
		done();
}

std::string Task::indent() {
	std::string spaces;
	for (int i=0; i<nesting; ++i)
		spaces += "| ";
	return spaces;
}

void Task::done() {
	assert(!isDone);
	assert(taskStack.back()==taskTreeNode);
	taskStack.pop_back();
	isDone = true;
	--currentNesting;
	assert(nesting==currentNesting);
	payload->setElapsedTime(timer.getElapsed());
	taskTreeNode->done();
}
} // namespace PCP_Project 
