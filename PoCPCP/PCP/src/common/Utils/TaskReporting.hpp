#include "Timing.hpp"
#include "common/Utils/DedupTree.hpp"

#include <memory>
#include <list>
#include <string>

#ifndef  __TaskReporting_HPP
#define  __TaskReporting_HPP

namespace PCP_Project {
/********************************************************/
/******************* Task reporting *********************/
/********************************************************/

class TaskDedupPayload;
class DedupTree;

class Task {
private:
	class TaskDedupPayload;
private:
	std::string name;
	bool isDone;
	Timer timer;
	int nesting; // nesting level of this task
	std::shared_ptr<DedupTree> taskTreeNode;
	TaskDedupPayload* payload;
	std::string indent();
private:
	static int currentNesting; // nesting level of next task to start
	static std::list<std::shared_ptr<DedupTree>> taskStack;
public:
	/** Report a new task starting */
	Task(std::string name);

	~Task();

	/** Report the task as done. If not called explicitly, will be done automatically when the Task object is destroyed. */
	void done();
};

/** Macro to start a new task.
 * Usage pattern is: { TASK("do stuff); stuff(); }.
 * Don't put multiple TASK() invocation in the same curlies-scope, as the identifier _currentTask will collide.
 */
#define TASK(_name)   Task _currentTask(_name) //TODO add debug level to switch this off

} // namespace PCP_Project 

#endif   // __TaskReporting_HPP
