#ifdef _MSC_VER // Visual C++
/* Workeraund, windows build doesnt generate mach_desch.h automatically
 * But Linux build does.
 * What we(Michael & Shaul) do is keep the windows file in repository
 * and use this ifdef macro to decide wether to use the constant file
 * (in the case of MS), or the automatically generated one (otherwise)
 */
#include <NTL/mach_desc_windows_const.h>
#else //Not Visual C++
#include <NTL/mach_desc_auto.h>
#endif //Not Visual C++
