#pragma once

namespace precice {
namespace utils {

/// Install signal handler to print timings and cleanup on simulation abort
/** When precice stops abruptly, e.g. an external solver crashes, the
    SolverInterfaceImpl destructor is never called. Since we still want
    to print the timings, we install the signal handler here.
*/
void installSignals();

/// To be called when everything is lost and we try to terminate gracefully.
/**
 * - Writes out the event timings
 * - Deletes the address files
 */
void terminationSignalHandler(int signal);

}
}
