#include "SignalHandler.hpp"
#include <csignal>
#include "utils/EventUtils.hpp"
#include <boost/filesystem.hpp>


#ifndef SIGXCPU
#define SIGXCPU 24 /* exceeded CPU time limit */
#endif

namespace precice {
namespace utils {

void installSignals()
{
  /* When precice stops abruptly, e.g. an external solver crashes, the
     SolverInterfaceImpl destructor is never called. Since we still want
     to print the timings, we install the signal handler here. */
  // Disable SIGSEGV handler, because we don't want to interfere with crash backtrace.
  // signal(SIGSEGV, precice::utils::terminationSignalHandler);
  signal(SIGABRT, terminationSignalHandler);
  signal(SIGTERM, terminationSignalHandler);
  // SIGXCPU is emitted when the job is killed due to walltime limit on SuperMUC
  signal(SIGXCPU, terminationSignalHandler);
  // signal(SIGINT,  precice::utils::terminationSignalHandler);
}

void terminationSignalHandler(int signal)
{
  // Print the events statistics
  precice::utils::EventRegistry::instance().signal_handler(signal);

  // Remove connection info files. This is just a guess and will only
  // work if the address directory is ".".
  boost::filesystem::remove_all("precice-run");
}

}
}
