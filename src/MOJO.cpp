#define BOOST_LOG_DYN_LINK

#include <fstream>
#include <iostream>
#include <sys/types.h>
#include <unistd.h>

pid_t getpid(void);
pid_t getppid(void);

#include <boost/unordered_map.hpp>
#include <boost/thread.hpp>
#include <boost/log/sources/global_logger_storage.hpp>
#include <boost/format.hpp>
#include <boost/timer/timer.hpp>
#include <boost/chrono/chrono.hpp>

#include "MOJOConfig.h"
#include "Utils.h"
#include "Logger.h"
#include "Config.h"
#include "GeneModel.h"
#include "DiscordantReadFinder.h"
#include "DiscordantClusterFinder.h"
#include "JunctionAligner.h"
#include "JunctionFilter.h"
#include "FusionCompiler.h"
#include "JunctionFilter.h"

using namespace std;
using namespace MOJO;

namespace MOJO {
	static src::channel_logger<> mainLogger(keywords::channel = "Main");
	Config a_tmp; //dummy call - to initialize the static variable;
	SystemCall c;
}

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace keywords = boost::log::keywords;
namespace sinks = boost::log::sinks;
namespace expr = boost::log::expressions;

static int pid = -1, sid = -1, pollingPauseMillisecs = 1000;
static double maxMemoryUsage = 0;

void PollMemory(bool pollOnce = false) 
{
#ifdef __APPLE__
	return;
#endif

	double sysmem = SystemInfo::GetTotalSystemMemory();
	while (1) {
		if ( !pollOnce )
			boost::this_thread::interruption_point();	// to allow for termination
		boost::this_thread::sleep(
			boost::posix_time::milliseconds(pollingPauseMillisecs));

		auto mem = SystemInfo::GetTotalMemoryUsage(sid, pid);
		if (mem.first > maxMemoryUsage) 
			maxMemoryUsage = mem.first;
		
		if ( (mem.first/sysmem) > 0.85 || pollOnce ) {
			double free = SystemInfo::GetTotalFreeMemory();
			BOOST_LOG(mainLogger) << "WARNING: MOJO memory usage at " 
				<< boost::format("%.2f") % (mem.first*100 / (sysmem)) << "% "
				<< "(MOJO-exec:" << boost::format("%.2f") % mem.second << "gb, "
				<< "MOJO-all:" << boost::format("%.2f") % mem.first << "gb, "
				<< "Free: " << boost::format("%.2f") % free << "gb)";
		}
		if (pollOnce)
			return;
	}
}
int main( int argc, char *argv[] )
{
	boost::timer::cpu_timer timer;
	boost::chrono::system_clock::time_point start = 
		boost::chrono::system_clock::now();

	pid = (int) getppid();
	sid = (int)getsid(pid);

	try {
		if (!Config::MOJORunConf.LoadConfiguration(argc, argv)) {
			std::cerr << "Error occurred while loading configuration file. " 
				<< " Exiting." << endl;
			exit(1);
		}
	}
	catch (std::exception &e) {
		std::cerr << "Exception occurred while loading configuration: " 
			<< e.what() << endl;
		exit(1);
	}
	Logger::Initializer(Config::MOJORunConf.outputDir + 
		Config::MOJORunConf.sampleName);
	boost::thread *memThread = new boost::thread(PollMemory, false);

	BOOST_LOG(mainLogger) << "MOJO v" << MOJO::version;
	BOOST_LOG(mainLogger) << "---------";
	BOOST_LOG(mainLogger) << "Config file: " << Config::MOJORunConf.configFile;
	BOOST_LOG(mainLogger) << "Sample: " << Config::MOJORunConf.sampleName;
	BOOST_LOG(mainLogger) << "Output: " << Config::MOJORunConf.workingDir;
	BOOST_LOG(mainLogger) << "Reference: " << Config::MOJORunConf.MOJOReferenceDir;

	if (!Config::MOJORunConf.PrerunFastqCheck()) {
		BOOST_LOG(mainLogger) << "Error: sample input not properly configured ";
		exit(1);
	}

	BOOST_LOG(mainLogger) << "MOJO configuration verified ";
	BOOST_LOG(mainLogger) << "Starting MOJO run with " 
		<< Config::MOJORunConf.maxCores << " cpus and " 
		<< Config::MOJORunConf.maxMem << "gb";

	ComputePerTask::MAX_MEM = Config::MOJORunConf.maxMem;
	ComputePerTask::MAX_CPU = Config::MOJORunConf.maxCores;

	PollMemory(true);
	GeneModel::gm.LoadGeneModel();
	PollMemory(true);
	DiscordantReadFinder::Run();
	PollMemory(true);
	auto clusters = DiscordantClusterFinder::LoadDiscordantClusters();
	PollMemory(true);

	JunctionAligner::Run(clusters);
	PollMemory(true);
	JunctionFilter::Run(clusters);
	PollMemory(true);
	FusionCompiler::Run(clusters);
	PollMemory(true);
	
	Config::MOJORunConf.FinalCleanup();

	boost::this_thread::sleep(boost::posix_time::milliseconds(10000));
	memThread->interrupt();
	BOOST_LOG(mainLogger) << "Run Successfully Completed! ";
	boost::chrono::duration<double> sec = 
		boost::chrono::system_clock::now() - start;
	int s = (int) sec.count();
	int hours = s / 3600;
	int mins = (s % 3600) / 60;
	int secs = (s - (hours * 3600) - (mins * 60));
	BOOST_LOG(mainLogger) << "Time Elapsed: "
		<< (boost::format("%02d:%02d:%02d") % hours % mins % secs).str()  
		<< " (approx. peak memory usage: "
		<< (boost::format("%.2f") % maxMemoryUsage).str() << "GB)";
	return 0;
}