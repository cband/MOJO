
#include "Logger.h"

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace sinks = boost::log::sinks;
namespace keywords = boost::log::keywords;
namespace expr = boost::log::expressions;
using namespace boost::posix_time;

typedef sinks::synchronous_sink< sinks::text_ostream_backend > text_sink;

namespace MOJO 
{
	const string Logger::DEFAULT_CHANNEL = "Main";

	void Logger::Initializer( string outPrefix ) 
	{			
		//// Status output link; output to Console and .Status file
		if (!outPrefix.empty() && outPrefix.length() > 1 
			&& outPrefix[outPrefix.length() - 1] == '.')
			outPrefix.erase(outPrefix.length() - 1);
		RegisterChannel(outPrefix + "." + Utils::GetTimeStamp() + ".log", "Main");
	}

	void Logger::RegisterChannel(string outputFileName, string channelName) 
	{
		//Create a sink for log file
		boost::shared_ptr< sinks::text_ostream_backend > backend = 
			boost::make_shared< sinks::text_ostream_backend >();
		backend->add_stream(boost::shared_ptr< std::ostream >
			(new std::ofstream(outputFileName)));
		if (channelName == "Main")
			backend->add_stream(boost::shared_ptr< std::ostream >
			(&std::clog, boost::empty_deleter()));

		backend->auto_flush(true);
		boost::shared_ptr< text_sink > sink(new text_sink(backend));
		sink->set_formatter(
			expr::stream << "["
			<< expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "%Y-%m-%d %H:%M:%S")
			<< "] " << expr::smessage
			);
		sink->set_filter(expr::attr< std::string >("Channel") == channelName);
		logging::core::get()->add_sink(sink);
		logging::add_common_attributes();

		//Create a sink for cmds executed in this channel
		boost::shared_ptr< sinks::text_ostream_backend > backend2 = 
			boost::make_shared< sinks::text_ostream_backend >();
		backend2->add_stream(boost::shared_ptr< std::ostream >
			(new std::ofstream(outputFileName + ".cmds")));
		backend2->auto_flush(true);
		boost::shared_ptr< text_sink > sink2(new text_sink(backend2));
		sink2->set_formatter(
			expr::stream << "["
			<< expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "%Y-%m-%d %H:%M:%S")
			<< "] " << endl << expr::smessage
			);
		sink2->set_filter(
			expr::attr< std::string >("Channel") == channelName + ".cmds");
		logging::core::get()->add_sink(sink2);
		logging::add_common_attributes();
	}
}



