#ifndef LOGGER_MOJO_H
#define LOGGER_MOJO_H

#define BOOST_LOG_DYN_LINK

#include <cstring>
#include <fstream>

#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/smart_ptr/make_shared_object.hpp>

#include <boost/locale/generator.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/expressions/formatters/date_time.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/sources/channel_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/expressions/formatter.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/formatter_parser.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/locale/generator.hpp>
#include <boost/utility/empty_deleter.hpp>

#include "Utils.h"

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace keywords = boost::log::keywords;

using namespace std;

namespace MOJO 
{
	class Logger 
	{
		public:
			static const string DEFAULT_CHANNEL;

			static void Initializer(string outputFilePrefix);

			static void RegisterChannel(string outputFileName, 
				string channelName);
	};
};

#endif