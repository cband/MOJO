#ifndef UTILS_H
#define UTILS_H

#pragma once

#define BOOST_LOG_DYN_LINK

#include <vector>
#include <cstring>
#include <sstream>
#include <iostream>

//File system api
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/thread.hpp>
#include <boost/unordered_map.hpp>
//Boost Logging api
#include <boost/lexical_cast.hpp>
#include <boost/log/core.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/sources/basic_logger.hpp>
#include <boost/log/sources/global_logger_storage.hpp>
#include <boost/log/sources/channel_feature.hpp>
#include <boost/log/sources/channel_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/formatter_parser.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/locale/generator.hpp>
#include <boost/utility/empty_deleter.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>        

//Includes to get current time
#include <boost/date_time/local_time/local_time.hpp>
#include <boost/date_time/local_time/local_date_time.hpp>
#include <boost/date_time/posix_time/ptime.hpp>
#include <boost/date_time/microsec_time_clock.hpp>
#include <boost/date_time/time.hpp>

using namespace std;
namespace src = boost::log::sources;
namespace keywords = boost::log::keywords;
namespace b_pt = boost::posix_time;
namespace b_lt = boost::local_time;

#define EXIT_ON_FAIL true
#define DONOT_EXIT_ON_FAIL false

namespace MOJO 
{
	class SystemCall 
	{
		public:
			string systemcall;
			string output;
			int exit_code;

			SystemCall() : exit_code(0) {}
			
			SystemCall(string cmd, string out, int code) : 
				systemcall(cmd), output(out), exit_code(code) {}
			
			friend ostream& operator<<(ostream &outStream, const SystemCall &call) 
			{
				//Remove newline character
				string s = call.output;
				if (!s.empty() && s.length() > 1 && s[s.length() - 1] == '\n')
					s.erase(s.length() - 1);

				outStream << "\tCommand: " << call.systemcall << endl;
				outStream << "\tOutput: " << s << endl;
				outStream << "\tExit Code: " << call.exit_code << endl;
				return outStream;
			}
	};

	class SystemInfo 
	{
		public:
			static int GetSystemCpuCount();
	
			static double GetTotalSystemMemory();
			
			static std::pair<double, double> GetTotalMemoryUsage(int sid, int pid = 0);
			
			static double GetTotalFreeMemory();
			
			SystemInfo();
	};

	class ComputePerTask 
	{
		public:
			static int MAX_CPU, MAX_MEM;
			int numSplits, numCoresPerSplit, numResidualCores;
			
			ComputePerTask() : numSplits(1), numCoresPerSplit(1),
				numResidualCores(0) {}
			
			ComputePerTask(int num_splits, int num_cores_per_split) :
				numSplits(num_splits), numCoresPerSplit(num_cores_per_split), 
				numResidualCores(0) {}
		
			static ComputePerTask CalculateComputePerTask(int max_splits,
				int min_cores_per_split, int max_mem_per_split);

			bool operator==(const ComputePerTask & cpt) const
			{
				return this->numSplits == cpt.numSplits &&
					this->numCoresPerSplit == cpt.numCoresPerSplit;
			};
	};

	class Utils
	{
		public:
			static char* StringCopy( const char* src );
			
			static char* ToLower( char* src );
			
			static char* ReverseString(char* str);
			
			static char* TrimString( string line, char chr=0 );
			
			//template<typename T = string> 
			//static vector<T> SplitToVector( string line, string chr ) ;

			static string ReverseComplement(string s);

			static vector<string> SplitToVector(string line, string chr);
			
			static SystemCall ExecuteCommand(const char* cmd, 
				string channelName = "",
				bool doPrintErrors = true, 
				bool exitOnFail = DONOT_EXIT_ON_FAIL);
			
			static SystemCall ExecuteCommand(char* cmd, string channelName = "",
				bool doPrintErrors = true, bool exitOnFail = DONOT_EXIT_ON_FAIL);
			
			static bool FileExists(string file, bool logMsg = false);
			
			static bool DeleteFiles(vector<string> files);
			
			static bool DeleteFile(string file);
			
			static int LineCount(string file);

			static vector<string> SortListOfFilesBySize(vector<string> files, 
				bool ascending = true);

			//system tools
			static string GetTimeStamp(string dtFormat = "%Y.%m.%d.%H:%M:%S");
	};
};
#endif
