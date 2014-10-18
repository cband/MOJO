
#include <sys/wait.h>
#include <sys/types.h>
#include <cstdio>

#include "Utils.h"

namespace MOJO 
{

	BOOST_LOG_INLINE_GLOBAL_LOGGER_CTOR_ARGS(logger, src::channel_logger_mt< >, 
		(keywords::channel = "Main"))
	int ComputePerTask::MAX_CPU = 0;
	int ComputePerTask::MAX_MEM = 0;

	char* Utils::StringCopy( const char* src ) 
	{
		auto line_len = strlen( src );
		char* linedata = new char[ line_len + 1 ];
		
		for ( size_t j = 0; j < line_len; j++ )
			linedata[ j ] = src[ j ];
		linedata [ line_len ] = '\0';

		return linedata;	////////WARNING:  Consumer is expected to delete this
	}

	char* Utils::ToLower( char* src ) 
	{
		for ( size_t i=0; i < strlen( src ); i++ )
			if (src[ i ] >= 0x41 && src[ i ] <= 0x5A)
				src[ i ] = src[ i ] + 0x20;
		return src;
	}

	char * Utils::ReverseString(char* str)
	{
		char *left = str;
		char *right = left + strlen(str) - 1;
		char tmp;

		while (left < right) {
			tmp = *left;
			*(left++) = *right;
			*(right--) = tmp;
		}
		return str;
	}

	char* Utils::TrimString( string line, char chr ) 
	{
		const char* ldata = line.c_str();
		int len = 1, idx = 0;
		for( string::size_type i = 0; i < line.length(); i++ )
			if ( ldata[i] != '\r' && ldata[i] != '\n' && ldata[i] != chr )  
				len++;
 
		char* buf = new char[ len + 1];
		for (string::size_type i = 0; i < line.length(); i++)
			if ( ldata[i] != '\r' && ldata[i] != '\n' && ldata[i] != chr )  
				buf[idx++] = ldata[i];
		buf[ idx ] = '\0';
		return buf;		////////WARNING:  Make sure the consumer deletes this after use
	}

//	template<typename T> 
//	vector<T> Utils::SplitToVector(string lineStr, string delimStr) {

	vector<string> Utils::SplitToVector(string lineStr, string delimStr) 
	{
		//Note: strtok_r treats consecutive delimitters as one; so cannot be used;
		vector<string> splits;
		int carat_ = 0;
		const char *ldata = lineStr.c_str(), *ddata = delimStr.c_str();
		int ldata_len = (int)strlen(ldata), ddata_len = (int)strlen(ddata);
		for( int i = 0; i < ldata_len; i++ ) {
			for (int j = 0; j < ddata_len; j++) {
				if (ldata[i] == ddata[j]) {
					if ( i == carat_ )
						splits.push_back("");
					else
						splits.push_back(std::string(ldata + carat_, (i - carat_)));
					carat_ = i+1;
				}
			}
		}
		if (strlen(ldata) == carat_ )
			splits.push_back("");
		else if (strlen(ldata) > carat_) 
			splits.push_back(std::string(ldata + carat_, (strlen(ldata) - carat_)));
		return splits;
	}

	string Utils::ReverseComplement(string s) 
	{
		if (s.length() == 0)
			return s;

		boost::unordered_map<char, char> comp({ { 'A', 'T' }, { 'T', 'A' },
		{ 'C', 'G' }, { 'G', 'C' }, { 'N', 'N' }, { 'a', 't' }, { 't', 'a' },
		{ 'c', 'g' }, { 'g', 'c' }, { 'n', 'n' } });

		const char* linedata = s.c_str();
		char* new_s = new char[s.length() + 1];

		for (int i = ((int) s.length()) - 1; i >= 0; i--) {
			new_s[s.length() - i - 1] = comp[linedata[i]];
		}
		new_s[s.length()] = '\0';

		string sequence(new_s);
		delete new_s;
		return sequence;
	}

	bool Utils::FileExists(string file, bool logMsg ) 
	{
		bool exists = boost::filesystem::exists(file);
		if (logMsg && !exists) {
			BOOST_LOG_CHANNEL(logger::get(), "Main")
				<< "File does not exist: " << file;
		}
		return exists;
	}

	bool Utils::DeleteFiles(vector<string> files) 
	{
		bool success = true;
		for (auto str : files) 
			if (!Utils::DeleteFile(str)) success = false;
		return success;
	}

	bool Utils::DeleteFile(string file) { 
		string cmdStr = ("rm -rf " + file);
		auto result = Utils::ExecuteCommand(cmdStr.c_str(), "Main", true);
		if (result.exit_code == 0)
			return true;
		return false;
	}

	int Utils::LineCount(string file) 
	{
		char buf[1000];
		sprintf(buf, "wc -l %s | cut -f 1 -d' '", file.c_str());
		SystemCall call = Utils::ExecuteCommand(buf);
		if (call.exit_code != 0 || call.output == "") {
			return 0;
		}
		return stoi(call.output);
	}

	SystemCall Utils::ExecuteCommand(const char* cmd, string channelName,
		bool doPrintErrors, bool exitOnFail) 
	{
		char *buf = Utils::StringCopy(cmd);
		SystemCall call = ExecuteCommand(buf, channelName, doPrintErrors);
		delete[] buf;
		return call;
	}

	SystemCall Utils::ExecuteCommand(char* cmd, string channelName, 
		bool doPrintErrors, bool exitOnFail ) 
	{
		/// SEE: http://stackoverflow.com/questions/15058876/how-to-capture-the-exit-code-and-stderr-of-the-command-that-is-run-in-c
		/// To capture stderr, use three pipes, see: http://jineshkj.wordpress.com/2006/12/22/how-to-capture-stdin-stdout-and-stderr-of-child-program/
		char bufcmd[50000];
		sprintf(bufcmd, 
			"exec bash -c \"set -e; set -o pipefail; (%s) 2>&1\"", cmd);

	#ifdef _WIN32
		FILE* pipe = _popen(bufcmd, "r");
	#else
		FILE* pipe = popen(bufcmd, "r");
	#endif
		if (!pipe) {
			BOOST_LOG_CHANNEL(logger::get(), "Main") 
				<< "Error: cannot execute shell command (popen failed)" 
				<< endl << "Command: " << bufcmd;
			exit(1);
		}
		char buffer[128];
		stringstream ss;
		while(!feof(pipe)) {
			if (fgets(buffer, 128, pipe) != NULL)
				ss << buffer;
		}
		int status = -1, exit_code = -1;
	#ifdef _WIN32
		_pclose(pipe);
		return SystemCall(std::string(bufcmd), ss.str(), exit_code);
	#else
		status = pclose(pipe);
		exit_code = WEXITSTATUS(status);
		SystemCall call(std::string(bufcmd), ss.str(), exit_code);
		if ( exit_code != 0 && doPrintErrors ) {
			//write to both "Thread" log and "Main" log
			BOOST_LOG_CHANNEL( logger::get(), "Main") 
				<< "Execution failed: "<< endl << call;
			if ( channelName != "Main")
				BOOST_LOG_CHANNEL( logger::get(), channelName) 
					<< "Execution failed: "<< endl << call;
		}
		else if ( channelName != "" ) {
			BOOST_LOG_CHANNEL( logger::get(), channelName + ".cmds") << endl<< call;
		}

		if ( exit_code != 0 && exitOnFail ) {
			BOOST_LOG_CHANNEL(logger::get(), "Main") << "Failure occurred. Exiting";
			exit(1);
		}

		//DEBUG code
		//if (WIFEXITED(status)) printf("normal termination, exit status = %d, %d\n", WEXITSTATUS(status), status);
		//else if (WIFSIGNALED(status)) printf("abnormal termination, exit status = %d%s, %d\n", WTERMSIG(status), status);
		//else if (WIFSTOPPED(status)) printf("execution stopped, signal number = %d\n", WSTOPSIG(status));

		return call;
	#endif
	}

	string Utils::GetTimeStamp(string dtFormat) {
		//Alternate implementation;
		//#include <time.h>
		//size_t strftime(char *s, size_t maxsize, const char *format, const struct tm *timeptr);
		//See: http://uw714doc.sco.com/en/man/html.3C/strftime.3C.html
	
		//http://stackoverflow.com/questions/7589170/boost-formatting-sub-second-precision-time-with-a-time-stamp
		std::stringstream   ss;
		b_pt::ptime  now = b_pt::second_clock::local_time();
		b_pt::ptime  utc = b_pt::second_clock::universal_time();
		ss << (now - utc);

		b_lt::time_zone_ptr zone(new b_lt::posix_time_zone(ss.str().c_str()));
		b_lt::local_date_time  ldt = b_lt::local_microsec_clock::local_time(zone);
		b_lt::local_time_facet* output_facet = new b_lt::local_time_facet();
		ss.imbue(std::locale(std::locale::classic(), output_facet));
		output_facet->format(dtFormat.c_str());
		ss.str("");
		ss << ldt;
		return ss.str();
	}

	vector<string> Utils::SortListOfFilesBySize(vector<string> files, 
		bool ascending)
	{
		vector<string> sortedFiles;
		std::map<int, string> fileSizes;
		for (auto file: files) {
			int sz = (int) (boost::filesystem::file_size(file));
			while (fileSizes.find(sz) != fileSizes.end())
				sz++;
			fileSizes[sz] = file;
		}
		if (ascending) {
			for (auto iter : fileSizes)
				sortedFiles.push_back(iter.second);
		}
		else {
			for (auto iter = fileSizes.rbegin(); iter != fileSizes.rend(); ++iter)
				sortedFiles.push_back((*iter).second);
		}
		return sortedFiles;
	}

	//	grep -i memtotal /proc/meminfo	// "MemTotal:       16535760 kB" in kb
	//	dmesg | grep Memory				// "Memory: 16525476k/16908288k available ...." in kb
	//	vmstat -s | grep 'total memory'	// "     16535760  total memory" in kb
	//	free | grep -i mem				// "Mem:      16535760   14168300....." in kb 
	double SystemInfo::GetTotalSystemMemory() {
		double mem = 0;
		try {
			boost::match_results<std::string::const_iterator> results;
			SystemCall call = Utils::ExecuteCommand(
				"grep -i memtotal /proc/meminfo", "", false);
			boost::regex_match(call.output, results,
				boost::regex("MemTotal:[\\s|\\t]+(\\d+)\\s.*"));
			try{
				mem = stof(results[1].str()) / (1024 * 1024);
			}
			catch (exception &e){
			}
			if (mem > 0) return mem;

			call = Utils::ExecuteCommand("dmesg | grep Memory", "", false);
			boost::regex_match(call.output, results, 
				boost::regex(".*Memory:\\s+(\\d+)\\w.*"));
			try{
				mem = stof(results[1].str()) / (1024 * 1024);
			}
			catch (exception &e){
			}
			if (mem > 0) return mem;

			call = Utils::ExecuteCommand("vmstat -s | grep 'total memory'", 
				"", false);
			boost::regex_match(call.output, results,
				boost::regex("[\\s|\\t]+(\\d+)[\\s|\\t]+.*"));
			try {
			mem = stof(results[1].str()) / (1024 * 1024);
			}
			catch (exception &e){
			}
			
			if (mem > 0) return mem;

			call = Utils::ExecuteCommand("free | grep -i mem", "", false);
			boost::regex_match(call.output, results,
				boost::regex("Mem:[\\s|\\t]+(\\d+)[\\s|\\t]+.*"));
			try{
				mem = stof(results[1].str()) / (1024 * 1024);
			}
			catch (exception &e){
			}
			if (mem > 0) return mem;
		}
		catch (std::exception const& e) {
			BOOST_LOG_CHANNEL(logger::get(), "Main") 
				<< "Warning: Unable to get system memory. "  << e.what() << endl;
		}
		return mem;
	}

	int SystemInfo::GetSystemCpuCount() {
		int cpus = 0;
		try {
			cpus = boost::thread::hardware_concurrency();
			if (cpus > 0) return cpus;

			boost::match_results<std::string::const_iterator> results;
			SystemCall call = Utils::ExecuteCommand(
				"grep --count processor /proc/cpuinfo", "", false);
			try{
			cpus = stoi(call.output);
			}
			catch (exception &e){
			}
			if (cpus > 0) return cpus;

			call = Utils::ExecuteCommand("nproc", "", false);
			try{
				cpus = stoi(call.output);
			}
			catch (exception &e){
			}
			if (cpus > 0) return cpus;
		}
		catch (std::exception &e) {
			BOOST_LOG_CHANNEL(logger::get(), "Main") << "Warning: Unable to get "
				<< "cpu count from the system: "  << e.what() << endl;
		}
		return cpus;
	}

	double SystemInfo::GetTotalFreeMemory() {
		double mem = 0;
		try {
			boost::match_results<std::string::const_iterator> results;
			SystemCall call = Utils::ExecuteCommand(
				"vmstat -s | grep 'free memory'", "", false);
			boost::regex_match(call.output, results, 
				boost::regex("[\\s|\\t]+(\\d+)[\\s|\\t]+.*"));
			try{
				mem = stof(results[1].str()) / (1024 * 1024);
			}
			catch (exception &e){
			}
			if (mem > 0) return mem;

			call = Utils::ExecuteCommand("free | grep -i mem", "", false);
			boost::regex re("Mem:[\\s|\\t]+\\d+[\\s|\\t]+\\d+[\\s|\\t]+(\\d+)[\\s|\\t]+.*");
			boost::regex_match(call.output, results, re);
			try{
				mem = stof(results[1].str()) / (1024 * 1024);
			}
			catch (exception &e){
			}
			if (mem > 0) return mem;
		}
		catch (std::exception const& e) {
			BOOST_LOG_CHANNEL(logger::get(), "Main") 
				<< "Warning: Unable to get free memory. " << e.what() << endl;
		}
		return mem;
	}

	std::pair<double, double> SystemInfo::GetTotalMemoryUsage(int sid, int pid){
		std::pair<double, double> mem (0, 0);

		boost::match_results<std::string::const_iterator> results;
		string padRe = "[\\s|\\t|\\n|\\r]*";
		boost::regex re(padRe + "(\\d+)" + padRe + "(\\d+)" + padRe + "(\\d+)" + padRe);
	
		SystemCall call = Utils::ExecuteCommand(("ps -o vsize=,pid=,ppid= --sid=" + 
			boost::lexical_cast<string>(sid)).c_str());
		if (call.exit_code == 0) {
			for (auto str : Utils::SplitToVector(call.output, "\n")) {
				if (str == "")
					continue;
				boost::regex_match(str, results, re);
				double mem_ = stof(results[1].str()) / (1024 * 1024);
				mem.first += mem_;

				if (stoi(results[2].str()) == pid || stoi(results[3].str()) == pid)
					mem.second += mem_;
			}
		}
		return mem;
	}

	ComputePerTask ComputePerTask::CalculateComputePerTask(int max_splits,
		int min_cores_per_split, int max_mem_per_split)
	{
		ComputePerTask cpt(1, MAX_CPU);
		cpt.numSplits = (int)(MAX_MEM / max_mem_per_split);
		if (cpt.numSplits > max_splits) cpt.numSplits = max_splits;

		cpt.numCoresPerSplit = max(1, (int)(MAX_CPU / cpt.numSplits));

		if (cpt.numCoresPerSplit < min_cores_per_split) {
			cpt.numSplits = max(1, (int)(MAX_CPU / min_cores_per_split));
			cpt.numCoresPerSplit = min(MAX_CPU, min_cores_per_split);
		}
		cpt.numResidualCores = MAX_CPU - (cpt.numCoresPerSplit * cpt.numSplits);
		return cpt;
	}

}