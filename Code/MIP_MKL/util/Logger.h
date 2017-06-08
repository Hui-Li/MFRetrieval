#ifndef LOGGER_H
#define LOGGER_H

#include "Base.h"
#include "Conf.h"
#include "FileUtil.h"

namespace Logger {

	string filePath;

	inline void open(string _filePath) {
		Logger::filePath = _filePath;
		FileUtil::createFolder(filePath);
	}

	inline void Log(string log, bool printToConsole=true) {
		if (printToConsole) {
			cout << log << endl;
		}
		if (Conf::log) {
			ofstream fileStream;
			fileStream.open(Logger::filePath.c_str(), ios_base::app);
			fileStream << log << endl;
			fileStream.close();
		}
	}

};

#endif //LOGGER_H
