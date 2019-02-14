/**
This file is added by Ian to control output of ORB-SLAM2
*/

#ifndef OUTPUT_H
#define OUTPUT_H

#include <string>

namespace ORB_SLAM2
{

class Output
{
public: 
	static Output& instance()
	{
		static Output* instance = new Output();
		return *instance;
	}
	
	void set(bool LC_on, std::string filename)
	{
		valueSet = true;
		loopClosingOn = LC_on;
		outfile = filename;
	}
	
	bool valueSet;
	bool loopClosingOn;
	std::string outfile;

private:
	Output() : valueSet(false) {}
};

}// namespace ORB_SLAM

#endif // OUTPUT_H

