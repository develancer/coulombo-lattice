// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for very large systems
#include "throwf.hpp"
#include "Parser.hpp"

Parser::Parser(void) { }

Parser& Parser::allowFlag(char shortName, const std::string& name)
{
	allowFlag(name);
	settings_.shortNames[shortName] = name;
	return *this;
}

Parser& Parser::allowFlag(const std::string& name)
{
	settings_.allowedFlags.insert(name);
	return *this;
}

Parser& Parser::allowValue(const std::string& name)
{
	settings_.allowedValues.insert(name);
	return *this;
}

ParseResults Parser::process(int argc, char** argv) const
{
	return ParseResults(settings_, argc, argv);
}

//----------------------------------------------------------------------

ParseResults::ParseResults(const ParseSettings& settings, int argc, char** argv)
{
	flags_.clear();
	values_.clear();
	args_.clear();

	for (int argi = 1; argi<argc; ++argi) {
		std::string s = argv[argi];
		if (s.size()>1 && s[0]=='-') {
			// flaga lub wartość
			if (s[1]=='-') {
				// długa flaga lub wartość
				std::string::size_type eq = s.find('=', 2);
				if (eq!=std::string::npos) {
					// wartość
					std::string name = s.substr(2, eq-2);
					std::string value = s.substr(eq+1);

					if (settings.allowedValues.count(name)) {
						values_[name] = value;
					}
					else if (settings.allowedFlags.count(name)) {
						throwfr("flag '%s' does not need a value", name.c_str());
					}
					else {
						throwfr("unknown value '%s'", name.c_str());
					}
				}
				else {
					// długa flaga
					std::string name = s.substr(2);

					if (settings.allowedFlags.count(name)) {
						flags_.insert(name);
					}
					else if (settings.allowedValues.count(name)) {
						throwfr("missing value for parameter '%s'", name.c_str());
					}
					else {
						throwfr("unknown flag '%s'", name.c_str());
					}
				}
			}
			else {
				// krótka flaga
				for (std::string::size_type i = 1; i<s.size(); ++i) {
					char shortName = s[i];

					std::map<char, std::string>::const_iterator it = settings.shortNames.find(shortName);
					if (it!=settings.shortNames.end()) {
						std::string name = it->second;
						flags_.insert(name);
					}
					else {
						throwfr("unknown flag '%c'", shortName);
					}
				}
			}
		}
		else {
			// argument
			args_.push_back(s);
		}
	}
}

bool ParseResults::hasFlag(const std::string& name) const
{
	return flags_.count(name);
}

bool ParseResults::hasValue(const std::string& name) const
{
	return values_.count(name);
}

bool ParseResults::hasValue(const std::string& name, std::string& value) const
{
	std::map<std::string, std::string>::const_iterator it = values_.find(name);
	if (it!=values_.end()) {
		value = it->second;
		return true;
	}
	return false;
}

const std::string& ParseResults::getValue(const std::string& name) const
{
	std::map<std::string, std::string>::const_iterator it = values_.find(name);
	if (it==values_.end()) throwfr("parameter '%s' not specified", name.c_str());
	return it->second;
}

int ParseResults::getArgCount(void) const
{
	return args_.size();
}

const std::string& ParseResults::getArg(int n) const
{
	return args_.at(n);
}
