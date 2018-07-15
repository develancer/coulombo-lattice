// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for very large systems
#ifndef COULOMBO_PARSER_HPP
#define COULOMBO_PARSER_HPP

#include <map>
#include <set>
#include <string>
#include <vector>

//----------------------------------------------------------------------

struct ParseSettings {
	std::map<char, std::string> shortNames;
	std::set<std::string> allowedFlags;
	std::set<std::string> allowedValues;
};

//----------------------------------------------------------------------

class ParseResults {
	std::set<std::string> flags_;

	std::map<std::string, std::string> values_;

	std::vector<std::string> args_;

public:
	ParseResults(void);

	ParseResults(const ParseSettings& settings, int argc, char** argv);

	bool hasFlag(const std::string& name) const;

	bool hasValue(const std::string& arg) const;

	bool hasValue(const std::string& arg, std::string& value) const;

	const std::string& getValue(const std::string& arg) const;

	int getArgCount(void) const;

	const std::string& getArg(int n) const;
};

//----------------------------------------------------------------------

class Parser {
	ParseSettings settings_;

public:
	Parser(void);

	Parser& allowFlag(char shortName, const std::string& name);

	Parser& allowFlag(const std::string& name);

	Parser& allowValue(const std::string& name);

	ParseResults process(int argc, char** argv) const;
};

//----------------------------------------------------------------------

#endif
