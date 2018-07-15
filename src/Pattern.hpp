// Coulombo Ⓒ 2018
// [Computer Physics Communications] Różański & Zieliński:
// Efficient computation of Coulomb and exchange integrals for very large systems
#ifndef COULOMBO_PATTERN_HPP
#define COULOMBO_PATTERN_HPP

#include <list>
#include <string>
#include <vector>

//----------------------------------------------------------------------

class Pattern {
public:
	Pattern(const std::string& description)
	{
		std::string::size_type start = 0, comma;
		while ((comma = description.find(',', start))!=std::string::npos) {
			addPattern(description.substr(start, comma-start));
			start = comma+1;
		}
		addPattern(description.substr(start));
	}

	bool match(int i0, int i1, int i2, int i3) const
	{
		for (const std::string& pattern : patterns) {
			if (matchPattern(pattern, i0, i1, i2, i3)) {
				return true;
			}
		}
		return false;
	}

private:
	std::list<std::string> patterns;

	static void checkCharacter(char c)
	{
		bool valid = (c=='*') || (std::isalnum(c) && c!='0');
		if (!valid) {
			throwfr("invalid character '%c' in pattern", c);
		}
	}

	static void checkPattern(const std::string& pattern)
	{
		if (pattern.empty()) {
			throwfr("invalid empty pattern");
		}
		else {
			for (char c : pattern) {
				checkCharacter(c);
			}
		}
	}

	static bool matchLetter(char c, int i, std::vector<int>& assignments)
	{
		assert(i>0);
		if (c=='*') {
			return true;
		}
		else if (c>='1' && c<='9') {
			return i==(c-'0');
		}
		else {
			if (assignments[c]) {
				return i==assignments[c];
			}
			else {
				assignments[c] = i;
				return true;
			}
		}
	}

	static bool matchPattern(const std::string& pattern, int i0, int i1, int i2, int i3)
	{
		std::vector<int> assignments(128, 0); // indexed by ASCII codes 0-127
		// TODO zmiana na va_list
		return matchLetter(pattern[0], i0, assignments)
				&& matchLetter(pattern[1], i1, assignments)
				&& matchLetter(pattern[2], i2, assignments)
				&& matchLetter(pattern[3], i3, assignments);
	}

	void addPattern(const std::string& pattern)
	{
		checkPattern(pattern);
		patterns.push_back(pattern);
	}
};

//----------------------------------------------------------------------

#endif
