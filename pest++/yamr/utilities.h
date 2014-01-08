#ifndef UTILITIES_H_
#define UTILITIES_H_

/** @file
 @brief Utility Functions
 
 This file contains a variety of utility functions for string and numeric operations.
*/

#include <string>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include "pest_error.h"

namespace pest_utils
{
	enum CASE_CONV{NO_CONV, TO_UPPER, TO_LOWER};

   /** @brief Sign of a number
 
   Returns the sign of a val (-1, 1 or 0).
  */
	inline int sign(double val)
	{
		if (val > 0) return 1;
		if (val < 0) return -1;
		return 0;
	}

	 /** @brief Sign of a number
 
	  
	/** @brief Splits a string and returns the sub-strings.
	
	String str is split into sub-strings using the characters specified 
	in delimiters and the sub-strings are are returned in the tokens constainer.
	If trimEmpty is specified as true, empty records are removed from
	tokens before it is returned
	*/
	template < class ContainerT >
	void tokenize(const std::string& str, ContainerT& tokens,
		const std::string& delimiters=" \t\n\r", const bool trimEmpty=true);


	/** @brief Convert a string to another type.
	
	String s is converted to the type of the return variable x.  If extra
	characters are left over after the conversion and failIfLeftoverChars is set
	to true, a PestConversionError exception is thrown
	*/
	template<typename T>
	inline void convert_ip(string const& s, T& x, bool failIfLeftoverChars = true)
	{
		istringstream i(s);
		char c;
		if (!(i >> x) || (failIfLeftoverChars && i.get(c))) {
			throw PestConversionError(s);
		}
	}

	/** @brief Convert a string to another type.
	
	String s is converted to the type T.  If extra
	characters are left over after the conversion and failIfLeftoverChars is set
	to true, a PestConversionError exception is thrown.

	The tepmplatized return type must be included in the fuction call.  The following example shows 
	how to call this function to convert a string to an integer
	convert_cp<int>("10")
	*/
	template<typename T>
	inline T convert_cp(std::string const& s,
		bool failIfLeftoverChars = true)
	{
		T x;
		convert_ip(s, x, failIfLeftoverChars);
		return x;
	}

	/** @brief Strip leading and/or trailing characters from a string
	
		The characters contained in arguement delimiters are stripped from 
		string s.  op can be specified as "front", "back" or "both" to control
		whether the characters are stripped from the begining, end or both sides 
		of string s.
	*/



	std::string& strip_ip(string &s, const string &op="both",
		const string &delimiters=" \t\n\r");

	/** @brief Strip leading and/or trailing characters from a string
	
		The characters contained in arguement delimiters are stripped from 
		string s and the updated string is returned without modifying s.  op can be specified as "front", "back" or "both" to control
		whether the characters are stripped from the begining, end or both sides 
		of string s.
	*/
	string strip_cp(const string &s, const string &op="both", 
		const string &delimiters=" \t\n\r");

	/** @brief Convert all the characters in a string to upper case
	
		All characters in string s are converted to upper case
	*/
	void upper_ip(string &s);

	/** @brief Convert all the characters in a string to upper case
	
		All characters in string s are converted to upper case and returned as 
		a new string without modifying the original string
	*/
	string upper_cp(const string &s);

		/** @brief Convert all the characters in a string to upper case.
	
		All characters in string s are converted to lower case.
	*/
	string upper(char *);

	void lower_ip(string &s);

	/** @brief Convert all the characters in a string to lower case.
	
		All characters in string s are converted to lower case and returned as 
		a new string without modifying the original string.
	*/
	string lower_cp(const string &s);

	/** @brief Return the base filename (filename without the "." extention).
	*/
	string get_base_filename(const string &s);

	/** @brief Converts a C++ string to a FORTRAN character array.

		Converts string "in" to a FORTRAN character array of length "length".
	*/
	void string_to_fortran_char(string in, char out[], int length, CASE_CONV conv_type=NO_CONV);


	string remove_file_ext(const string &filename, size_t max_len=string::npos);
	/** @brief Given a combined path and filename return just the filename.

		Given path the combined path and filname complete_path, return just the filename.
	*/

	string get_filename(const string &complete_path);

	/** @brief Given a combined path and filename return just the pathname.

		Given path the combined path and filname complete_path, return just the pathname.
	*/
	string get_pathname(const string &complete_path);


	template <class keyType, class dataType>
	vector<keyType> get_map_keys(const map<keyType,dataType> &my_map);

class String2CharPtr
{
public:
	String2CharPtr(const std::string &str);
	char *get_char_ptr();
	~String2CharPtr() {}
private:
	std::vector<char> my_str;
};


class StringvecFortranCharArray
{
public:
	StringvecFortranCharArray(const vector<string> in, int length, CASE_CONV conv_type=NO_CONV );
	char *get_prt();
	~StringvecFortranCharArray();
private:
	char *fort_array;

};

template <class type>
class CompareItemInSet
{
	public:
		CompareItemInSet(const std::set<type> &_set_ref) : set_ref(_set_ref){};
		bool operator()(const type &item) {return set_ref.count(item) > 0;}
		~CompareItemInSet(){}
	private:
		const std::set<type> &set_ref;
};

void copyfile(const string &from_file, const string &to_file);

std::string fortran_str_2_string(char *fstr, int str_len);

std::vector<std::string> fortran_str_array_2_vec(char *fstr, int str_len, int fstr_len);
}  // end namespace pest_utils
#endif /* UTILITIES_H_ */
