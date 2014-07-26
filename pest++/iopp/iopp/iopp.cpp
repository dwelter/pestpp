// class libraries to replicate the fortran pest IO functionailty


#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <string>
#include <regex>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <exception>
#include <io.h>
#include "iopp.h"

using namespace std;

string string2upper(string& str)
{
	for(string::iterator s = str.begin();s != str.end();++s)
	{
		*s = toupper((unsigned char)*s);
	}
	return str;
}

string string2lower(string& str)
{
	for(string::iterator s = str.begin();s != str.end();++s)
	{
		*s = tolower((unsigned char)*s);
	}
	return str;
}

string::size_type find_index(const string strg,const char target)
{
	string::size_type i;
	bool found = false;
	for (i = 0;i != strg.size();i++)
	{
		if (strg[i] == target)
		{
			found = true;
			break;
		}
	}
	if (found == false)
	{
		throw runtime_error(" fixed parameter misssing ']'" + strg);
	}
	return i;
}
	
pair<int,int> parse_indices(const string strg, const char delim)
{
	vector<string> elems;
	stringstream ss(strg);
	string item;
	while(getline(ss,item,delim))
	{
		elems.push_back(item);
	}
	if (elems.size() != 2)
	{
		throw runtime_error("couldn't parse column indices for instruction :" + strg);
	}

	int start,end;
	
	try
	{
		stringstream convert(elems[0]);
		convert >> start;
		if (convert.fail())
		{
			throw runtime_error("could not cast start column index : " + elems[0] +" for instruction : "+strg);
		}
	}
	catch (exception e)
	{
		throw runtime_error("could not cast start column index : " + elems[0] +" for instruction : "+strg);
	}
	
	try
	{
		stringstream convert(elems[1]);
		convert >> end;
		if (convert.fail())
		{
			throw runtime_error("could not cast end column index : " + elems[1] +" for instruction : "+strg);
		}
	}
	catch (exception e)
	{
		throw runtime_error("could not cast end column index : " + elems[1] +" for instruction : "+strg);
	}
	pair<int,int> idx (start,end);
	return idx;
}

void rewind_file(ifstream *f, streampos line_start,int num_chars)
{
	f->clear();
	f->seekg(line_start);
	if (f->fail())
	{
		throw RewindError("could not reset g pointer to position: " + line_start);
	}
	f->clear();
	char c;
	for (int i=0;i<num_chars;i++)
	{
		f->clear(); 
		c = f->get();
		if (f->fail())
		{
			f->clear();
			f->seekg(line_start);
			string line;
			getline(*f,line);
			throw RewindError("could not read character from line "+line);
		}

	}
	f->clear();
}

int get_nonnumeric(string str)
{
	//read the string in reverse - avoids base,sign,radix
	char c;
	for (int i=str.size()-1;i>=0;i--)
	{
		c = str[i];
		if (isdigit(c))
		{
			return i;
		}
	}
	return -1;
}

double text_to_num(string text,bool strict)
{
	double val;
	stringstream converter;
	converter << text;
	converter >> val;
	if ((converter.fail()))
	{
		throw Text2NumParseError("Could not convert string to numeric: "+text);
	}
	if (strict)
	{
		int nnpos = get_nonnumeric(text);
		if (nnpos == -1)
		{
			throw Text2NumParseError("nonnumeric characters found in string "+text);
		}
	}
	return val;
}

string read_line(ifstream* out)
{
	string line;
	out->clear();
	getline(*out,line);
	if (out->fail())
	{
		throw ReadLineError();
	}
	out->clear();
	return line;
}

string read_line(ifstream* out,streampos seekpoint)
{
	string line;
	out->clear();
	out->seekg(seekpoint);
	if (out->fail())
	{
		throw ReadLineError(" unable to seek to streampos");
	}
	out->clear();
	getline(*out,line);
	if (out->fail())
	{
		throw ReadLineError(" unable to read line ");
	}
	out->clear();
	return line;
}

string escape_regex(string marker_string)
{
	//vector<string> meta_chars = {"\\","*","+","?","|","{","}","[","]","^","$",".","#"};
	string meta_chars = "\\*+?|{}[]^$.#";
	string escaped_string = "";
	size_t found;
	for (auto& c : marker_string)
	{
		found = meta_chars.find(c);
		if (found != string::npos)
		{
			escaped_string += "\\";
		}
		escaped_string += c;
	}
	return escaped_string;
}



Instruction::Instruction (string i_string,char m_delim)
{
	ins_string = i_string;
	marker_delim = m_delim;
	mtype = Instruction::marker_type::secondary;
}

void Instruction::execute(ifstream* out,int* lpos)
{
	switch (itype)
	{
	case instruction_type::fixedLine:
		execute_fixed(out,lpos);
		break;
	case instruction_type::marker:
		execute_marker(out,lpos);
		break;
	case instruction_type::whitespace:
		execute_whitespace(out,lpos);
		break;
	default:
		throw runtime_error("unknown instruction_type: "+itype);
	}
}

double Instruction::read(ifstream* out,int* lpos)
{
	double val;
	switch (itype)
	{
	case instruction_type::fixedObs:
		val = read_fixedObs(out,lpos);
		break;
	case instruction_type::nonFixedObs:
		val = read_nonFixedObs(out,lpos);
		break;
	default:
		throw runtime_error(" unrecognized instruction type "+itype);
	}

	return val;
}

void Instruction::parse()
{
	char first_char = ins_string[0]; 
		
	//observation instructions
	if ((first_char == '[') || (first_char == '(') || (first_char ==  '!'))
	{
		if (first_char ==  '!')
		{
			itype = instruction_type::nonFixedObs;
			parse_nonFixedObs();
		}
		else
		{
			itype = instruction_type::fixedObs;
			parse_fixedObs();
		}
	}	
	
        
	//other instructions
	else 
	{
		if((first_char ==  'l') || (first_char == 't'))
		{
			itype = instruction_type::fixedLine;
			parse_fixed();
		}
		else if (first_char == marker_delim)
		{
			itype = instruction_type::marker;
			parse_marker();
		}
		else if (toupper(first_char) == 'W')
		{
			itype = instruction_type::whitespace;
			parse_whitespace();
		}	
		else
		{
			throw runtime_error(" unrecognized instruction: "+ ins_string);
		}
	}
	return;
}

void Instruction::parse_whitespace()
{
}

void Instruction::execute_whitespace(ifstream* out, int* lpos)
{
	start_point = out->tellg();
	string line,subline;
	try
	{
		line = read_line(out);
	}
	catch (ReadLineError)
	{
		throw WhitespaceExecutionError(" reading line");
	}
	//read to a whitespace
	char c;
	int i,ii;
	for (i=0;i<line.size();i++)
	{
		c = line[i];		
		if (isspace(c))
		{
			break;
		}
	}
	//read to the end of the whitespace
	for (ii=i;i<line.size();ii++)
	{
		c = line[ii];
		if (!isspace(c))
		{
			break;
		}
	}
	*lpos += ii;

	try
	{
		rewind_file(out,start_point,ii);
	}
	catch (RewindError& e)
	{
		string message = " error rewinding file ";
		throw WhitespaceExecutionError(message.append(e.what()));
	}
	return;
}

void Instruction::parse_fixed()
{
	(ins_string[0] == 'l') ? rtype = read_type::line:rtype = read_type::tab;
	string rstring = ins_string;
	rstring.erase(0,1);
	stringstream convert(rstring);
	convert >> number;
	if (number < 1) throw invalid_argument("could not cast to integer in instruction : " + ins_string);
	return;
}

void Instruction::execute_fixed(ifstream* out,int* lpos)
{
	start_point = out->tellg();
	streampos curr_point;
	long line_length;
	string line;
	if (rtype == Instruction::read_type::line)
	{

		for (int i=0;i<number;i++)
		{
			start_point = out->tellg();
			line = read_line(out);
			curr_point == out->tellg();	
		}
		out->seekg(start_point);
		if (out->fail())
		{
			throw runtime_error(" could not set g-pointer to start of line in fixed-line read");
		}
	}
	else if (rtype == Instruction::read_type::tab)
	{
		line = read_line(out);
		try
		{
			rewind_file(out,start_point,number-*lpos);
		}
		catch (RewindError& e)
		{
			throw runtime_error("could not reset file g-pointer to line position: "+number);
		}
		//getline(*out,line);
		int i=0;

	}
	return;		
}

void Instruction::parse_marker()
{
	marker_string = ins_string;
	//we know the first char is the marker delim
	marker_string.erase(0,1);
	int idx = find_index(marker_string,marker_delim);
	//make sure the the marker delim is the last char
	if (idx != marker_string.size()-1)
	{
		throw MarkerError(" parsing marker: marker not bracketed by marker delim"+ins_string);
	}
	marker_string.erase(idx,1);
	return;
}

void Instruction::execute_marker(ifstream* out,int* lpos)
{
	//build a regex for the marker - fast
	regex rmarker;
	try
	{
		rmarker.assign(escape_regex(marker_string));
	}
	catch (...)
	{
		throw MarkerError(" building regex from expression: "+marker_string);
	}
	smatch mresults; 
	string line,subline;
	start_point = out->tellg();
	if (mtype == marker_type::primary)
	{
		
		while (true)
		{
			try
			{
				line = read_line(out);
			}
			catch (ReadLineError)
			{
				break;
			}
			//check for marker string in line
			if (regex_search(line,mresults,rmarker))
			{	
				//seek to the end of the marker
				int mstart = line.find(marker_string);
				try
				{
					rewind_file(out,start_point,mstart + marker_string.size());
				}
				catch (RewindError& e)
				{
					string message = " error rewinding file ";
					throw MarkerError(message.append(e.what()));
				}		
				
				*lpos = mstart + marker_string.size();
				return;
			}
			start_point = out->tellg();
			*lpos = 0;
		}
		throw PrimaryMarkerEOFError(marker_string);
	}
	//this is a secondary marker
	else
	{
		//get the current line
		try
		{
			line = read_line(out);
			if (regex_search(line,mresults,rmarker))
			{	
				//seek to the end of the marker
				int mstart = line.find(marker_string);
				try
				{
					rewind_file(out,start_point,mstart+marker_string.size());
				}
				catch (RewindError& e)
				{
					string message = " error rewinding file ";
					throw MarkerError(message.append(e.what()));
				}		
				*lpos += mstart + marker_string.size();
				return;
			}
			else
			{
				throw SecondaryMarkerNotFound(marker_string);
			}	
		}
		catch (ReadLineError)
		{
			throw SecondaryMarkerReadError("unable to read line from file while searching for secondary marker: " + marker_string);
		}
	}

	
}

double Instruction::read_nonFixedObs(ifstream* out,int* lpos)
{
	start_point = out->tellg();
	string line,subline;
	double dval = -1.0E+10;
	line = read_line(out);
	stringstream sline(line);
	sline >> subline;
	if (sline.fail())
	{
		throw NonFixedObsReadError(" from line: " + line + " for instruction: " + ins_string);
	}
	try
	{
		//read the double with out strict - might be trailing marker characters
		dval = text_to_num(subline,false);
	}
	catch (Text2NumParseError)
	{
		rewind_file(out,start_point,line.find(subline)+subline.size());
		throw NonFixedObsParseError(" casting string to double " + subline);
	}
	int nnpos = get_nonnumeric(subline);
	int lstart = line.find(subline);
	int lend = lstart + subline.size();
	//if nonnumeric characters are detected, then read the file forward to the last numeric character
	if (nnpos != -1)
	{
		lend = lstart + nnpos+1;	
	}
	try
	{
		rewind_file(out,start_point,lend);
	}
	catch (RewindError& e)
	{
		string message = " error rewinding file ";
		throw MarkerError(message.append(e.what()));
	}		
	*lpos += lend;
	//getline(*out,line);
	return dval;
}

void Instruction::parse_nonFixedObs()
{
	ins_string.erase(0,1);
	string::size_type i = find_index(ins_string,'!');
	obs_name = ins_string.substr(0,i);
	return;
}

double Instruction::read_fixedObs(ifstream* out,int* lpos)
{
	start_point = out->tellg();
	string line;
	double dval = -1.0E+10;
	line = read_line(out);
	if (ftype == fixed_type::strict)
	{
		string sval;
		try
		{
			sval = line.substr(start-*lpos-1,end-start+1);
		}
		catch (...)
		{
			int istart = start - *lpos -1;
			int iend = end - start + 1;
			throw FixedObsReadError(" forming substring from line " + line + "from positions: " + to_string(istart) + " " + to_string(iend));
		}
		try
		{
			dval = text_to_num(sval,true);
		}
		catch (Text2NumParseError)
		{
			throw FixedObsReadError(" casting string to double " + sval);
		}	
		int i;
		for (i=sval.size()-1;i>=0;i--)
		{
			if (!isspace(sval[i]))
			{
				break;
			}
		}
		*lpos += start + i;
		try
		{
			rewind_file(out,start_point,*lpos);
		}
		catch (RewindError& e)
		{
			string message = " error rewinding file ";
			throw MarkerError(message.append(e.what()));
		}
		


		
	}
	else
	{
		//split the line on white spaces and built a vector of start index and length of each entry
		string subline;
		int s,c = 0;
		vector<string> entries;
		vector<int> entry_start,entry_end;
		stringstream sline(line);
		while (sline)
		{
			sline >> subline;
			if (!sline)
				break;
			s = line.find(subline,c);
			entries.push_back(subline);
			entry_start.push_back(s);
			entry_end.push_back(s + subline.size());
			c = s;
		}
		
		//find an entry that meets the requirements of the instruction
		int i;
		bool found = false;
		for (i=0;i<entries.size();i++)
		{
			//if the starting index is bracketed
			s = entry_start[i];
			int e = entry_end[i];
			if ((entry_start[i] <= (start - *lpos)) && ((entry_end[i]) >= (start - *lpos)))
			{
				found = true;
				break;
			}
			//if the ending index is bracketed
			if ((entry_start[i] <= (end - *lpos)) && ((entry_end[i]) >= (end - *lpos)))
			{
				found = true;
				break;
			}
			//if the entire entry is bracketed
			if ((entry_start[i] >= (start - *lpos)) && (entry_end[i] <= (end - *lpos)))
			{
				found = true;
				break;
			}
		}
		if (!found)
		{
			throw SemiFixedObsReadError(" couldn't find token meeting index requirements on line: " + line);
		}
		try
		{
			dval = text_to_num(entries[i],true);
		}
		catch (Text2NumParseError& e)
		{
			throw SemiFixedObsReadError(" casting string to double " + entries[i] + e.what());
		}
		
		try
		{
			rewind_file(out,start_point,entry_end[i]);
		}
		catch (RewindError& e)
		{
			string message = " error rewinding file ";
			throw MarkerError(message.append(e.what()));
		}		
		*lpos = *lpos + entry_end[i];
	}
	return dval;
}

void Instruction::parse_fixedObs()
{
	char end_char;
	if (ins_string[0] == '['){
		end_char = ']';
		ftype = fixed_type::strict;
	}
	else
	{
		end_char = ')';
		ftype = fixed_type::flexible;
	}
	ins_string.erase(0,1);
	string::size_type i = find_index(ins_string,end_char);
	//get the obs name
	obs_name = ins_string.substr(0,i);
	ins_string.erase(0,i+1);
	pair<int,int> idx = parse_indices(ins_string,':');
	start = idx.first;
	end = idx.second;
}



InstructionLine::InstructionLine(vector<Instruction> instructs)
{
	lpos = 0;
	instructions = instructs;

	if (instructions[0].itype == Instruction::instruction_type::fixedLine)
	{
		markerFirst = false;
	}
	else if (instructions[0].itype == Instruction::instruction_type::marker)
	{
		markerFirst = true;
	}
	else
	{
		throw InstructionLineError("first instruction must be marker of fixed-line read");
	}
}

unordered_map<string,double> InstructionLine::execute(ifstream* out)
{
	vector<double> ovals;
	double oval;
	string line;
	lpos = 0;
	Instruction i;
	size_t icount = 0;
	streampos line_start;
	unordered_map<string,double> obs_map; 
	pair<string,double> opair;
	while (icount < instructions.size())
	{
		i = instructions[icount];
		if ((i.itype == Instruction::instruction_type::fixedObs) ||
			(i.itype == Instruction::instruction_type::nonFixedObs))
		{
			line_start = out->tellg();
			try
			{
				oval = i.read(out,&lpos);
			}
			//if a non-fixed obs failed to cast, check for a secondary marker right behind it
			catch (NonFixedObsParseError& e)
			{
				if ((icount+1<instructions.size()) && (instructions[icount+1].itype == Instruction::instruction_type::marker))
				{
					throw InstructionLineError("ifstream could not handle trailing marker");
				}
				else if(i.get_name() != "dum")
				{
					throw InstructionLineError("Nonfixed-obs parse error and no trailing secondary marker ");
				}
			}
			catch (runtime_error& e)
			{
				throw InstructionLineError("Error reading instruction "+i.get_name()+" "+e.what());
			}
			if (string2lower(i.get_name()) != "dum")
			{
				ovals.push_back(oval);
				opair.first = string2upper(i.get_name());
				opair.second = oval;
				obs_map.insert(opair);

			}
			icount++;
		}
		else
		{
			try
			{
				i.execute(out,&lpos);
				icount++;
			}
			
			//if a secondary marker was not found, start looking for a new line
			catch (SecondaryMarkerNotFound& e)
			{
				if (!markerFirst)
				{
					throw InstructionLineError("Secondary Marker not found, but not a primary marker line");
				}
				else
				{
					//cout << e.what() << endl;
					icount = 0;
					lpos = 0;
				}
			}
			catch (runtime_error& e)
			{
				getline(*out,line);
				throw InstructionLineError(" executing instruction " + i.get_name()+" "+e.what());
			}
		}
	}
	//read any remaining junk on the line
	line = read_line(out,line_start);
	return obs_map;
}

vector<string> InstructionLine::get_observation_names()
{
	string oname;
	vector<string> obs_names;
	for (auto& i : instructions)
	{
		oname = i.get_name();
		if (oname != "")
		{
			obs_names.push_back(oname);
		}
	}
	return obs_names;
}


InstructionFile::InstructionFile(string ins_filename)
{
	instruction_filename = ins_filename;
	ins.open(ins_filename);
	if (!ins.is_open())
	{
		throw InstructionFileError("could not open instruction file "+ins_filename);
	}
	string line,ftype,mstring;
	getline(ins,line);
	istringstream iss(line);

	iss >> ftype >> mstring;
	if (string2upper(ftype) != "PIF")
	{
		throw InstructionFileError("instruction file doesn't start with 'PIF'");
	}
	if (mstring.size() != 1)
	{
		throw InstructionFileError("marker delimiter must be a single character, not " + mstring);
	}
	marker_delim = mstring[0];
	build_instruction_set();
}

void InstructionFile::build_instruction_set()
{
	string iline;
	while (!ins.eof())
	{
		iline = get_instruction_line();
		parse_instruction_line(iline);
	}
	return;

}

string InstructionFile::get_instruction_line()
{
	string iline,iline2;
	//read atleast one line
	if (!getline(ins,iline)) return iline;
	//ncheck if the next starts with '&'
	while (true)
	{
		if (ins.peek() == '&')
		{
			if (!getline(ins,iline2))
			{
				throw InstructionFileError("Unable to read line following continuation character on line " + iline);
			}
			iline += iline2.substr(1,iline.size());
		}
		else break;
	}
	return iline;
}

void InstructionFile::parse_instruction_line(string ins_line)
{
	vector<string> ins_strings;
	string itoken;
	//check for marker instructions
	int i = ins_line.find(marker_delim,0);
	if (i >= 0)
	{
		vector<pair<int,int>> indices = get_marker_indices(ins_line);
		string subline;
		size_t start = 0,end;
		for (auto &idx : indices)
		{
			end = idx.first;
			//first process tokens to the left of this starting index
			if (start != end)
			{
				subline = ins_line.substr(start,end-start);
				istringstream iss(subline);
				while (iss >> itoken)
				{
					ins_strings.push_back(itoken);
				}

			}
			//set the entire marker substring as a instruction string
			start = idx.first;
			end = idx.second;
			ins_strings.push_back(ins_line.substr(start,end-start));

			//now reset the start index to the end 
			start = end;
		}
		//process any tokens right of the last index
		if (indices[indices.size()-1].second != ins_line.size())
		{
			start = indices[indices.size()-1].second;
			end = ins_line.size();
			subline = ins_line.substr(start,end);
			istringstream iss(subline);
			while (iss >> itoken)
			{
				ins_strings.push_back(itoken);
			}
		}
	}
	else
	{
		//process each token of the instruction line
		istringstream iss(ins_line);
		while (iss >> itoken)
		{
			ins_strings.push_back(itoken);
		}
	}

	vector<Instruction> instructs;
	for ( auto &ins : ins_strings)
	{
		Instruction i(ins,marker_delim);
		i.parse();
		instructs.push_back(i);
	}
	InstructionLine iline(instructs);
	instruction_lines.push_back(iline);
	return;
}

vector<pair<int,int>> InstructionFile::get_marker_indices(string line)
{
	size_t key = -1;
	size_t value = -1;
	vector<pair<int,int>> indices;
	for (string::size_type i=0;i < line.size();i++)
	{
		if (line[i] == marker_delim)
		{						
			if (key == -1)
			{
				key = i;
			}
			else
			{				
				if (i == (key+1))
				{					
					cout << " only one space between marker delimiters in file" << instruction_filename << "on line " << line << endl;
				}
				
				pair<int,int> p(key,i+1);
				indices.push_back(p);
				key = -1;
			}
		}
	}
	if (key != -1)
	{		
		cout << "unbalanced marker on line " << line << endl;
	}
	return indices;
}

unordered_map<string,double> InstructionFile::read(ifstream* out)
{
	unordered_map<string,double> obs_map,line_map;
	for (auto &i : instruction_lines)
	{
		line_map = i.execute(out);
		obs_map.insert(line_map.begin(),line_map.end());
	}
	return obs_map;
}

vector<string> InstructionFile::get_observation_names()
{
	vector<string> obs_names,line_names;
	for (auto& iline : instruction_lines)
	{
		line_names = iline.get_observation_names();
		obs_names.insert(obs_names.end(),line_names.begin(),line_names.end());
	}
	return obs_names;
}


InstructionFiles::InstructionFiles(vector<string> ins_files,vector<string>out_files,vector<string>onames)
{
	instruction_filenames = ins_files;
	output_filenames = out_files;
	obs_names = onames;
}

vector<double> InstructionFiles::readins()
{
	size_t i;
	vector<double> obs_vals(obs_names.size());
	unordered_map<string,double> ifile_map;
	vector<string> ifile_onames;
	unordered_map<string,double>::iterator miter;
	//should be initialize to false
	vector<bool> obs_found(obs_names.size());

	for (i=0;i<instruction_filenames.size();i++)
	{
		InstructionFile ins(instruction_filenames[i]);
		ifile_onames = ins.get_observation_names();
		ifstream out(output_filenames[i]);
		if (out.fail())
		{
			throw InstructionFileError(" could not open model output file " + output_filenames[i]);
		}
		try
		{
			ifile_map = ins.read(&out);
		}
		catch (runtime_error& e)
		{
			cout << "unable to read instruction file" << instruction_filenames[i] << e.what() << endl;
			throw e;
		}

		for (int ii=0;ii<obs_names.size();ii++)
		{
			miter = ifile_map.find(obs_names[ii]);
			if (miter != ifile_map.end())
			{
				obs_vals[ii] = miter->second;
				obs_found[ii] = true;
			}
		}

	}
	//check that all requested obs have been read
	bool missing=false;
	string missing_obs = "";
	for (i=0;i<obs_names.size();i++)
	{
		if (!obs_found[i])
		{
			missing_obs += " "+obs_names[i];
			missing = true;
		}
	}
	if (missing)
	{
		throw InstructionFileError(" observation values not found in model output files: "+missing_obs);
	}

	return obs_vals;
}


FixedWidthValue::FixedWidthValue(bool isDblPres, bool frcRad,int wdth)
{
	isDoublePrecision = isDblPres;
	forceRadix = frcRad;
	text_width = wdth;
	base = 'E';
	needSign = true;
	needEsign = true;
	needScientific = true;
	needRadix = true;
	isRace = false;
	(isDoublePrecision) ? max_exponent=278:max_exponent=38;
	(isDoublePrecision) ? max_sig=16:max_sig=8;	
	(isDoublePrecision) ? max_width=23:max_width=16;
	width = min(text_width,max_width);
}

string FixedWidthValue::get_templatefile_string(double val)
{	
	value = val;
	set_precision_components();
	//check for overflow
	if (abs(exponent) > max_exponent)
	{
		//cout << "value over flow " << value;		
		throw ValueOverflowError(max_str);
	}	

	//try to squeeze some digits using clever tricks that john developed
	maximize_precision();

	// build the output significand
	prepare_output_significand();

	//get the output string
	string value_str = as_string();
	
	size_t s = value_str.size(); 
	if (s != text_width)
		throw FixedWidthError("internal error - string representation not the right number of characters"+value_str);

	return value_str;
}	

void FixedWidthValue::maximize_precision()
{
	update();
	// try to drop scientific - the base, exponent sign and exponent
	drop_scientific();
	update();
	

	//try to gain a digit by reducing the exponent
	if ((needScientific) && (forceRadix)) 
	{
		reduce_exponent();
		update();
	}
	
	//if we can, try to drop the radix
	if (!forceRadix)
	{
		drop_radix();
		update();
	}
	
	//if we still don't have any digits, this space is too narrow
	if (sig_digits == 0)
		throw FixedWidthError("not enought spaces to represent value "+ max_str);
	
}

void FixedWidthValue::prepare_output_significand()
{
	
	// reset output_significand
	output_significand.clear();

	//if the significant is too short, pad with zeros if needed - only on the right
	if (sig_digits > significand.size())
	{
		output_significand = pad_zeros(significand,false);						
	}

	//truncate if needed
	else if (sig_digits < significand.size())
	{				
		for (vector<int>::size_type i=0;i<sig_digits;i++)
		{
			output_significand.push_back(significand[i]);
		}		
		//round if nessecary
		if (sig_digits < significand.size())
		{
			if (significand[sig_digits] >= 5)
			{
				//convert to integer
				long long trunc_d = digits_to_number(output_significand);
				//add one
				trunc_d += 1;
				//convert back to int vector
				output_significand = number_to_digits(trunc_d);	
				//if the rounding triggered a digit increase (if the rightmost value if 9)			
				while (output_significand.size() > sig_digits)
				{
					shift_left();
					output_significand.pop_back();
				}
				
				//replace any missing zeros, could be left or right, depending on if abs(value) < 1.0
				output_significand = pad_zeros(output_significand,needEsign);						
			}
		}
	}
	else
	{
		output_significand = significand;
	}
	return;
}

void FixedWidthValue::drop_radix()
{
	//if we need to squeeze this value into a narrower space
	if ((width < max_width) && (sig_digits < width))
	{		
		//hypothesis: we can drop the radix
		needRadix = false;
		update();
		//see if we are in a position to gain any digits by reducing the exponent
		reduce_exponent();		
		while (radix_pos < sig_digits)
		{										
			//if we dropped the exponent to zero, we don't need scientific but we need the radix
			if (exponent == 0)
			{
				needScientific = false;
				needRadix = true;	
				update();
				break;
			}				
			shift_right();													
		}
		//check for race condition - happens when exponent gains a digit: -9 -> -10 or -99 -> -100
		if (radix_pos > sig_digits)
		{			
			isRace = true;
			exponent++;
			radix_pos--;
			update();
		}
		//if we successfully moved the radix to the last position, try to move the exponent to gain a digit...		
		if (needRadix) reduce_exponent();		
		update();
	}
}

void FixedWidthValue::drop_scientific()
{	
	//if we need to squeeze the value into a narrower space
	if (width < max_width)
	{		
		if ((exponent == 0) || (value == 0.0)) needScientific = false;
		//if this is a really small value	
		else if ((needEsign) && (abs(exponent) < 3))
		{
			needScientific = false;
			update();
			while (exponent != 0)	
				shift_left();				
		}
		else 
		{						
			//hypothesis: we can drop scientific
			needScientific = false;	
			//...and the radix
			if(!forceRadix) needRadix = false;
			update();
			// abs(values) greater than 1.0
			if ((!needEsign) && (sig_digits >= exponent))
			{
				while (exponent != 0)	
					shift_right();	
			}
			//otherwise, we need to undo our hypothesis and keep scientific
			else
			{
				needScientific = true;
				needRadix = true;
				update();
			}
		}
		//if we didn't shift the radix all the way to the end of the significand, we still need it
		if (radix_pos < sig_digits) needRadix = true;
	}
	return;
}

void FixedWidthValue::reduce_exponent()
{	
	if (sig_digits < max_sig)
	{
		if ((exponent == 100) || ((exponent > 100) && (exponent - sig_digits -1 < 100)))
		{
			while (exponent >= 100)
			{
				shift_right();			
			}
		}	

		if ((exponent == 10) || ((exponent > 10) && (exponent - sig_digits -1 < 10)))
		{
			while (exponent >= 10)
			{
				shift_right();	
			}
		}
	}
}

string FixedWidthValue::as_string()
{
	
	stringstream value_str;	
	int position = 0;		
	if (needSign)
		value_str << sign;
	
	if ((needRadix) && (position == radix_pos))
	{
		value_str << '.';
		position++;
	}

	for (vector<int>::size_type i=0;i<output_significand.size();i++)
	{
		if ((needRadix) && (position == radix_pos)) value_str << '.';			
		value_str << output_significand[i];		
		position++;
	}
	if ((needRadix) && (position == radix_pos))
		value_str << '.';	
		
	if (needScientific)
	{
		//append the base	
		value_str << base;

		//if the value needs a sign on the exponent
		if (needEsign) value_str << eSign;
		
		//append the exponent - pretty hackish, but needed to control the digits of the exponent explicitly	
		char buf[5];		
		int pos_exp = abs(exponent);	
		//switch on the number of exponent digits needed 
		switch (exp_digits)
		{
			case (1):
			{			
				sprintf_s(buf,"%01d",pos_exp);
				for (int i=0;i<exp_digits;i++)
				{		
					value_str << buf[i];
				}
				break;
			}
			case (2):
			{			
				sprintf_s(buf,"%02d",pos_exp);
				for (int i=0;i<exp_digits;i++)
				{		
					value_str << buf[i];
				}
				break;
			}
			case (3):
			{		
				sprintf_s(buf,"%03d",pos_exp);
				for (int i=0;i<exp_digits;i++)
				{		
					value_str << buf[i];
				}
				break;
			}	
		}
	}
	
	string v_string = value_str.str();

	//fill any remaining characters - on the left with spaces
	while (v_string.size() < text_width)
	{
		v_string.insert(0," ");
	}
	return v_string;
}

double FixedWidthValue::as_double()
{
	string value_str = as_string();
	stringstream ss;
	ss << value_str;
	double val;
	ss >> val;
	return val;
}

void FixedWidthValue::update()
{
	//update current exp_digits and sig_digits using current state
	exp_digits = 3;
	if ((abs(exponent) < 100) && (exp_digits > 2))
	{
		exp_digits--;
	}

	if ((abs(exponent) < 10) && (exp_digits > 1))
	{
		exp_digits--;
	}
	//count up the number of non-numeric characters
	int nn = 0;	
	if (needSign) nn++;	
	if (needRadix) nn++;
	if (isRace) nn++;
	if (needScientific)
	{
		//the base 'E'
		nn++;
		//the exp sign 
		if (needEsign) nn++;
			
		nn += exp_digits;
	}
	(nn > width) ? sig_digits=0:sig_digits=width-nn;
	if ((needScientific) && (sig_digits > max_sig)) sig_digits = max_sig;
	return;
}

vector<int> FixedWidthValue::pad_zeros(vector<int> vi,bool onLeft)
{
	if (!onLeft)
	{
		while (vi.size() < sig_digits)
		{
			vi.push_back(0);
		}
	}
	else
	{
		while (vi.size() < sig_digits)
		{
			vector<int>::iterator it = vi.begin();
			vi.insert(it,0);
		}
	}
	return vi;
}

void FixedWidthValue::shift_left()
{
	if (radix_pos == 0)
	{
		vector<int>::iterator it = significand.begin();
		significand.insert(it,0);
	}
	else
		radix_pos--;
	exponent++;
	(exponent<0) ? needEsign=true:needEsign==false;
	update();
}

void FixedWidthValue::shift_right()
{
	if (radix_pos > significand.size()) significand.push_back(0);
	radix_pos++;
	exponent--;		
	(exponent<0) ? needEsign=true:needEsign==false;
	update();
}

void FixedWidthValue::set_precision_components()
{
	needSign = false;	
	if (value < 0.0)
	{
		sign = '-';
		needSign = true;		
	}
	
	needEsign = false;
	if ((abs(value) < 1.0) && (value != 0.0))	
	{
		eSign = '-';
		needEsign = true;
	}
	
	//print the string representation maximum precision
	char buf[50];	
	sprintf_s(buf,"%#23.16E",value);
	max_str = buf;
	
	//strip off any whitespace
	max_str.erase(remove_if(max_str.begin(),max_str.end(),isspace),max_str.end());

	size_t max_size = max_str.size();

	//set exponent
	string str_exp = max_str.substr(max_str.size()-3,3);	
	exponent = atoi(str_exp.c_str());
	if (needEsign) exponent *= -1;
	exp_digits = 3;
	
	//get the radix position
	radix_pos = max_str.find('.');
	//if negative, move the radix left one
	if (needSign) radix_pos--;
	if (radix_pos == max_str.size()-1)
	{
		//cout << "no radix found in max_str representation " << max_str;
		throw FixedWidthError("no radix found in max_str representation " + max_str);
	}

	//find the base	
	size_t e_idx = max_str.find('E');
	if (e_idx == max_str.size()-1)
	{
		//cout << "No 'E' found in max_str representation " << max_str;
		throw FixedWidthError("no 'E' found in max_str representation " + max_str);;
	}
	string significand_str = max_str.substr(0,e_idx);
	
	//remove the leading sign
	if (needSign) significand_str.erase(0,1);

	//remove the radix
	significand_str = significand_str.erase(radix_pos,1);	

	//populate significand	
	for (vector<int>::size_type i=0;i<significand_str.size();i++)	
	{		
		int si = significand_str[i] - '0';
		if ((si < 0) || (si > 9))
		{
			//cout << "error casting significand component to intergers" << si;
			throw FixedWidthError("error casting significand component to intergers" + significand_str[i]);
		}
		significand.push_back(si);			
	}
	//set the radix to position 0
	while (radix_pos > 0)
		shift_left();

	output_significand = significand;
	update();
	return;
}

vector<int> FixedWidthValue::number_to_digits(long long val)
{
	vector<int> iv;		
	do
	{
		long long digit = val % 10;
		iv.push_back((int)digit);
		val  /= 10;
	}while (val > 0);
	reverse(iv.begin(),iv.end());	
	return iv;
}

long long FixedWidthValue::digits_to_number(vector<int> iv)
{	
	long long val = 0 ;
	long long s;
	double e;	
	for (vector<int>::size_type i=0;i<iv.size();i++)
	{
		s = iv[i];
		e = iv.size() - i - 1;
		val += (s * pow(10,e));
	}
	return val;
}

TemplateParameter::TemplateParameter(string nm,double val, int start, int end,int lnum,bool isDbl, bool frcRad)
{
	
	name = nm;
	value = val;
	start_idx = start;
	end_idx = end;
	line_num = lnum;
	int width = end_idx - start_idx;
	isDoublePrecision = isDbl;
	forceRadix = frcRad;

}

bool TemplateParameter::operator==(const TemplateParameter &other) const
{
	if ((start_idx == other.start_idx) && (end_idx == other.end_idx) && (name == other.name))
	{
		return true;
	}
	return false;
}

void TemplateParameter::write_value(string &line)
{
	
	FixedWidthValue fwv(isDoublePrecision,forceRadix,end_idx - start_idx);
	string val_string = "";
	try
	{
		val_string = fwv.get_templatefile_string(value);	
	}
	catch (exception &e)
	{
		throw TemplateParameterError("error generating string representation of value for parameter "+name+" : "+e.what());
	}
	try
	{
		string::size_type istr = start_idx;
		string::size_type iend1 = end_idx;
		string::size_type iend2 = line.size();
		string start = line.substr(0,istr);
		string end = line.substr(iend1,iend2 - iend1);
		start.append(val_string);
		start.append(end);		
		line = start;
	}
	catch (exception &e)
	{
		throw TemplateParameterError("could not write value string "+val_string+" to line "+line+" : "+e.what());
	}
}




void TemplateFile::write_inputfile()
{
	ofstream ofs(input_filename);
	if (!ofs.is_open())
	{		
		throw TemplateFileError("unable to open model input file: "+input_filename);				
	}
	ifstream ifs(template_filename);
	if (!ifs.is_open())
	{
		throw TemplateFileError("unable to open model input file: "+template_filename);						
	}	
	string smarker;
	string tpl_typ;
	string line;
	//read the template file header
	getline(ifs,line);
	//read-write remaining lines
	int lnum = 2;
	while (getline(ifs,line))
	{				
		vector<int>::iterator it = find(line_numbers.begin(),line_numbers.end(),lnum);
		if (it != line_numbers.end())
		{
			line = build_input_line(line,parameter_line_index[lnum]);
		}		
		ofs << line << "\n";
		lnum++;
	}
}

string TemplateFile::build_input_line(string line,vector<TemplateParameter> line_parameters)
{		
	for (vector<TemplateParameter>::iterator tp=line_parameters.begin();tp != line_parameters.end();++tp)
	{		
		try
		{
			tp->write_value(line);
		}
		catch (exception &e)
		{
			throw TemplateFileError("error writing template file "+template_filename+" on line "+line+" : "+e.what());			
		}
	}
	return line;
}

vector<string> TemplateFile::get_parameter_names() 
{
	//returns a non-unique list
	vector<string> parameter_names;
	if (parameter_line_index.size() == 0)
	{
		return parameter_names;
	}
	vector<TemplateParameter>::iterator tp;
	for (map<int,vector<TemplateParameter>>::iterator tpv = parameter_line_index.begin();tpv != parameter_line_index.end();++tpv)
		for (tp = tpv->second.begin(); tp != tpv->second.end();++tp)
			parameter_names.push_back(string2upper(tp->get_name()));
	return parameter_names;
}

void TemplateFile::set_parameter_values(map<string,double> parameter_map)
{	
	for (map<int,vector<TemplateParameter>>::iterator tpv = parameter_line_index.begin();tpv!=parameter_line_index.end();++tpv)
	{
		for (vector<TemplateParameter>::iterator tp = tpv->second.begin();tp != tpv->second.end();++tp)
		{
			map<string,double>::iterator pm = parameter_map.find(string2upper(tp->get_name()));
			if (pm != parameter_map.end())
			{
				tp->set_value(pm->second);
			}
		}
	}	
}

void TemplateFile::read_templatefile()
{
	ifstream f(template_filename);
	if (!f.is_open())
	{
		throw TemplateFileError("unable to open template file: "+template_filename);
	}	
	string smarker;
	f >> template_filetype >> smarker;
	if (smarker.size() != 1)
	{		
		throw TemplateFileError("template file "+template_filename+" marker must be a single character, not "+smarker);	
	}
	marker = char(smarker[0]);

	if ((string2upper(template_filetype) != "PTF") && (string2upper(template_filetype) != "JTF"))			
	{
		throw TemplateFileError("template file "+template_filename+" must start with 'PTF' or 'JTF', not "+template_filetype);
	}
	//build a regex for the marker - fast
	regex rmarker(smarker);
	smatch mresults; 
	//read the remaining lines and error check markers
	string line;
	int lnum = 1;
	while (getline(f,line))
	{
		//check for marker in line
		if (regex_search(line,mresults,rmarker))
		{						
			//get marker indices on this line
			vector<pair<int,int>> indices = get_marker_indices(marker,line);	
			//get the parameter names within the markers
			vector<string> line_parameter_names = get_line_parameters(indices,line);
			vector<TemplateParameter> line_parameters;
			//for each parameter on this line, create TemplateParameter objects		
			for (vector<string>::size_type i=0;i<line_parameter_names.size();i++)
			{
				TemplateParameter tp(line_parameter_names[i],-1.0e+10,indices[i].first,indices[i].second,lnum,isDouble,forceRadix);				
				line_parameters.push_back(tp);
			}
			parameter_line_index.insert(pair<int,vector<TemplateParameter>>(lnum,line_parameters));
			line_numbers.push_back(lnum);
		}
		lnum++;
	}
	
}

vector<string> TemplateFile::get_line_parameters(vector<pair<int,int>> indices,string line)
{
	vector<string> line_parameters;	
	int npos;
	for (vector<pair<int,int>>::iterator index=indices.begin();index != indices.end();++index)	
	{
	    npos= index->second - index->first;
		string pname = line.substr(index->first,npos);			
		reduce_parameter_name(pname);
		line_parameters.push_back(pname);
	}
	if (line_parameters.size() == 0)
	{		
		throw TemplateFileError(" Internal Error - no valid parameter names found in " + template_filename + "on line " + line);
	}
	return line_parameters;
}

void TemplateFile::reduce_parameter_name(string &pname)
{
	pname = pname.erase(pname.size()-1,1);
	pname = pname.erase(0,1);
	for (size_t i=pname.size();i>0;i--)
	{		
		if ((pname[i-1] == ' ') || (pname[i-1] == '	'))
		{
			pname = pname.erase(i-1,1);
		}
	}
}

vector<pair<int,int>> TemplateFile::get_marker_indices(char marker,string line)
{
	warn_width = 5;
	size_t key = -1;
	size_t value = -1;
	vector<pair<int,int>> indices;
	for (string::size_type i=0;i < line.size();i++)
	{
		if (line[i] == marker)
		{						
			if (key == -1)
			{
				key = i;
			}
			else
			{				
				if (i == (key+1))
				{					
					throw TemplateFileError(" only one space between markers in file" + template_filename + "on line " + line);
				}
				else if (i-key+1<=warn_width)
				{
					cout << "warning - narrow parameter space in template file " << template_filename+" on line " << line << endl;
				}
				pair<int,int> p(key,i+1);
				indices.push_back(p);
				key = -1;
			}
		}
	}
	if (key != -1)
	{		
		throw TemplateFileError("unbalanced marker on line" + line);
	}
	return indices;
}



TemplateFiles::TemplateFiles(bool isDouble, bool forceRadix,vector<string> tpl_filenames,vector<string> ipt_filenames,vector<string> pnames)
{
	if (tpl_filenames.size() != ipt_filenames.size())
	{
		throw TemplateFileError(" number of template files not equal to number of input file");
	}
	template_filenames = tpl_filenames;
	input_filenames = ipt_filenames;
	parameter_names = pnames;
	for (vector<string>::size_type i=0;i<tpl_filenames.size();i++)
	{
		TemplateFile tpl(tpl_filenames[i],ipt_filenames[i],isDouble,forceRadix);
		tpl.read_templatefile();
		template_files.push_back(tpl);
	}
	//this one will throw errors...
	check_parameter_names();
	
}

void TemplateFiles::check_parameter_names()
{
	vector<bool> parameter_found;
	for (vector<string>::size_type i=0;i< parameter_names.size();i++)
	{
		parameter_found.push_back(false);
	}
	for (vector<string>::size_type i=0;i< parameter_names.size();i++)
	{
		string pname = parameter_names[i];
		for (vector<TemplateFile>::iterator tpl=template_files.begin();tpl != template_files.end();++tpl)
		{
			vector<string> tpl_names = tpl->get_parameter_names();
			if(find(tpl_names.begin(),tpl_names.end(),pname) != tpl_names.end())
			{
				parameter_found[i] = true;
			}
		}
	}
	bool missing = false;
	string missing_str = "";
	for (vector<string>::size_type i=0;i< parameter_names.size();i++)
	{
		if (!parameter_found[i])
		{
			cout << "\n\nError: parameter not found in template files: " << parameter_names[i] <<endl;
			missing_str += "parameter not found in template files: " + parameter_names[i] + " \n";
			missing = true;
		}
	}

	if (missing)
	{
		throw TemplateFileError(" missing parameters :"+missing_str);
	}
}

int TemplateFiles::writtpl(vector<double> parameter_values)
{	
	map <string,double> parameter_map;
	for (vector<double>::size_type i=0;i<parameter_values.size();i++)
	{
		pair <string,double> p(parameter_names[i],parameter_values[i]);
		parameter_map.insert(p);
	}
	for (vector<TemplateFile>::iterator tpl=template_files.begin();tpl != template_files.end();++tpl)
	{
		tpl->set_parameter_values(parameter_map);
		tpl->write_inputfile();
	}
	return 0;
}
