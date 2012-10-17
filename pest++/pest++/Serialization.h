/*  
    © Copyright 2012, David Welter
    
    This file is part of PEST++.
   
    PEST++ is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PEST++ is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PEST++.  If not, see<http://www.gnu.org/licenses/>.
*/

#ifndef SERIALIZE_H_
#define SERIALIZE_H_

#include <vector>

class Transformable;
class Parameters;
class Observations;

class Serialization
{
public:
	static std::vector<char> serialize(const Transformable &tr_data);
	static std::vector<char> serialize(const std::vector<const Transformable*> tr_vec);
	static std::vector<char> serialize(const std::vector<Transformable*> &tr_vec);
	static std::vector<char> serialize(const Parameters &pars, const Observations &obs);
	static std::vector<char> serialize(const std::vector<std::string> &string_vec);
	static unsigned unserialize(const std::vector<char> ser_data, Transformable &tr_data, unsigned ser_data_loc=0);
	static unsigned unserialize(const std::vector<char>, std::vector<Transformable*> &tr_vec);
    static unsigned unserialize(const std::vector<char> data, Parameters &pars, Observations &obs);
	static unsigned unserialize(const std::vector<char> ser_data, std::vector<std::string> &string_vec);
private:
};


#endif /* SERIALIZE_H_ */




