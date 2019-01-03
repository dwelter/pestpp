/*
© Copyright 2017, Chas Egan

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

#include "picosha2.h"

using std::string;

namespace hmacsha2 {

	string xor(const string& key, const string& pad) {
		string answer = pad;
		int n = (key.size() < pad.size()) ? key.size() : pad.size();
		for (int i = 0; i < n; i++)
			answer[i] = key[i] ^ pad[i];
		return answer;
	}


	string hmac(const string& message, const string& key) {
		//Following https://en.wikipedia.org/wiki/HMAC we will calculate a Hashed Message Authentication Code (HMAC)
		//for the given message using the supplied key. The process is H(outerKey + H(innerKey + message)), where
		//H(x) is a suitable hash fucntion. We will use a SHA-256 as the hash function. InnerKey and OuterKey are 
		//different, but both are derived from the key supplied in the method call. The HMAC provides a secure 
		//way to verify the authenticity of the message because perpetrators are unable to modify the message and
		//provide a matching HMAC without knowing the key.

		//Hashing the key to make sure it is 64 bytes long, and then calculating inner_key and outer_key. Note 
		//SHA-256 has a 64 byte block size.
		string hashedKey = picosha2::hash256_hex_string(key);
		string innerKey = xor(hashedKey, "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"); //64 "0x5C" values
		string outerKey = xor(hashedKey, "6666666666666666666666666666666666666666666666666666666666666666"); //64 "0x36" values

		//H(outer_key + H(inner_key + message))
		string innerHash = picosha2::hash256_hex_string(innerKey + message);
		string outerHash = picosha2::hash256_hex_string(outerKey + innerHash);
		return outerHash;
	}
}