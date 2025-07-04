/*
Multi-Component Mesh Generator
Copyright (C) 2025 Okinawa Institute of Science and Technology, Japan.

Developer: Weiliang Chen  & Jules Lallouette

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef MULTICOMPMESHER_UTILITY_H
#define MULTICOMPMESHER_UTILITY_H
#include <string>
#include <vector>

inline bool append_str_before_suffix(std::string& str, const std::string& append_str) {
    size_t start_pos = str.rfind(".");
    if (start_pos == std::string::npos)
        return false;
    str.insert(start_pos, append_str);
    return true;
}

inline bool ends_with(std::string const& value, std::string const& ending) {
    if (ending.size() > value.size())
        return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

// https://thispointer.com/c-how-to-read-a-file-line-by-line-into-a-vector/
inline bool getFileContent(std::string fileName, std::vector<std::string>& vecOfStrs) {
    // Open the File
    std::ifstream in(fileName.c_str());

    // Check if object is valid
    if (!in) {
        std::cerr << "Cannot open the File : " << fileName << std::endl;
        return false;
    }

    std::string str;
    while (std::getline(in, str)) {
        if (str.size() > 0 && str.find("#") != 0) {
            vecOfStrs.push_back(str);
        }
    }
    in.close();
    return true;
}
#endif
