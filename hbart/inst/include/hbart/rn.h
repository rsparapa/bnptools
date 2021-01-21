// rn.h: Random number generator virtual class for using BART in R package.
// Copyright (C) 2012-2016 Matthew T. Pratola, Robert E. McCulloch and Hugh A. Chipman
//
// This file is part of BART.
//
// BART is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// BART is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Author contact information
// Matthew T. Pratola: mpratola@gmail.com
// Robert E. McCulloch: robert.e.mculloch@gmail.com
// Hugh A. Chipman: hughchipman@gmail.com


#ifndef GUARD_rn
#define GUARD_rn

//pure virtual base class for random numbers
class rn
{
public:
   virtual double normal() = 0; //standard normal
   virtual double uniform() = 0; //uniform(0,1)
   virtual double chi_square() = 0; //chi-square
   virtual void set_df(int df) = 0; //set df for chi-square
   virtual ~rn() {}
};

#endif
