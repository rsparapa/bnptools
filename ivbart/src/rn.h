// Copyright (C) 2013-2015 Matthew T. Pratola, Robert E. McCulloch and Hugh A. Chipman - All Rights Reserved.
// You may NOT use, distribute or modify this code without the explicit written permission of all
// the above authors.
// mpratola@gmail.com
// robert.e.mculloch@gmail.com
// hughchipman@gmail.com
//
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
   virtual double gamma() = 0; 
   virtual double exp() = 0; 
   virtual void set_alpha(double alpha) = 0; //set df for chi-square
   virtual ~rn() {}
};

#endif
