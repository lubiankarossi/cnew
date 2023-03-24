#include "storage.h"

#ifdef _WIN32_
#endif
#ifdef _LINUX_
#endif

//_________________________________________________________________________________________________________

storage::~storage()
{
	if (!(Yield.empty()) ){
		Yield.clear();
		count.clear();
	}
}

//_________________________________________________________________________________________________________

double storage::GetYield(int Zaid)
{
	if (Yield.find(Zaid) != Yield.end()) {
		return (Yield[Zaid]/count[Zaid]); 
	}
	else{
		return 0.0;
	}
}

//_________________________________________________________________________________________________________

void storage::AddYield(int Zaid, double M)
{
	if (Yield.find(Zaid) != Yield.end()) {
		(Yield[Zaid]) = double(M); 
		(count[Zaid]) += 1.0; 
	}
	else{
		Yield.insert(map<int,double>::value_type(Zaid, double(M) ));
		count.insert(map<int,double>::value_type(Zaid,       1.0 ));
	}
}

//_________________________________________________________________________________________________________

void storage::AddYield(int Zaid, double M, double w)
{
	if (Yield.find(Zaid) != Yield.end()) {
		(Yield[Zaid]) = w * double(M); 
		(count[Zaid]) += w; 
	}
	else{
		Yield.insert(map<int,double>::value_type(Zaid, w * double(M) ));
		count.insert(map<int,double>::value_type(Zaid,       w ));
	}
}

//_________________________________________________________________________________________________________

void storage::SumYield(int Zaid, double M)
{
	if (Yield.find(Zaid) != Yield.end()) {
		(Yield[Zaid]) += double(M); 
		(count[Zaid]) += 1.0; 
	}
	else{
		Yield.insert(map<int,double>::value_type(Zaid, double(M) ));
		count.insert(map<int,double>::value_type(Zaid,       1.0 ));
	}
}

//_________________________________________________________________________________________________________

void storage::SumYield(int Zaid, double M, double w)
{
	if (Yield.find(Zaid) != Yield.end()) {
		(Yield[Zaid]) += w * double(M); 
		(count[Zaid]) += w; 
	}
	else{
		Yield.insert(map<int,double>::value_type(Zaid, w * double(M) ));
		count.insert(map<int,double>::value_type(Zaid,           w   ));
	}
}

//_________________________________________________________________________________________________________

void storage::reset(){

	if (!(Yield.empty()) && !(count.empty()) ){
		Yield.clear();
		count.clear();
	}
}

//_________________________________________________________________________________________________________

bool storage::exist(int Zaid)
{
	if (Yield.find(Zaid) != Yield.end()) {
		return true;
	}
	else{
		return false;
	}
}

//_________________________________________________________________________________________________________


