#ifndef __STORAGE_HH
#define __STORAGE_HH

#ifdef _WIN32_
#endif
#ifdef _LINUX_
#endif

#include <iostream>
#include <map>
#include <cmath>

#include <string>
#include <fstream>


using namespace std;

/// storage é uma classe, não muito especializada, criada para guardar valores 
//_________________________________________________________________________________________________________

class storage {

public:

    typedef std::map<int, double>::iterator iterator;
    typedef std::map<int, double>::const_iterator const_iterator;

    iterator       begin()       { return Yield.begin(); }
    const_iterator begin() const { return Yield.begin(); }
    iterator         end()       { return Yield.end(); }
    const_iterator   end() const { return Yield.end(); }


/// Construtor Padrão.
//_________________________________________________________________________________________________________

	storage(){}

/// Destrutor Padrão.
//_________________________________________________________________________________________________________

	~storage();

/// Limpa Conteiners e zera variaveis.
//_________________________________________________________________________________________________________

	void reset();

/// Guarda.
//_________________________________________________________________________________________________________
	
	void AddYield(int Zaid, double M);

	void AddYield(int Zaid, double M, double w);

        void SumYield(int Zaid, double M);

        void SumYield(int Zaid, double M, double w);

        double GetYield(int Zaid);

        bool exist(int Zaid);

/// Operadores.
//_________________________________________________________________________________________________________

        storage& operator=( storage& stor){

                reset();

		//map<int, double>::iterator imap;
                storage::iterator imap;

		for (imap = stor.begin(); imap != stor.end(); imap++ ){
                       Yield.insert(map<int,double>::value_type(imap->first, imap->second/stor.count[imap->first] ));
                       count.insert(map<int,double>::value_type(imap->first,                                  1.0 ));
		}

                return *this; 

        }


	friend ostream& operator<<(ostream& out, storage& stor){
		
		out.setf(ios::scientific, ios::floatfield);
		

		out << "     Prodction [Nuclide(at./barn/cm)]    \n\n";
		out << "		ZAID		Yield	     \n\n";

		out.adjustfield;

		map<int, double>::iterator imap;
		
		for (imap = stor.Yield.begin(); imap != stor.Yield.end(); imap++ ){
			out 	<< "		"<< imap->first 
				<< "		"<< imap->second/stor.count[imap->first] 
				<< "\n";
		}

		return out;
	}

/// Sobrecarga de operador << para saida de disco.
//_________________________________________________________________________________________________________

	friend ofstream& operator<<(ofstream& out , storage& stor){

		out.setf(ios::scientific, ios::floatfield);
		
		out << "\n      \n\n";


		out << "\n\n";

		out << "     Prodction [Nuclide(at./barn/cm)]    \n\n";
		out << "		ZAID		Yield	     \n\n";

		out.adjustfield;

		map<int, double>::iterator imap;
		
		for (imap = stor.Yield.begin(); imap != stor.Yield.end(); imap++ ){
			out 	<< "		"<< imap->first 
				<< "		"<< imap->second/stor.count[imap->first] 
				<< "\n";
		}

		return out;
	}


private:

	/// Mapa para armazenar Yield (ZZAAA)
	map<int, double> Yield;
	map<int, double> count;

};
#endif

