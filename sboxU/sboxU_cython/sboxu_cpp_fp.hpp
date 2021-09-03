/* Time-stamp: <2021-07-01 16:50:12 mjoly>
 *
 * LICENSE
 */ 

#ifndef _SBOXU_CPP_FP_H_
#define _SBOXU_CPP_FP_H_

#include "sboxu_cpp.hpp"

Integer power(Integer p, BinWord n);
std::vector<Integer> add(std::vector<Integer> a, std::vector<Integer> b, Integer p);
Integer sub(std::vector<Integer> a, std::vector<Integer> b, Integer p);

class CppFptGenerator
{
	private :
		Integer p ;
		Integer t ;
		std::vector<Integer> vect ;
		bool not_finished ;

	public :
		CppFptGenerator()
		{}

		CppFptGenerator(Integer p, BinWord t)
		{
			this->p = p ;
			this->t = t ;
			this->vect = std::vector<Integer>(t,0) ;
			this->not_finished = true ;
		}

		std::vector<Integer> getVector()
		{
			return this->vect ;
		}


		bool is_not_finished()
		{
			return not_finished ;
		}

		void next()
		{
			unsigned int i = this->vect.size() ;
			do
			{
				i -- ;
				this->vect[i] +=1 ;
				if (this->vect[i] == this->p) 
				{
					this->vect[i] = 0 ;
					if (i == 0) this->not_finished = false ;
				}
			}
			while (i > 0 && this->vect[i] == 0) ;
		}

};

class CppFptFunction
{
	public:
		Integer p ;
		BinWord t ;
		BinWord u ;
		std::vector<Integer> powers ;
		std::vector<std::vector<Integer>> indexes ;
		std::vector<std::vector<Integer>> values ;

		CppFptFunction()
		{
			this->p = 2 ;
			this->t = 8 ;
			this->u = 8 ;
		}

		CppFptFunction(Integer p, BinWord t, BinWord u, std::vector<std::vector<Integer>> indexes,  std::vector<std::vector<Integer>> values)
		{
			this->p = p ;
			this->t = t ;
			this->u = u ;
			this->powers = std::vector<Integer>(t,1) ;
			for (unsigned int i = t ; i > 1 ; i --)
			{
				this->powers[i-2] = p * this->powers[i-1] ;
			}
			this->indexes = indexes ;
			this->values = values ;

		}


		std::vector<Integer> eval(std::vector<Integer> index)
		{
			Integer i = std::inner_product(index.begin(), index.end(), this->powers.begin(), 0) ;
			return this->values[i] ;
		}

		std::vector<Integer> fpt_ddt_row(Integer a);

		std::vector<std::vector<Integer>> fpt_ddt();

		std::map<Integer,Integer> fpt_differential_spectrum_fast(const unsigned int n_threads);
		
};

void fpt_ddt_rows_count(std::map<Integer,Integer>& result, CppFptFunction S, const BinWord a_min, const BinWord a_max);


#endif /* _SBOXU_CPP_FP_H_ */