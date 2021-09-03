#include "sboxu_cpp_fp.hpp"

Integer power(Integer p, BinWord n)
{
	if (n == 0) 
	{
		return 1 ;
	}
	else 
	{
		if (n % 2 == 0) return power(p*p,n/2) ; else return p * power(p*p,(n-1)/2) ;
	}
}
std::vector<Integer> add(std::vector<Integer> a, std::vector<Integer> b, Integer p)
{
	BinWord n = a.size(), m = b.size();
	std::vector<Integer> c(n,0) ;
	if (n != m)
	{ std::cout << "Size Error add : a " << n << ", b " << m ; }
	else
	{
		for (unsigned int i = 0 ; i < n; i++)
		{
			 c[i] = ( a[i] + b[i] ) % p ;
		}	
	}
	return c ;
}

Integer sub(std::vector<Integer> a, std::vector<Integer> b, Integer p)
{
	Integer k = 0 ;
	BinWord n = a.size(), m = b.size();
	if (n != m)
	{ std::cout << "Size Error add : a " << n << ", b " << m ; }
	else
	{
		Integer factor = 1;
		for (int i = n-1; i >= 0; i--)
		{
			k += ( (a[i] + p - b[i]) % p ) * factor ;
			factor *= p ;
		}	
	}
	return k ;
}

std::vector<Integer> CppFptFunction::fpt_ddt_row(Integer a)	
{
	Integer n = this->indexes.size() ;
	std::vector<Integer> row(n,0) ;
	for( unsigned int b = 0;b < n; b++ )
	{
		Integer k = sub(this->eval(add(this->indexes[a],this->indexes[b],this->p)),this->values[b],this->p) ;
		row[k] ++ ;
	}
	return row ;
}

std::vector<std::vector<Integer>> CppFptFunction::fpt_ddt()	
{
	Integer n = this->indexes.size() ;
	std::vector<std::vector<Integer>> ddt ;
	for (unsigned int a = 1 ; a < n; a ++)
	{
		ddt.push_back( this->fpt_ddt_row(a) ) ;
	}
	return ddt ;
}

void fpt_ddt_rows_count(std::map<Integer,Integer>& result, CppFptFunction S, const BinWord a_min, const BinWord a_max)
{
	for (unsigned int a=a_min; a<a_max; a++)
	{
	std::vector<Integer> row(S.fpt_ddt_row(a)); 
	for (auto &v : row)
		result[v] ++;
	}
}

std::map<Integer,Integer> CppFptFunction::fpt_differential_spectrum_fast( const unsigned int n_threads)
{
	std::map<Integer,Integer> count;
	BinWord n = this->values.size() ;
	if (n_threads == 1)
	{
		// small S-Box
		fpt_ddt_rows_count(std::ref(count),*this, 1, n);
	}
	else
	{
		std::vector<std::thread> threads;
		std::vector<std::map<Integer,Integer> > local_counts(n_threads);
		unsigned int slice_size = n/n_threads;
		for (unsigned int i=0; i<n_threads; i++)
		{
			unsigned int
			lower_bound = i*slice_size,
			upper_bound = (i+1)*slice_size;
			if (lower_bound == 0)
			lower_bound = 1;
			if (upper_bound > n)
			upper_bound = n;
			threads.push_back(std::thread(std::ref(fpt_ddt_rows_count), std::ref(local_counts[i]), *this, lower_bound, upper_bound));
		}
		for (unsigned int i=0; i<n_threads; i++)
		{
			threads[i].join();
			for (auto &entry : local_counts[i])
				count[entry.first] += entry.second ;
		}
	}
	return count;
}
