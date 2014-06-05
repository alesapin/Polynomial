#include <iostream>
#include <vector>
#include <algorithm>
#include <complex>
#include <cstdarg>
#include <cstdlib>
#include <fftw3.h>	
#include <utility>
#include <initializer_list>
#include <type_traits>
#include <stdexcept> 
#include <typeinfo>
#include <fstream>
#define REAL 0
#define IMAG 1
typedef fftw_complex complex;
template< typename T,class Enable = void>
class Polynom;
template<typename T>
class Polynom<T, typename std::enable_if<std::is_integral<T>::value>::type>{
public:
	Polynom(const Polynom<T>&); //theta(n)
	Polynom(int size,T *data);
	Polynom(std::initializer_list<T> l):coefficients(l){}
	Polynom(Polynom<T>&&);
	
	Polynom<T>& operator=(Polynom<T>&&);
	Polynom<T>& operator=(const Polynom<T>&); //theta(n)
	
	Polynom<T> operator+(Polynom<T>); //theta(n)
	Polynom<T> operator+(); //O(1)
	Polynom<T>& operator+=(const Polynom<T>& o){return *this=*this+o;}; 

	Polynom<T> operator-(); //theta(n)
	Polynom<T> operator-(Polynom<T>); //theta(n)
	Polynom<T>& operator-=(const Polynom<T>& o){return *this=*this-o;};

	Polynom<T> operator*(const Polynom<T>&);
	Polynom<T>& operator*=(const Polynom<T>& o){return *this=*this*o;};
	Polynom<T> operator*(const T&);	//theta(n)

	Polynom<T> operator/(const Polynom<T>& o) const {return euclidian(o).first;};
	Polynom<T> operator/(const T&)const; //theta(n)
	Polynom<T>& operator/=(const Polynom<T>& o){return *this=*this/o;};
	
	Polynom<T>& operator%=(const Polynom<T>& o){return *this=*this%o;};
	Polynom<T> operator%(const Polynom<T>& o) const{return euclidian(o).second;};

	bool operator==(const Polynom<T>&); //theta(n)
	bool operator!=(const Polynom<T>&); 
	bool operator<(const Polynom<T>&);	//O(1)
	bool operator>(const Polynom<T>&);

	void reduce(){reduce(gcd());};
	int border() const {return coefficients.size();};	
	T operator()(const T&); //theta(n)
	Polynom<T> derivate(std::size_t); //theta(n)
	
	template <typename U>
	friend std::ostream& operator<<(std::ostream& os,const Polynom<U>&); //theta(n)
	
	template <typename U,typename std::enable_if<std::is_integral<U>::value,int>::type>
	friend Polynom<U> gcd(const Polynom<U>& first, const Polynom<U>& second);

	template <typename U,typename std::enable_if<std::is_integral<U>::value,int>::type>
	friend Polynom<U> lcm(const Polynom<U>& first, const Polynom<U>& second);
	
private:
	std::vector<T> coefficients;
	static const char *super[];

	void align(Polynom<T>&,Polynom<T>&);
	void toPowerOfTwo();
	complex* multDFT(const complex*,int,const complex*,int);
	void trimZeros();

	Polynom<T> mult(int,T)const;
	std::pair<Polynom<T>,Polynom<T>> euclidian(const Polynom<T>&) const;
	Polynom(); 
	std::string sup(int num) const;
	void reduce(T num){for(T& x:coefficients)x/=num;};
	T gcd(T a,T b){return b==0?a:gcd(b,a%b);};
	T gcd(){
		T result=coefficients[0];
		for(int i=1;i<coefficients.size();++i) result = gcd(result,coefficients[i]);
		return result;
	}
};
