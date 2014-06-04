#include "PolynomInt.h"

template<typename T>
class Polynom<T, typename std::enable_if<!std::is_integral<T>::value>::type>{
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

	Polynom<T> operator/(const Polynom<T>& o) const{return euclidian(o).first;}; //O((n-m)*m)	
	Polynom<T> operator/(const T&); //theta(n)
	Polynom<T>& operator/=(const Polynom<T>& o){return *this=*this/o;};

	Polynom<T> operator%(const Polynom<T>& o) const{return euclidian(o).second;}; //O((n-m)*m)	
	Polynom<T>& operator%=(const Polynom<T>& o){return *this=*this%o;};
	Polynom<T> operator()(T); //theta(n)
	Polynom<T> derivate(std::size_t); //theta(n)
	int border() const {return coefficients.size();};	
	
	bool operator==(const Polynom<T>&); //theta(n)
	bool operator!=(const Polynom<T>&); 
	bool operator<(const Polynom<T>&);	//O(1)
	bool operator>(const Polynom<T>&);

	template <typename U>
	friend std::ostream& operator<<(std::ostream& os,const Polynom<U>&); //theta(n)

private:
	std::vector<T> coefficients;
	static const char *super[];

	std::pair<Polynom<T>,Polynom<T>> euclidian(const Polynom<T>&) const;
	void align(Polynom<T>&,Polynom<T>&);
	void trimZeros();
	Polynom(); //theta(deg)
	Polynom<T> mult(int,T)const;
	std::string sup(int num) const;

};

