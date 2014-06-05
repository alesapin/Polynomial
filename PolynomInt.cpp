#include "PolynomInt.h"
template<typename T>
const char* Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
super[] = {"\xe2\x81\xb0", "\xc2\xb9", "\xc2\xb2","\xc2\xb3",
 	"\xe2\x81\xb4", "\xe2\x81\xb5", "\xe2\x81\xb6","\xe2\x81\xb7",
  	"\xe2\x81\xb8", "\xe2\x81\xb9"};

template<typename T> Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
Polynom(){
		coefficients.resize(1);
}

template<typename T> Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
Polynom(Polynom<T>&& other){
	coefficients=std::vector<T>();
	coefficients=other.coefficients;
}

template<typename T> Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
Polynom(const Polynom<T>& other){
	coefficients.resize(other.border());
	for(int i=0;i<other.border();++i){
		coefficients[i]=other.coefficients[i];
	}
}

template<typename T> Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
Polynom(int size,T *data){
	coefficients.resize(size);
	for(int i=0;i<size;++i){
		coefficients[i]=data[i];
	}
}
template<typename T> 
T Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
operator()(const T& n){
	T result=coefficients[border()-1];
	for(int i=border()-2;i>=0;--i){
		result=result*n+coefficients[i];
	}
	return result;
}
template<typename T> 
Polynom<T>& Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
operator=(const Polynom<T>& a){
	if(a.coefficients.size()==0){
		coefficients=std::vector<T>();
		coefficients.resize(1);
		return *this;
	} 
	coefficients.resize(a.border());
	for(int i=0;i<a.coefficients.size();++i){
			coefficients[i]=a.coefficients[i];
	}
	return *this;
} 
template<typename T> 
Polynom<T>& Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
operator=(Polynom<T>&& other){
	coefficients=std::vector<T>();
	coefficients=other.coefficients;
	return *this;
}

template<typename U> 
std::ostream& operator<<(std::ostream& os,const Polynom<U>& a){
	int i=a.border()-1;
	if(i==0){
		os<<a.coefficients[i];
		return os;
	}
	else if(i==1){ 
		if(a.coefficients[i]==-1) os<<"-";
		else if(a.coefficients[i]!=1) os<<a.coefficients[i];
		os<<'x';
		if(a.coefficients[i-1]>0) os<<" + ";
		if(a.coefficients[i-1]<0) os<<" - ";
		os<<abs(a.coefficients[i]);
		return os;
	}
	if(a.coefficients[i]==-1) os<<"-";
	else if(a.coefficients[i]!=1) os<<a.coefficients[i];
	if(typeid(os)==typeid(std::ofstream)){
		os<<"x^"<<i;
		--i;
		for(;i>1;--i){
			if(a.coefficients[i]>0) os<<" + ";
			else if(a.coefficients[i]<0) os<<" - ";
			if(abs(a.coefficients[i])!=1 && a.coefficients[i]!=0){
			 	os<<abs(a.coefficients[i]);
				os<<"x^"<<i;
			}
		}
	}else{
		os<<'x'<<a.sup(i);
		--i;
		for(;i>1;--i){
			if(a.coefficients[i]>0) os<<" + ";
			else if(a.coefficients[i]<0) os<<" - ";
			if(abs(a.coefficients[i])!=1 && a.coefficients[i]!=0){
			 	os<<abs(a.coefficients[i]);
				os<<"x"<<a.sup(i);
			}
		}
	}
	if(i>0){
		if(a.coefficients[i]>0) os<<" + ";
		else if(a.coefficients[i]<0) os<<" - ";
		if(abs(a.coefficients[i]!=1)) os<<abs(a.coefficients[i--]);
		os<<'x';
	}
	if(a.coefficients[i]>0) os<<" + "<<a.coefficients[i];
	else if(a.coefficients[i]<0) os<<" - "<<abs(a.coefficients[i]);
	return os;
}

template<typename T> 
void Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
align(Polynom<T>& f,Polynom<T>& s){
	int bigger=std::max(f.border(),s.border());
	f.coefficients.resize(bigger);
	s.coefficients.resize(bigger);
}

template<typename T> 
Polynom<T> Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
operator+(Polynom<T> other){
	align(*this,other);
	for(int i=0;i<other.border();++i){
		other.coefficients[i]+=coefficients[i];
	}
	trimZeros();
	other.trimZeros();
	return other;
}

template<typename T> 
Polynom<T> Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
operator+(){
	return *this;
}

template<typename T> 
Polynom<T> Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
operator*(const T& a){
	Polynom<T> result(*this);
	for(T& x:result.coefficients){
		x*=a;
	}
	return result;
}

template<typename T> 
Polynom<T> Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
operator-(Polynom<T> other){
	return *this+(other*-1);
}

template<typename T> 
Polynom<T> Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
operator-(){
	Polynom<T> result(*this);
	return result*-1;
}

template<typename T> 
complex* Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
multDFT(const complex* f,int fs,const complex* s,int ss){
	complex* result=new complex[fs+ss];
	for(int i=0;i<fs+ss;++i){
		result[i][REAL]=f[i][REAL]*s[i][REAL]-f[i][IMAG]*s[i][IMAG];
		result[i][IMAG]=f[i][REAL]*s[i][IMAG]+f[i][IMAG]*s[i][REAL];
	}
	return result;
}

template<typename T> 
void Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
toPowerOfTwo(){
	for(int i=coefficients.size();i & (i-1);++i){
		coefficients.push_back(0);
	}
}

template<typename T> 
Polynom<T> Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
operator*(const Polynom<T>& second){
	int resultSize=border()+second.border()-1;
	Polynom<T> result(second);
	result.coefficients.resize(resultSize);
	result.toPowerOfTwo();
	align(*this,result);
	std::vector<double> f(coefficients.begin(),coefficients.end());
	std::vector<double> s(result.coefficients.begin(),result.coefficients.end());
	double* one=&f[0];
	double* two=&s[0];
	complex* result1=new complex[border()];
	complex* result2=new complex[border()];
	fftw_plan plan1 = fftw_plan_dft_r2c_1d(border(),one,result1,FFTW_ESTIMATE);
    fftw_plan plan2 = fftw_plan_dft_r2c_1d(border(),two,result2,FFTW_ESTIMATE);
    fftw_execute(plan1);
	fftw_execute(plan2);
	complex* res=multDFT(result1,border(),result2,border());
	fftw_plan plan3 = fftw_plan_dft_c2r_1d(border(),res,two,FFTW_ESTIMATE);
	fftw_execute(plan3);
	for(int i=0;i<border();++i) result.coefficients[i]=floor(two[i]/border()+0.5);
	result.coefficients.resize(resultSize);
	fftw_destroy_plan(plan1);
	fftw_destroy_plan(plan2);
    fftw_destroy_plan(plan3);
    delete [] result1;
    delete [] result2;
    this->trimZeros();
	return result;
}

template<typename T> 
void Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
trimZeros(){
	int i=border()-1;
	while(coefficients.size()>1 && !coefficients[i]){
		coefficients.pop_back();
		--i;
	}
}

template<typename T> 
Polynom<T> Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
operator/(const T& num) const{
	Polynom<T> result(*this);
	for(T& x:result.coefficients){
		if(x%num==0)x/=num;
		else throw std::invalid_argument("Wrong division");
	}
	return result;
}

template<typename T> 
std::pair<Polynom<T>,Polynom<T>> Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
euclidian(const Polynom<T>& other) const{
	if(other.border()==1) return std::pair<Polynom<T>,Polynom<T>>(*this/other.coefficients[0],Polynom<T>{0});
	Polynom<T> r;
	r.coefficients=std::vector<T>(coefficients.begin(),coefficients.end());
	Polynom<T> sec;
	sec.coefficients=std::vector<T>(other.coefficients.begin(),other.coefficients.end());
	Polynom<T> q;
	T h=r.coefficients[0];
	q.coefficients.resize(border()-1);
	int d=sec.border();
	T c=sec.coefficients[d-1];
	int rsize=r.border();
	while(r.border()>=d){
		T s=r.coefficients[r.border()-1]/c;
		int deg=r.border()-d;
		q.coefficients[deg]=s;
		r=r-sec.mult(deg,s);
		if(rsize==r.border() && r.border()!=1) throw std::invalid_argument("Wrong division");
		else rsize=r.border();
	}
	q.trimZeros();
	r.trimZeros();
	return std::pair<Polynom<T>,Polynom<T>>(q,r);
}

template<typename T> 
bool Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
operator==(const Polynom<T>& other){
	if(border()!=other.border()) return false;
	for(int i=0;i<border();++i){
		if(coefficients[i]!=other.coefficients[i]) return false;
	}
	return true;
}
template<typename U,typename std::enable_if<std::is_integral<U>::value,int>::type=0> 
Polynom<U> gcd(const Polynom<U>& first, const Polynom<U>& second){
	int f=first.border();
	int s=second.border();
	if(f>s || (f==s && first.coefficients[f-1]<second.coefficients[s-1])) return gcd(second,first);
	Polynom<U> u(first);
	Polynom<U> v(second);
	Polynom<U> t;
	try{
		while(v.border()>0 && v.coefficients[0]){
			t=u;
			u=v;
			v=t%u;
			if(v.border()>1) v.reduce();
		}
		return u;
	}catch(const std::invalid_argument& a){
		return Polynom<U>({1});
	}
}
template<typename U,typename std::enable_if<std::is_integral<U>::value,int>::type=0> 
Polynom<U> lcm(const Polynom<U>& first, const Polynom<U>& second){
	return (first/gcd(first,second))*second;
}

template<typename T> 
Polynom<T> Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
mult(int deg,T num)const {
	Polynom<T> result;
	result.coefficients.resize(border()+deg);
	int j=0;
	for(int i=deg;i<border()+deg;++i,++j){
		result.coefficients[i]=coefficients[j]*num;
	}
	return result;
}

template<typename T> 
bool Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
operator!=(const Polynom<T>& other){
	return !(*this==other);
} 
template<typename T>
bool Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
operator<(const Polynom<T>& other){
	return border()<other.border();
}

template<typename T>
bool Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
operator>(const Polynom<T>& other){
	return border()>other.border();
}

template<typename T>
Polynom<T> Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
derivate(std::size_t num){
	if(border()-num<=0) return Polynom<T>();
	Polynom<T> result;
	result.coefficients.resize(border()-num);
	for(int i=border()-num-1;i>=0;--i){
		result.coefficients[i]=coefficients[i+num];
		for(int j=0;j<num;++j){
			result.coefficients[i]*=(i+num-j);
		}
	}
	return result;
}

template<typename T>
std::string Polynom<T,typename std::enable_if<std::is_integral<T>::value>::type>::
sup(int num) const{
    std::string result;
    do{ 
        result.insert(0,super[num%10]);
    }while(num/=10);
    return result;
}
int main(){
	Polynom<long> a{-3,5,-6,1};
	Polynom<long> c{1,5,4,67,8,3,4,5,654,24};
	Polynom<long> b=a;
	std::cout<<(b==a)<<'\n';
	b=b*2;
	std::cout<<a(6)<<'\n';
	return 0;
}