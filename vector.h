#include <math.h>
#include <iostream>

class Vector3
{
public:
	double x,y,z;
	Vector3( double x=0.0, double y=0.0, double z=0.0 ){
		Vector3::x = x;
		Vector3::y = y;
		Vector3::z = z;
	}
	double mag(void){
		return sqrt(x*x + y*y + z*z);
	}
	double normalize( void ){
		double m = mag();
		if(m==0.0){
			x=1.0;
			y=0.0;
			z=0.0;
		}else{
			*this/=m;
		}
		return m;
	}
	void print(const char* p_str){
		std::cout << p_str << "\n";
		std::cout << "x:" << x << " y:" << y << " z:" << z << "\n";
	}
	Vector3& operator=(const Vector3& a){
		x = a.x;
		y = a.y;
		z = a.z;
		return *this;
	}
	Vector3 operator+(void){
		return *this;
	}
	Vector3 operator-(void){
		Vector3 l_r(-x,-y,-z);
		return l_r;
	}
	Vector3& operator/=(const double a){
		x/=a;
		y/=a;
		z/=a;
		return *this;
	}
	Vector3& operator+=(const Vector3& a){
		x += a.x;
		y += a.y;
		z += a.z;
		return *this;
	}
	Vector3& operator-=(const Vector3& a){
		x -= a.x;
		y -= a.y;
		z -= a.z;
		return *this;
	}
	friend Vector3 operator+(const Vector3& a, const Vector3& b);
	friend Vector3 operator-(const Vector3& a, const Vector3& b);
	friend double  operator*(const Vector3& a, const Vector3& b);
	friend Vector3 operator*(double a, const Vector3& b);
	friend Vector3 operator*(const Vector3& a, double b);
	friend Vector3 operator^(const Vector3& a, const Vector3& b);
};

Vector3 operator+(const Vector3& a,const Vector3& b)
{
	Vector3 l_r(a.x+b.x,a.y+b.y,a.z+b.z);
	return l_r;
}

Vector3 operator-(const Vector3& a, const Vector3& b)
{
	Vector3 l_r(a.x-b.x,a.y-b.y,a.z-b.z);
	return l_r;
}

double operator*(const Vector3& a, const Vector3& b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

Vector3 operator*(double a, const Vector3& b)
{
	Vector3 l_r(b.x*a,b.y*a,b.z*a);
	return l_r;
}

Vector3 operator*(const Vector3& a, double b)
{
	Vector3 l_r(a.x*b,a.y*b,a.z*b);
	return l_r;
}

Vector3 operator^(const Vector3& a, const Vector3& b)
{
	Vector3 l_r(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);
	return l_r;
}
