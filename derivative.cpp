#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

double h = 0.01;
const double adjustment = 1.3;

double f(double x)
{
	return ( ( 1.0 + x ) * exp(-x*x) );
}

double f(int i)
{
	double x = static_cast<double>(i)*h;
	return ( ( 1.0 + x ) * exp(-x*x) );
}


int main()
{
	const int n_pts = 20;

	//vector<double> f;
	//for (int i = 0; i < n_pts; ++i)
	//	f.push_back(func(double(i)*h));

	//for (int i = 0; i < 10; ++i)
	//	cout << i << "   " << static_cast<double>(i)*h << "   " << f(i) << endl;

	// function value
	cout << "f(0) = " << f(0) << endl;
	// first derivative
	cout << "f'(0) = " << (-7129*f(0)+22680*f(1)-45360*f(2)+70560*f(3)-79380*f(4)+63504*f(5)-35280*f(6)+12960*f(7)-2835*f(8)+280*f(9))/(2520*1.0*h) << endl;
	h *= adjustment;

	// second derivative
	cout << "2nd deriv" << endl;
	for (int i = 0; i < 10; ++i)
		cout << i << "   " << f(i) << endl;

	cout << "f''(0) = " << (32575*f(0)-165924*f(1)+422568*f(2)-704368*f(3)+818874*f(4)-667800*f(5)+375704*f(6)-139248*f(7)+30663*f(8)-3044*f(9))/(5040*1.0*h*h) << endl;
	h *= adjustment;

	// third derivative
	cout << "f'''(0) = " << (-180920*f(0)+1145259*f(1)-3375594*f(2)+6095796*f(3)-7392546*f(4)+6185970*f(5)-3540894*f(6)+1328724*f(7)-295326*f(8)+29531*f(9))/(15120*1.0*h*h*h) << endl;
	h *= adjustment;

	// fourth derivative
	cout << "f''''(0) = " << (-1062*f(-7)+20137*f(-6)-186930*f(-5)+1150452*f(-4)-5435782*f(-3)+22035519*f(-2)-53626914*f(-1)+72089160*f(+0)-53626914*f(+1)+22035519*f(+2)-5435782*f(+3)+1150452*f(+4)-186930*f(+5)+20137*f(+6)-1062*f(+7))/(4989600*1.0*h*h*h*h) << endl;
	h *= adjustment;

	// fifth derivative
	cout << "f'''''(0) = " << (-518*f(-7)+8301*f(-6)-62710*f(-5)+295244*f(-4)-944862*f(-3)+1819681*f(-2)-1718382*f(-1)+0*f(+0)+1718382*f(+1)-1819681*f(+2)+944862*f(+3)-295244*f(+4)+62710*f(+5)-8301*f(+6)+518*f(+7))/(181440*1.0*h*h*h*h*h) << endl;
	h *= adjustment;

	// sixth derivative
	cout << "f''''''(0) = " << (148*f(-7)-2767*f(-6)+25084*f(-5)-147622*f(-4)+629908*f(-3)-1819681*f(-2)+3436764*f(-1)-4243668*f(+0)+3436764*f(+1)-1819681*f(+2)+629908*f(+3)-147622*f(+4)+25084*f(+5)-2767*f(+6)+148*f(+7))/(60480*1.0*h*h*h*h*h*h) << endl;
	h *= adjustment;

	// seventh derivative
	cout << "f'''''''(0) = " << (311*f(-7)-4848*f(-6)+34975*f(-5)-151232*f(-4)+405219*f(-3)-655792*f(-2)+552891*f(-1)+0*f(0)-552891*f(1)+655792*f(2)-405219*f(3)+151232*f(4)-34975*f(5)+4848*f(6)-311*f(7))/(17280*h*h*h*h*h*h*h) << endl;
	h *= adjustment;

	// eighth derivative
	cout << "f''''''''(0) = " << (17311*f(-8)-351616*f(-7)+3434760*f(-6)-21445760*f(-5)+95023460*f(-4)-302537088*f(-3)+689491768*f(-2)-1126894720*f(-1)+1326523770*f(0)-1126894720*f(1)+689491768*f(2)-302537088*f(3)+95023460*f(4)-21445760*f(5)+3434760*f(6)-351616*f(7)+17311*f(8))/(3628800*1.0*h*h*h*h*h*h*h*h) << endl;

	cout << h << endl;

	return 0;	
}
