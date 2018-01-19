//////#include <fstream>
//////#include <iostream>
//////using namespace std;
//////
//////int main() {
//////
//////	char data[100];
//////
//////	// open a file in write mode.
//////	ofstream outfile;
//////	outfile.open("afile.dat");
//////
//////	cout << "Writing to the file" << endl;
//////	cout << "Enter your name: ";
//////	cin.getline(data, 100);
//////
//////	// write inputted data into the file.
//////	outfile << data << endl;
//////
//////	cout << "Enter your age: ";
//////	cin >> data;
//////	cin.ignore();
//////
//////	// again write inputted data into the file.
//////	outfile << data << endl;
//////
//////	// close the opened file.
//////	outfile.close();
//////
//////	// open a file in read mode.
//////	ifstream infile;
//////	infile.open("afile.dat");
//////
//////	cout << "Reading from the file" << endl;
//////	infile >> data;
//////
//////	// write the data at the screen.
//////	cout << data << endl;
//////
//////	// again read the data from the file and display it.
//////	infile >> data;
//////	cout << data << endl;
//////
//////	// close the opened file.
//////	infile.close();
//////	system("pause");
//////
//////	return 0;
//////}
////
//////#include <iostream>
//////using namespace std;
//////
//////int main()
//////{
//////	char character1 = 'a';
//////	char *character2 = "phrase";
//////	char *character3[] = { "array os strings","w","dfdf","ewewe" }	;
//////
//////	cout << character1 << endl;
//////	cout << character2 << endl;
//////	cout << *(*(character3+2)+2) << endl;
//////	cout << character3[3] << endl;
//////
//////
//////	system("pause");
//////	return 0;
//////}
////
////
////#include <iostream>
////#include <string>
////
////using namespace std;
////int main()
////{
////	char ans = ' ';
////	string ques = "Are you 18 years old or above? ";
////
////	do
////	{
////		cout << ques;
////		cin >> ans;
////
////		if (ans == 'y' || ans == 'Y' || ans == 'N' || ans == 'n')
////		{
////			cout << "VALID LETTER" << endl;
////		}
////		else
////		{
////			cout << "INVALID LETTER" << endl;
////		}
////
////
////	} while (ans == 'y' || ans == 'Y' || ans == 'N' || ans == 'n');
////
////	system("pause");
////}
////
////#define _CRT_SECURE_NO_WARNINGS
////
//#include <iostream>
//#include <fstream>
//#include <string>
//#include <iomanip>
//
//using namespace std;
////public functions
//
//
//double function1(double &a);
//int main()
//{
//	double pa=0;
//	double w = 20;
//	double po[3] = { 1, 2, 3 };
//
//	cout << *(po+1) << endl;
//
//	pa = w / 3;
//	cout << "a = " << fixed << setprecision(4) << pa << endl;
//
//	double xos;
//	double da = 2.3;
//	xos = function1(da);
//	cout <<"x = " <<xos << endl;
//	cout << "da = " << da << endl;
//	
//	system("pause");
//	return 1;
//}
//
//double function1(double &a)
//{
//	double y;
//	a = 2 * a;
//
//	return a;
//}
//
//
