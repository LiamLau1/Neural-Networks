#include<iostream>
#include<cmath>
#include<cstdlib>

using namespace std;

int main(){
	//RANDOM training sets for XOR-- two inputs and one output

	cout << "topology: 2 4 1" <<endl;
	for(int i =8000; i >= 0; --i){
		//2000 passes
		int n1 = (int)(2.0 * rand()/ double(RAND_MAX));
		int n2 = (int)(2.0 * rand()/ double(RAND_MAX));
		//int n3 = (int)(2.0 * rand()/ double(RAND_MAX));
		int t = (n1 ^ n2) /*^ n3*/; //should be 0 or 1 XOR or AND
		cout << "in: " << n1 << ".0 " << n2 << ".0 " /*<< n3 << ".0 "*/ <<endl;
		cout << "out: " << t << ".0" <<endl;
	}
	return 0;
}
