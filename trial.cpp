#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include<sstream>

using namespace std;

class WeightsData
{
public:
    WeightsData(const string filename);
    bool isEof(void) { return m_WeightsDataFile.eof(); }
    void getTopology(vector<unsigned> &topology);

    // Returns the number of input values read from the file:
    unsigned getNextWeights(vector<double> &inputVals);
    //unsigned getTargetOutputs(vector<double> &targetOutputVals);

private:
    ifstream m_WeightsDataFile;
};

void WeightsData::getTopology(vector<unsigned> &topology)
{
    string line;
    string label;

    getline(m_WeightsDataFile, line);
    stringstream ss(line);
    ss >> label;
	unsigned n;
	ss >> n;
	topology.push_back(n);
    /*while (!ss.eof()) {
        unsigned n;
        ss >> n;
        topology.push_back(n);
    }*/

    return;
}


WeightsData::WeightsData(const string filename)
{
    m_WeightsDataFile.open(filename.c_str());
}

unsigned WeightsData::getNextWeights(vector<double> &weightVals)
{
    weightVals.clear();

    string line;
    getline(m_WeightsDataFile, line);
    stringstream ss(line);

    string label;
    ss>> label;
    if (label.compare("in:") == 0) {
        double oneValue;
        while (ss >> oneValue) {
            weightVals.push_back(oneValue);
        }
    }

    return weightVals.size();
}

int main(){
	fstream& GotoLine(fstream& file, unsigned int num);
	WeightsData WData("/tmp/WeightsData.txt");
	vector<unsigned> topology;
	WData.getTopology(topology);
	if(topology[0] != 2){
		cout << topology[1];
	}
	/*string line9;
	string weight;
	GotoLine(m_WeightsDataFile,9);
	getline(m_WeightsDataFile, line9);
	stringstream ss(line9);
	ss >> weight;
	cout << weight <<endl;
	*/
	return 0;
}


fstream& GotoLine(fstream& file, unsigned int num)
{
    file.seekg(ios::beg);
    for(unsigned int i=0; i < num - 1; ++i)
    {
        file.ignore(numeric_limits<streamsize>::max(),'\n');
    }
    return file;
}
