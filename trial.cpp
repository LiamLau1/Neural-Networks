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
	fstream& GotoLine(fstream& file, unsigned int num);
    double getNextWeights(vector<double> &weightVals);
    //unsigned getTargetOutputs(vector<double> &targetOutputVals);

private:
    fstream m_WeightsDataFile;
};

void WeightsData::getTopology(vector<unsigned> &topology)
{
    string line;

    getline(m_WeightsDataFile, line);
    stringstream ss(line);
    while (!ss.eof()) {
        unsigned n;
        ss >> n;
        topology.push_back(n);
		//fill the topology vector with the number of neurons per layer written in weights.txt
    }

    return;
}


WeightsData::WeightsData(const string filename)
{
    m_WeightsDataFile.open(filename.c_str());
}


double WeightsData::getNextWeights(vector<double> &weightVals)
{
	//countLines(m_WeightsDataFile);
	weightVals.clear();
    string line;
	GotoLine(m_WeightsDataFile, 2);
    while (!m_WeightsDataFile.eof()) {
        getline(m_WeightsDataFile,line);
        double oneValue;
        stringstream ss(line);
        ss >> oneValue;
        weightVals.push_back(oneValue);
    }
    weightVals.pop_back(); //remove the repeated element


    return weightVals.size();
}

fstream& WeightsData::GotoLine(fstream& file, unsigned int num)
{
    file.seekg(ios::beg);
    for(unsigned int i=0; i < num - 1; ++i)
    {
        file.ignore(numeric_limits<streamsize>::max(),'\n');
    }
    return file;
}

int main(){
	WeightsData WData("/tmp/weights.txt");
	vector<unsigned> topology;
	vector<double> weights;
	WData.getTopology(topology);
	WData.getNextWeights(weights);
	return 0;
}
