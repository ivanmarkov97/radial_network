#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include "kmeans.h"

using namespace std;

ofstream colFile, resFile;
ofstream func, rosenFile;
ofstream colPlotFile;
vector<ofstream> colocls;

vector<Cluster> clusters;

double fRand(double min, double max){
	double f = (double) rand() / RAND_MAX;
	return min + f * (max - min);
}

bool sort_pair_func(pair<double, int> p1, pair<double, int> p2){
	return p1.first < p2.first;
}

double rosenFunction(double x1, double x2){
	return sin(x1) + cos(x2);
	//return (1 - x1)*(1 - x1) + 0.001 * pow(x2 - x1*x1, 2.0);
}

void multi(double **Mat1, int N1, int M1, double **Mat2, int N2, int M2, double **Res){
	if (M1 != N2){
		cout << "Cant multi matrix" << endl;
		//return NULL;
	}

	//double **Res = new double*[N1];
	//for (int i = 0; i < N1; i++)
	//	Res[i] = new double[M2]();

	cout << "M1 " << endl;
	for(int i = 0; i < N1; i++){
		for (int j = 0; j < M1; j++){
			cout << Mat1[i][j] << " ";
		}
		cout << endl;
	}

	cout << "M1 " << endl;
	for(int i = 0; i < N2; i++){
		for (int j = 0; j < M2; j++){
			cout << Mat2[i][j] << " ";
		}
		cout << endl;
	}

	cout << endl;
	for (int i = 0; i < N1; i++){
		for (int j = 0; j < M2; j++){
			cout << "Res [" << i << "][" << j << "] = ";
			Res[i][j] = 0.0;
			for (int k = 0; k < M1; k++){
				Res[i][j] += Mat1[i][k] * Mat2[k][j];
				cout << Mat1[i][k] << " * " << Mat2[k][j];
				cout << " + ";
			}
			cout << endl;
		}
	}
	//return Res;
}


void transpose(double **Mat, int N, int M, double **Res){
	//cout << N << "x" << M << endl;
	//double **Res = new double*[M];
	//for (int i = 0; i < M; i++)
	//	Res[i] = new double[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
			Res[j][i] = Mat[i][j];
	//return Res;
}


int GetMinor(double **src, double **dest, int row, int col, int order){
    // indicate which col and row is being copied to dest
    int colCount=0,rowCount=0;

    for(int i = 0; i < order; i++ ){
    	if( i != row ){
    		colCount = 0;
    		for(int j = 0; j < order; j++){
    			if( j != col){
    				dest[rowCount][colCount] = src[i][j];
    				colCount++;
    			}
    		}
    		rowCount++;
    	}
    }
    return 1;
}


double CalcDeterminant(double **A, int N){
	int i, j, k, r;
	double c, M, max, s, det=1;
	double **a;

	a = new double *[N];
	for(int i = 0; i < N; i++)
		a[i] = new double [N]();

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			a[i][j] = A[i][j];

	for (k = 0; k < N; k++){
		max = fabs(a[k][k]); r = k;
		for (i = k + 1; i < N; i++){
			if (fabs(a[i][k]) > max){ max = fabs(a[i][k]); r = i; }
			if (r != k) det = -det;
			for(j = 0; j < N; j++){
				c = a[k][j];
				a[k][j] = a[r][j];
				a[r][j] = c;
			}
			for(i = k + 1; i < N; i++)
				for(M = a[i][k] / a[k][k], j = k; j<N; j++)
					a[i][j] -= M * a[k][j];
		}
	}

	for (i = 0; i < N; i++)
		det *= a[i][i];

	for (i = 0; i < N; i++)
		delete[] a[i];
	delete [] a;

	return det;

}


void MatrixInversion(double **A, int order, double **Y){
    double det = 1.0/CalcDeterminant(A,order);

    double *temp = new double[(order-1)*(order-1)];
    double **minor = new double*[order-1];
    for(int i = 0;i < order-1; i++)
    	minor[i] = temp + (i * (order-1));

    for(int j=0; j<order; j++){
    	for(int i=0; i<order; i++){
    	// get the co-factor (matrix) of A(j,i)
    	GetMinor(A,minor,j,i,order);
    	Y[i][j] = det*CalcDeterminant(minor,order-1);
    	if( (i+j)%2 == 1)
    	Y[i][j] = -Y[i][j];
    	}
    }
    delete [] temp;
    delete [] minor;
}


void invert(double **Mat, int N, double **Res){
	//double **Res = new double*[N];
	//for (int i = 0; i < N; i++)
	//	Res[i] = new double[N];

	MatrixInversion(Mat, N, Res);
	//return Res;
}


class SurfacePoint{
	private:
		double x1, x2, value;
	public:
		SurfacePoint(double _x1, double _x2) { x1 = _x1, x2 = _x2, value = rosenFunction(_x1, _x2); };
		double getValue(){ return rosenFunction(x1, x2); };
		double getX1(){return x1;};
		double getX2(){return x2;};
};

vector<SurfacePoint> surface;


class RadialNeuron{
	private:
		int id;
	public:
		RadialNeuron( double c1, double c2, double s1, double s2, int _id): C1(c1), C2(c2), S1(s1), S2(s2), id(_id) {}; 

		double C1, C2, S1, S2;

		int getId(){return id;};

		double radFunction(SurfacePoint point){
			return exp( - (pow(point.getX1() - C1, 2.0) / (2.0 * S1 * S1) + 
				 		   pow(point.getX2() - C2, 2.0) / (2.0 * S2 * S2) 
				 		   )
					  );
		};
};


class RadialNetwork{
	private:
		int neuronNum;
		double **outputDense;
		double *W_local;
		vector<RadialNeuron> neurons;
	public:
		RadialNetwork(int _neuronNum) {
			for (int i = 0; i < _neuronNum; i++)
				neurons.push_back(RadialNeuron(fRand(-4, 4), fRand(-4, 4), 1.0, 1.0, i));
			neuronNum = neurons.size();
			W_local = new double[neuronNum]();
		}

		RadialNetwork(vector<Cluster> clusters) {
			for (int i = 0; i < clusters.size(); i++){
				double x1 = clusters[i].getCentralValue(0);
				double x2 = clusters[i].getCentralValue(1);
				neurons.push_back(RadialNeuron(x1, x2, 1.0, 1.0, i));
			}

			neuronNum = neurons.size();
			W_local = new double[neuronNum]();
		}

		RadialNeuron getNeuronById(int id){
			for (int i = 0; i < neurons.size(); i++){
				if (neurons[i].getId() == id){
					return neurons[i];
				}
			}
			return neurons[0];
		}

		double getOutput(SurfacePoint point){
			double sum = 0.0;
			for (int i = 0; i < neurons.size(); i++){
				sum += W_local[i] * neurons[i].radFunction(point);
			}
			return sum;
		};

		double getError(vector<SurfacePoint> surface){
			double sum = 0.0;
			for (int k = 0; k < surface.size(); k++)
				sum += pow(getOutput(surface[k]) - surface[k].getValue(), 2.0);
			return sum * 0.5;
		};

		double *getW() {return W_local;};
		vector<RadialNeuron> getNeurons() {return neurons;};

		void setW(double **W){
			//cout << "neuron num " << neuronNum << endl;
			for (int i = 0; i < neuronNum; i++)
				W_local[i] = W[i][0];
		}

		RadialNeuron getNearestNeuron(SurfacePoint point){
			double dist = sqrt( pow(neurons[0].C1 - point.getX1(), 2.0) + pow(neurons[0].C2 - point.getX2(), 2.0) );
			int neuronIndex = 0;

			for (int i = 1; i < neurons.size(); i++){
				double d = sqrt( pow(neurons[i].C1 - point.getX1(), 2.0) + pow(neurons[i].C2 - point.getX2(), 2.0) );
				if (d < dist) {
					dist = d;
					neuronIndex = i;
				}
			}

			//printf("min dist %lf\n", dist);
			return neurons[neuronIndex];
		}

		
		double neuronDistance(RadialNeuron n1, RadialNeuron n2){
			return sqrt( pow(n1.C1 - n2.C1, 2.0) + pow(n1.C2 - n2.C2, 2.0) );
		}

		double computeGauseDistance(RadialNeuron neuron, int neighborNum=1){
			vector<pair<double, int>> dist_indx_pair;
			for (int i = 0; i < neurons.size(); i++){
				if (neuron.getId() != neurons[i].getId()){
					double d = neuronDistance(neuron, neurons[i]);
					dist_indx_pair.push_back(make_pair(d, i));
				}
			}
			for (int i = 0; i < dist_indx_pair.size(); i++){
				//cout << dist_indx_pair[i].first << " " << dist_indx_pair[i].second << endl;
			}
			sort(dist_indx_pair.begin(), dist_indx_pair.end(), sort_pair_func);
			//cout << "AFTER SORT"<<endl;
			for (int i = 0; i < dist_indx_pair.size(); i++){
			//	cout << dist_indx_pair[i].first << " " << dist_indx_pair[i].second << endl;
			}
			double avDist = 0;
			for (int i = 0; i < neighborNum; i++){
				double ad = neuronDistance(getNeuronById(dist_indx_pair[i].second), neuron);
				RadialNeuron n = getNeuronById(dist_indx_pair[i].second);
				//cout << "DIST FROM (" << n.C1 << " " << n.C2 << ") ";
				//cout << "TO (" << neuron.C1 << " " << neuron.C2 << ") ";
				//cout << "IS " << ad << endl;
				avDist += ad;
			}
			avDist /= neighborNum;
			//cout << "AVDIST " << avDist << endl;
			return avDist;
		}

		void setNeuron(RadialNeuron neuron){
			for (int i = 0; i < neurons.size(); i++){
				if (neurons[i].getId() == neuron.getId()){
					neurons[i] = neuron;
					break;
				}
			}
		}

		void fit(vector<SurfacePoint> surface) {
			double nu = 0.2;
			int t = 0;
			int c = 0;
			double err1 = getError(surface);
			double err2 = getError(surface);
			cout << "error " << err1 << endl;

			//while(c != -1){
			for(int i = 0; i < surface.size(); i++){
				//cout << "ALL NEURONS" << endl;
				for (int i = 0; i < neurons.size(); i++){
				//	cout << neurons[i].getId() << " (" << neurons[i].C1 << " " << neurons[i].C2 << ")" << endl;
				}
				RadialNeuron neuron = getNearestNeuron(surface[i]);
				//cout << "NEAR NEUR 2 DOT " << surface[i].getX1() << ":" << surface[i].getX2() << " is " << neuron.getId();
				//cout << " (" << neuron.C1 << " " << neuron.C2 << ")" << endl;
				neuron.C1 = neuron.C1 + nu * (surface[i].getX1() - neuron.C1);
				neuron.C2 = neuron.C2 + nu * (surface[i].getX2() - neuron.C2);
				//cout << "AFTER MOVE " << " (" << neuron.C1 << " " << neuron.C2 << ")" << endl;

				setNeuron(neuron);
				//cout << "TOTAL" << endl;
				for (int i = 0; i < neurons.size(); i++){
				//	cout << neurons[i].getId() << " (" << neurons[i].C1 << " " << neurons[i].C2 << ")" << endl;
				}


				for (int i = 0; i < clusters.size(); i++){
					double x1 = clusters[i].getCentralValue(0);
					double x2 = clusters[i].getCentralValue(1);
					for (int j = 0; j < clusters[i].getTotalPoints(); j++){
						rosenFile << clusters[i].getPoint(j).getValue(0) << " ";
						rosenFile << clusters[i].getPoint(j).getValue(1) << " ";
						rosenFile << clusters[i].getID() + 1 << endl;
					}
					//rosenFile << x1 << " " << x2 << " "  << 0 << endl;
					//cout << x1 << " " << x2 << endl;
				}
				for (int i = 0; i < neurons.size(); i++){
					rosenFile << neurons[i].C1 << " " << neurons[i].C2 << " " << 0 << endl;
				}
				rosenFile << "\n";
			}
			

			//cout << endl;
			//cout << "PAIR" << endl;
			for (int i = 0; i < neurons.size(); i++){
				double si = computeGauseDistance(neurons[i], 3);
				neurons[i].S1 = si / 1;
				neurons[i].S2 = si / 1;
			}


			int GN = surface.size(), GM = neurons.size();

			//cout << "G" << endl;
			double **G = new double*[GN];
			for (int i = 0; i < GN; i++)
				G[i] = new double[GM];

			for (int i = 0; i < GN; i++){
				for (int j = 0; j < GM; j++){
					//cout << "NEUR" << neurons[j].getId() << " (" << neurons[j].C1 << " " << neurons[j].C2 << ")" << endl;
					//cout << "Point " << surface[i].getX1() << " " << surface[i].getX2() << endl;
					G[i][j] = neurons[j].radFunction(surface[i]);
					//cin >> c;
				}
			}

			//cout << "D" << endl;
			double **D = new double*[GN];
			for (int i = 0; i < GN; i++)
				D[i] = new double[1];
			for (int i = 0; i < GN; i++){
				//cout << "D [" << i << "] " << surface[i].getValue() << endl;
				D[i][0] = surface[i].getValue();
				//cin >> c;
			}

			// GW = D
			
			cout << "G" << endl;
			for(int i = 0; i < GN; i++){
				for(int j = 0; j < GM; j++){
					cout << G[i][j] << " ";
				}
				cout << endl;
			}

			//cout << "D" << endl;
			/*for(int i = 0; i < GN; i++){
				for(int j = 0; j < 1; j++){
					cout << D[i][j] << " ";
				}
				cout << endl;
			}*/

			cout << "GT" << endl;
			double **GT = new double*[GM];
			for (int i = 0; i < GM; i++)
				GT[i] = new double[GN];
			transpose(G, GN, GM, GT);
			for(int i = 0; i < GM; i++){
				for(int j = 0; j < GN; j++){
					cout << GT[i][j] << " ";
				}
				cout << endl;
			}


			cout << "GTG" << endl;
			double **GTG = new double*[GM];
			for (int i = 0; i < GM; i++)
				GTG[i] = new double[GM];
			multi(GT, GM, GN, G, GN, GM, GTG);
			for(int i = 0; i < GM; i++){
				for(int j = 0; j < GM; j++){
					cout << GTG[i][j] << " ";
				}
				cout << endl;
			}
			
			cout << "GTG_INV" << endl;
			double **GTG_INV = new double*[GM];
			for (int i = 0; i < GM; i++)
				GTG_INV[i] = new double[GM];
			invert(GTG, GM, GTG_INV);
			for(int i = 0; i < GM; i++){
				for(int j = 0; j < GM; j++){
					cout << GTG_INV[i][j] << " ";
				}
				cout << endl;
			}
			
			cout << "G_plus" << endl;
			double **G_plus = new double*[GM];
			for (int i = 0; i < GM; i++)
				G_plus[i] = new double[GN];
			multi(GTG_INV, GM, GM, GT, GM, GN, G_plus);
			for(int i = 0; i < GM; i++){
				for(int j = 0; j < GN; j++){
					cout << G_plus[i][j] << " ";
				}
				cout << endl;
			}


			double **W = new double*[GM];
			for (int i = 0; i < GM; i++)
				W[i] = new double[1];
			multi(G_plus, GM, GN, D, GN, 1, W);
			cout << "W" << endl;
			for(int i = 0; i < GM; i++){
				for(int j = 0; j < 1; j++){
					cout << W[i][j] << " ";
				}
				cout << endl;
			}


			cout << "SET W" << endl;
			setW(W);

			err2 = getError(surface);
			cout << "new error " << err2 << endl;
			cout << "old error " << err1 << endl;
			cout << getError(surface) << endl;
			t++;
			cin >> c;
			
			cout << "write to file" << endl;
			for (int k = 0; k < neurons.size(); k++){
				for (double i = -6; i < 3; i += 0.5){
					for(double j = -3; j < 3; j += 0.5){
						colocls[k] << i << " " << j << " " << W_local[k] * neurons[k].radFunction(SurfacePoint{i, j})  << 0 << endl;
					}
				}
			}

			for (double i = -6; i < 3; i += 0.5){
				for(double j = -3; j < 3; j += 0.5){
					func << i << " " << j << " " << rosenFunction(i, j) << endl;
 					//resFile << i << " " << j << " " << getOutput(SurfacePoint{i, j})  << endl;
				}
			}
			func << "\n";
			for (double i = -6; i < 3; i += 0.5){
				for(double j = -3; j < 3; j += 0.5){
					//func << i << " " << j << " " << rosenFunction(i, j) << endl;
 					resFile << i << " " << j << " " << getOutput(SurfacePoint{i, j})  << endl;
				}
			}

			//cout << "W" <<endl;
			for (int i = 0; i < neurons.size(); i++){
				cout << W_local[i] << endl;
			}

			colPlotFile << "set hidden3d\nset dgrid 20,20 qnorm 4\nsplot ";
			for (int i = 0; i < colocls.size(); i++){
				colPlotFile << "'colocl" + to_string(i) + ".txt' with lines notitle, ";
			}
			colPlotFile << "\npause -1";
		}
};


int main(){

	srand(time(NULL));

	colPlotFile.open("plot_colocols.gp");
	rosenFile.open("rosen_func.txt");
	//colFile.open("colocl.txt");
	resFile.open("res.txt");
	func.open("f.txt");

	int L = 10;

	vector<Point> points;
	int id = 0;
	for (double x1 = -6.; x1 < 3; x1 += 1){
		for (double x2 = -6; x2 < 3; x2 += 1){
			vector<double> values;
			values.push_back(x1);
			values.push_back(x2);
			points.push_back(Point(id, values));
			id++;
		}
	}

	random_shuffle(points.begin(), points.end());

	for (int i = 0; i < points.size(); i++){
		double x1 = points[i].getValue(0);
		double x2 = points[i].getValue(1);
		surface.push_back(SurfacePoint(x1, x2));
	}

	cout <<"TOTAL POINTS" << endl;
	cout << surface.size() << endl;

	cin >> L;

	KMeans kmeans(L, surface.size(), 2, 10000);
	kmeans.run(points);

	for (int i = 0; i < L; i++){
		string fileName = "colocl" + to_string(i) + ".txt";
		cout << "File name " << fileName << endl;
		ofstream file;
		file.open(fileName);
		colocls.push_back(move(file));
	}

	clusters = kmeans.getClusters();

	for (int i = 0; i < clusters.size(); i++){
		double x1 = clusters[i].getCentralValue(0);
		double x2 = clusters[i].getCentralValue(1);
		for (int j = 0; j < clusters[i].getTotalPoints(); j++){
			rosenFile << clusters[i].getPoint(j).getValue(0) << " ";
			rosenFile << clusters[i].getPoint(j).getValue(1) << " ";
			rosenFile << clusters[i].getID() + 1 << endl;
		}
		rosenFile << x1 << " " << x2 << " "  << 0 << endl;
		cout << x1 << " " << x2 << endl;
	}
	rosenFile << "\n";

	RadialNetwork radialNetwork(clusters);
	radialNetwork.fit(surface);


	/*cout << "AFTER FIT"<<endl;
	auto neurons = radialNetwork.getNeurons();
	for (int i = 0; i < neurons.size(); i++){
		cout << neurons[i].C1 << " " << neurons[i].C2 << endl;
	}*/

	return 0;
}
