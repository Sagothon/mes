#include <iostream>;
using namespace std;

class global_data {   
	friend class FEM_GRID; 
	friend class element;
	friend class SOE;

int ne;    //iloœæ elementów
int	nh;     //iloœæ wi¹zañ
double L;   //d³ugoœæ
double k;   //wspó³czynnik przewodzenia ciep³a
double a;   //wspólczynnik wymiany konwekcyjnej
int tn;  //temperatura otoczenia
int q;      //gêstoœæ strumienia ciep³a

int kon;   //konwekcja
int zet;

public:
	global_data() {

		ne = 10;
		nh = ne+1;
		L = 5;
		a = 10;
		tn = 400;
		q = -150;

		kon = 0;
		zet = 10;
	
	}

	void read_data() {
 
	cout<<"podaj liczbe elementów"<<endl;
	cin>>ne;
	nh = ne+1;
	
	cout<<"podaj dlugosc"<<endl;
	cin>>L;

	cout<<"podaj temperature otoczenia"<<endl;
	cin>>tn;

	cout<<"podaj na którym wêŸle jest konwekcja"<<endl;
	cin>>kon;

	if (kon==0)
		zet=ne;
	else 
		kon=ne;

}
};

class element {      
public:
	int id1, id2;
	double dL;	//d³ugoœæ elementu
	double C;	//wspó³czynnik C
	double S;	//pole przekroju
	double k;	//wspolczynnik przewodzenia ciepla
	double H[2][2]; //macierz H
	double P[2];	//macierz P

	element() {};
	element(global_data* g, int id1, int id2, double S, double k) {

		this->id1 = id1;
		this->id2 = id2;
		this->k = k;
		this->S = S;
		dL = g->L/g->ne;
		C = (S*k)/dL;

		H[0][0] = C; 
		H[0][1] = -C;  
		H[1][0] = -C;  
		H[1][1] = C; 

		P[0] = 0;
		P[1] = 0;
	};
};

class node {           
public:
	int BC;
	int x;
	int t;

	node() {
		BC = 0;
	};
};

class FEM_GRID {
	friend class SOE;
	element* table_of_elements;
	node* table_of_nodes;

public:
	FEM_GRID(global_data* g) {

		table_of_elements = new element[g->ne];					//tworzenie odpowiedniej iloœci elementów

		for (int i=0; i<g->ne; i++) {
			double k,S;
			
		cout<<"element nr "<<i<<"\n podaj powierzchnie: ";        
		cin>>S;
        cout<<endl;      

        cout<<"podaj wspolczynnik przewodzenia: ";
		cin>>k;
        cout<<endl;  

		table_of_elements[i] = element(g, i, i+1 ,S, k);
		
	}


		table_of_nodes = new node[g->nh];				//tworzenie odpowiedniej iloœci wêz³ów

		table_of_nodes[g->kon].BC = 1;
		table_of_nodes[g->zet].BC = 2;
		
		if (g->kon==0)
		{
			table_of_elements[g->kon].P[0] = g->q* table_of_elements[g->kon].S;

			table_of_elements[g->zet-1].H[1][1] += table_of_elements[g->zet-2].S*g->a;					//wprowadzenie odpowiednich zmian w tablicach tam gdzie warunki brzegowe
			table_of_elements[g->zet-1].P[1] = -(g->a*g->tn*table_of_elements[g->zet-1].S);
		}
		
		else if (g->zet==0)
		{
			table_of_elements[g->kon-1].P[1] = g->q* table_of_elements[g->kon-1].S;

			table_of_elements[g->zet].H[0][0] += table_of_elements[g->zet].S*g->a;    
			table_of_elements[g->zet].P[0] = -(g->a*g->tn*table_of_elements[g->zet].S);
		}
		

		for ( int j=0;j<g->ne;j++) { 
			cout<<"lokalna macierz H elementu "<<j<<endl;
		for (int x=0; x<2; x++)
        {
            cout<<endl;
            for (int y=0; y<2; y++)                    //wypisywanie macierzy H
            {
				cout<<table_of_elements[j].H[x][y]<<" ";
            }
			
        }
		cout<<endl<<endl<<"lokalny wektor P elementu "<<j<<endl<<"[ ";
		for (int i=0;i<2;i++)													//wypisywanie wektorów P
			{
				cout<<table_of_elements[j].P[i]<<" ";
			}
		cout<<"]"<<endl<<endl;
	}
	}
	
};

class SOE {
public:
	double** global_H;
	double* global_P;
	int* t;

	SOE(global_data* g, FEM_GRID* f) {
		
		build_GH_GP (g,f);
		solve (global_H, global_P, g->nh);

	};

void build_GH_GP (global_data* g, FEM_GRID* f) {

	global_H = new double* [g->nh];    //utworzenie macierzy globalnej
	for (int i=0; i<g->nh; i++)
	global_H[i]=new double[g->nh];

	for (int i=0; i<g->nh; i++)
       for(int j=0; j<g->nh; j++)
			global_H[i][j]=0;                //wyzerowanie macierzy globalnej


	global_P = new double [g->nh];         //utworzenie wektora globalnego P
	for (int i=0; i<g->nh; i++)
		global_P[i]=0;                                 // zerowanie wektora globalnego P;

	for (int i=0;i<g->ne;i++) {

		global_P[f->table_of_elements[i].id1] += f->table_of_elements[i].P[0];               //uk³adanie wektorów P w globalny
		global_P[f->table_of_elements[i].id2] += f->table_of_elements[i].P[1];
	}

	for (int i=0; i<g->ne; i++)
	{
		global_H[f->table_of_elements[i].id1][f->table_of_elements[i].id1]+= f->table_of_elements[i].H[0][0];
		global_H[f->table_of_elements[i].id1][f->table_of_elements[i].id2]+= f->table_of_elements[i].H[0][1];          //uk³adanie macierzy H w macierz globaln¹
		global_H[f->table_of_elements[i].id2][f->table_of_elements[i].id1]+= f->table_of_elements[i].H[1][0];
		global_H[f->table_of_elements[i].id2][f->table_of_elements[i].id2]+= f->table_of_elements[i].H[1][1];
	}

	        
	cout<<endl<<endl<<"Macierz globalna H:"<<endl;

	for (int i=0; i<g->nh; i++)
		{
			cout<<endl;
			for(int j=0; j<g->nh; j++)                //wypisanie tablicy globalnej
				cout<<global_H[i][j]<<"  ";
		}

	cout<<endl;
	cout<<endl<<"Wektor globalny P: [ ";

	for (int i=0; i<g->nh; i++)          //wypisanie globalnego wektora P
		{
			cout<<global_P[i]<<" ";
		}
		cout<<"]"<<endl<<endl;
	}

void solve (double** H, double* P, int nh) {

	for (int i = 0; i < nh; i++)                                      //rozwiazanie ukladu rownan
            {
				double x = H[i][i];
				P[i] /= x;
                for (int j = 0; j < nh; j++)                            
                {
					H[i][j] /= x;
                }
                for (int j = 0; j < nh; j++)
                {
                    if (j != i)                                         
                    {
						x = H[j][i];
                        for (int k = 0; k < nh; k++)
                        {
							H[j][k] -= (H[i][k] * x);
                        }
						P[j] -= (P[i] * x);
                    }
                }
            }
	}

void print(global_data* g, double* P) {
	cout<<"Wyniki to: ";
for (int i=0; i<g->nh; i++)                    //wypisywanie wyników
{
cout<<"t"<<(i+1)<<": "<<-P[i]<<"   ";
}
};
};


int main() {
 
	cout<<"Hello Allice! \n";

	global_data* g = new global_data();
	//g->read_data();
	FEM_GRID* f = new FEM_GRID(g);
	SOE* s = new SOE (g,f);
	s->print(g,s->global_P);

	system("PAUSE");
	return 0;
}

	
