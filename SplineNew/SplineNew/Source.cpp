#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <errno.h>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;


// интерполяционный лагранжев сплайн 3 порядка
int count_of_elements; // количество конечных элементов 
vector<double> grid;
vector<double> weights; // массив весов
vector<double> q_vector; // массив q

vector<vector<double>> operator+(const vector<vector<double>>& matrix, const vector<vector<double>>& matrix2) {
	vector<vector<double>> c;
	int size = matrix.size();
	c.resize(size);
	for (auto& vect : c) {
		vect.resize(size);
	}
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			c[i][j] = matrix[i][j] + matrix2[i][j];
		}
	}
	return c;
}

vector<vector<double>> operator*(const double& a, const vector<vector<double>>& matrix) {
	vector<vector<double>> c;
	int size = matrix.size();
	c.resize(size);
	for (auto& vect : c) {
		vect.resize(size);
	}
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			c[i][j] = matrix[i][j] * a;
		}
	}
	return c;
}

vector<double> operator*(const vector<vector<double>>& matrix, const vector<double>& vect) {
	vector<double> c;
	int size = matrix.size();
	c.resize(size);

	for (int i = 0; i < size; ++i) {
		double sum = 0;
		for (int j = 0; j < size; ++j) {
			sum += matrix[i][j] * vect[j];
		}
		c[i] = sum;
	}
	return c;
}

//         
vector<double> operator%(const vector<vector<double>>& matrix, const vector<double>& vect)
{
	vector<double> c;
	int size = matrix.size();
	c.resize(size);

	for (int i = 0; i < size; ++i) {
		double sum = 0;
		for (int j = 0; j < size; ++j) {
			sum += matrix[j][i] * vect[j];
		}
		c[i] = sum;
	}
	return c;
}

vector<double> operator* (const double& c, const vector<double>& vect)
{
	vector<double> result = vect;
	for (int i = 0; i < result.size(); i++)
	{
		result[i] = vect[i] * c;
	}
	return result;
}

//      
double  operator* (const vector<double>& vector1, const vector<double>& vector2)
{
	double result = 0;
	for (int i = 0; i < vector1.size(); i++)
	{
		result += (vector1[i] * vector2[i]);
	}
	return result;
}

vector<double> operator+(const vector<double>& vect1, const vector<double>& vect2) {
	vector<double> c;
	int size = vect1.size();
	c.resize(size);

	for (int i = 0; i < size; ++i) {
		c[i] = vect1[i] + vect2[i];
	}
	return c;
}

vector<double> operator-(const vector<double>& vect1, const vector<double>& vect2) {
	vector<double> c;
	int size = vect1.size();
	c.resize(size);

	for (int i = 0; i < size; ++i) {
		c[i] = vect1[i] - vect2[i];
	}
	return c;
}

void Read_matrix(string matrix_filename, vector<vector<double>>& matrix_A, int Size)
{
	fstream file_matrix;
	file_matrix.open(matrix_filename);
	if (file_matrix.is_open())
	{
		for (int i = 0; i < Size; i++)
		{
			for (int j = 0; j < Size; j++)
			{
				file_matrix >> matrix_A[i][j];
			}
		}
	}
	file_matrix.close();
}
void Read_vector(string vector_filename, vector<double>& vector, int Size)
{
	fstream file_vector;
	file_vector.open(vector_filename);
	if (file_vector.is_open())
	{
		for (int i = 0; i < Size; i++)
		{
			file_vector >> vector[i];
		}
	}
	file_vector.close();
}

double ScalVectOnVect(vector<double> vector1, vector<double> vector2, int Size)
{
	int N = vector1.size();
	double sum = 0;
	for (int i = 0; i < N; i++)
	{
		sum += (vector1[i] * vector2[i]);
	}
	return sum;
}

vector<double> LOS(vector<vector<double>> matrix_A, vector<double> vector_B, int Size)  //    A
{
	double epsilon = 1e-13;
	vector<double> r0(Size);
	vector<double> r_next(Size);
	vector<double> z0(Size);
	vector<double> z_next(Size);
	vector<double> x_next(Size);
	vector<double> p0(Size);
	vector<double> p_next(Size);
	vector<double> x0(Size);
	vector<double> result = x0;

	double normF;
	double sum = 0;
	double alpha = 0;
	int iter_count = 0;

	int max_iter = 10000;


	double normRkkvad = 0;
	r0 = vector_B - matrix_A * result;
	z0 = r0;
	p0 = matrix_A * z0;

	//    
	FILE* out;
	out = fopen("LOSout.txt", "wt");
	fprintf(out, "Iter  x                           y                           z                               NormRk\n");
	fprintf(out, "%d;\t%2.18le;\t%2.18le;\t%2.18le;\t%2.18le;\n", iter_count, x0[0], x0[1], x0[2], normRkkvad);
	//    

	do
	{
		iter_count++;
		//alpha nach
		double alphaChisl = 0;
		double alphaZnam = 0;
		alphaChisl = ScalVectOnVect(p0, r0, Size);
		alphaZnam = ScalVectOnVect(p0, p0, Size);

		alpha = alphaChisl / alphaZnam;
		//alpha kon

		x_next = result + alpha * z0;

		//rk
		r_next = r0 - alpha * p0;

		normRkkvad = ScalVectOnVect(r_next, r_next, Size);


		//betta nach
		double bettaChisl = 0;
		double bettaZnam = 0;
		double betta = 0;
		vector<double> temp;
		temp = matrix_A * r_next;

		bettaZnam = ScalVectOnVect(p0, p0, Size);
		bettaChisl = ScalVectOnVect(p0, temp, Size);

		betta = -(bettaChisl / bettaZnam);
		//betta kon

		//zk
		z_next = r_next + betta * z0;


		p_next = matrix_A * r_next + betta * p0;


		result = x_next;
		r0 = r_next;
		z0 = z_next;
		p0 = p_next;
		fprintf(out, "%d;\t%2.18le;\t%2.18le;\t%2.18le;\t%2.18le;\n", iter_count, result[0], result[1], result[2], normRkkvad);
	} while (normRkkvad > epsilon && iter_count < max_iter);
	cout << "iter_count LOS: " << iter_count << "\n";
	fclose(out);

	return result;
}


int read_grid(vector<double>& grid)
{
	fstream file_grid;										// файл с сеткой
	file_grid.open("grid.txt");

	if (file_grid.is_open())
	{
		file_grid >> count_of_elements;

		grid.resize(count_of_elements * 3 + 1);
		weights.resize(count_of_elements * 3 + 1);

		for (int i = 0; i < count_of_elements * 3 + 1; i++)
		{
			file_grid >> grid[i];
			weights[i] = 1;
		}
	}
	else
	{
		file_grid.close();
		perror("ERROR: open file");
		return (1);
	}
	file_grid.close();
	return(0);
}

double Ksi(double x, int i, double h) // страшная буква 
{
	return (x - grid[i*3]) / h;
}

double f_func(double x)
{
	return x*x*x*x;
}

// 38 в методе по сплайнам
double psi_0(double ksi)
{	
	return (-9. / 2.) * (ksi - (1. / 3.)) * (ksi - (2. / 3.)) * (ksi - 1.);
}
double psi_1(double ksi)
{
	return (27. / 2.) * ksi * (ksi - (2. / 3.)) * (ksi - 1.);
}
double psi_2(double ksi)
{
	return (27. / 2.) * ksi * (ksi - (1. / 3.)) * (ksi - 1.);
}
double psi_3(double ksi)
{
	return (9. / 2.) * ksi * (ksi - (1. / 3.)) * (ksi - (2. / 3.));
}

void local_matrix_A(vector<vector<double>>& local_A, int node_number)
{
	// обнуляем локальную матрицу
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			local_A[i][j] = 0;
		}
	}


	/*double h = grid[3*node_number] - grid[node_number];*/
	double h = 3.0;
	int first_local = 3 * node_number;
	int last_local = 3 * node_number + 4;

	for (int j = first_local; j < last_local; j++) // проходимся по локальным узлам конечного элемента 3n...3n+4
	{
		double ksi = (grid[j] - grid[3*node_number]) / h;

		local_A[0][0] += weights[j] * psi_0(ksi) * psi_0(ksi); // первый элемент на бумаге посчитай
		local_A[0][1] += weights[j] * psi_0(ksi) * psi_1(ksi);
		local_A[0][2] += weights[j] * psi_0(ksi) * psi_2(ksi);
		local_A[0][3] += weights[j] * psi_0(ksi) * psi_3(ksi);

		local_A[1][0] += weights[j] * psi_1(ksi) * psi_0(ksi);
		local_A[1][1] += weights[j] * psi_1(ksi) * psi_1(ksi);
		local_A[1][2] += weights[j] * psi_1(ksi) * psi_2(ksi);
		local_A[1][3] += weights[j] * psi_1(ksi) * psi_3(ksi);

		local_A[2][0] += weights[j] * psi_2(ksi) * psi_0(ksi);
		local_A[2][1] += weights[j] * psi_2(ksi) * psi_1(ksi);
		local_A[2][2] += weights[j] * psi_2(ksi) * psi_2(ksi);
		local_A[2][3] += weights[j] * psi_2(ksi) * psi_3(ksi);

		local_A[3][0] += weights[j] * psi_3(ksi) * psi_0(ksi);
		local_A[3][1] += weights[j] * psi_3(ksi) * psi_1(ksi);
		local_A[3][2] += weights[j] * psi_3(ksi) * psi_2(ksi);
		local_A[3][3] += weights[j] * psi_3(ksi) * psi_3(ksi);
	}

	cout << "\n\nLocal_Matrix_A" << "\n\n";
	for (int i = 0; i < local_A.size(); i++)
	{
		for (int j = 0; j < local_A.size(); j++)
		{
			cout << local_A[i][j] << ' ';
		}
		cout << "\n";
	}

}

void local_vector_B(vector<double>& vector_B, double node_number) 
{
	//double h = grid[node_number + 3] - grid[node_number];
	double h = 3.0;

	// обнуление локального вектора В
	for (int i = 0; i < 4; i++)
	{
		vector_B[i] = 0;
	}

	int first_local = 3 * node_number;
	int last_local = 3 * node_number + 4;

	for (int i = first_local; i < last_local; i++)
	{
		double ksi = (grid[i] - grid[3*node_number]) / h;
		double f = f_func(grid[i]);
		
		vector_B[0] += weights[i] * psi_0(ksi) * f;
		vector_B[1] += weights[i] * psi_1(ksi) * f;
		vector_B[2] += weights[i] * psi_2(ksi) * f;
		vector_B[3] += weights[i] * psi_3(ksi) * f;
	}
	cout << "\n\nlocal_B\n\n";
	for (int i = 0; i < vector_B.size(); i++)
	{
		cout << vector_B[i] << '\n';
	}
	cout << "\n\n";
}

void global_matrix_A( vector<vector<double>> &local_matrix_A, int node_number) 
{


}

void global_vector_B(vector<double> &vector_B, int node_number) 
{


}





double Spline(double x, vector<double> & q_vector)
{
	double ksi;
	double h;

	double q0;
	double q1;
	double q2;
	double q3;

	double phi0;
	double phi1;
	double phi2;
	double phi3;

	double spline = 0;


	for (int i = 0; i < count_of_elements; i++)
	{
		//Получение значения сплайна в точке заключается в поиске нужного интервала 
		//и вычисления линейной комбинации локальных базисных функций//

		if (x >= grid[i * 3] && x <= grid[3 * i + 3]) // если дохожу до последнего  grid[3*2 +3 = 9] 
		{
			h = grid[i * 3 + 3] - grid[i * 3];
			ksi = Ksi(x, i, h);
			
			q0 = q_vector[3*i];
			q1 = q_vector[3*i + 1];
			q2 = q_vector[3*i + 2];
			q3 = q_vector[3*i + 3];


			phi0 = psi_0(ksi);
			phi1 = psi_1(ksi);
			phi2 = psi_2(ksi);
			phi3 = psi_3(ksi);

			spline = (phi0 * q0) + (phi1 * q1) + (phi2 * q2) + (phi3 * q3);
			return spline;
		}
	}
}


void main()
{

	vector<vector<double>> local_A(4);
	for (int i = 0; i < 4; i++)
	{
		local_A[i].resize(4);
	}
	vector<double> local_B(4);

	vector<vector<double>> global_A;
	vector<double> global_B;



	setlocale(LC_ALL, "Russian");
	read_grid(grid);
	// выделить память под глобальную матрицу и глобальный вектор
	global_B.resize((4 * count_of_elements) - 1);
	global_A.resize((4 * count_of_elements) - 1);
	for (int i = 0; i < global_A.size(); i++)
	{
		global_A[i].resize((4 * count_of_elements - 1));
	}



	// положить в глобальную матрицу и в глобальный веектор B
	for (int i = 0; i < count_of_elements; i++)
	{
		local_matrix_A(local_A, i);
		local_vector_B(local_B, i);
		// положить в глобальную матрицу и в глобальный веектор B

		// в глобальную матрицу А положили
		int n = i;

		global_A[3 * n][3 * n] += local_A[0][0];
		global_A[3 * n][3 * n + 1] += local_A[0][1];
		global_A[3 * n][3 * n + 2] += local_A[0][2];
		global_A[3 * n][3 * n + 3] += local_A[0][3];

		global_A[3 * n + 1][3 * n] += local_A[1][0];
		global_A[3 * n + 1][3 * n + 1] += local_A[1][1];
		global_A[3 * n + 1][3 * n + 2] += local_A[1][2];
		global_A[3 * n + 1][3 * n + 3] += local_A[1][3];

		global_A[3 * n + 2][3 * n] = local_A[2][0];
		global_A[3 * n + 2][3 * n + 1] += local_A[2][1];
		global_A[3 * n + 2][3 * n + 2] += local_A[2][2];
		global_A[3 * n + 2][3 * n + 3] += local_A[2][3];

		global_A[3 * n + 3][3 * n] += local_A[3][0];
		global_A[3 * n + 3][3 * n + 1] += local_A[3][1];
		global_A[3 * n + 3][3 * n + 2] += local_A[3][2];
		global_A[3 * n + 3][3 * n + 3] += local_A[3][3];

		// в глобальный вектор В положили
		global_B[3 * n] += local_B[0];
		global_B[3 * n + 1] += local_B[1];
		global_B[3 * n + 2] += local_B[2];
		global_B[3 * n + 3] += local_B[3];
	}

	// чисто отладка
	cout << "\n\nGlobal_Matrix_A\n\n";
	for (int i = 0; i < global_A.size(); i++)
	{
		for (int j = 0; j < global_A.size(); j++)
		{
			cout << global_A[i][j] << ' ';
		}
		cout << '\n';
	}

	cout << "\n\nGlobal_B\n\n";
	for (int i = 0; i < global_B.size(); i++)
	{
		cout << global_B[i] << '\n';
	}
	cout << "\n\n";
	// чисто отладка


	q_vector = LOS(global_A, global_B, 7);

	// по идее теперь у меня есть все чтобы вычислять значение сплайна в точке
	double Sp2 = Spline(4, q_vector);


	double x = 0;
	double step = 0.5;
	printf("x               f(x)	               Spline(x)	            |f(x) - spline(x)|\n");
	for (int i = 0; i <= grid[count_of_elements * 3] / step; i++)
	{

		printf("%2.4le  %2.18le  %2.18le  %2.18le\n", x, f_func(x), Spline(x, q_vector), abs(f_func(x) - Spline(x, q_vector)));
		x += step;
	}

}

//Стр 212