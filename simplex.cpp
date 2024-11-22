#include <iostream>
#include <array>

using namespace std;

#define RESET   "\033[0m"
#define GREEN  "\x1B[32m"

void print_simplex_tableau(int n, int m, double (&a)[3][5], double *b, double *c, int *xb, double *z, double *c_z, double *theta, double fmax, int in = -1, int out = -1)
{
	bool is_result = (in == -1 && out == -1);
	////////////////////////////////////////
	cout << "\t\tc_j\t";
	for (int j = 0; j < n; j++)
		cout << c[j] << '\t';
	cout << "\t\n";
	////////////////////////////////////////
	cout << "i\tc_b\tx_b\t";
	for (int j = 0; j < n; j++)
		cout << "x_" << j + 1 << '\t';
	cout << "b_i\ttheta_i\n";
	////////////////////////////////////////
	for (int i = 0; i < m; i++)
	{
		cout << i + 1 << '\t' << c[xb[i]] << "\tx_" << xb[i] << '\t';
		for (int j = 0; j < n; j++)
			cout << a[i][j] << '\t';
		cout << (is_result ? GREEN : RESET) << b[i] << RESET << '\t';
		cout << (!is_result ? to_string(theta[i]) : "");
		cout << ((i == out) ? "\t-> out" : "") << '\n';
	}
	////////////////////////////////////////
	cout << "\t\tz_j\t";
	for (int j = 0; j < n; j++)
		cout << z[j] << '\t';
	cout << (is_result ? GREEN : RESET) << fmax << RESET << '\n';
	////////////////////////////////////////
	cout << "\t\tc_j-z_j\t";
	for (int j = 0; j < n; j++)
		cout << c[j] - z[j] << '\t';
	cout << '\n';
	////////////////////////////////////////
	if (!is_result)
	{
		cout << "\t\t\t";
		for (int j = 0; j < in; j++)
			cout << '\t';
		cout << "^ in\n";
	}
	cout << "--------------------------------------------------------------------------------\n"; 
}

int main()
{
	/// Initialize all the inputs
	int n = 5, m = 3;

	double a[3][5] = {
		{5, 2, 1, 0, 0},
		{2, 3, 0, 1, 0},
		{4, 2, 0, 0, 1}},
	b[m] = {150,
			100,
			80},
	c[n] = {12, 8, 0, 0, 0};
	int xb[m] = {2,
				 3,
				 4};

	/// Declare all the outputs
	double z[n + 1], c_z[n], theta[m];
	double fmax;

	/// Check if any b is negative
	for (int i = 0; i < m; i++)
	{
		if (b[i] < 0)
		{
			for (int j = 0; j < n; j++)
				a[i][j] = -a[i][j];
			b[i] = -b[i];
		}
	}

	while (true)
	{

	/// Compute the z_j and c_j - z_j rows
		for (int j = 0; j < n; j++)
		{
			z[j] = 0;
			for (int i = 0; i < m; i++)
				z[j] += c[xb[i]] * a[i][j],
					c_z[j] = c[j] - z[j];
		}

	/// Determine the next incoming variable
		int in = -1; // < n

		for (int j = 0; j < n; j++)
			if (c_z[j] > 0 && (in == -1 || (c_z[j] > c_z[in])))
				in = j;

		if (in == -1)
			break;

	/// Determine the next outgoing variable
		int out = -1; // < m

		for (int i = 0; i < m; i++)
			theta[i] = b[i] / a[i][in];

		for (int i = 0; i < m; i++)
			if (theta[i] > 0 && (out == -1 || theta[i] < theta[out]))
				out = i;

		if (out == -1)
			break;

	/// Print the currently computed tableau
		fmax = 0;
		for (int i = 0; i < m; i++)
			fmax += c[xb[i]] * b[i];
		
		print_simplex_tableau(n, m, a, b, c, xb, z, c_z, theta, fmax, in, out);

	/// Scale the outgoing variable row
		double pivot_element = a[out][in];

		for (int j = 0; j < n; j++)
			a[out][j] *= 1 / pivot_element;
		b[out] *= 1 / pivot_element;

	/// Eliminate all the rows except the outgoing
		for (int i = 0; i < m; i++) if (i != out) {
			double factor = a[i][in] / a[out][in];

			for (int j = 0; j < n; j++)
				a[i][j] -= factor * a[out][j];
			b[i] -= factor * b[out];
		}

		xb[out] = in;
	}

	/// Print the final result
	fmax = 0;
	for (int i = 0; i < m; i++)
		fmax += c[xb[i]] * b[i];
	
	print_simplex_tableau(n, m, a, b, c, xb, z, c_z, theta, fmax);
}

