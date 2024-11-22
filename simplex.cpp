#include <iostream>
#include <vector>

using namespace std;

#define RESET "\033[0m"
#define GREEN "\x1B[32m"

template <typename T>
using matrix = vector<vector<T>>;

class simplex_tableau
{
	/// Initialize all the inputs
private:
	size_t n, m;
	matrix<double> a;
	vector<double> b, c;
	vector<size_t> xb;

	vector<double> z, c_z, theta;

	size_t in, out;

public:
	simplex_tableau(size_t _n, size_t _m,
					const matrix<double> &_a,
					const vector<double> &_b,
					const vector<double> &_c,
					const vector<size_t> &_xb)
		: n{_n}, m{_m}, a{_a}, b{_b}, c{_c}, xb{_xb},
		  z(n), c_z(n), theta(m)
	{
		/// Check if any b is negative
		for (int i = 0; i < m; i++)
			if (b[i] < 0)
			{
				for (int j = 0; j < n; j++)
					a[i][j] = -a[i][j];
				b[i] = -b[i];
			}
	}

	size_t find_pivot_column()
	{
		if (in == size_t(-1))
			return -1;

		/// Compute the z_j and c_j - z_j rows
		for (int j = 0; j < n; j++)
		{
			z[j] = 0;
			for (int i = 0; i < m; i++)
				z[j] += c[xb[i]] * a[i][j],
					c_z[j] = c[j] - z[j];
		}

		/// Determine the next incoming variable
		in = -1; // < n

		for (int j = 0; j < n; j++)
			if (c_z[j] > 0 && (in == -1 || (c_z[j] > c_z[in])))
				in = j;

		return in;
	}

	size_t find_pivot_row()
	{
		if (out == size_t(-1))
			return -1;

		/// Compute the theta column
		for (int i = 0; i < m; i++)
			theta[i] = b[i] / a[i][in];

		/// Determine the next outgoing variable
		out = -1; // < m

		for (int i = 0; i < m; i++)
			if (theta[i] > 0 && (out == -1 || theta[i] < theta[out]))
				out = i;

		return out;
	}

	bool check_if_optimal()
	{
		return !(find_pivot_column() != size_t(-1) && find_pivot_row() != size_t(-1));
	}

	bool is_optimal() const
	{
		return (in == size_t(-1) || out == size_t(-1));
	}

	void pivot_variables()
	{
		/// Scale the outgoing variable row
		double pivot_element = a[out][in];

		for (int j = 0; j < n; j++)
			a[out][j] *= 1 / pivot_element;
		b[out] *= 1 / pivot_element;

		/// Eliminate all the rows except the outgoing
		for (int i = 0; i < m; i++)
			if (i != out)
			{
				double factor = a[i][in] / a[out][in];

				for (int j = 0; j < n; j++)
					a[i][j] -= factor * a[out][j];
				b[i] -= factor * b[out];
			}

		xb[out] = in;
	}

	void print_tableau() const
	{
		cout << (*this);
	}

	double fmax() const
	{
		double fmax = 0;

		for (int i = 0; i < m; i++)
			fmax += c[xb[i]] * b[i];

		return fmax;
	}

	vector<double> coefficients() const
	{
		vector<double> coefficients(n);

		for (int i = 0; i < m; i++)
			coefficients[xb[i]] = b[i];

		return coefficients;
	}

	friend ostream &operator<<(ostream &os, const simplex_tableau &t)
	{
		size_t n = t.n, m = t.m;
		const matrix<double> &a = t.a;
		const vector<double> &b = t.b, &c = t.c;
		const vector<size_t> &xb = t.xb;

		vector<double> z = t.z, c_z = t.c_z, theta = t.theta;
		bool optimal = t.is_optimal();

		size_t in = t.in, out = t.out;

		///----------------------------------------///
		os << "\t\tc_j\t";
		for (int j = 0; j < n; j++)
			os << c[j] << '\t';
		os << "\t\n";
		///----------------------------------------///
		os << "i\tc_b\tx_b\t";
		for (int j = 0; j < n; j++)
			os << "x_" << j + 1 << '\t';
		os << "b_i\t" << ((!optimal) ? "theta_i\n" : "\n");
		///----------------------------------------///
		for (int i = 0; i < m; i++)
		{
			os << i + 1 << '\t' << c[xb[i]] << "\tx_" << xb[i] << '\t';
			for (int j = 0; j < n; j++)
				os << a[i][j] << '\t';
			os << (optimal ? GREEN : RESET) << b[i] << RESET << '\t';
			if (!optimal)
				os << theta[i];
			os << ((!optimal && i == out) ? "\t-> out" : "") << '\n';
		}
		///----------------------------------------///
		os << "\t\tz_j\t";
		for (int j = 0; j < n; j++)
			os << z[j] << '\t';
		os << (optimal ? GREEN : RESET) << t.fmax() << RESET << '\n';
		///----------------------------------------///
		os << "\t\tc_j-z_j\t";
		for (int j = 0; j < n; j++)
			os << c[j] - z[j] << '\t';
		os << '\n';
		///----------------------------------------///
		if (!optimal)
		{
			os << "\t\t\t";
			for (int j = 0; j < in; j++)
				os << '\t';
			os << "^ in\n";
		}
		return os << "--------------------------------------------------------------------------------\n";
	}
};

int main()
{
	/// Initialize all the inputs
	matrix<double> a({{5, 2, 1, 0, 0},
					  {2, 3, 0, 1, 0},
					  {4, 2, 0, 0, 1}});
	vector<double> b({150,
					  100,
					  80});
	vector<double> c({12, 8, 0, 0, 0});
	vector<size_t> xb({2,
					   3,
					   4});

	simplex_tableau solver(5, 3, a, b, c, xb);

	while (!solver.check_if_optimal())
	{
		solver.print_tableau();
		solver.pivot_variables();
	}

	solver.print_tableau();

	cout << "fmax: " << solver.fmax() << endl;
	cout << "coefficients:";
	for (auto c : solver.coefficients())
		cout << ' ' << c;
	cout << endl;
}
