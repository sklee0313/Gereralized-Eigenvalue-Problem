// C++ program for Implementation of
// Reverse Cuthill Mckee Algorithm

#include <bits/stdc++.h>
#include <Eigen/Sparse>

// using namespace std;

std::vector<double> globalDegree;
int iter_findIndex;

int findIndex(std::vector<std::pair<int, double>> &a, const int &x)
{
	// for (int i = 0; i < a.size(); i++)
	// 	if (a[i].first == x)
	// 		return i;
	iter_findIndex = 0; // counter
	for (std::vector<std::pair<int, double>>::iterator it = a.begin(); it != a.end(); it++, iter_findIndex++)
		if (it->first == x)
			return iter_findIndex;
	return -1;
}

bool compareDegree(int &i, int &j)
{
	return ::globalDegree[i] < ::globalDegree[j];
}

template <typename T>
std::ostream &operator<<(std::ostream &out, std::vector<T> const &v)
{
	for (int i = 0; i < v.size(); i++)
		out << v[i] << ' ';
	return out;
}

class ReorderingSSM
{
private:
	// vector<vector<double> > _matrix;
	Eigen::SparseMatrix<double> _matrix;

public:
	// Constructor and Destructor
	ReorderingSSM(Eigen::SparseMatrix<double> m)
	{
		_matrix = m;
	}

	ReorderingSSM() {}
	~ReorderingSSM() {}

	// class methods

	// Function to generate degree of all the nodes
	std::vector<double> degreeGenerator()
	{

		std::vector<double> degrees;
		// int i = 0; i < _matrix.size(); i++
		for (int i = 0; i < _matrix.outerSize(); ++i)
		{
			double count = 0;
			// int j = 0; j < _matrix[0].size(); j++
			for (Eigen::SparseMatrix<double>::InnerIterator it(_matrix, i); it; ++it)
			{
				count += 1; //(it.value()); //_matrix[i][j];
			}

			degrees.push_back(count);
		}

		return degrees;
	}

	// Implementation of Cuthill-Mckee algorithm
	std::vector<int> CuthillMckee()
	{
		std::vector<double> degrees = degreeGenerator();
		//	std::cout << "degrees for RCM is calculated" << std::endl;
		::globalDegree = degrees;

		std::queue<int> Q;
		std::vector<int> R;
		std::vector<std::pair<int, double>> notVisited;

		int size_degrees = degrees.size();
		for (int i = 0; i < size_degrees; i++)
			notVisited.push_back(std::make_pair(i, degrees[i]));

		// Vector notVisited helps in running BFS
		// even when there are dijoind graphs
		while (notVisited.size())
		{

			int minNodeIndex = 0;

			// find minimum in notvisted
			int size_notVisited = notVisited.size();
			for (int i = 0; i < size_notVisited; i++)
				if (notVisited[i].second < notVisited[minNodeIndex].second)
					minNodeIndex = i;

			Q.push(notVisited[minNodeIndex].first);

			notVisited.erase(notVisited.begin() + findIndex(notVisited,
															notVisited[Q.front()].first));

			// Simple BFS
			std::vector<int> toSort;
			int isfindIndex;
			while (!Q.empty())
			{
				toSort.clear();

				for (SparseMatrix<double>::InnerIterator it(_matrix, Q.front()); it; ++it)
				{
					isfindIndex = findIndex(notVisited, it.row());
					if (it.row() != Q.front() && isfindIndex != -1) //
					{
						toSort.push_back(it.row());
						notVisited.erase(notVisited.begin() + isfindIndex);
					}
				}

				sort(toSort.begin(), toSort.end(), compareDegree);

				int size_toSort = toSort.size();
				for (int i = 0; i < size_toSort; i++)
					Q.push(toSort[i]);

				R.push_back(Q.front());
				Q.pop();
			}
		}

		return R;
	}

	// Implementation of reverse Cuthill-Mckee algorithm
	Eigen::VectorXi ReverseCuthillMckee()
	{

		std::vector<int> cuthill = CuthillMckee();

		int n = cuthill.size();

		if (n % 2 == 0)
			n -= 1;

		n = n / 2;

		for (int i = 0; i <= n; i++)
		{
			int j = cuthill[cuthill.size() - 1 - i];
			cuthill[cuthill.size() - 1 - i] = cuthill[i];
			cuthill[i] = j;
		}
		return Eigen::VectorXi::Map(cuthill.data(), static_cast<int>(cuthill.size()));
	}
};
