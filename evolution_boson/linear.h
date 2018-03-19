#ifndef LINEAR_H_
#define LINEAR_H_

#include <vector>
#include <stdlib.h>
#include <assert.h>
#include <bbhutil.h>
#include <sdf.h>
#include <memory.h>
#include <math.h>

struct LinearSystemEntry
{
	int VariableId;
	double Coefficient;
};

enum LinearSystemResult
{
	LINEAR_SYSTEM_RESULT_SUCCESS,
	LINEAR_SYSTEM_RESULT_FAILURE,
	LINEAR_SYSTEM_RESULT_NO_CONVERGENCE
};

char debug_name[] = "debug";

class LinearSystem
{
public:
	std::vector<std::vector<LinearSystemEntry> > Entries;

	LinearSystem(const int variables)
	{
		Entries.resize(variables);
	}

	void Clear()
	{
		for (auto& row : Entries)
			row.clear();
	}

	void Add(const int row, const int column, const double value)
	{
		Entries[row].push_back({ column, value });
	}
};


class GaussSolver
{
private:
	double* gaussSolverData;
	size_t gaussSolverDataSize;

public:
	GaussSolver() : gaussSolverData(nullptr), gaussSolverDataSize(0)
	{
	}

	~GaussSolver()
	{
		if (gaussSolverData != nullptr)
			delete [] gaussSolverData;
	}

	LinearSystemResult SolveGaussSiedel(const LinearSystem& system, const double input[], double output[], const double epsilon, const size_t max_iterations)
	{
		double* guess = new double[system.Entries.size()];
		memset(guess, 0, sizeof(double) * system.Entries.size());

		auto result = SolveGaussSiedel(system, input, guess, output, epsilon, max_iterations);

		delete [] guess;
		return result;
	}

	LinearSystemResult SolveGaussSiedel(const LinearSystem& system, const double input[], const double guess[], double output[], const double epsilon, const size_t max_iterations)
	{
		int n = system.Entries.size();

		if (gaussSolverDataSize != system.Entries.size())
		{
			if (gaussSolverData != nullptr)
				delete [] gaussSolverData;
			gaussSolverData = new double[system.Entries.size()];
			gaussSolverDataSize = system.Entries.size();
		}

		memcpy(gaussSolverData, guess, sizeof(double) * n);

		bool bConvergence = false;
		double newValue = 0.0;
		double diff = 0.0;
		double error = 0.0;

		size_t iteration = 0;
		for (; !bConvergence && (iteration < max_iterations); ++iteration)
		{
			bConvergence = true;
			error = 0.0;

			for (int i = 0; i < n; ++i)
			{
				double sigma = 0.0;

				auto& expr = system.Entries[i];
				const LinearSystemEntry* diag = nullptr;

				for (auto& term : expr)
				{
					if (i != term.VariableId)
						sigma += term.Coefficient * gaussSolverData[term.VariableId];
					else
						diag = &term;
				}

				assert(diag != nullptr);

				newValue = (input[i] - sigma) / diag->Coefficient;
				diff = newValue - gaussSolverData[i];
				error = fabs(diff);
				bConvergence = bConvergence && (sqrt(error) < epsilon);

				gaussSolverData[i] = newValue;
			}
		}

		if (iteration == max_iterations)
		{
			std::cout << "No convergence!" << std::endl;
			return LINEAR_SYSTEM_RESULT_NO_CONVERGENCE;
		}

		memcpy(output, gaussSolverData, sizeof(double) * n);

		return LINEAR_SYSTEM_RESULT_SUCCESS;
	}
};

#endif
