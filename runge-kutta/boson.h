/*
 * boson.h
 *
 *  Created on: Oct 21, 2016
 *      Author: phil
 */

#ifndef BOSON_H_
#define BOSON_H_

#include <iostream>
#include <vector>
#include <sdf.h>
#include <bbhutil.h>
#include <functional>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iomanip>
#include <gmp.h>

#define BOSON_SYSTEM_DIMENSION 4
#define DEFAULT_LAMBDA 100.0
#define DEFAULT_MASS 1.0
#define DEFAULT_SIGMA 0.05
#define DEFAULT_METRIC_THRESHOLD 100.0
#define DEFAULT_FIELD_THRESHOLD 2.0
#define DEFAULT_FIELD_DERIVATIVE_THRESHOLD 100.0
#define DEFAULT_KAPPA 8.0 * M_PI

template <size_t size, typename T>
struct vec
{
public:
	T v[size];

	inline double& operator[](const size_t i)
	{
		return v[i];
	}
};

template <size_t size, typename T>
using d_func = std::function<vec<size, T>(const T, vec<size, T>)>;
using bosonvec = vec<BOSON_SYSTEM_DIMENSION, double>;
using bosonfunc = std::function<bosonvec(const double, bosonvec, const double)>;

enum BosonVars
{
	VAR_A0 = 0,
	VAR_ALPHA0 = 1,
	VAR_PHI0 = 2,
	VAR_PHI1 = 3
};

enum BosonStarType
{
	BOSON_STAR_MINI = 0,
	BOSON_STAR_MASSIVE = 1,
	BOSON_STAR_SOLITONIC = 2
};

inline double sqr(const double x) { return x * x; }
template <size_t n>
inline double ipow(const double x) { return x * ipow<n - 1>(x); }
template <>
inline double ipow<1>(const double x) { return x; }
static double Power(const double x, const int i) { return (i == 1) ? x : x * Power(x, i - 1); }

template <size_t size, typename T>
inline vec<size, T> add(vec<size, T> x1, vec<size, T> x2)
{
	vec<size, T> v;
	for (size_t i = 0; i < size; ++i)
		v[i] = x1[i] + x2[i];
	return v;
}

template <size_t size, typename T>
inline vec<size, T> mult(const double scalar, vec<size, T> x)
{
	vec<size, T> v;
	for (size_t i = 0; i < size; ++i)
			v[i] = scalar * x[i];
	return v;
}

template <size_t size, typename T>
struct IntegrateParams
{
	vec<size, T> initial_value;
	T domain_start;
	T domain_end;
	T step_size;
	d_func<size, T> func;
	std::vector<vec<size, T>>* result;
	vec<size, T> thresholds;
	size_t divergence_variable;
};

using BosonIntegrateParams = IntegrateParams<BOSON_SYSTEM_DIMENSION, double>;

enum Direction
{
	DIRECTION_POSITIVE,
	DIRECTION_NEGATIVE
};

struct IntegrateResult
{
	bool bThresholdReached;
	Direction divergenceDirection;
	size_t lastIndex;
};

struct BosonSearchParams
{
	double initial_phi_value;
	double domain_start;
	double domain_end;
	double step_size;
	bosonfunc func;
	double bracket_begin;
	double bracket_end;
	double requested_final_bracket_size;
	size_t requested_energy_state;
	std::vector<bosonvec>* result;
};

template <size_t size, typename T>
IntegrateResult integrate(IntegrateParams<size, T>& params)
{
	IntegrateResult returnValue;
	returnValue.bThresholdReached = false;

	auto& result = *params.result;
	double t0 = params.domain_start;
	double h = params.step_size;

	// Resize the output list to correct size
	size_t step_count = (size_t)((params.domain_end - params.domain_start) / params.step_size);
	result.resize(step_count);
	memset(result.data(), 0, sizeof(vec<size, T>) * step_count);

	returnValue.lastIndex = step_count - 1;

	// Set the initial values
	result[0] = params.initial_value;

	// Integrate outwards with K4
	for (size_t i = 0; i < step_count - 1 && !returnValue.bThresholdReached; ++i)
	{
		double t = t0 + i * h;
		vec<size, T> k1 = params.func(t, result[i]);
		vec<size, T> k2 = params.func(t + h / 2, add(result[i], mult(h / 2, k1)));
		vec<size, T> k3 = params.func(t + h / 2, add(result[i], mult(h / 2, k2)));
		vec<size, T> k4 = params.func(t + h, add(result[i], mult(h, k3)));

		result[i + 1] = add(result[i], mult(h / 6, add(k1, add(mult(2, k2), add(mult(2, k3), k4)))));

		// Check if threshold has been reached, break if true
		for (size_t j = 0; j < size; ++j)
		{
			if (std::abs(result[i + 1][j]) > params.thresholds[j])
			{
				std::cout << "Threshold reached! Breaking integration..." << std::endl;

				returnValue.bThresholdReached = true;
				if (result[i + 1][params.divergence_variable] > result[i][params.divergence_variable])
					returnValue.divergenceDirection = DIRECTION_POSITIVE;
				else
					returnValue.divergenceDirection = DIRECTION_NEGATIVE;

				returnValue.lastIndex = i;
				memset(&result[i + 1], 0, sizeof(vec<size, T>));

				break;
			}
		}
	}

	return returnValue;
}

double boson_search(BosonSearchParams& params);

bosonvec func_miniboson(const double r, bosonvec y, const double omega, const double mass);
bosonvec func_miniboson(const double r, bosonvec y, const double omega, const double mass,
		const double kappa);
bosonvec func_massiveboson(const double r, bosonvec y, const double omega,
		const double mass, const double lambda);
bosonvec func_massiveboson(const double r, bosonvec y, const double omega,
		const double mass, const double lambda, const double kappa);
bosonvec func_solitonicboson(const double r, bosonvec y, const double omega,
		const double mass, const double sigma);
bosonvec func_solitonicboson(const double r, bosonvec y, const double omega,
		const double mass, const double sigma, const double kappa);

size_t boson_get_zeroes(std::vector<bosonvec>& bosonStar, const size_t lastIndex);

#endif /* BOSON_H_ */
