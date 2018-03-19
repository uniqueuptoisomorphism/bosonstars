/*
 * boson_hp.h
 *
 *  Created on: Oct 21, 2016
 *      Author: phil
 */

#ifndef BOSON_HP_H_
#define BOSON_HP_H_

#include "expr.h"
#include "boson.h"

template <size_t size>
struct vec<size, mpf_t>
{
public:
	mpf_t v[size];

	vec()
	{
		for (size_t i = 0; i < size; ++i)
			mpf_init(v[i]);
	}

	~vec()
	{
		for (size_t i = 0; i < size; ++i)
			mpf_clear(v[i]);
	}

	inline mpf_t& operator[](const size_t i)
	{
		return v[i];
	}
};

using bosonvec_hp = vec<BOSON_SYSTEM_DIMENSION, mpf_t>;

template <size_t size>
void set_vec_mpf_d(vec<size, mpf_t>* dest, vec<size, double>& src)
{
	for (size_t i = 0; i < size; ++i)
		mpf_set_d((*dest)[i], src[i]);
}

template <size_t size>
void set_vec_mpf(vec<size, mpf_t>* dest, vec<size, mpf_t>& src)
{
	for (size_t i = 0; i < size; ++i)
		mpf_set((*dest)[i], src[i]);
}

template <size_t size>
struct exprvec
{
public:
	ExprPtr v[size];

	inline ExprPtr& operator[](const size_t i)
	{
		return v[i];
	}

	void release()
	{
		for (size_t i = 0; i < size; ++i)
			v[i].p->release();
	}
};

template <size_t size>
using hp_func = std::function<void(vec<size + 1, mpf_t>&, vec<size, mpf_t>*)>;

template <size_t size>
struct IntegrateParams<size, mpf_t>
{
	vec<size, mpf_t> initial_value;
	mpf_t domain_start;
	mpf_t domain_end;
	mpf_t step_size;
	hp_func<size> func;
	std::vector<vec<size, mpf_t>>* result;
	vec<size, mpf_t> thresholds;
	size_t divergence_variable;
};

using BosonIntegrateParamsHP = IntegrateParams<BOSON_SYSTEM_DIMENSION, mpf_t>;
using bosonexprvec = exprvec<BOSON_SYSTEM_DIMENSION>;
using hp_gen_func = std::function<bosonexprvec(const mpf_t)>;

template <size_t size>
void init_integrate_params(IntegrateParams<size, mpf_t>* params)
{
	mpf_init(params->domain_start);
	mpf_init(params->domain_end);
	mpf_init(params->step_size);
}

template <size_t size>
void clear_integrate_params(IntegrateParams<size, mpf_t>* params)
{
	mpf_clear(params->domain_start);
	mpf_clear(params->domain_end);
	mpf_clear(params->step_size);
}

struct BosonSearchParamsHP
{
	mpf_t initial_phi_value;
	mpf_t domain_start;
	mpf_t domain_end;
	mpf_t step_size;
	hp_gen_func func;
	mpf_t bracket_begin;
	mpf_t bracket_end;
	size_t suggested_bisection_count;
	size_t bisection_count_max;
	size_t requested_energy_state;
	std::vector<bosonvec_hp>* result;
};

void init_boson_search_params(BosonSearchParamsHP* params);
void clear_boson_search_params(BosonSearchParamsHP* params);

template <size_t size>
void print_vars(vec<size, mpf_t>& vars)
{
	for (size_t i = 0; i < size; ++i)
		std::cout << mpf_get_d(vars[i]) << "\t";
	std::cout << std::endl;
}

template <size_t size>
IntegrateResult integrate_hp(IntegrateParams<size, mpf_t>& params)
{
	IntegrateResult returnValue;
	returnValue.bThresholdReached = false;

	auto& result = *params.result;
	mpf_t t0, h, t, h_half, tmp, tmp2;
	mpf_init(t0);
	mpf_init(h);
	mpf_init(h_half);
	mpf_init(tmp);
	mpf_init(tmp2);
	mpf_init(t);

	mpf_set(t0, params.domain_start);
	mpf_set(h, params.step_size);
	mpf_div_ui(h_half, h, 2);

	// Resize the output list to correct size
	mpf_sub(tmp, params.domain_end, params.domain_start);
	mpf_div(tmp, tmp, params.step_size);
	size_t step_count = (size_t)mpf_get_d(tmp);
	result.resize(step_count);
	returnValue.lastIndex = step_count - 1;

	// Set the initial values
	set_vec_mpf(&result[0], params.initial_value);

	vec<size, mpf_t> k1;
	vec<size, mpf_t> k2;
	vec<size, mpf_t> k3;
	vec<size, mpf_t> k4;

	vec<size + 1, mpf_t> vars;

	// Note that vars[0] is the independent variable t
	// Integrate outwards with K4
	for (size_t i = 0; i < step_count - 1 && !returnValue.bThresholdReached; ++i)
	{
		mpf_mul_ui(tmp, h, i);
		mpf_add(t, t0, tmp);
		mpf_set(vars[0], t);

		for (size_t j = 0; j < size; ++j)
			mpf_set(vars[j + 1], result[i][j]);
		params.func(vars, &k1);

		mpf_add(vars[0], h_half, t);
		for (size_t j = 0; j < size; ++j)
		{
			mpf_mul(vars[j + 1], h_half, k1[j]);
			mpf_add(vars[j + 1], result[i][j], vars[j + 1]);
		}
		params.func(vars, &k2);

		for (size_t j = 0; j < size; ++j)
		{
			mpf_mul(vars[j + 1], h_half, k2[j]);
			mpf_add(vars[j + 1], result[i][j], vars[j + 1]);
		}
		params.func(vars, &k3);

		mpf_add(vars[0], h, t);
		for (size_t j = 0; j < size; ++j)
		{
			mpf_mul(vars[j + 1], h, k3[j]);
			mpf_add(vars[j + 1], result[i][j], vars[j + 1]);
		}
		params.func(vars, &k4);

		for (size_t j = 0; j < size; ++j)
		{
			mpf_mul_ui(tmp, k3[j], 2);
			mpf_add(tmp, k4[j], tmp);
			mpf_mul_ui(tmp2, k2[j], 2);
			mpf_add(tmp, tmp2, tmp);
			mpf_add(tmp, k1[j], tmp);
			mpf_mul(tmp, h, tmp);
			mpf_div_ui(tmp, tmp, 6);
			mpf_add(tmp, result[i][j], tmp);
			mpf_set(result[i + 1][j], tmp);
		}

		// Check if threshold has been reached, break if true
		for (size_t j = 0; j < size; ++j)
		{
			mpf_abs(tmp, result[i + 1][j]);
			if (mpf_cmp(tmp, params.thresholds[j]) > 0)
			{
				std::cout << "Threshold reached! Variables: [";
				for (size_t k = 0; k < size; ++k)
					std::cout << mpf_get_d(result[i + 1][k]) << ", ";
				std::cout << "]" << std::endl;

				returnValue.bThresholdReached = true;
				if (mpf_cmp(result[i + 1][params.divergence_variable], result[i][params.divergence_variable]))
					returnValue.divergenceDirection = DIRECTION_POSITIVE;
				else
					returnValue.divergenceDirection = DIRECTION_NEGATIVE;

				returnValue.lastIndex = i;
				break;
			}
		}
	}

	return returnValue;
}

template <size_t size>
void func_boson_hp(vec<size + 1, mpf_t>& vars,
	vec<size, mpf_t>* output, exprvec<size>& exprs)
{
	for (size_t i = 0; i < size; ++i)
	{
		mpf_t* result = exprs[i].p->eval(vars.v);
		mpf_set((*output)[i], *result);
	}
}

void boson_search_hp(BosonSearchParamsHP& params, mpf_t* omegaOut);

bosonexprvec gen_func_miniboson_hp(const mpf_t omega_t, const mpf_t mass_t);
bosonexprvec gen_func_massiveboson_hp(const mpf_t omega_t, const mpf_t mass_t, const mpf_t lambda_t);
bosonexprvec gen_func_solitonicboson_hp(const mpf_t omega_t, const mpf_t mass_t, const mpf_t sigma_t);

bosonexprvec gen_func_miniboson_hp(const mpf_t omega_t, const mpf_t mass_t, const mpf_t kappa_t);
bosonexprvec gen_func_massiveboson_hp(const mpf_t omega_t, const mpf_t mass_t, const mpf_t lambda_t, const mpf_t kappa_t);
bosonexprvec gen_func_solitonicboson_hp(const mpf_t omega_t, const mpf_t mass_t, const mpf_t sigma_t, const mpf_t kappa_t);

size_t boson_get_zeroes_hp(std::vector<bosonvec_hp>& bosonStar, const size_t lastIndex);

#endif /* BOSON_HP_H_ */
