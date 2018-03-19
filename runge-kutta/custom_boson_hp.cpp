/*
 * custom_boson_hp.cpp
 *
 *  Created on: Jan 9, 2017
 *      Author: phil
 */

#include "custom_boson_hp.h"

using namespace std;

void init_custom_boson_search_params(BosonSearchParamsCustomHP* params)
{
	mpf_init(params->domain_start);
	mpf_init(params->domain_end);
	mpf_init(params->step_size);
}

void clear_custom_boson_search_params(BosonSearchParamsCustomHP* params)
{
	mpf_clear(params->domain_start);
	mpf_clear(params->domain_end);
	mpf_clear(params->step_size);
}

void boson_search_custom_hp(BosonSearchParamsCustomHP& params, ICustomBosonSearchHP* searchInterface)
{
	auto& data = *params.result;

	mpf_t init_phi;
	mpf_t omega;

	mpf_init(init_phi);
	mpf_init(omega);

	BosonIntegrateParamsHP integrateParams;
	init_integrate_params(&integrateParams);
	mpf_set(integrateParams.domain_start, params.domain_start);
	mpf_set(integrateParams.domain_end, params.domain_end);

	mpf_set_ui(integrateParams.initial_value[0], 1);
	mpf_set_ui(integrateParams.initial_value[1], 1);
	mpf_set_ui(integrateParams.initial_value[3], 0);

	mpf_set(integrateParams.step_size, params.step_size);
	integrateParams.result = params.result;
	integrateParams.divergence_variable = VAR_PHI0;

	IntegrateResult result;
	result.lastIndex = 0;

	while (!searchInterface->GetNextParameters(&omega, &init_phi))
	{
		mpf_set(integrateParams.initial_value[2], init_phi);
		bosonvec thresholds = { DEFAULT_METRIC_THRESHOLD, DEFAULT_METRIC_THRESHOLD,
			DEFAULT_FIELD_THRESHOLD * mpf_get_d(init_phi), DEFAULT_FIELD_DERIVATIVE_THRESHOLD };
		set_vec_mpf_d(&integrateParams.thresholds, thresholds);

		bosonexprvec boson_expr = params.func(omega);
		integrateParams.func = [&boson_expr](vec<BOSON_SYSTEM_DIMENSION + 1, mpf_t>& vars, bosonvec_hp* output) {
			func_boson_hp<BOSON_SYSTEM_DIMENSION>(vars, output, boson_expr);
		};

		result = integrate_hp(integrateParams);

		searchInterface->ProcessResult(result, data);

		boson_expr.release();
	}

	clear_integrate_params(&integrateParams);

	mpf_clear(omega);
	mpf_clear(init_phi);

	// Clear to zero everything after final index
	for (size_t i = result.lastIndex + 1; i < params.result->size(); ++i)
		for (size_t j = 0; j < BOSON_SYSTEM_DIMENSION; ++j)
			mpf_set_ui((*params.result)[i][j], 0);
}

ManualRangeBosonSearchHP::ManualRangeBosonSearchHP(const ManualRangeBosonSearchParamsHP& params)
	: search_variable(params.search_variable),
	  suggested_bisection_count(params.suggested_bisection_count),
	  bisection_count_max(params.bisection_count_max),
	  current_bisection(0U),
	  b_inverted_search(params.b_inverted_search),
	  b_terminate(false)
{
	mpf_init(lower);
	mpf_init(upper);

	mpf_init(bisection_var);
	mpf_init(other_var);

	mpf_set_d(lower, params.lower_bound);
	mpf_set_d(upper, params.upper_bound);
	mpf_add(bisection_var, upper, lower);
	mpf_div_ui(bisection_var, bisection_var, 2UL);

	mpf_set_d(other_var, params.fix_variable_value);
}

ManualRangeBosonSearchHP::~ManualRangeBosonSearchHP()
{
	mpf_clear(lower);
	mpf_clear(upper);

	mpf_clear(bisection_var);
	mpf_clear(other_var);
}

bool ManualRangeBosonSearchHP::GetNextParameters(mpf_t* omegaOut, mpf_t* initPhiOut)
{
	switch (search_variable)
	{
	case MANUAL_RANGE_SEARCH_VAR_OMEGA:
		mpf_set(*omegaOut, bisection_var);
		mpf_set(*initPhiOut, other_var);
		break;

	case MANUAL_RANGE_SEARCH_VAR_PHI0:
		mpf_set(*omegaOut, other_var);
		mpf_set(*initPhiOut, bisection_var);
		break;
	}

	cout << endl;
	cout << "Iteration: " << current_bisection + 1 << " / " << suggested_bisection_count << endl;
	cout << "Current bracket: [" << std::setprecision(17) << mpf_get_d(lower) << ", " << mpf_get_d(upper) << "]" << endl;

	return b_terminate;
}

void ManualRangeBosonSearchHP::ProcessResult(IntegrateResult& result, std::vector<bosonvec_hp>& resultData)
{
	++current_bisection;

	if (result.bThresholdReached)
	{
		if (current_bisection >= bisection_count_max)
			b_terminate = true;
	}
	else
	{
		if (current_bisection >= suggested_bisection_count)
			b_terminate = true;
	}

	cout << "LAST INDEX VALUE: " << mpf_get_d(resultData[result.lastIndex - 1][VAR_PHI0]) << endl;

	if (!b_inverted_search)
	{
		if (mpf_cmp_ui(resultData[result.lastIndex - 1][VAR_PHI0], 0UL) >= 0)
			mpf_set(lower, bisection_var);
		else
			mpf_set(upper, bisection_var);
	}
	else
	{
		if (mpf_cmp_ui(resultData[result.lastIndex - 1][VAR_PHI0], 0UL) >= 0)
			mpf_set(upper, bisection_var);
		else
			mpf_set(lower, bisection_var);
	}

	mpf_add(bisection_var, upper, lower);
	mpf_div_ui(bisection_var, bisection_var, 2UL);
	cout << "CURRENT BISECT " << mpf_get_d(bisection_var) << endl;
}

const mpf_t* ManualRangeBosonSearchHP::GetBisectionVar() const
{
	return &bisection_var;
}

const mpf_t* ManualRangeBosonSearchHP::GetOtherVar() const
{
	return &other_var;
}

