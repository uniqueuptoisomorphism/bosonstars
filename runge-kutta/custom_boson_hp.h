/*
 * custom_boson_hp.h
 *
 *  Created on: Jan 9, 2017
 *      Author: phil
 */

#ifndef CUSTOM_BOSON_HP_H_
#define CUSTOM_BOSON_HP_H_

#include "boson_hp.h"

struct BosonSearchParamsCustomHP
{
	mpf_t domain_start;
	mpf_t domain_end;
	mpf_t step_size;
	hp_gen_func func;
	std::vector<bosonvec_hp>* result;
};

void init_custom_boson_search_params(BosonSearchParamsCustomHP* params);
void clear_custom_boson_search_params(BosonSearchParamsCustomHP* params);

class ICustomBosonSearchHP
{
public:
	virtual ~ICustomBosonSearchHP() { }

	/* Gets the next parameters to run the integrator on.
	 * Return value: boolean flag for whether search should terminate.
	 */
	virtual bool GetNextParameters(mpf_t* omegaOut, mpf_t* initPhiOut) = 0;

	/* Process the results of integration, called after GetNextOmega and the
	 * integrator have been run.
	 */
	virtual void ProcessResult(IntegrateResult& result, std::vector<bosonvec_hp>& resultData) = 0;
};

void boson_search_custom_hp(BosonSearchParamsCustomHP& params, ICustomBosonSearchHP* searchInterface);

enum ManualRangeSearchVariable
{
	MANUAL_RANGE_SEARCH_VAR_OMEGA = 0,
	MANUAL_RANGE_SEARCH_VAR_PHI0 = 1
};

struct ManualRangeBosonSearchParamsHP
{
	double upper_bound;
	double lower_bound;
	double fix_variable_value;
	ManualRangeSearchVariable search_variable;
	size_t suggested_bisection_count;
	size_t bisection_count_max;
	bool b_inverted_search;
};

class ManualRangeBosonSearchHP : public ICustomBosonSearchHP
{
private:
	ManualRangeSearchVariable search_variable;
	mpf_t lower;
	mpf_t upper;
	size_t suggested_bisection_count;
	size_t bisection_count_max;
	size_t current_bisection;
	bool b_inverted_search;
	bool b_terminate;

	mpf_t bisection_var;
	mpf_t other_var;

public:
	ManualRangeBosonSearchHP(const ManualRangeBosonSearchParamsHP& params);

	virtual ~ManualRangeBosonSearchHP();

	/* Gets the next parameters to run the integrator on.
	 * Return value: boolean flag for whether search should terminate.
	 */
	virtual bool GetNextParameters(mpf_t* omegaOut, mpf_t* initPhiOut);

	/* Process the results of integration, called after GetNextOmega and the
	 * integrator have been run.
	 */
	virtual void ProcessResult(IntegrateResult& result, std::vector<bosonvec_hp>& resultData);

	const mpf_t* GetBisectionVar() const;
	const mpf_t* GetOtherVar() const;
};

#endif /* CUSTOM_BOSON_HP_H_ */
