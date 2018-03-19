#include "boson_hp.h"

using namespace std;

bosonexprvec gen_func_miniboson_hp_old(const double omega, const double mass)
{
	ExprPtr r = new ExprVar(0);
	ExprPtr A0 = new ExprVar(1);
	ExprPtr alpha0 = new ExprVar(2);
	ExprPtr phi0 = new ExprVar(3);
	ExprPtr phi1 = new ExprVar(4);

	ExprPtr	resultA0 = A0 / (2 * r) - ipow<3>(A0) / (2 * r) +
		(sqr(mass) / 2) * r * ipow<3>(A0) * sqr(phi0) +
		(r * sqr(omega) * ipow<3>(A0) * sqr(phi0)) / (2 * sqr(alpha0)) +
		(r / 2) * A0 * sqr(phi1);

	ExprPtr resultAlpha0 = -(1 / (
		2 * r * alpha0)) *(sqr(alpha0) - sqr(A0) * sqr(alpha0) -
		sqr(r) * sqr(omega) * sqr(A0) * sqr(phi0) +
		sqr(mass) * sqr(r) * sqr(A0) * sqr(alpha0) * sqr(phi0) -
		sqr(r) * sqr(alpha0) * sqr(phi1));

	ExprPtr resultPhi0 =  phi1;

	ExprPtr resultPhi1 =  -(1 / (r * sqr(alpha0))) * (r * sqr(omega) * sqr(A0) * phi0 -
		r * sqr(A0) * sqr(alpha0) * phi0 + sqr(alpha0) * phi1 +
		sqr(A0) * sqr(alpha0) * phi1 -
		sqr(mass) * sqr(r) * sqr(A0) * sqr(alpha0) * sqr(phi0) * phi1);

	return { resultA0, resultAlpha0, resultPhi0, resultPhi1 };
}

bosonexprvec gen_func_miniboson_hp(const mpf_t omega_t, const mpf_t mass_t)
{
	ExprPtr r = new ExprVar(0);
	ExprPtr A0 = new ExprVar(1);
	ExprPtr alpha0 = new ExprVar(2);
	ExprPtr phi0 = new ExprVar(3);
	ExprPtr phi1 = new ExprVar(4);

	ExprPtr omega = new ExprDataRootT<mpf_t>(omega_t);
	ExprPtr mass = new ExprDataRootT<mpf_t>(mass_t);

	// At r = 0
	ExprPtr resultA0_0 = new ExprDataRootT<unsigned long int>(0);
	ExprPtr resultAlpha0_0 = new ExprDataRootT<unsigned long int>(0);
	ExprPtr resultPhi0_0 = phi1;
	ExprPtr resultPhi1_0 = -((sqr(omega) * sqr(A0) * phi0 - sqr(A0) * sqr(alpha0) * phi0) / (3 * sqr(alpha0)));

	// For r > 0
	ExprPtr	resultA0 = A0 / (2 * r) - ipow<3>(A0) / (2 * r) +
		(sqr(mass) * r * ipow<3>(A0) * sqr(phi0)) / 2 +
		(r * sqr(omega) * ipow<3>(A0) * sqr(phi0)) / (2 * sqr(alpha0)) +
		(r * A0 * sqr(phi1)) / 2;

	ExprPtr resultAlpha0 = -(1 / (
		2 * r * alpha0)) *(sqr(alpha0) - sqr(A0) * sqr(alpha0) -
		sqr(r) * sqr(omega) * sqr(A0) * sqr(phi0) +
		sqr(mass) * sqr(r) * sqr(A0) * sqr(alpha0) * sqr(phi0) -
		sqr(r) * sqr(alpha0) * sqr(phi1));

	ExprPtr resultPhi0 =  phi1;

	ExprPtr resultPhi1 =  -(1 / (r * sqr(alpha0))) * (r * sqr(omega) * sqr(A0) * phi0 -
		r * sqr(A0) * sqr(alpha0) * phi0 + sqr(alpha0) * phi1 +
		sqr(A0) * sqr(alpha0) * phi1 -
		sqr(mass) * sqr(r) * sqr(A0) * sqr(alpha0) * sqr(phi0) * phi1);

	ExprPtr cmpA0 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0, resultA0, resultA0_0, resultA0_0);
	ExprPtr cmpAlpha0 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0, resultAlpha0, resultAlpha0_0, resultAlpha0_0);
	ExprPtr cmpPhi0 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0, resultPhi0, resultPhi0_0, resultPhi0_0);
	ExprPtr cmpPhi1 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0, resultPhi1, resultPhi1_0, resultPhi1_0);

	return { cmpA0, cmpAlpha0, cmpPhi0, cmpPhi1 };
}

bosonexprvec gen_func_miniboson_hp(const mpf_t omega_t, const mpf_t mass_t, const mpf_t kappa_t)
{
	ExprPtr r = new ExprVar(0);
	ExprPtr A0 = new ExprVar(1);
	ExprPtr alpha0 = new ExprVar(2);
	ExprPtr phi0 = new ExprVar(3);
	ExprPtr phi1 = new ExprVar(4);

	ExprPtr omega = new ExprDataRootT<mpf_t>(omega_t);
	ExprPtr m = new ExprDataRootT<mpf_t>(mass_t);
	ExprPtr kappa = new ExprDataRootT<mpf_t>(kappa_t);

	ExprPtr resultA0_0 = 0;
	ExprPtr resultAlpha0_0 = 0;
	ExprPtr resultPhi0_0 = phi1;
	ExprPtr resultPhi1_0 = -(Power(omega,2)*Power(A0,2)*phi0 - \
	Power(m,2)*Power(A0,2)*Power(alpha0,2)*phi0)/(3.*Power(alpha0,2));

	ExprPtr resultA0 = A0/(2.*r) - Power(A0,3)/(2.*r) + \
	(Power(m,2)*r*kappa*Power(A0,3)*Power(phi0,2))/2. + \
	(r*kappa*Power(omega,2)*Power(A0,3)*Power(phi0,2))/(2.*Power(alpha0,2)\
	) + (r*kappa*A0*Power(phi1,2))/2.;
	ExprPtr resultAlpha0 = -(Power(alpha0,2) - \
	Power(A0,2)*Power(alpha0,2) - \
	Power(r,2)*kappa*Power(omega,2)*Power(A0,2)*Power(phi0,2) + \
	Power(m,2)*Power(r,2)*kappa*Power(A0,2)*Power(alpha0,2)*Power(phi0,2) \
	- Power(r,2)*kappa*Power(alpha0,2)*Power(phi1,2))/(2.*r*alpha0);
	ExprPtr resultPhi0 = phi1;
	ExprPtr resultPhi1 = -((r*Power(omega,2)*Power(A0,2)*phi0 - \
	Power(m,2)*r*Power(A0,2)*Power(alpha0,2)*phi0 + Power(alpha0,2)*phi1 \
	+ Power(A0,2)*Power(alpha0,2)*phi1 - \
	Power(m,2)*Power(r,2)*kappa*Power(A0,2)*Power(alpha0,2)*Power(phi0,2)*\
	phi1)/(r*Power(alpha0,2)));

	ExprPtr cmpA0 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0UL, resultA0, resultA0_0, resultA0_0);
	ExprPtr cmpAlpha0 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0UL, resultAlpha0, resultAlpha0_0, resultAlpha0_0);
	ExprPtr cmpPhi0 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0UL, resultPhi0, resultPhi0_0, resultPhi0_0);
	ExprPtr cmpPhi1 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0UL, resultPhi1, resultPhi1_0, resultPhi1_0);

	return { cmpA0, cmpAlpha0, cmpPhi0, cmpPhi1 };
}

bosonexprvec gen_func_massiveboson_hp(const mpf_t omega_t, const mpf_t mass_t, const mpf_t lambda_t)
{
	ExprPtr r = new ExprVar(0);
	ExprPtr A0 = new ExprVar(1);
	ExprPtr alpha0 = new ExprVar(2);
	ExprPtr phi0 = new ExprVar(3);
	ExprPtr phi1 = new ExprVar(4);

	ExprPtr omega = new ExprDataRootT<mpf_t>(omega_t);
	ExprPtr mass = new ExprDataRootT<mpf_t>(mass_t);
	ExprPtr lambda = new ExprDataRootT<mpf_t>(lambda_t);

	ExprPtr resultA0_0 = new ExprDataRootT<unsigned long int>(0);
	ExprPtr resultAlpha0_0 = new ExprDataRootT<unsigned long int>(0);
	ExprPtr resultPhi0_0 = phi1;
	ExprPtr resultPhi1_0 = sqr(A0) * phi0 * (- sqr(omega) + sqr(mass) * sqr(alpha0) + lambda * sqr(alpha0) * sqr(phi0)) / (3 * sqr(alpha0));

	ExprPtr resultA0 = - (A0 / 4) * (-(2 / r) + sqr(A0) * (2 / r +
		r * sqr(phi0) * (-2 * sqr(mass) -
		(2 * sqr(omega)) / sqr(alpha0) - lambda * sqr(phi0))) -
		2 * r * sqr(phi1));

	ExprPtr resultAlpha0 = (1 / (4 * r * alpha0)) *
		(2 * sqr(r) * sqr(omega) * sqr(A0) * sqr(phi0)
		- sqr(alpha0) * (2 + sqr(A0) * (-2 +
		sqr(r) * sqr(phi0) * (2 * sqr(mass) + lambda * sqr(phi0))) -
		2 * sqr(r) * sqr(phi1)));

	ExprPtr resultPhi0 = phi1;

	ExprPtr resultPhi1 = -(phi1 / r) +
		sqr(A0) * ((sqr(mass) - sqr(omega) / sqr(alpha0)) * phi0 + lambda * sqr(phi0) * phi0 - phi1 / r +
		sqr(mass) * r * sqr(phi0) * phi1 +
		(r / 2) * lambda * sqr(phi0) * sqr(phi0) * phi1);

	ExprPtr cmpA0 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0, resultA0, resultA0_0, resultA0_0);
	ExprPtr cmpAlpha0 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0, resultAlpha0, resultAlpha0_0, resultAlpha0_0);
	ExprPtr cmpPhi0 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0, resultPhi0, resultPhi0_0, resultPhi0_0);
	ExprPtr cmpPhi1 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0, resultPhi1, resultPhi1_0, resultPhi1_0);

	return { cmpA0, cmpAlpha0, cmpPhi0, cmpPhi1 };
}

bosonexprvec gen_func_massiveboson_hp(const mpf_t omega_t, const mpf_t mass_t, const mpf_t lambda_t, const mpf_t kappa_t)
{
	ExprPtr r = new ExprVar(0);
	ExprPtr A0 = new ExprVar(1);
	ExprPtr alpha0 = new ExprVar(2);
	ExprPtr phi0 = new ExprVar(3);
	ExprPtr phi1 = new ExprVar(4);

	ExprPtr omega = new ExprDataRootT<mpf_t>(omega_t);
	ExprPtr m = new ExprDataRootT<mpf_t>(mass_t);
	ExprPtr kappa = new ExprDataRootT<mpf_t>(kappa_t);
	ExprPtr lambda = new ExprDataRootT<mpf_t>(lambda_t);

	ExprPtr resultA0_0 = 0;
	ExprPtr resultAlpha0_0 = 0;
	ExprPtr resultPhi0_0 = phi1;
	ExprPtr resultPhi1_0 = -(Power(omega,2)*Power(A0,2)*phi0 - \
	Power(m,2)*Power(A0,2)*Power(alpha0,2)*phi0 - \
	lambda*Power(A0,2)*Power(alpha0,2)*Power(phi0,3))/(3.*Power(alpha0,2));

	ExprPtr resultA0 = A0/(2.*r) - Power(A0,3)/(2.*r) + \
	(Power(m,2)*r*kappa*Power(A0,3)*Power(phi0,2))/2. + \
	(r*kappa*Power(omega,2)*Power(A0,3)*Power(phi0,2))/(2.*Power(alpha0,2)\
	) + (r*kappa*lambda*Power(A0,3)*Power(phi0,4))/4. + \
	(r*kappa*A0*Power(phi1,2))/2.;
	ExprPtr resultAlpha0 = -(2*Power(alpha0,2) - \
	2*Power(A0,2)*Power(alpha0,2) - \
	2*Power(r,2)*kappa*Power(omega,2)*Power(A0,2)*Power(phi0,2) + \
	2*Power(m,2)*Power(r,2)*kappa*Power(A0,2)*Power(alpha0,2)*Power(phi0,\
	2) + Power(r,2)*kappa*lambda*Power(A0,2)*Power(alpha0,2)*Power(phi0,4)\
	 - 2*Power(r,2)*kappa*Power(alpha0,2)*Power(phi1,2))/(4.*r*alpha0);
	ExprPtr resultPhi0 = phi1;
	ExprPtr resultPhi1 = -(2*r*Power(omega,2)*Power(A0,2)*phi0 - \
	2*Power(m,2)*r*Power(A0,2)*Power(alpha0,2)*phi0 - \
	2*r*lambda*Power(A0,2)*Power(alpha0,2)*Power(phi0,3) + \
	2*Power(alpha0,2)*phi1 + 2*Power(A0,2)*Power(alpha0,2)*phi1 - \
	2*Power(m,2)*Power(r,2)*kappa*Power(A0,2)*Power(alpha0,2)*Power(phi0,\
	2)*phi1 - \
	Power(r,2)*kappa*lambda*Power(A0,2)*Power(alpha0,2)*Power(phi0,4)*\
	phi1)/(2.*r*Power(alpha0,2));

	ExprPtr cmpA0 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0UL, resultA0, resultA0_0, resultA0_0);
	ExprPtr cmpAlpha0 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0UL, resultAlpha0, resultAlpha0_0, resultAlpha0_0);
	ExprPtr cmpPhi0 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0UL, resultPhi0, resultPhi0_0, resultPhi0_0);
	ExprPtr cmpPhi1 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0UL, resultPhi1, resultPhi1_0, resultPhi1_0);

	return { cmpA0, cmpAlpha0, cmpPhi0, cmpPhi1 };
}

bosonexprvec gen_func_solitonicboson_hp(const mpf_t omega_t, const mpf_t mass_t, const mpf_t sigma_t)
{
	ExprPtr r = new ExprVar(0);
	ExprPtr A0 = new ExprVar(1);
	ExprPtr alpha0 = new ExprVar(2);
	ExprPtr phi0 = new ExprVar(3);
	ExprPtr phi1 = new ExprVar(4);

	ExprPtr omega = new ExprDataRootT<mpf_t>(omega_t);
	ExprPtr mass = new ExprDataRootT<mpf_t>(mass_t);
	ExprPtr sigma = new ExprDataRootT<mpf_t>(sigma_t);

	ExprPtr resultA0_0 = new ExprDataRootT<unsigned long int>(0);
	ExprPtr resultAlpha0_0 = new ExprDataRootT<unsigned long int>(0);
	ExprPtr resultPhi0_0 = phi1;
	ExprPtr resultPhi1_0 = -(1 / (3 * ipow<4>(sigma) * sqr(alpha0))) * (ipow<4>(sigma) * sqr(omega) * sqr(A0) * phi0 -
		sqr(mass) * ipow<4>(sigma) * sqr(A0) * sqr(alpha0) * phi0 +
		4 * sqr(mass) * sqr(sigma) * sqr(A0) * sqr(alpha0) * ipow<3>(phi0) -
		3 * sqr(mass) * sqr(A0) * sqr(alpha0) * ipow<5>(phi0));

	ExprPtr resultA0 = A0 / (2 * r) - ipow<3>(A0) / (2 * r) +
		(sqr(mass) / 2) * r * ipow<3>(A0) * sqr(phi0) +
		(r * sqr(omega) * ipow<3>(A0) * sqr(phi0)) / (2 * sqr(alpha0)) -
		(sqr(mass) * r * ipow<3>(A0) * ipow<4>(phi0)) / sqr(sigma) +
		(sqr(mass) * r * ipow<3>(A0) * ipow<6>(phi0)) / (2 * ipow<4>(sigma)) +
		(r / 2) * A0 * sqr(phi1);

	ExprPtr resultAlpha0 = -(1.0 / (2.0 * r * ipow<4>(sigma) * alpha0)) * (ipow<4>(sigma) * sqr(alpha0) -
		ipow<4>(sigma) * sqr(A0) * sqr(alpha0) -
		sqr(r) * ipow<4>(sigma) * sqr(omega) * sqr(A0) * sqr(phi0) +
		sqr(mass) * sqr(r) * ipow<4>(sigma) * sqr(A0) * sqr(alpha0) * sqr(phi0) -
		2.0 * sqr(mass) * sqr(r) * sqr(sigma) * sqr(A0) * sqr(alpha0) * ipow<4>(phi0) +
		sqr(mass) * sqr(r) * sqr(A0) * sqr(alpha0) * ipow<6>(phi0) -
		sqr(r) * ipow<4>(sigma) * sqr(alpha0) * sqr(phi1));

	ExprPtr resultPhi0 = phi1;

	ExprPtr resultPhi1 = sqr(mass) * sqr(A0) * phi0 - (sqr(omega) * sqr(A0) * phi0) / sqr(alpha0) -
		(4.0 * sqr(mass) * sqr(A0) * ipow<3>(phi0)) / sqr(sigma) +
		(3.0 * sqr(mass) * sqr(A0) * ipow<5>(phi0)) / ipow<4>(sigma) - (2.0 * phi1)/r +
		(1.0 / alpha0) * (alpha0 / r - (sqr(A0) * alpha0) / r +
		sqr(mass) * r * sqr(A0) * alpha0 * sqr(phi0) -
		(2.0 * sqr(mass) * r * sqr(A0) * alpha0 * ipow<4>(phi0)) / sqr(sigma) +
		(sqr(mass) * r * sqr(A0) * alpha0 * ipow<6>(phi0)) / ipow<4>(sigma)) * phi1;

	ExprPtr cmpA0 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0UL, resultA0, resultA0_0, resultA0_0);
	ExprPtr cmpAlpha0 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0UL, resultAlpha0, resultAlpha0_0, resultAlpha0_0);
	ExprPtr cmpPhi0 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0UL, resultPhi0, resultPhi0_0, resultPhi0_0);
	ExprPtr cmpPhi1 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0UL, resultPhi1, resultPhi1_0, resultPhi1_0);

	return { cmpA0, cmpAlpha0, cmpPhi0, cmpPhi1 };
}

bosonexprvec gen_func_solitonicboson_hp(const mpf_t omega_t, const mpf_t mass_t, const mpf_t sigma_t, const mpf_t kappa_t)
{
	ExprPtr r = new ExprVar(0);
	ExprPtr A0 = new ExprVar(1);
	ExprPtr alpha0 = new ExprVar(2);
	ExprPtr phi0 = new ExprVar(3);
	ExprPtr phi1 = new ExprVar(4);

	ExprPtr omega = new ExprDataRootT<mpf_t>(omega_t);
	ExprPtr m = new ExprDataRootT<mpf_t>(mass_t);
	ExprPtr sigma = new ExprDataRootT<mpf_t>(sigma_t);
	ExprPtr kappa = new ExprDataRootT<mpf_t>(kappa_t);

	ExprPtr resultA0_0 = 0;
	ExprPtr resultAlpha0_0 = 0;
	ExprPtr resultPhi0_0 = phi1;
	ExprPtr resultPhi1_0 = \
	(Power(A0,2)*phi0*(-(Power(sigma,4)*Power(omega,2)) + \
	Power(m,2)*Power(sigma,4)*Power(alpha0,2) - \
	4*Power(m,2)*Power(sigma,2)*Power(alpha0,2)*Power(phi0,2) + \
	3*Power(m,2)*Power(alpha0,2)*Power(phi0,4)))/(3.*Power(sigma,4)*Power(\
	alpha0,2));

	ExprPtr resultA0 = A0/(2.*r) - Power(A0,3)/(2.*r) + \
	(Power(m,2)*r*kappa*Power(A0,3)*Power(phi0,2))/2. + \
	(r*kappa*Power(omega,2)*Power(A0,3)*Power(phi0,2))/(2.*Power(alpha0,2)\
	) - (Power(m,2)*r*kappa*Power(A0,3)*Power(phi0,4))/Power(sigma,2) + \
	(Power(m,2)*r*kappa*Power(A0,3)*Power(phi0,6))/(2.*Power(sigma,4)) + \
	(r*kappa*A0*Power(phi1,2))/2.;
	ExprPtr resultAlpha0 = -(Power(sigma,4)*Power(alpha0,2) - \
	Power(sigma,4)*Power(A0,2)*Power(alpha0,2) - \
	Power(r,2)*kappa*Power(sigma,4)*Power(omega,2)*Power(A0,2)*Power(phi0,\
	2) + Power(m,2)*Power(r,2)*kappa*Power(sigma,4)*Power(A0,2)*Power(\
	alpha0,2)*Power(phi0,2) - \
	2*Power(m,2)*Power(r,2)*kappa*Power(sigma,2)*Power(A0,2)*Power(alpha0,\
	2)*Power(phi0,4) + \
	Power(m,2)*Power(r,2)*kappa*Power(A0,2)*Power(alpha0,2)*Power(phi0,6) \
	- Power(r,2)*kappa*Power(sigma,4)*Power(alpha0,2)*Power(phi1,2))/(2.*\
	r*Power(sigma,4)*alpha0);
	ExprPtr resultPhi0 = phi1;
	ExprPtr resultPhi1 = Power(m,2)*Power(A0,2)*phi0 - \
	(Power(omega,2)*Power(A0,2)*phi0)/Power(alpha0,2) - \
	(4*Power(m,2)*Power(A0,2)*Power(phi0,3))/Power(sigma,2) + \
	(3*Power(m,2)*Power(A0,2)*Power(phi0,5))/Power(sigma,4) - (2*phi1)/r \
	+ ((alpha0/r - (Power(A0,2)*alpha0)/r + \
	Power(m,2)*r*kappa*Power(A0,2)*alpha0*Power(phi0,2) - \
	(2*Power(m,2)*r*kappa*Power(A0,2)*alpha0*Power(phi0,4))/Power(sigma,2)\
	 + (Power(m,2)*r*kappa*Power(A0,2)*alpha0*Power(phi0,6))/Power(sigma,\
	4))*phi1)/alpha0;

	ExprPtr cmpA0 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0UL, resultA0, resultA0_0, resultA0_0);
	ExprPtr cmpAlpha0 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0UL, resultAlpha0, resultAlpha0_0, resultAlpha0_0);
	ExprPtr cmpPhi0 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0UL, resultPhi0, resultPhi0_0, resultPhi0_0);
	ExprPtr cmpPhi1 = new ExprVarConditional<unsigned long int, mpf_cmp_ui>(0, 0UL, resultPhi1, resultPhi1_0, resultPhi1_0);

	return { cmpA0, cmpAlpha0, cmpPhi0, cmpPhi1 };
}

size_t boson_get_zeroes_hp(std::vector<bosonvec_hp>& bosonStar, const size_t lastIndex)
{
	mpf_t tmp;
	mpf_init(tmp);
	size_t zeroes = 0;
	for (size_t i = 1; i <= lastIndex; ++i)
	{
		if (mpf_cmp_ui(bosonStar[i][VAR_PHI0], 0) == 0)
			++zeroes;
		else
		{
			mpf_mul(tmp, bosonStar[i][VAR_PHI0], bosonStar[i - 1][VAR_PHI0]);
			if (mpf_cmp_ui(tmp, 0) < 0)
				++zeroes;
		}

	}
	mpf_clear(tmp);
	return zeroes;
}

void boson_search_hp(BosonSearchParamsHP& params, mpf_t* omegaOut)
{
	auto& data = *params.result;

	mpf_t upper;
	mpf_t lower;
	mpf_t tmp;
	mpf_t bisect;

	mpf_init(upper);
	mpf_init(lower);
	mpf_init(tmp);
	mpf_init(bisect);
	mpf_set(upper, params.bracket_end);
	mpf_set(lower, params.bracket_begin);

	BosonIntegrateParamsHP integrateParams;
	init_integrate_params(&integrateParams);
	mpf_set(integrateParams.domain_start, params.domain_start);
	mpf_set(integrateParams.domain_end, params.domain_end);

	mpf_set_ui(integrateParams.initial_value[0], 1);
	mpf_set_ui(integrateParams.initial_value[1], 1);
	mpf_set(integrateParams.initial_value[2], params.initial_phi_value);
	mpf_set_ui(integrateParams.initial_value[3], 0);

	integrateParams.result = params.result;

	mpf_set(integrateParams.step_size, params.step_size);
	bosonvec thresholds = { DEFAULT_METRIC_THRESHOLD, DEFAULT_METRIC_THRESHOLD,
			DEFAULT_FIELD_THRESHOLD * mpf_get_d(params.initial_phi_value), DEFAULT_FIELD_DERIVATIVE_THRESHOLD };
	set_vec_mpf_d(&integrateParams.thresholds, thresholds);
	integrateParams.divergence_variable = VAR_PHI0;

	bool bThresholdReached = false;
	size_t lastIndex = 0;

	for (size_t iterations = 0; (iterations < params.suggested_bisection_count || bThresholdReached) &&
		(iterations < params.bisection_count_max); ++iterations)
	{
		cout << "Iteration: " << iterations + 1 << " / " << params.suggested_bisection_count << endl;
		cout << "Current bracket: [" << std::setprecision(17) << mpf_get_d(lower) << ", " << mpf_get_d(upper) << "]" << endl;

		mpf_add(bisect, upper, lower);
		mpf_div_ui(bisect, bisect, 2);

		bosonexprvec boson_expr = params.func(bisect);
		integrateParams.func = [&boson_expr](vec<BOSON_SYSTEM_DIMENSION + 1, mpf_t>& vars, bosonvec_hp* output) {
			func_boson_hp<BOSON_SYSTEM_DIMENSION>(vars, output, boson_expr);
		};

		auto result = integrate_hp(integrateParams);
		bThresholdReached = result.bThresholdReached;
		lastIndex = result.lastIndex;
		size_t zeroes = boson_get_zeroes_hp(*integrateParams.result, result.lastIndex);

		mpf_mul_ui(tmp, params.step_size, lastIndex);
		cout << "Final i = " << result.lastIndex << ", t = " << mpf_get_d(tmp) << endl;
		cout << "Variables [i - 1] = [";
		for (size_t i = 0; i < BOSON_SYSTEM_DIMENSION; ++i)
			cout << mpf_get_d(data[result.lastIndex - 1][i]) << ", ";
		cout << "]" << endl << endl;

		if (zeroes > params.requested_energy_state)
			mpf_set(upper, bisect);
		else if (zeroes < params.requested_energy_state)
			mpf_set(lower, bisect);
		else
		{
			if (result.bThresholdReached)
			{
				if (result.divergenceDirection == DIRECTION_POSITIVE)
					mpf_set(lower, bisect);
				else
					mpf_set(upper, bisect);
			}
			else
			{
				if (mpf_cmp_ui(data[result.lastIndex][VAR_PHI0], 0) >= 0)
					mpf_set(lower, bisect);
				else
					mpf_set(upper, bisect);
			}
		}

		boson_expr.release();
	}

	clear_integrate_params(&integrateParams);

	mpf_clear(lower);
	mpf_clear(upper);
	mpf_clear(tmp);
	mpf_clear(bisect);

	mpf_set(*omegaOut, bisect);

	// Clear to zero everything after final index
	for (size_t i = lastIndex + 1; i < params.result->size(); ++i)
		for (size_t j = 0; j < BOSON_SYSTEM_DIMENSION; ++j)
			mpf_set_ui((*params.result)[i][j], 0);
}

void init_boson_search_params(BosonSearchParamsHP* params)
{
	mpf_init(params->initial_phi_value);
	mpf_init(params->domain_start);
	mpf_init(params->domain_end);
	mpf_init(params->step_size);
	mpf_init(params->bracket_begin);
	mpf_init(params->bracket_end);
}

void clear_boson_search_params(BosonSearchParamsHP* params)
{
	mpf_clear(params->initial_phi_value);
	mpf_clear(params->domain_start);
	mpf_clear(params->domain_end);
	mpf_clear(params->step_size);
	mpf_clear(params->bracket_begin);
	mpf_clear(params->bracket_end);
}
