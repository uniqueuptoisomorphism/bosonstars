#include "boson.h"

using namespace std;

bosonvec func_miniboson_old(const double r, bosonvec y, const double omega)
{
	bosonvec v;

	double A0 = y[VAR_A0];
	double alpha0 = y[VAR_ALPHA0];
	double phi0 = y[VAR_PHI0];
	double phi1 = y[VAR_PHI1];

	if (r == 0)
	{
		v[VAR_A0] = 0;
		v[VAR_ALPHA0] = 0;
		v[VAR_PHI0] = phi1;
		v[VAR_PHI1] = -((sqr(omega) * sqr(A0) * phi0 - sqr(A0) * sqr(alpha0) * phi0) / (3.0 * sqr(alpha0)));
	}
	else
	{
		v[VAR_A0] = 1.0 / (2.0 * r * sqr(alpha0)) *
			A0 * (sqr(alpha0) - sqr(A0) * sqr(alpha0) +
			sqr(r) * sqr(omega) * sqr(A0) * sqr(phi0) +
			sqr(r) * sqr(A0) * sqr(alpha0) * sqr(phi0) +
			sqr(r) * sqr(alpha0) * sqr(phi1));

		v[VAR_ALPHA0] = 1.0 / (2.0 * r * alpha0) *
			(-sqr(alpha0) + sqr(A0) * sqr(alpha0) +
			sqr(r) * sqr(omega) * sqr(A0) * sqr(phi0) -
			sqr(r) * sqr(A0) * sqr(alpha0) * sqr(phi0) +
			sqr(r) * sqr(alpha0) * sqr(phi1));

		v[VAR_PHI0] = phi1;

		v[VAR_PHI1] = 1.0 / (r * sqr(alpha0)) *
			(-r * sqr(omega) * sqr(A0) * phi0 +
			r * sqr(A0) * sqr(alpha0) * phi0 - sqr(alpha0) * phi1 -
			sqr(A0) * sqr(alpha0) * phi1 +
			sqr(r) * sqr(A0) * sqr(alpha0) * sqr(phi0) * phi1);
	}

	return v;
}

bosonvec func_miniboson(const double r, bosonvec y, const double omega, const double mass)
{
	bosonvec v;

	double A0 = y[VAR_A0];
	double alpha0 = y[VAR_ALPHA0];
	double phi0 = y[VAR_PHI0];
	double phi1 = y[VAR_PHI1];

	if (r == 0)
	{
		v[VAR_A0] = 0;
		v[VAR_ALPHA0] = 0;
		v[VAR_PHI0] = phi1;
		v[VAR_PHI1] = -((sqr(omega) * sqr(A0) * phi0 - sqr(A0) * sqr(alpha0) * phi0) / (3.0 * sqr(alpha0)));
	}
	else
	{
		v[VAR_A0] = A0 / (2.0 * r) - ipow<3>(A0) / (2.0 * r) +
			1.0 / 2.0 * sqr(mass) * r * ipow<3>(A0) * sqr(phi0) +
			(r * sqr(omega) * ipow<3>(A0) * sqr(phi0)) / (2.0 * sqr(alpha0)) +
			1.0 / 2.0 * r * A0 * sqr(phi1);

		v[VAR_ALPHA0] = -(1.0 / (
			2.0 * r * alpha0)) *(sqr(alpha0) - sqr(A0) * sqr(alpha0) -
			sqr(r) * sqr(omega) * sqr(A0) * sqr(phi0) +
			sqr(mass) * sqr(r) * sqr(A0) * sqr(alpha0) * sqr(phi0) -
			sqr(r) * sqr(alpha0) * sqr(phi1));

		v[VAR_PHI0] =  phi1;

		v[VAR_PHI1] =  -(1.0 / (r * sqr(alpha0))) * (r * sqr(omega) * sqr(A0) * phi0 -
			r * sqr(A0) * sqr(alpha0) * phi0 + sqr(alpha0) * phi1 +
			sqr(A0) * sqr(alpha0) * phi1 -
			sqr(mass) * sqr(r) * sqr(A0) * sqr(alpha0) * sqr(phi0) * phi1);
	}

	return v;
}

/*bosonvec func_miniboson(const double r, bosonvec y, const double omega, const double mass,
		const double kappa)
{
	bosonvec v;

	double A0 = y[VAR_A0];
	double alpha0 = y[VAR_ALPHA0];
	double phi0 = y[VAR_PHI0];
	double phi1 = y[VAR_PHI1];

	if (r == 0)
	{
		v[VAR_A0] = 0;
		v[VAR_ALPHA0] = 0;
		v[VAR_PHI0] = phi1;
		v[VAR_PHI1] = -((sqr(omega) * sqr(A0) * phi0 - sqr(A0) * sqr(alpha0) * phi0) / (3.0 * sqr(alpha0)));
	}
	else
	{
		v[VAR_A0] = (1.0 / (2.0 * r * sqr(A0))) * A0 *
			(sqr(r) * kappa * sqr(omega) * sqr(A0) * sqr(phi0) +
			sqr(alpha0) * (1 + sqr(A0) * (-1.0 + sqr(mass) * sqr(r) * kappa * sqr(phi0)) +
			sqr(r) * kappa * sqr(phi1)));

		v[VAR_ALPHA0] = (sqr(r) * kappa * sqr(omega) * sqr(A0) * sqr(phi0) +
			sqr(alpha0) * (-1.0 + sqr(A0) * (1.0 - sqr(mass) * sqr(r) * kappa * sqr(phi0)) +
			sqr(r) * kappa * sqr(phi1))) / (2.0 * r * A0);

		v[VAR_PHI0] = phi1;

		v[VAR_PHI1] = -(phi1 / r) +
			sqr(A0) * (phi0 - (sqr(omega) * phi0) / sqr(alpha0) - phi1 / r +
			sqr(mass) * r * kappa * sqr(phi0) * phi1);
	}

	return v;
}*/

bosonvec func_miniboson(const double r, bosonvec y, const double omega, const double m,
		const double kappa)
{
	bosonvec v;

	double A0 = y[VAR_A0];
	double alpha0 = y[VAR_ALPHA0];
	double phi0 = y[VAR_PHI0];
	double phi1 = y[VAR_PHI1];

	if (r == 0)
	{
		v[VAR_A0] = 0;
		v[VAR_ALPHA0] = 0;
		v[VAR_PHI0] = phi1;
		v[VAR_PHI1] = -(Power(omega,2)*Power(A0,2)*phi0 - \
		Power(m,2)*Power(A0,2)*Power(alpha0,2)*phi0)/(3.*Power(alpha0,2));
	}
	else
	{
		v[VAR_A0] = A0/(2.*r) - Power(A0,3)/(2.*r) + \
		(Power(m,2)*r*kappa*Power(A0,3)*Power(phi0,2))/2. + \
		(r*kappa*Power(omega,2)*Power(A0,3)*Power(phi0,2))/(2.*Power(alpha0,2)\
		) + (r*kappa*A0*Power(phi1,2))/2.;
		v[VAR_ALPHA0] = -(Power(alpha0,2) - Power(A0,2)*Power(alpha0,2) - \
		Power(r,2)*kappa*Power(omega,2)*Power(A0,2)*Power(phi0,2) + \
		Power(m,2)*Power(r,2)*kappa*Power(A0,2)*Power(alpha0,2)*Power(phi0,2) \
		- Power(r,2)*kappa*Power(alpha0,2)*Power(phi1,2))/(2.*r*alpha0);
		v[VAR_PHI0] = phi1;
		v[VAR_PHI1] = -((r*Power(omega,2)*Power(A0,2)*phi0 - \
		Power(m,2)*r*Power(A0,2)*Power(alpha0,2)*phi0 + Power(alpha0,2)*phi1 \
		+ Power(A0,2)*Power(alpha0,2)*phi1 - \
		Power(m,2)*Power(r,2)*kappa*Power(A0,2)*Power(alpha0,2)*Power(phi0,2)*\
		phi1)/(r*Power(alpha0,2)));
	}

	return v;
}

bosonvec func_massiveboson(const double r, bosonvec y, const double omega,
		const double mass, const double lambda)
{
	bosonvec v;

	double A0 = y[VAR_A0];
	double alpha0 = y[VAR_ALPHA0];
	double phi0 = y[VAR_PHI0];
	double phi0d = y[VAR_PHI1];

	if (r == 0)
	{
		v[VAR_A0] = 0;
		v[VAR_ALPHA0] = 0;
		v[VAR_PHI0] = phi0d;
		v[VAR_PHI1] =  sqr(A0) * phi0 * (- sqr(omega) + sqr(mass) * sqr(alpha0) + lambda * sqr(alpha0) * sqr(phi0)) / (3.0 * sqr(alpha0));
	}
	else
	{
		v[VAR_A0] = -(1.0/4.0) * A0 * (-(2.0/r) + sqr(A0) * (2.0/r +
			r * sqr(phi0) * (-2.0 * sqr(mass) -
			(2.0 * sqr(omega)) / sqr(alpha0) - lambda * sqr(phi0))) -
			2.0 * r * sqr(phi0d));

		v[VAR_ALPHA0] = (1.0/(4.0 * r * alpha0)) *
			(2.0 * sqr(r) * sqr(omega) * sqr(A0) * sqr(phi0)
			- sqr(alpha0) * (2.0 + sqr(A0) * (-2.0 +
			sqr(r) * sqr(phi0) * (2.0 * sqr(mass) + lambda * sqr(phi0))) -
			2.0 * sqr(r) * sqr(phi0d)));

		v[VAR_PHI0] = phi0d;

		v[VAR_PHI1] = -(phi0d / r) +
			sqr(A0) * ((sqr(mass) - sqr(omega) / sqr(alpha0)) * phi0 + lambda * sqr(phi0) * phi0 - phi0d / r +
			sqr(mass) * r * sqr(phi0) * phi0d +
			1.0 / 2.0 * r * lambda * sqr(phi0) * sqr(phi0) * phi0d);
	}

	return v;
}

bosonvec func_massiveboson(const double r, bosonvec y, const double omega,
		const double m, const double lambda, const double kappa)
{
	bosonvec v;

	double A0 = y[VAR_A0];
	double alpha0 = y[VAR_ALPHA0];
	double phi0 = y[VAR_PHI0];
	double phi1 = y[VAR_PHI1];

	if (r == 0)
	{
		v[VAR_A0] = 0;
		v[VAR_ALPHA0] = 0;
		v[VAR_PHI0] = phi1;
		v[VAR_PHI1] = -(Power(omega,2)*Power(A0,2)*phi0 - \
		Power(m,2)*Power(A0,2)*Power(alpha0,2)*phi0 - \
		lambda*Power(A0,2)*Power(alpha0,2)*Power(phi0,3))/(3.*Power(alpha0,2));
	}
	else
	{
		v[VAR_A0] = A0/(2.*r) - Power(A0,3)/(2.*r) + \
		(Power(m,2)*r*kappa*Power(A0,3)*Power(phi0,2))/2. + \
		(r*kappa*Power(omega,2)*Power(A0,3)*Power(phi0,2))/(2.*Power(alpha0,2)\
		) + (r*kappa*lambda*Power(A0,3)*Power(phi0,4))/4. + \
		(r*kappa*A0*Power(phi1,2))/2.;
		v[VAR_ALPHA0] = -(2*Power(alpha0,2) - 2*Power(A0,2)*Power(alpha0,2) - \
		2*Power(r,2)*kappa*Power(omega,2)*Power(A0,2)*Power(phi0,2) + \
		2*Power(m,2)*Power(r,2)*kappa*Power(A0,2)*Power(alpha0,2)*Power(phi0,\
		2) + Power(r,2)*kappa*lambda*Power(A0,2)*Power(alpha0,2)*Power(phi0,4)\
		 - 2*Power(r,2)*kappa*Power(alpha0,2)*Power(phi1,2))/(4.*r*alpha0);
		v[VAR_PHI0] = phi1;
		v[VAR_PHI1] = -(2*r*Power(omega,2)*Power(A0,2)*phi0 - \
		2*Power(m,2)*r*Power(A0,2)*Power(alpha0,2)*phi0 - \
		2*r*lambda*Power(A0,2)*Power(alpha0,2)*Power(phi0,3) + \
		2*Power(alpha0,2)*phi1 + 2*Power(A0,2)*Power(alpha0,2)*phi1 - \
		2*Power(m,2)*Power(r,2)*kappa*Power(A0,2)*Power(alpha0,2)*Power(phi0,\
		2)*phi1 - \
		Power(r,2)*kappa*lambda*Power(A0,2)*Power(alpha0,2)*Power(phi0,4)*\
		phi1)/(2.*r*Power(alpha0,2));
	}


	return v;
}

bosonvec func_solitonicboson(const double r, bosonvec y, const double omega,
		const double mass, const double sigma)
{
	bosonvec v;

	double A0 = y[VAR_A0];
	double alpha0 = y[VAR_ALPHA0];
	double phi0 = y[VAR_PHI0];
	double phi1 = y[VAR_PHI1];

	if (r == 0.0)
	{
		v[VAR_A0] = 0;
		v[VAR_ALPHA0] = 0;
		v[VAR_PHI0] = phi1;
		v[VAR_PHI1] = -(1.0 / ( 3.0 * ipow<4>(sigma) * sqr(alpha0))) * (ipow<4>(sigma) * sqr(omega) * sqr(A0) * phi0 -
			sqr(mass) * ipow<4>(sigma) * sqr(A0) * sqr(alpha0) * phi0 +
			4.0 * sqr(mass) * sqr(sigma) * sqr(A0) * sqr(alpha0) * ipow<3>(phi0) -
			3.0 * sqr(mass) * sqr(A0) * sqr(alpha0) * ipow<5>(phi0));
	}
	else
	{
		v[VAR_A0] = A0 /(2.0 * r) - ipow<3>(A0)/(2.0 * r) +
			1.0 / 2.0 * sqr(mass) * r * ipow<3>(A0) *sqr(phi0) +
			(r * sqr(omega) * ipow<3>(A0) * sqr(phi0))/(2 * sqr(alpha0)) -
			(sqr(mass) * r * ipow<3>(A0) * ipow<4>(phi0))/sqr(sigma) +
			(sqr(mass) * r * ipow<3>(A0) * ipow<6>(phi0))/(2.0 * ipow<4>(sigma)) +
			1.0 / 2.0 * r * A0 * sqr(phi1);

		v[VAR_ALPHA0] = -(1.0 / (2.0 * r * ipow<4>(sigma) * alpha0)) * (ipow<4>(sigma) * sqr(alpha0) -
			ipow<4>(sigma) * sqr(A0) * sqr(alpha0) -
			sqr(r) * ipow<4>(sigma) * sqr(omega) * sqr(A0) * sqr(phi0) +
			sqr(mass) * sqr(r) * ipow<4>(sigma) * sqr(A0) * sqr(alpha0) * sqr(phi0) -
			2.0 * sqr(mass) * sqr(r) * sqr(sigma) * sqr(A0) * sqr(alpha0) * ipow<4>(phi0) +
			sqr(mass) * sqr(r) * sqr(A0) * sqr(alpha0) * ipow<6>(phi0) -
			sqr(r) * ipow<4>(sigma) * sqr(alpha0) * sqr(phi1));

		v[VAR_PHI0] = phi1;

		v[VAR_PHI1] = sqr(mass) * sqr(A0) * phi0 - (sqr(omega) * sqr(A0) * phi0) / sqr(alpha0) -
			(4.0 * sqr(mass) * sqr(A0) * ipow<3>(phi0)) / sqr(sigma) +
			(3.0 * sqr(mass) * sqr(A0) * ipow<5>(phi0)) / ipow<4>(sigma) - (2.0 * phi1)/r +
			(1.0 / alpha0) * (alpha0 / r - (sqr(A0) * alpha0) / r +
			sqr(mass) * r * sqr(A0) * alpha0 * sqr(phi0) -
			(2.0 * sqr(mass) * r * sqr(A0) * alpha0 * ipow<4>(phi0)) / sqr(sigma) +
			(sqr(mass) * r * sqr(A0) * alpha0 * ipow<6>(phi0)) / ipow<4>(sigma)) * phi1;
	}

	return v;
}

bosonvec func_solitonicboson(const double r, bosonvec y, const double omega,
		const double m, const double sigma, const double kappa)
{
	bosonvec v;

	double A0 = y[VAR_A0];
	double alpha0 = y[VAR_ALPHA0];
	double phi0 = y[VAR_PHI0];
	double phi1 = y[VAR_PHI1];

	if (r == 0)
	{
		v[VAR_A0] = 0;
		v[VAR_ALPHA0] = 0;
		v[VAR_PHI0] = phi1;
		v[VAR_PHI1] = (Power(A0,2)*phi0*(-(Power(sigma,4)*Power(omega,2)) + \
		Power(m,2)*Power(sigma,4)*Power(alpha0,2) - \
		4*Power(m,2)*Power(sigma,2)*Power(alpha0,2)*Power(phi0,2) + \
		3*Power(m,2)*Power(alpha0,2)*Power(phi0,4)))/(3.*Power(sigma,4)*Power(\
		alpha0,2));
	}
	else
	{
		v[VAR_A0] = A0/(2.*r) - Power(A0,3)/(2.*r) + \
		(Power(m,2)*r*kappa*Power(A0,3)*Power(phi0,2))/2. + \
		(r*kappa*Power(omega,2)*Power(A0,3)*Power(phi0,2))/(2.*Power(alpha0,2)\
		) - (Power(m,2)*r*kappa*Power(A0,3)*Power(phi0,4))/Power(sigma,2) + \
		(Power(m,2)*r*kappa*Power(A0,3)*Power(phi0,6))/(2.*Power(sigma,4)) + \
		(r*kappa*A0*Power(phi1,2))/2.;
		v[VAR_ALPHA0] = -(Power(sigma,4)*Power(alpha0,2) - \
		Power(sigma,4)*Power(A0,2)*Power(alpha0,2) - \
		Power(r,2)*kappa*Power(sigma,4)*Power(omega,2)*Power(A0,2)*Power(phi0,\
		2) + Power(m,2)*Power(r,2)*kappa*Power(sigma,4)*Power(A0,2)*Power(\
		alpha0,2)*Power(phi0,2) - \
		2*Power(m,2)*Power(r,2)*kappa*Power(sigma,2)*Power(A0,2)*Power(alpha0,\
		2)*Power(phi0,4) + \
		Power(m,2)*Power(r,2)*kappa*Power(A0,2)*Power(alpha0,2)*Power(phi0,6) \
		- Power(r,2)*kappa*Power(sigma,4)*Power(alpha0,2)*Power(phi1,2))/(2.*\
		r*Power(sigma,4)*alpha0);
		v[VAR_PHI0] = phi1;
		v[VAR_PHI1] = Power(m,2)*Power(A0,2)*phi0 - \
		(Power(omega,2)*Power(A0,2)*phi0)/Power(alpha0,2) - \
		(4*Power(m,2)*Power(A0,2)*Power(phi0,3))/Power(sigma,2) + \
		(3*Power(m,2)*Power(A0,2)*Power(phi0,5))/Power(sigma,4) - (2*phi1)/r \
		+ ((alpha0/r - (Power(A0,2)*alpha0)/r + \
		Power(m,2)*r*kappa*Power(A0,2)*alpha0*Power(phi0,2) - \
		(2*Power(m,2)*r*kappa*Power(A0,2)*alpha0*Power(phi0,4))/Power(sigma,2)\
		 + (Power(m,2)*r*kappa*Power(A0,2)*alpha0*Power(phi0,6))/Power(sigma,\
		4))*phi1)/alpha0;
	}

	return v;
}

size_t boson_get_zeroes(vector<bosonvec>& bosonStar, const size_t lastIndex)
{
	size_t zeroes = 0;
	for (size_t i = 1; i <= lastIndex; ++i)
	{
		if (bosonStar[i][VAR_PHI0] == 0.0)
			++zeroes;
		else if (bosonStar[i][VAR_PHI0] * bosonStar[i - 1][VAR_PHI0] < 0)
			++zeroes;
	}
	return zeroes;
}

double boson_search(BosonSearchParams& params)
{
	auto& data = *params.result;

	double upper = params.bracket_end;
	double lower = params.bracket_begin;

	BosonIntegrateParams integrateParams;
	integrateParams.domain_start = params.domain_start;
	integrateParams.domain_end = params.domain_end;
	integrateParams.initial_value = { 1.0, 1.0, params.initial_phi_value, 0.0 };
	integrateParams.result = params.result;
	integrateParams.step_size = params.step_size;
	integrateParams.thresholds = mult<BOSON_SYSTEM_DIMENSION, double>(params.initial_phi_value,
			{ DEFAULT_METRIC_THRESHOLD, DEFAULT_METRIC_THRESHOLD,
			DEFAULT_FIELD_THRESHOLD, DEFAULT_FIELD_DERIVATIVE_THRESHOLD });
	integrateParams.divergence_variable = VAR_PHI0;

	double bisect = 0.0;

	do
	{
		cout << "Current bracket: [" << std::setprecision(17) << lower << ", " << upper << "]" << endl;

		bisect = (upper + lower) / 2.0;
		auto base_func = params.func;
		integrateParams.func = [&bisect, &base_func](const double r, bosonvec y) {
			return base_func(r, y, bisect);
		};

		auto result = integrate(integrateParams);
		size_t zeroes = boson_get_zeroes(*integrateParams.result, result.lastIndex);

		cout << "Final i = " << result.lastIndex << endl;
		cout << "phi[i - 1] = " << data[result.lastIndex - 1][VAR_PHI0] << endl;

		if (zeroes > params.requested_energy_state)
			upper = bisect;
		else if (zeroes < params.requested_energy_state)
			lower = bisect;
		else
		{
			if (result.bThresholdReached)
			{
				if (result.divergenceDirection == DIRECTION_POSITIVE)
					lower = bisect;
				else
					upper = bisect;
			}
			else
			{
				if (data[result.lastIndex][VAR_PHI0] >= 0.0)
					lower = bisect;
				else
					upper = bisect;
			}
		}
	}
	while (upper - lower > params.requested_final_bracket_size);

	return bisect;
}

