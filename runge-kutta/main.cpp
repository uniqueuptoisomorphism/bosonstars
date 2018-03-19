#include "boson.h"
#include "expr.h"
#include "boson_hp.h"
#include "custom_boson_hp.h"
#include <bbhutil.h>
#include <string>
#include <fstream>
#include <array>
#include <algorithm>

#define RADIUS_MASS_PERCENTAGE 0.99

using namespace std;

char out_name[][50] =
{
	"A0", "alpha0", "phi0", "phi1"
};

const string mass_v_radius_output_file = "mvr.dat";

string out_name_files[] =
{
	"A0.dat", "alpha0.dat", "phi0.dat", "phi1.dat"
};

string out_name_files_bin[] =
{
	"A0.bin", "alpha0.bin", "phi0.bin", "phi1.bin"
};

double mass_value = DEFAULT_MASS;
double sigma_value = DEFAULT_SIGMA;
double lambda_value = DEFAULT_LAMBDA;
double kappa_value = DEFAULT_KAPPA;
double bracket_begin_value = 0.2;
double bracket_end_value = 20.0;
double step_size_value = 0.01;
double domain_end_value = 40.0;
double init_phi_value = 0.1;
double radius_v_mass_phi_start = 0.05;
double radius_v_mass_phi_end = 0.06;
double radius_v_mass_phi_step = 0.001;
double radius_v_mass_adaptive_threshold = 0.1;
double find_transition_init_phi_lower = 0.053;
double find_transition_init_phi_upper = 0.054;
double find_transition_potential_peak_location = 7.0;
double omega_value = 1.0;
int find_transition_potential_peak_iterations = 10;
int program_mode = 0;
int output_mode = 0;
int requested_energy_state_value = 0;
int suggested_bisection_count_value = 60;
int bisection_count_max_value = 300;
BosonStarType boson_star_type = BOSON_STAR_MINI;
int use_high_precision = 1;
int precision = 512;
double manual_lower_bound = 0.0;
double manual_upper_bound = 0.0;
int search_variable = 0;
double fix_variable_value = 0;
int b_inverted_search = 0;

enum ProgramMode
{
	PROGRAM_MODE_DEFAULT = 0,
	PROGRAM_MODE_MASS_V_RADIUS = 1,
	PROGRAM_MODE_FIND_TRANSITION = 2,
	PROGRAM_MODE_MASS_V_RADIUS_ADAPTIVE = 3,
	PROGRAM_MODE_SPECIFY_OMEGA = 4,
	PROGRAM_MODE_MANUAL_SEARCH = 5
};

enum OutputMode
{
	OUTPUT_MODE_SDF = 0,
	OUTPUT_MODE_RAW_DATA_TXT = 1,
	OUTPUT_MODE_RAW_DATA_BIN = 2
};

typedef array<double, 2> point;

template <size_t size, typename T>
void clone(vector<vec<size, T>>& src, const size_t index, double out[])
{
	size_t length = src.size();
	for (size_t i = 0; i < length; ++i)
		out[i] = src[i][index];
}

template <size_t size, typename T>
void output_result_sdf(vector<vec<size, T>>& result, const double step_size)
{
	// Output result
	size_t step_count = result.size();

	double bbox[] = { 0.0, (step_count - 1) * step_size };
	double* u = new double[step_count];
	for (size_t i = 0; i < size; ++i)
	{
		clone<size>(result, i, u);
		int count = step_count;
		gft_out_bbox(out_name[i], 0.0, &count, 1, bbox, u);
	}
	delete [] u;
}

template <size_t size, typename T>
void output_result_mathematica(vector<vec<size, T>>& result, const double step_size)
{
	for (size_t i = 0; i < size; ++i)
	{
		ofstream file;
		file.open(out_name_files[i]);
		for (auto& p : result)
			file << p[i] << " ";
		file.close();
	}
}

template <size_t size, typename T>
void output_result_binary(vector<vec<size, T>>& result, const double step_size)
{
	double* buf = new double[result.size()];
	for (size_t i = 0; i < size; ++i)
	{
		ofstream file(out_name_files_bin[i], ios::out | ios::binary);
		for (size_t j = 0, len = result.size(); j < len; ++j)
			buf[j] = result[j][i];
		file.write((char*)buf, result.size() * sizeof(double));
		file.close();
	}
	delete [] buf;
}

template <size_t size, typename T>
void output_result(vector<vec<size, T>>& result, const double step_size)
{
	if (output_mode == OUTPUT_MODE_RAW_DATA_BIN)
		output_result_binary(result, step_size);
	else if (output_mode == OUTPUT_MODE_RAW_DATA_TXT)
		output_result_mathematica(result, step_size);
	else
		output_result_sdf(result, step_size);
}

struct BosonStarProperties
{
	double TotalMass;
	double Radius;
	double Compactness;
};

template <typename T>
void analyze_properties(vector<vec<BOSON_SYSTEM_DIMENSION, T>>& result, const double step_size,
	BosonStarProperties* properties_out)
{
	size_t step_count = result.size();

	double* u = new double[step_count];
	clone<BOSON_SYSTEM_DIMENSION>(result, VAR_A0, u);
	for (size_t i = 0; i < step_count; ++i)
	{
		double r = i * step_size;
		u[i] = - r * (1.0 / (u[i] * u[i]) - 1) / 2.0;
	}

	double total_mass = u[step_count - 1];
	size_t radius_int = 0;
	for (; radius_int < step_count; ++radius_int)
		if (u[radius_int] >= total_mass * RADIUS_MASS_PERCENTAGE)
			break;
	double radius = radius_int * step_size;
	cout << "Radius: " << radius << endl;
	cout << "Compactness (R/M) : " << radius / total_mass << endl << endl;

	properties_out->TotalMass = total_mass;
	properties_out->Radius = radius;
	properties_out->Compactness = radius / total_mass;

	delete [] u;
}

template <typename T>
void output_mass(vector<vec<BOSON_SYSTEM_DIMENSION, T>>& result, const double step_size)
{
	// Output result
	size_t step_count = result.size();

	// Compute mass
	double bbox[] = { 0.0, (step_count - 1) * step_size };
	double* u = new double[step_count];
	clone<BOSON_SYSTEM_DIMENSION>(result, VAR_A0, u);
	for (size_t i = 0; i < step_count; ++i)
	{
		double r = i * step_size;
		u[i] = - r * (1.0 / (u[i] * u[i]) - 1) / 2.0;
	}

	int count = step_count;
	char mass_name[] = "mass";
	gft_out_bbox(mass_name, 0.0, &count, 1, bbox, u);

	// Compute radius
	cout << endl;
	double total_mass = u[step_count - 1];
	cout << "Total mass: " << total_mass << endl;
	size_t radius_int = 0;
	for (; radius_int < step_count; ++radius_int)
		if (u[radius_int] >= total_mass * RADIUS_MASS_PERCENTAGE)
			break;
	double radius = radius_int * step_size;
	cout << "Radius: " << radius << endl;
	cout << "Compactness (R/M) : " << radius / total_mass << endl << endl;

	// Compute mass density
	double* d = new double[step_count];
	for (size_t i = 1; i < step_count - 1; ++i)
		d[i] = (u[i + 1] - u[i - 1]) / (2.0 * step_size);
	d[0] = (u[1] - u[0]) / step_size;
	d[step_count - 1] = (u[step_count - 1] - u[step_count - 2]) / step_size;

	char mass_density_name[] = "mass_density";
	gft_out_bbox(mass_density_name, 0.0, &count, 1, bbox, d);

	double* potential = new double[step_count];
	double* phi = new double[step_count];
	clone<BOSON_SYSTEM_DIMENSION>(result, VAR_PHI0, phi);

	for (size_t i = 0; i < step_count; ++i)
		potential[i] = sqr(mass_value) * sqr(phi[i]) * sqr(1.0 - sqr(phi[i]) / sqr(sigma_value));

	char potential_name[] = "potential";
	gft_out_bbox(potential_name, 0.0, &count, 1, bbox, potential);

	double* potential_r = new double[step_count];
	for (size_t i = 0; i < step_count; ++i)
		potential_r[i] = phi[i] * ((sqr(mass_value) * (-3.0 * sqr(phi[i]) + sqr(sigma_value)) * (- sqr(phi[i]) + sqr(sigma_value))) / pow(sigma_value, 4));

	char potential_r_name[] = "potential_r";
	gft_out_bbox(potential_r_name, 0.0, &count, 1, bbox, potential_r);

	delete [] potential_r;
	delete [] phi;
	delete [] potential;

	delete [] d;
	delete [] u;
}

template <typename T>
void compute_potential(vector<vec<BOSON_SYSTEM_DIMENSION, T>>& result, vector<double>* output)
{
	size_t step_count = result.size();
	auto& out = *output;
	out.resize(result.size());
	clone<BOSON_SYSTEM_DIMENSION>(result, VAR_PHI0, out.data());

	for (size_t i = 0; i < step_count; ++i)
		out[i] = sqr(mass_value) * sqr(out[i]) * sqr(1.0 - sqr(out[i]) / sqr(sigma_value));
}

void bisection_search(BosonStarType type, BosonSearchParams& params)
{
	switch (type)
	{
	case BOSON_STAR_MINI:
		params.func = [](const double r, bosonvec y, const double omega) {
			return func_miniboson(r, y, omega, mass_value, kappa_value);
		};
		break;
	case BOSON_STAR_MASSIVE:
		params.func = [](const double r, bosonvec y, const double omega) {
			return func_massiveboson(r, y, omega, mass_value, lambda_value, kappa_value);
		};
		break;
	case BOSON_STAR_SOLITONIC:
		params.func = [](const double r, bosonvec y, const double omega) {
			return func_solitonicboson(r, y, omega, mass_value, sigma_value, kappa_value);
		};
		break;
	default:
		params.func = [](const double r, bosonvec y, const double omega) {
			return func_miniboson(r, y, omega, mass_value, kappa_value);
		};
		break;
	}

	double omega = boson_search(params);
	cout << "Final omega = " << omega << endl;
}

void bisection_search(BosonStarType type)
{
	vector<bosonvec> result;
	BosonSearchParams params;
	params.bracket_begin = bracket_begin_value;
	params.bracket_end = bracket_end_value;
	params.domain_start = 0;
	params.domain_end = domain_end_value;
	params.initial_phi_value = init_phi_value;
	params.requested_energy_state = requested_energy_state_value;
	params.requested_final_bracket_size = 2.0E-15;
	params.step_size = step_size_value;
	params.result = &result;

	bisection_search(type, params);

	output_result(result, params.step_size);
	output_mass(result, params.step_size);
}

template <size_t size, typename T>
void convert_to_double(vector<vec<size, T>>* dest, vector<vec<size, mpf_t>>& src)
{
	dest->resize(src.size());
	for (size_t i = 0; i < src.size(); ++i)
		for (size_t j = 0; j < size; ++j)
			(*dest)[i][j] = mpf_get_d(src[i][j]);
}

void bisection_search_hp(BosonStarType type, BosonSearchParamsHP& params)
{
	mpf_t mass;
	mpf_t omega;
	mpf_t lambda;
	mpf_t sigma;
	mpf_t kappa;
	mpf_init(mass);
	mpf_init(omega);
	mpf_init(sigma);
	mpf_init(lambda);
	mpf_init(kappa);
	mpf_set_d(mass, mass_value);
	mpf_set_d(lambda, lambda_value);
	mpf_set_d(sigma, sigma_value);
	mpf_set_d(kappa, kappa_value);

	switch (type)
	{
	case BOSON_STAR_MINI:
		params.func = [&mass, &kappa](const mpf_t omega) {
			return gen_func_miniboson_hp(omega, mass, kappa);
		};
		break;

	case BOSON_STAR_MASSIVE:
		params.func = [&mass, &lambda, &kappa](const mpf_t omega) {
			return gen_func_massiveboson_hp(omega, mass, lambda, kappa);
		};
		break;

	case BOSON_STAR_SOLITONIC:
		params.func = [&mass, &sigma, &kappa](const mpf_t omega) {
			return gen_func_solitonicboson_hp(omega, mass, sigma, kappa);
		};
		break;

	default:
		params.func = [&mass, &kappa](const mpf_t omega) {
			return gen_func_miniboson_hp(omega, mass, kappa);
		};
		break;
	}

	boson_search_hp(params, &omega);
	cout << "Final omega = " << mpf_get_d(omega) << endl;

	mpf_clear(omega);
	mpf_clear(mass);
	mpf_clear(lambda);
	mpf_clear(sigma);
}

void bisection_search_hp(BosonStarType type, double initial_phi_value, vector<bosonvec_hp>& result)
{
	BosonSearchParamsHP params;
	init_boson_search_params(&params);
	mpf_set_d(params.bracket_begin, bracket_begin_value);
	mpf_set_d(params.bracket_end, bracket_end_value);
	mpf_set_ui(params.domain_start, 0);
	mpf_set_ui(params.domain_end, domain_end_value);
	mpf_set_d(params.initial_phi_value, initial_phi_value);
	params.requested_energy_state = requested_energy_state_value;
	params.suggested_bisection_count = suggested_bisection_count_value;
	params.bisection_count_max = bisection_count_max_value;
	mpf_set_d(params.step_size, step_size_value);
	params.result = &result;

	bisection_search_hp(type, params);

	clear_boson_search_params(&params);
}

void bisection_search_hp(BosonStarType type)
{
	vector<bosonvec_hp> result;

	bisection_search_hp(type, init_phi_value, result);

	vector<bosonvec> result_sp;
	convert_to_double(&result_sp, result);
	output_result(result_sp, step_size_value);
	output_mass(result_sp, step_size_value);
}

// Function to load all the program parameters from a parameter file (argv[1])
void loadParameters(int argc, const char* argv[])
{
	if (argc < 2)
	{
		cout << "Missing parameter file - defaults will be used" << endl;
	}
	else
	{
		auto pFile = argv[1];
		get_real_param(pFile, "lambda", &lambda_value, 1);
		get_real_param(pFile, "mass", &mass_value, 1);
		get_real_param(pFile, "sigma", &sigma_value, 1);
		get_real_param(pFile, "kappa", &kappa_value, 1);
		get_real_param(pFile, "bracket_begin", &bracket_begin_value, 1);
		get_real_param(pFile, "bracket_end", &bracket_end_value, 1);
		get_real_param(pFile, "step_size", &step_size_value, 1);
		get_real_param(pFile, "domain_end", &domain_end_value, 1);
		get_real_param(pFile, "init_phi", &init_phi_value, 1);
		get_int_param(pFile, "requested_energy_state", &requested_energy_state_value, 1);
		get_int_param(pFile, "suggested_bisection_count", &suggested_bisection_count_value, 1);
		get_int_param(pFile, "bisection_count_max", &bisection_count_max_value, 1);

		int temp;
		get_int_param(pFile, "boson_star_type", &temp, 1);
		boson_star_type = (BosonStarType)temp;

		get_int_param(pFile, "use_high_precision", &use_high_precision, 1);
		get_int_param(pFile, "precision", &precision, 1);
		get_int_param(pFile, "program_mode", &program_mode, 1);
		get_int_param(pFile, "output_mode", &output_mode, 1);
		get_real_param(pFile, "radius_v_mass_phi_start", &radius_v_mass_phi_start, 1);
		get_real_param(pFile, "radius_v_mass_phi_end", &radius_v_mass_phi_end, 1);
		get_real_param(pFile, "radius_v_mass_phi_step", &radius_v_mass_phi_step, 1);
		get_real_param(pFile, "radius_v_mass_adaptive_threshold", &radius_v_mass_adaptive_threshold, 1);
		get_real_param(pFile, "find_transition_init_phi_lower", &find_transition_init_phi_lower, 1);
		get_real_param(pFile, "find_transition_init_phi_upper", &find_transition_init_phi_upper, 1);
		get_real_param(pFile, "find_transition_potential_peak_location", &find_transition_potential_peak_location, 1);
		get_int_param(pFile, "find_transition_potential_peak_iterations", &find_transition_potential_peak_iterations, 1);
		get_real_param(pFile, "omega", &omega_value, 1);

		get_real_param(pFile, "manual_lower_bound", &manual_lower_bound, 1);
		get_real_param(pFile, "manual_upper_bound", &manual_upper_bound, 1);
		get_int_param(pFile, "search_variable", &search_variable, 1);
		get_real_param(pFile, "fix_variable_value", &fix_variable_value, 1);
		get_int_param(pFile, "b_inverted_search", &b_inverted_search, 1);
	}

	cout << "lambda = " << lambda_value << endl;
	cout << "mass = " << mass_value << endl;
	cout << "sigma = " << sigma_value << endl;
	cout << "kappa = " << kappa_value << endl;
	cout << "bracket_begin = " << bracket_begin_value << endl;
	cout << "bracket_end = " << bracket_end_value << endl;
	cout << "step_size = " << step_size_value << endl;
	cout << "domain_end = " << domain_end_value << endl;
	cout << "init_phi = " << init_phi_value << endl;
	cout << "requested_energy_state = " << requested_energy_state_value << endl;
	cout << "suggested_bisection_count = " << suggested_bisection_count_value << endl;
	cout << "bisection_count_max = " << bisection_count_max_value << endl;
	cout << "boson_star_type = " << boson_star_type << endl;
	cout << "use_high_precision = " << use_high_precision << endl;
	cout << "precision = " << precision << endl;
	cout << "program_mode = " << program_mode << endl;
	cout << "output_mode = " << output_mode << endl;

	switch (program_mode)
	{
	case PROGRAM_MODE_MASS_V_RADIUS:
		cout << "radius_v_mass_phi_start = " << radius_v_mass_phi_start << endl;
		cout << "radius_v_mass_phi_end = " << radius_v_mass_phi_end << endl;
		cout << "radius_v_mass_phi_step = " << radius_v_mass_phi_end << endl;
		cout << "radius_v_mass_adaptive_threshold = " << radius_v_mass_adaptive_threshold << endl;
		break;
	case PROGRAM_MODE_FIND_TRANSITION:
		cout << "find_transition_init_phi_lower = " << find_transition_init_phi_lower << endl;
		cout << "find_transition_init_phi_upper = " << find_transition_init_phi_upper << endl;
		cout << "find_transition_potential_peak_location = " << find_transition_potential_peak_location << endl;
		cout << "find_transition_potential_peak_iterations = " << find_transition_potential_peak_iterations << endl;
		break;
	case PROGRAM_MODE_SPECIFY_OMEGA:
		cout << "omega = " << omega_value << endl;
		break;
	case PROGRAM_MODE_MANUAL_SEARCH:
		cout << "manual_lower_bound = " << manual_lower_bound << endl;
		cout << "manual_upper_bound = " << manual_upper_bound << endl;
		cout << "search_variable = " << search_variable << endl;
		cout << "fix_variable_value = " << fix_variable_value << endl;
		cout << "b_inverted_search = " << b_inverted_search << endl;
		break;
	}

	cout << endl;

	mpf_set_default_prec(precision);
	cout << "Precision: " << mpf_get_default_prec() << endl;
}

void output_points_mathematica_data(vector<point>& points, const string& output_file)
{
	ofstream file;
	file.open(output_file);
	for (auto p : points)
		file << p[0] << " " << p[1] << endl;
	file.close();
}

void mass_v_radius_hp(BosonStarType type)
{
	vector<point> points;

	for (double current_phi = radius_v_mass_phi_start;
			current_phi <= radius_v_mass_phi_end;
			current_phi += radius_v_mass_phi_step)
	{
		vector<bosonvec_hp> result;

		cout << "START PHI = " << current_phi << endl << endl;

		bisection_search_hp(type, current_phi, result);

		vector<bosonvec> result_sp;
		convert_to_double(&result_sp, result);

		BosonStarProperties props;
		analyze_properties(result_sp, step_size_value, &props);
		point p = { props.Radius, props.TotalMass };
		points.push_back(p);
	}

	cout << "Outputing result to file..." << endl;
	output_points_mathematica_data(points, mass_v_radius_output_file);
}

void critical_transition_find_hp()
{
	double upper = find_transition_init_phi_upper;
	double lower = find_transition_init_phi_lower;

	vector<double> potential;
	vector<bosonvec> result_sp;

	int comparison_index = static_cast<size_t>(find_transition_potential_peak_location / step_size_value);

	for (int i = 0; i < find_transition_potential_peak_iterations; ++i)
	{
		vector<bosonvec_hp> result;

		double value = (upper + lower) / 2.0;
		cout << "ITERATION = " << i + 1 << endl;
		cout << "CURRENT BRACKET: [" << lower << ", " << upper << "]" << endl << endl;

		bisection_search_hp(BOSON_STAR_SOLITONIC, value, result);
		convert_to_double(&result_sp, result);
		compute_potential(result_sp, &potential);

		// Find index of maximum of potential
		int maximum_index = distance(potential.begin(), max_element(potential.begin(), potential.end()));

		if (maximum_index > comparison_index)
			upper = value;
		else
			lower = value;
	}

	cout << "CURRENT BRACKET: [" << lower << ", " << upper << "]" << endl << endl;

	// Output result
	cout << "Outputing result to file... " << endl;
	output_result(result_sp, step_size_value);
	output_mass(result_sp, step_size_value);
}

void run_in_specify_frequency_mode()
{
	// Set input parameters
	mpf_t omega;
	mpf_t mass;
	mpf_t lambda;
	mpf_t sigma;
	mpf_t kappa;
	mpf_init(mass);
	mpf_init(omega);
	mpf_init(sigma);
	mpf_init(lambda);
	mpf_init(kappa);
	mpf_init(omega);
	mpf_init(mass);
	mpf_set_d(omega, omega_value);
	mpf_set_d(mass, mass_value);
	mpf_set_d(mass, mass_value);
	mpf_set_d(lambda, lambda_value);
	mpf_set_d(sigma, sigma_value);
	mpf_set_d(kappa, kappa_value);

	// Set the initial values for the four variables { A0, alpha0, phi0, phi0' }
	bosonvec initialvalues = { 1.0, 1.0, init_phi_value, 0.0 };
	// Set the thresholds for the four variables (integration stops when these are reached)
	bosonvec thresholds = { DEFAULT_METRIC_THRESHOLD, DEFAULT_METRIC_THRESHOLD,
			DEFAULT_FIELD_THRESHOLD * init_phi_value , DEFAULT_FIELD_DERIVATIVE_THRESHOLD };

	// Choose the correct generator for the derivative function used in integration
	bosonexprvec boson_expr = { nullptr, nullptr, nullptr, nullptr };
	switch (boson_star_type)
	{
	case BOSON_STAR_MINI:
		boson_expr = gen_func_miniboson_hp(omega, mass);
		break;

	case BOSON_STAR_MASSIVE:
		boson_expr = gen_func_massiveboson_hp(omega, mass, lambda, kappa);
		break;

	case BOSON_STAR_SOLITONIC:
		boson_expr = gen_func_solitonicboson_hp(omega, mass, sigma, kappa);
		break;

	default:
		boson_expr = gen_func_miniboson_hp(omega, mass, kappa);
		break;
	}

	BosonIntegrateParamsHP params;
	init_integrate_params(&params);
	set_vec_mpf_d(&params.initial_value, initialvalues);
	mpf_set_d(params.step_size, step_size_value);
	mpf_set_ui(params.domain_start, 0);
	mpf_set_ui(params.domain_end, domain_end_value);
	set_vec_mpf_d(&params.thresholds, thresholds);
	params.divergence_variable = VAR_PHI0;
	params.func = [&boson_expr](vec<BOSON_SYSTEM_DIMENSION + 1, mpf_t>& vars, bosonvec_hp* output) {
		func_boson_hp<BOSON_SYSTEM_DIMENSION>(vars, output, boson_expr);
	};

	// Actually integrate the system
	vector<bosonvec_hp> result;
	params.result = &result;
	cout << "Integrating system..." << endl;
	IntegrateResult output = integrate_hp(params);
	cout << "Integration finished... writing output..." << endl;

	if (output.bThresholdReached)
	{
		cout << endl << "THRESHOLD REACHED!" << endl;
		cout << (output.divergenceDirection == DIRECTION_POSITIVE ? "POSITIVE DIRECTION" : "NEGATIVE DIRECTION") << endl;
		cout << "POSITION: " << output.lastIndex + 1 << endl << endl;
	}

	// Output the result to file
	vector<bosonvec> result_sp;
	convert_to_double(&result_sp, result);
	// Remove the divergence from the result
	if (output.bThresholdReached)
		for (size_t j = 0; j < BOSON_SYSTEM_DIMENSION; ++j)
			result_sp[output.lastIndex + 1][j] = 0.0;
	output_result(result_sp, mpf_get_d(params.step_size));

	// Clean up
	boson_expr.release();
	clear_integrate_params(&params);
	mpf_clear(omega);
	mpf_clear(mass);
	mpf_clear(omega);
	mpf_clear(mass);
	mpf_clear(lambda);
	mpf_clear(sigma);

	cout << "Success!" << endl;
}

void run_in_manual_search_mode()
{
	mpf_t mass;
	mpf_t omega;
	mpf_t lambda;
	mpf_t sigma;
	mpf_t kappa;
	mpf_init(mass);
	mpf_init(omega);
	mpf_init(sigma);
	mpf_init(lambda);
	mpf_init(kappa);
	mpf_set_d(mass, mass_value);
	mpf_set_d(lambda, lambda_value);
	mpf_set_d(sigma, sigma_value);
	mpf_set_d(kappa, kappa_value);

	vector<bosonvec_hp> result;

	ManualRangeBosonSearchParamsHP man_range_params;
	man_range_params.lower_bound = manual_lower_bound;
	man_range_params.upper_bound = manual_upper_bound;
	man_range_params.bisection_count_max = bisection_count_max_value;
	man_range_params.suggested_bisection_count = suggested_bisection_count_value;
	man_range_params.b_inverted_search = (bool)b_inverted_search;
	man_range_params.search_variable = (ManualRangeSearchVariable)search_variable;
	man_range_params.fix_variable_value = fix_variable_value;

	BosonSearchParamsCustomHP search_params;
	init_custom_boson_search_params(&search_params);
	mpf_set_ui(search_params.domain_start, 0);
	mpf_set_ui(search_params.domain_end, domain_end_value);
	mpf_set_d(search_params.step_size, step_size_value);
	search_params.result = &result;

	switch (boson_star_type)
	{
	case BOSON_STAR_MINI:
		search_params.func = [&mass, &kappa](const mpf_t omega) {
			return gen_func_miniboson_hp(omega, mass, kappa);
		};
		break;

	case BOSON_STAR_MASSIVE:
		search_params.func = [&mass, &lambda, &kappa](const mpf_t omega) {
			return gen_func_massiveboson_hp(omega, mass, lambda, kappa);
		};
		break;

	case BOSON_STAR_SOLITONIC:
		search_params.func = [&mass, &sigma, &kappa](const mpf_t omega) {
			return gen_func_solitonicboson_hp(omega, mass, sigma, kappa);
		};
		break;

	default:
		search_params.func = [&mass, &kappa](const mpf_t omega) {
			return gen_func_miniboson_hp(omega, mass, kappa);
		};
		break;
	}

	ManualRangeBosonSearchHP search(man_range_params);
	boson_search_custom_hp(search_params, &search);

	cout << "Final Bisection Variable Value = " << mpf_get_d(*(search.GetBisectionVar())) << endl;
	cout << "Final Fix Variable Value = " << mpf_get_d(*(search.GetOtherVar())) << endl;

	vector<bosonvec> result_sp;
	convert_to_double(&result_sp, result);
	output_result(result_sp, step_size_value);
	output_mass(result_sp, step_size_value);

	clear_custom_boson_search_params(&search_params);

	mpf_clear(omega);
	mpf_clear(mass);
	mpf_clear(lambda);
	mpf_clear(sigma);
}

int main(int argc, const char* argv[])
{
	loadParameters(argc, argv);

	switch(program_mode)
	{
	case PROGRAM_MODE_DEFAULT:
		if (use_high_precision)
			bisection_search_hp(boson_star_type);
		else
			bisection_search(boson_star_type);
		break;
	case PROGRAM_MODE_MASS_V_RADIUS:
		cout << "PROGRAM MODE: OUTPUT MASS V. RADIUS CURVE" << endl << endl;
		mass_v_radius_hp(boson_star_type);
		break;
	case PROGRAM_MODE_FIND_TRANSITION:
		cout << "PROGRAM MODE: FIND CRITICAL TRANSITION" << endl << endl;
		critical_transition_find_hp();
		break;
	case PROGRAM_MODE_MASS_V_RADIUS_ADAPTIVE:
		cout << "PROGRAM MODE: OUTPUT MASS V. RADIUS CURVE ADAPTIVE" << endl << endl;
		cout << "Mode not implemented" << endl;
		break;
	case PROGRAM_MODE_SPECIFY_OMEGA:
		cout << "PROGRAM MODE: COMPUTE BOSON STAR FOR SPECIFIC FREQUENCY VALUE" << endl << endl;
		run_in_specify_frequency_mode();
		break;
	case PROGRAM_MODE_MANUAL_SEARCH:
		cout << "PROGRAM MODE: MANUAL SEARCH" << endl << endl;
		run_in_manual_search_mode();
		break;
	}

	return 0;
}
