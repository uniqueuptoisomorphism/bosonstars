/*
 * main.cpp
 */

#include <iostream>
#include <stdlib.h>
#include <sdf.h>
#include <bbhutil.h>
#include <math.h>
#include <memory.h>
#include <string>
#include <vector>
#include <functional>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <fstream>
#include "linear.h"

#define VAR_COUNT 12
#define NAME_MAX_LENGTH 128
#define TAG_MAX_LENGTH 32
#define OUTPUT_COUNT 12
#define RUN_DIAGONAL_CHECK
#define INITIAL_DATA_RELAXATION_COUNT 50
#define ERROR_TOLERANCE 0.00000001
#define MAX_MATRIX_ITERATIONS 100
#define SPARSE_MASK_COUNT 6
// #define OUTPUT_RELAXATIONS
#define SPARSE_MASK_SIZE VAR_COUNT + 1
#define SPARSE_MASK_RESIDUAL_INDEX VAR_COUNT

using namespace std;
using namespace Eigen;

static string var_names[] = { "phi1", "xi1", "pi1", "phi2", "xi2", "pi2", "phi3", "xi3", "pi3", "psi", "beta", "alpha" };
static string psi_residual_name = "psi_residual";

static const string var_names_debug[] = { "phi1_debug", "xi1_debug", "pi1_debug", "phi2_debug", "xi2_debug",
		"pi2_debug", "phi3_debug", "xi3_debug", "pi3_debug", "psi_debug", "beta_debug", "alpha_debug" };
static const string bin_input_file_name = "phi0.bin";
//static char mass_aspect_name[] = "mass_aspect";
//static char field_name[] = "field";
static const double kappa = 8.0 * M_PI;

static const bool var_output_masks[SPARSE_MASK_COUNT][VAR_COUNT + 1] = {
	{ true, 	true, 	true, 	true, 	true, 	true, 	true, 	true, 	true, 	true, 	true, 	true, 	true },
	{ true, 	false, 	false, 	true, 	false, 	false, 	true, 	false, 	false, 	true, 	true, 	true, 	true },
	{ true, 	false, 	false, 	true, 	false, 	false, 	true, 	false, 	false, 	false, 	false, 	false, 	true },
	{ true, 	false, 	false, 	true, 	false, 	false, 	false, 	false, 	false, 	false, 	false, 	false, 	true },
	{ true, 	false, 	false, 	false, 	false, 	false, 	false, 	false, 	false, 	false, 	false, 	false, 	true },
	{ false, 	false, 	false, 	false, 	false, 	false, 	false, 	false, 	false, 	false, 	false, 	false, 	true }
};

static bool var_output_mask[VAR_COUNT + 1];

enum Variables_debug
{
	VAR_PHI1,
	VAR_XI1,
	VAR_PI1,
	VAR_PHI2,
	VAR_XI2,
	VAR_PI2,
	VAR_PHI3,
	VAR_XI3,
	VAR_PI3,
	VAR_PSI,
	VAR_BETA,
	VAR_ALPHA,
};

enum InitDataType
{
	DATA_TYPE_NOTHING = 0,

	DATA_TYPE_WAVE_PACKET1 = 1,
	DATA_TYPE_LEFT_MOVING_WAVE_PACKET1 = 2,
	DATA_TYPE_RIGHT_MOVING_WAVE_PACKET1 = 3,

	DATA_TYPE_WAVE_PACKET2 = 4,
	DATA_TYPE_LEFT_MOVING_WAVE_PACKET2 = 5,
	DATA_TYPE_RIGHT_MOVING_WAVE_PACKET2 = 6,

	DATA_TYPE_WAVE_PACKET3 = 7,
	DATA_TYPE_LEFT_MOVING_WAVE_PACKET3 = 8,
	DATA_TYPE_RIGHT_MOVING_WAVE_PACKET3 = 9,

	DATA_TYPE_PERTURBATION_TEST = 10,

	DATA_TYPE_PERTURBATION_GAUSSIAN = 11,

	DATA_TYPE_PERTURBATION_FROM_FILE = 12
};

void ReadBinFile(const string& fileName, vector<double>* output, bool* bSuccess)
{
	ifstream file(fileName, ios::in | ios::binary);
	if (file)
	{
		*bSuccess = true;
		file.seekg(0, std::ios::end);
		size_t len = file.tellg();
		file.seekg(0, std::ios::beg);
		output->resize(len / sizeof(double));
		file.read((char*)output->data(), len);
	}
	else
		*bSuccess = false;
}

class Application
{
private:
    int gridsize;
    int Nt;
    int var_count;
    int b_use_minkowski;
    int init_data_type;
    int output_time_skip;
    int sparse_output_level;
    double lambda;
    double phi_amp;
    double phi_delta;
    double phi_r0;
    double phi_perturb_amp;
	double phi_perturb_delta;
	double phi_perturb_r0;
    double radius;
    double event_horizon_radius;
    double mass;
    double init_data_omega;
    char* pFile;
    char output_tag[TAG_MAX_LENGTH];

    double bbox[2];

    double dr;
    double dt;

    double* phi_init;
    double* var_buffer;
    double* var_buffer_new;
    double* mass_aspect;
    double* field;

    int time;

    GaussSolver solver;
    size_t max_iterations_per_step = 100;

    inline double sqr(const double d)
    {
    	return d * d;
    }

    double Power(const double d, const int i)
    {
    	if (i < 0)
    		return 1.0 / Power(d, -i);
    	else if (i == 0)
    		return 1.0;
    	else
    		return d * Power(d, i - 1);
    }

	inline double xi1(const int i) 						{ return var_buffer_new[VAR_XI1 * gridsize + i]; }
	inline double xi1(const int i, const int t) 		{ return var_buffer[VAR_XI1 * gridsize + i]; }
	inline double pi1(const int i) 						{ return var_buffer_new[VAR_PI1 * gridsize + i]; }
	inline double pi1(const int i, const int t) 		{ return var_buffer[VAR_PI1 * gridsize + i]; }
	inline double phi1(const int i)						{ return var_buffer_new[VAR_PHI1 * gridsize + i]; }
	inline double phi1(const int i, const int t)		{ return var_buffer[VAR_PHI1 * gridsize + i]; }
	inline double xi2(const int i) 						{ return var_buffer_new[VAR_XI2 * gridsize + i]; }
	inline double xi2(const int i, const int t) 		{ return var_buffer[VAR_XI2 * gridsize + i]; }
	inline double pi2(const int i) 						{ return var_buffer_new[VAR_PI2 * gridsize + i]; }
	inline double pi2(const int i, const int t) 		{ return var_buffer[VAR_PI2 * gridsize + i]; }
	inline double phi2(const int i)						{ return var_buffer_new[VAR_PHI2 * gridsize + i]; }
	inline double phi2(const int i, const int t)		{ return var_buffer[VAR_PHI2 * gridsize + i]; }
	inline double xi3(const int i) 						{ return var_buffer_new[VAR_XI3 * gridsize + i]; }
	inline double xi3(const int i, const int t) 		{ return var_buffer[VAR_XI3 * gridsize + i]; }
	inline double pi3(const int i) 						{ return var_buffer_new[VAR_PI3 * gridsize + i]; }
	inline double pi3(const int i, const int t) 		{ return var_buffer[VAR_PI3 * gridsize + i]; }
	inline double phi3(const int i)						{ return var_buffer_new[VAR_PHI3 * gridsize + i]; }
	inline double phi3(const int i, const int t)		{ return var_buffer[VAR_PHI3 * gridsize + i]; }
	inline double psi(const int i) 						{ return var_buffer_new[VAR_PSI * gridsize + i]; }
	inline double psi(const int i, const int t) 		{ return var_buffer[VAR_PSI * gridsize + i]; }
	inline double beta(const int i) 					{ return var_buffer_new[VAR_BETA * gridsize + i]; }
	inline double beta(const int i, const int t) 		{ return var_buffer[VAR_BETA * gridsize + i]; }
	inline double alpha(const int i) 					{ return var_buffer_new[VAR_ALPHA * gridsize + i]; }
	inline double alpha(const int i, const int t) 		{ return var_buffer[VAR_ALPHA * gridsize + i]; }
    inline double r(const int i) 						{ return bbox[0] + i * dr; }
    inline double r(const int i, const int t) 			{ return r(i); }
    inline double gamma_static(const int i) 			{ return (r(i) / (r(i) + 2.0 * mass)); }
    inline double beta_static(const int i) 				{ return (2.0 * mass) / (r(i) + 2.0 * mass); }
    inline double V(const double phiSqr) 				{ return mass * phiSqr; }
    inline double DV(const double phiSqr)				{ return mass; }
    inline double V(const double phiSqr, const int t) 	{ return V(phiSqr); }
    inline double DV(const double phiSqr, const int t)	{ return DV(phiSqr); }

	void InitGaussian(double vec[], const double amp, const double delta, const double x0)
	{
	    double dr = (bbox[1] - bbox[0]) / (gridsize - 1);

	    for (int i = 0; i < gridsize; ++i)
	    {
	        double r = i * dr + bbox[0];
	        vec[i] = amp * exp(-((r - x0) * (r - x0) / (delta * delta)));
	    }
	}

	void ComputeSpatialDerivative(double dest_vec[], const double src_vec[])
	{
	    double dr = (bbox[1] - bbox[0]) / (gridsize - 1);

	    for (int i = 1; i < gridsize - 1; ++i)
	        dest_vec[i] = (src_vec[i + 1] - src_vec[i - 1]) / (2.0 * dr);
	    dest_vec[0] = (src_vec[1] - src_vec[0]) / dr;
	    dest_vec[gridsize - 1] = (src_vec[gridsize - 1] - src_vec[gridsize - 2]) / dr;
	}

	void LoadParameters()
	{
	    get_int_param(pFile, "Nr", &gridsize, 1);
	    get_int_param(pFile, "Nt", &Nt, 1);
	    get_int_param(pFile, "b_use_minkowski", &b_use_minkowski, 1);
	    get_int_param(pFile, "init_data_type", &init_data_type, 1);
	    get_int_param(pFile, "output_time_skip", &output_time_skip, 1);
	    get_int_param(pFile, "sparse_output_level", &sparse_output_level, 1);
	    get_real_param(pFile, "lambda", &lambda, 1);
	    get_real_param(pFile, "phi_amp", &phi_amp, 1);
	    get_real_param(pFile, "phi_delta", &phi_delta, 1);
	    get_real_param(pFile, "phi_x0", &phi_r0, 1);
	    get_real_param(pFile, "phi_perturb_amp", &phi_perturb_amp, 1);
	    get_real_param(pFile, "phi_perturb_delta", &phi_perturb_delta, 1);
	    get_real_param(pFile, "phi_perturb_x0", &phi_perturb_r0, 1);
	    get_real_param(pFile, "radius", &radius, 1);
	    get_real_param(pFile, "mass", &mass, 1);
	    get_real_param(pFile, "init_data_omega", &init_data_omega, 1);

	    char* output_tag_ptr = nullptr;
	    get_str_param(pFile, "output_tag", &output_tag_ptr, 1);

	    for (auto& str : var_names)
	    	str.append(output_tag_ptr);
	    psi_residual_name.append(output_tag_ptr);

	    sparse_output_level = std::min(std::max(sparse_output_level, 0), SPARSE_MASK_COUNT - 1);
	    memcpy(var_output_mask, var_output_masks[sparse_output_level], SPARSE_MASK_SIZE);
	}

	void InitVariables()
	{
		event_horizon_radius = 2.0;

		bbox[0] = 			0.0;
		// bbox[0] = 			event_horizon_radius + 1.0;
		bbox[1] = 			radius;
	    dr = 				(bbox[1] - bbox[0]) / (gridsize - 1);
	    dt = 				lambda * dr;

	    phi_init =					new double[gridsize];
	    var_buffer = 				new double[gridsize * var_count];
	    var_buffer_new = 			new double[gridsize * var_count];
	    mass_aspect = 				new double[gridsize];
	    field =						new double[gridsize];

	    double* phi1_ptr =	&var_buffer[gridsize * VAR_PHI1];
	    double* phi2_ptr = 	&var_buffer[gridsize * VAR_PHI2];
	    double* phi3_ptr =	&var_buffer[gridsize * VAR_PHI3];
	    double* xi1_ptr =	&var_buffer[gridsize * VAR_XI1];
	    double* pi1_ptr =	&var_buffer[gridsize * VAR_PI1];
	    double* xi2_ptr =	&var_buffer[gridsize * VAR_XI2];
		double* pi2_ptr =	&var_buffer[gridsize * VAR_PI2];
		double* xi3_ptr = 	&var_buffer[gridsize * VAR_XI3];
		double* pi3_ptr = 	&var_buffer[gridsize * VAR_PI3];
	    double* psi_ptr =	&var_buffer[gridsize * VAR_PSI];
	    double* beta_ptr =	&var_buffer[gridsize * VAR_BETA];
	    double* alpha_ptr =	&var_buffer[gridsize * VAR_ALPHA];

	    memset(var_buffer, 0, sizeof(double) * gridsize * var_count);

	    switch (init_data_type)
	    {
	    case DATA_TYPE_NOTHING:
	    	break;

	    case DATA_TYPE_WAVE_PACKET1:
	    case DATA_TYPE_RIGHT_MOVING_WAVE_PACKET1:
	    case DATA_TYPE_LEFT_MOVING_WAVE_PACKET1:
		    InitGaussian(phi1_ptr, phi_amp, phi_delta, phi_r0);
		    ComputeSpatialDerivative(xi1_ptr, phi1_ptr);
		    if (init_data_type == DATA_TYPE_LEFT_MOVING_WAVE_PACKET1)
			    for (int i = 0; i < gridsize; ++i)
			    	pi1_ptr[i] = xi1_ptr[i];
		    else if (init_data_type == DATA_TYPE_RIGHT_MOVING_WAVE_PACKET1)
				for (int i = 0; i < gridsize; ++i)
					pi1_ptr[i] = -xi1_ptr[i];
	    	break;

	    case DATA_TYPE_WAVE_PACKET2:
	    case DATA_TYPE_RIGHT_MOVING_WAVE_PACKET2:
	    case DATA_TYPE_LEFT_MOVING_WAVE_PACKET2:
		    InitGaussian(phi2_ptr, phi_amp, phi_delta, phi_r0);
		    ComputeSpatialDerivative(xi2_ptr, phi2_ptr);
		    if (init_data_type == DATA_TYPE_LEFT_MOVING_WAVE_PACKET2)
			    for (int i = 0; i < gridsize; ++i)
			    	pi2_ptr[i] = xi2_ptr[i];
		    else if (init_data_type == DATA_TYPE_RIGHT_MOVING_WAVE_PACKET2)
				for (int i = 0; i < gridsize; ++i)
					pi2_ptr[i] = -xi2_ptr[i];
	    	break;

	    case DATA_TYPE_WAVE_PACKET3:
		case DATA_TYPE_RIGHT_MOVING_WAVE_PACKET3:
		case DATA_TYPE_LEFT_MOVING_WAVE_PACKET3:
			InitGaussian(phi3_ptr, phi_amp, phi_delta, phi_r0);
			ComputeSpatialDerivative(xi3_ptr, phi3_ptr);
			if (init_data_type == DATA_TYPE_LEFT_MOVING_WAVE_PACKET3)
				for (int i = 0; i < gridsize; ++i)
					pi3_ptr[i] = xi3_ptr[i];
			else if (init_data_type == DATA_TYPE_RIGHT_MOVING_WAVE_PACKET3)
				for (int i = 0; i < gridsize; ++i)
					pi3_ptr[i] = -xi3_ptr[i];
			break;

		case DATA_TYPE_PERTURBATION_TEST:
			InitGaussian(phi1_ptr, phi_amp, phi_delta, phi_r0);
			ComputeSpatialDerivative(xi1_ptr, phi1_ptr);
			InitGaussian(phi3_ptr, phi_perturb_amp, phi_perturb_delta, phi_perturb_r0);
			ComputeSpatialDerivative(xi3_ptr, phi3_ptr);
			for (int i = 0; i < gridsize; ++i)
				pi3_ptr[i] = xi3_ptr[i];
			break;

		case DATA_TYPE_PERTURBATION_GAUSSIAN:
			InitGaussian(phi3_ptr, phi_perturb_amp, phi_perturb_delta, phi_perturb_r0);
			ComputeSpatialDerivative(xi3_ptr, phi3_ptr);
			for (int i = 0; i < gridsize; ++i)
				pi3_ptr[i] = xi3_ptr[i];

			InitGaussian(phi1_ptr, phi_amp, phi_delta, phi_r0);
			ComputeSpatialDerivative(xi1_ptr, phi1_ptr);
			for (int i = 0; i < gridsize; ++i)
				pi2_ptr[i] = phi1_ptr[i] * init_data_omega;
			break;

		case DATA_TYPE_PERTURBATION_FROM_FILE:
			InitGaussian(phi3_ptr, phi_perturb_amp, phi_perturb_delta, phi_perturb_r0);
			ComputeSpatialDerivative(xi3_ptr, phi3_ptr);
			for (int i = 0; i < gridsize; ++i)
				pi3_ptr[i] = xi3_ptr[i];

			vector<double> data;
			bool bSuccess;
			ReadBinFile(bin_input_file_name, &data, &bSuccess);
			if (bSuccess)
			{
				size_t copyLen = std::min(data.size(), (size_t)gridsize);
				memcpy(phi1_ptr, data.data(), copyLen * sizeof(double));
				ComputeSpatialDerivative(xi1_ptr, phi1_ptr);
				for (int i = 0; i < gridsize; ++i)
					pi2_ptr[i] = phi1_ptr[i] * init_data_omega;
			}

			break;
	    }

	    for (int i = 0; i < gridsize; ++i)
	    {
	    	psi_ptr[i] = 1.0;
			beta_ptr[i] = 0.0;
			alpha_ptr[i] = 1.0;
	    }

	    memcpy(var_buffer_new, var_buffer, sizeof(double) * gridsize * var_count);

	    if (!b_use_minkowski)
	    	RelaxInitialData();

	    // ComputeResiduals();
	}

	// For debug purposes
	/*void ComputeResiduals()
	{
		VectorXd res(gridsize);
		ComputeResidualPsi(res);
		ComputeResidualBeta(res);
		ComputeResidualAlpha(res);
	}*/

	bool IsResidualSmall(double vec[], int size)
	{
		double largest = 0.0;
		for (int i = 0; i < size; ++i)
			if (largest < fabs(vec[i]))
				largest = fabs(vec[i]);

		return largest < ERROR_TOLERANCE;
	}

	bool IsResidualSmall(VectorXd vec, int size)
	{
		double largest = 0.0;
		for (int i = 0; i < size; ++i)
		{
			double d = vec(i);
			if (largest < fabs(d))
				largest = fabs(d);
		}

		return largest < ERROR_TOLERANCE;
	}

	void ComputeResidualPsi(VectorXd& residual)
	{
		int i = 0;
		residual(gridsize * 0 + i) = (-3*psi(i) + 4*psi(1 + i) - psi(2 + \
	i))/(2.*dr);

		for (i = 1; i < gridsize - 1; ++i)
		{
			residual(gridsize * 0 + i) = (psi(-1 + i) - 2*psi(i) + psi(1 + \
	i))/Power(dr,2) + (Power(psi(i),5)*Power((-beta(-1 + i) + beta(1 + \
	i))/(2.*dr) - beta(i)/r(i),2))/(12.*Power(alpha(i),2)) + (-psi(-1 + \
	i) + psi(1 + i))/(dr*r(i)) + (kappa*psi(i)*(Power(pi1(i),2) + \
	Power(pi2(i),2) + Power(pi3(i),2) + \
	Power(psi(i),4)*V(Power(phi1(i),2) + Power(phi2(i),2)) + \
	Power(xi1(i),2) + Power(xi2(i),2) + Power(xi3(i),2)))/4.;
		}

		residual(gridsize * 0 + i) = -1 + psi(i);
	}

	void IteratePsi(bool* bResidualSmall)
	{
		double* psi_ptr =	&var_buffer_new[gridsize * VAR_PSI];

		VectorXd res(gridsize);
		VectorXd delta(gridsize);

		ComputeResidualPsi(res);

		*bResidualSmall = IsResidualSmall(res, gridsize);

		if (!*bResidualSmall)
		{
			SparseMatrix<double> system(gridsize, gridsize);
			vector<Triplet<double> > coeffs;

			{
				int i = 0;
				coeffs.push_back(Triplet<double>(gridsize * 0 + i, gridsize * 0 + i \
			+ 0, -3/(2.*dr)));
				coeffs.push_back(Triplet<double>(gridsize * 0 + i, gridsize * 0 + i \
			+ 1, 2/dr));
				coeffs.push_back(Triplet<double>(gridsize * 0 + i, gridsize * 0 + i \
			+ 2, -1/(2.*dr)));

				for (i = 1; i < gridsize - 1; ++i)
				{
					coeffs.push_back(Triplet<double>(gridsize * 0 + i, gridsize * 0 + i \
			+ -1, (1 - dr/r(i))/Power(dr,2)));
					coeffs.push_back(Triplet<double>(gridsize * 0 + i, gridsize * 0 + i \
			+ 0, -2/Power(dr,2) + (5*Power(psi(i),4)*Power(2*dr*beta(i) + \
			(beta(-1 + i) - beta(1 + \
			i))*r(i),2))/(48.*Power(dr,2)*Power(alpha(i),2)*Power(r(i),2)) + \
			kappa*Power(psi(i),4)*V(Power(phi1(i),2) + Power(phi2(i),2)) + \
			(kappa*(Power(pi1(i),2) + Power(pi2(i),2) + Power(pi3(i),2) + \
			Power(psi(i),4)*V(Power(phi1(i),2) + Power(phi2(i),2)) + \
			Power(xi1(i),2) + Power(xi2(i),2) + Power(xi3(i),2)))/4.));
					coeffs.push_back(Triplet<double>(gridsize * 0 + i, gridsize * 0 + i \
			+ 1, (dr + r(i))/(Power(dr,2)*r(i))));
				}

				coeffs.push_back(Triplet<double>(gridsize * 0 + i, gridsize * 0 + i \
			+ 0, 1));
			}

			system.setFromTriplets(coeffs.begin(), coeffs.end());
			SparseLU<SparseMatrix<double> > solver(system);
			delta = solver.solve(res);

			VectorXd t = system * delta;
			VectorXd t2 = t - res;

			for (int i = 0; i < gridsize; ++i)
				psi_ptr[i] -= delta(i);
		}
	}

	void ComputeResidualBeta(VectorXd& residual)
	{
		int i = 0;
		residual(gridsize * 0 + i) = beta(i);

		for (i = 1; i < gridsize - 1; ++i)
		{
			residual(gridsize * 0 + i) = (beta(-1 + i) - 2*beta(i) + beta(1 + \
	i))/Power(dr,2) + (-(-alpha(-1 + i) + alpha(1 + i))/(2.*dr*alpha(i)) \
	+ (3*(-psi(-1 + i) + psi(1 + i)))/(dr*psi(i)) + 2/r(i))*((-beta(-1 + \
	i) + beta(1 + i))/(2.*dr) - beta(i)/r(i)) + \
	(3*kappa*alpha(i)*(pi1(i)*xi1(i) + pi2(i)*xi2(i) + \
	pi3(i)*xi3(i)))/Power(psi(i),2);
		}

		residual(gridsize * 0 + i) = beta(i);
	}

	void IterateBeta(bool* bResidualSmall)
	{
		double* beta_ptr =	&var_buffer_new[gridsize * VAR_BETA];

		VectorXd res(gridsize);
		VectorXd delta(gridsize);

		ComputeResidualBeta(res);

		*bResidualSmall = IsResidualSmall(res, gridsize);

		if (!*bResidualSmall)
		{
			SparseMatrix<double> system(gridsize, gridsize);
			vector<Triplet<double> > coeffs;

			{
				int i = 0;
				coeffs.push_back(Triplet<double>(gridsize * 0 + i, gridsize * 0 + i \
			+ 0, 1));

				for (i = 1; i < gridsize - 1; ++i)
				{
					coeffs.push_back(Triplet<double>(gridsize * 0 + i, gridsize * 0 + i \
			+ -1, (4 - alpha(-1 + i)/alpha(i) + alpha(1 + i)/alpha(i) + (6*psi(-1 \
			+ i))/psi(i) - (6*psi(1 + i))/psi(i) - (4*dr)/r(i))/(4.*Power(dr,2))));
					coeffs.push_back(Triplet<double>(gridsize * 0 + i, gridsize * 0 + i \
			+ 0, (-4 - (4*Power(dr,2))/Power(r(i),2) + (dr*((-alpha(-1 + i) + \
			alpha(1 + i))*psi(i) + 6*alpha(i)*(psi(-1 + i) - psi(1 + \
			i))))/(alpha(i)*psi(i)*r(i)))/(2.*Power(dr,2))));
					coeffs.push_back(Triplet<double>(gridsize * 0 + i, gridsize * 0 + i \
			+ 1, ((alpha(-1 + i) - alpha(1 + i))*psi(i)*r(i) + \
			alpha(i)*(6*(-psi(-1 + i) + psi(1 + i))*r(i) + 4*psi(i)*(dr + \
			r(i))))/(4.*Power(dr,2)*alpha(i)*psi(i)*r(i))));
				}

				coeffs.push_back(Triplet<double>(gridsize * 0 + i, gridsize * 0 + i \
			+ 0, 1));
			}

			system.setFromTriplets(coeffs.begin(), coeffs.end());
			SparseLU<SparseMatrix<double> > solver(system);
			delta = solver.solve(res);

			for (int i = 0; i < gridsize; ++i)
				beta_ptr[i] -= delta(i);
		}
	}

	void ComputeResidualAlpha(VectorXd& residual)
	{
		int i = 0;
		residual(gridsize * 0 + i) = (-3*alpha(i) + 4*alpha(1 + i) - alpha(2 \
	+ i))/(2.*dr);

		for (i = 1; i < gridsize - 1; ++i)
		{
			residual(gridsize * 0 + i) = (alpha(-1 + i) - 2*alpha(i) + alpha(1 \
	+ i))/Power(dr,2) - 2*kappa*alpha(i)*(Power(pi1(i),2) + \
	Power(pi2(i),2) + Power(pi3(i),2)) + ((-alpha(-1 + i) + alpha(1 + \
	i))*((-psi(-1 + i) + psi(1 + i))/(dr*psi(i)) + 2/r(i)))/(2.*dr) - \
	(2*Power(psi(i),4)*Power((-beta(-1 + i) + beta(1 + i))/(2.*dr) - \
	beta(i)/r(i),2))/(3.*alpha(i)) + \
	kappa*alpha(i)*Power(psi(i),4)*V(Power(phi1(i),2) + \
	Power(phi2(i),2));
		}

		residual(gridsize * 0 + i) = -1 + alpha(i);
	}

	void IterateAlpha(bool* bResidualSmall)
	{
		double* alpha_ptr =	&var_buffer_new[gridsize * VAR_ALPHA];

		VectorXd res(gridsize);
		VectorXd delta(gridsize);

		ComputeResidualAlpha(res);

		*bResidualSmall = IsResidualSmall(res, gridsize);

		if (!*bResidualSmall)
		{
			SparseMatrix<double> system(gridsize, gridsize);
			vector<Triplet<double> > coeffs;

			{
				int i = 0;
				coeffs.push_back(Triplet<double>(gridsize * 0 + i, gridsize * 0 + i \
			+ 0, -3/(2.*dr)));
				coeffs.push_back(Triplet<double>(gridsize * 0 + i, gridsize * 0 + i \
			+ 1, 2/dr));
				coeffs.push_back(Triplet<double>(gridsize * 0 + i, gridsize * 0 + i \
			+ 2, -1/(2.*dr)));

				for (i = 1; i < gridsize - 1; ++i)
				{
					coeffs.push_back(Triplet<double>(gridsize * 0 + i, gridsize * 0 + i \
			+ -1, (2 + psi(-1 + i)/psi(i) - psi(1 + i)/psi(i) - \
			(2*dr)/r(i))/(2.*Power(dr,2))));
					coeffs.push_back(Triplet<double>(gridsize * 0 + i, gridsize * 0 + i \
			+ 0, -2/Power(dr,2) - 2*kappa*(Power(pi1(i),2) + Power(pi2(i),2) + \
			Power(pi3(i),2)) + (Power(psi(i),4)*Power(2*dr*beta(i) + (beta(-1 + \
			i) - beta(1 + \
			i))*r(i),2))/(6.*Power(dr,2)*Power(alpha(i),2)*Power(r(i),2)) + \
			kappa*Power(psi(i),4)*V(Power(phi1(i),2) + Power(phi2(i),2))));
					coeffs.push_back(Triplet<double>(gridsize * 0 + i, gridsize * 0 + i \
			+ 1, ((-psi(-1 + i) + psi(1 + i))*r(i) + 2*psi(i)*(dr + \
			r(i)))/(2.*Power(dr,2)*psi(i)*r(i))));
				}

				coeffs.push_back(Triplet<double>(gridsize * 0 + i, gridsize * 0 + i \
			+ 0, 1));
			}

			system.setFromTriplets(coeffs.begin(), coeffs.end());
			SparseLU<SparseMatrix<double> > solver(system);
			delta = solver.solve(res);

			for (int i = 0; i < gridsize; ++i)
				alpha_ptr[i] -= delta(i);
		}
	}

	void ComputeResidualPsiEvolution(double residual[])
	{
		residual[0] = 0.;
		for (int i = 1; i < gridsize - 1; ++i)
		{
			/*residual[i] = ((-beta(-1 + i) + beta(1 + i))*psi(i))/(12.*dr) + (beta(i)*(-psi(-1 + i) + psi(1 + i)))/(2.*dr) - (psi(i) - psi(i,-1))/dt +
					   (beta(i)*psi(i))/(3.*r(i));*/
			residual[i] = 0.5 * (((-beta(-1 + i) + beta(1 + i))*psi(i))/(12.*dr) + (beta(i)*(-psi(-1 + i) + psi(1 + i)))/(2.*dr) + (beta(i)*psi(i))/(3.*r(i)))
							+ 0.5 * (((-beta(-1 + i, -1) + beta(1 + i, -1))*psi(i, -1))/(12.*dr) + (beta(i, -1)*(-psi(-1 + i, -1) + psi(1 + i, -1)))/(2.*dr) + (beta(i, -1)*psi(i, -1))/(3.*r(i, -1)))
							- (psi(i) - psi(i,-1))/dt;
		}
		residual[gridsize - 1] = 0.;
	}

	void RelaxInitialData()
	{
		char name_psi[100] = "psi_relax";
		char name_beta[100] = "beta_relax";
		char name_alpha[100] = "alpha_relax";

		gft_out_bbox(name_psi, 0, &gridsize, 1, bbox, &var_buffer_new[gridsize * VAR_PSI]);
		gft_out_bbox(name_beta, 0, &gridsize, 1, bbox, &var_buffer_new[gridsize * VAR_BETA]);
		gft_out_bbox(name_alpha, 0, &gridsize, 1, bbox, &var_buffer_new[gridsize * VAR_ALPHA]);

		cout << "Relaxing initial data..." << endl;

		bool bConvergence = false;
		for (size_t i = 0; i < INITIAL_DATA_RELAXATION_COUNT && !bConvergence; ++i)
		{
			cout << "iteration " << i + 1 << " / " << INITIAL_DATA_RELAXATION_COUNT << endl;

			bool bResidualSmall = false;
			bConvergence = true;

			IteratePsi(&bResidualSmall);
			gft_out_bbox(name_psi, i + 1, &gridsize, 1, bbox, &var_buffer_new[gridsize * VAR_PSI]);
			bConvergence = bConvergence && bResidualSmall;

			IterateBeta(&bResidualSmall);
			gft_out_bbox(name_beta, i + 1, &gridsize, 1, bbox, &var_buffer_new[gridsize * VAR_BETA]);
			bConvergence = bConvergence && bResidualSmall;

			IterateAlpha(&bResidualSmall);
			gft_out_bbox(name_alpha, i + 1, &gridsize, 1, bbox, &var_buffer_new[gridsize * VAR_ALPHA]);
			bConvergence = bConvergence && bResidualSmall;
		}

		memcpy(var_buffer, var_buffer_new, gridsize * VAR_COUNT * sizeof(double));

		cout << endl;
	}


	void UpdatePhi(bool* bConvergence)
	{
		*bConvergence = true;

		double* phi1_ptr =	&var_buffer_new[gridsize * VAR_PHI1];
		double* phi2_ptr =	&var_buffer_new[gridsize * VAR_PHI2];
		double* phi3_ptr = 	&var_buffer_new[gridsize * VAR_PHI3];

		for (int i = 0; i < gridsize; ++i)
		{
			double temp = phi1(i, -1) + ((pi1(i) * alpha(i) / (psi(i) * psi(i)) + beta(i) * xi1(i)) +
					(pi1(i, -1) * alpha(i, -1) / (psi(i, - 1) * psi(i, - 1)) + beta(i, -1) * xi1(i, -1))) * dt / 2.;
			if (fabs(temp - phi1_ptr[i]) > ERROR_TOLERANCE)
				*bConvergence = false;
			phi1_ptr[i] = temp;
		}

		for (int i = 0; i < gridsize; ++i)
		{
			double temp = phi2(i, -1) + ((pi2(i) * alpha(i) / (psi(i) * psi(i)) + beta(i) * xi2(i)) +
					(pi2(i, -1) * alpha(i, -1) / (psi(i, - 1) * psi(i, - 1)) + beta(i, -1) * xi2(i, -1))) * dt / 2.;
			if (fabs(temp - phi2_ptr[i]) > ERROR_TOLERANCE)
				*bConvergence = false;
			phi2_ptr[i] = temp;
		}

		for (int i = 0; i < gridsize; ++i)
		{
			double temp = phi3(i, -1) + ((pi3(i) * alpha(i) / (psi(i) * psi(i)) + beta(i) * xi3(i)) +
					(pi3(i, -1) * alpha(i, -1) / (psi(i, - 1) * psi(i, - 1)) + beta(i, -1) * xi3(i, -1))) * dt / 2.;
			if (fabs(temp - phi3_ptr[i]) > ERROR_TOLERANCE)
				*bConvergence = false;
			phi3_ptr[i] = temp;
		}
	}

	void Evolve()
	{
		const size_t iteration_count = 100;

	    double* xi1_ptr =	&var_buffer_new[gridsize * VAR_XI1];
	    double* pi1_ptr =	&var_buffer_new[gridsize * VAR_PI1];
	    double* xi2_ptr =	&var_buffer_new[gridsize * VAR_XI2];
	    double* pi2_ptr =	&var_buffer_new[gridsize * VAR_PI2];
	    double* xi3_ptr =	&var_buffer_new[gridsize * VAR_XI3];
	    double* pi3_ptr = 	&var_buffer_new[gridsize * VAR_PI3];

	    bool bConvergence = false;

		for (size_t iterations = 0; iterations < iteration_count && !bConvergence; ++iterations)
		{
			bConvergence = true;

			// Hyperbolic evolution
			// KG Equation
			// Left boundary
			double error = 0.0;
			int i = 0;
			double xi1_temp = 0.0;
			double xi2_temp = 0.0;
			double xi3_temp = 0.0;
			double pi1_temp = (4*pi1(1 + i) - pi1(2 + i) - 3*pi1(i,-1) + 4*pi1(1 + i,-1) - pi1(2 + i,-1))/3.;
			double pi2_temp = (4*pi2(1 + i) - pi2(2 + i) - 3*pi2(i,-1) + 4*pi2(1 + i,-1) - pi2(2 + i,-1))/3.;
			double pi3_temp = (4*pi3(1 + i) - pi3(2 + i) - 3*pi3(i,-1) + 4*pi3(1 + i,-1) - pi3(2 + i,-1))/3.;
			error = fabs(xi1_temp - xi1_ptr[i]) + fabs(pi1_temp - pi1_ptr[i]) + fabs(xi2_temp - xi2_ptr[i]) + fabs(pi2_temp - pi2_ptr[i])
					+ fabs(xi3_temp - xi3_ptr[i]) + fabs(pi3_temp - pi3_ptr[i]);
			bConvergence = bConvergence && (error < ERROR_TOLERANCE);
			xi1_ptr[i] = xi1_temp;
			pi1_ptr[i] = pi1_temp;
			xi2_ptr[i] = xi2_temp;
			pi2_ptr[i] = pi2_temp;
			xi3_ptr[i] = xi3_temp;
			pi3_ptr[i] = pi3_temp;

			// Middle
			for (i = 1; i < gridsize - 1; ++i)
			{
				pi1_temp = (dr*dt*psi(i)*r(i)*((12*pi1(i,-1))/dt - 6*alpha(i)*DV(Power(phi1(i),2) + Power(phi2(i),2))*phi1(i)*Power(psi(i),2) -
					       (3*beta(-1 + i)*pi1(-1 + i)*Power(psi(-1 + i),4)*Power(r(-1 + i),2))/(dr*Power(psi(i),4)*Power(r(i),2)) + (3*beta(1 + i)*pi1(1 + i)*Power(psi(1 + i),4)*Power(r(1 + i),2))/(dr*Power(psi(i),4)*Power(r(i),2)) -
					       (3*alpha(-1 + i)*Power(psi(-1 + i),2)*Power(r(-1 + i),2)*xi1(-1 + i))/(dr*Power(psi(i),4)*Power(r(i),2)) + (3*alpha(1 + i)*Power(psi(1 + i),2)*Power(r(1 + i),2)*xi1(1 + i))/(dr*Power(psi(i),4)*Power(r(i),2)) +
					       (-2*beta(1 + i,-1)*pi1(i,-1)*Power(psi(i,-1),4)*Power(r(i,-1),2) - 6*dr*alpha(i,-1)*DV(Power(phi1(i),2) + Power(phi2(i),2),-1)*phi1(i,-1)*Power(psi(i,-1),6)*Power(r(i,-1),2) -
					          4*beta(i,-1)*pi1(i,-1)*Power(psi(i,-1),3)*r(i,-1)*(2*dr*psi(i,-1) + 3*(-psi(-1 + i,-1) + psi(1 + i,-1))*r(i,-1)) +
					          beta(-1 + i,-1)*(-3*pi1(-1 + i,-1)*Power(psi(-1 + i,-1),4)*Power(r(-1 + i,-1),2) + 2*pi1(i,-1)*Power(psi(i,-1),4)*Power(r(i,-1),2)) +
					          3*beta(1 + i,-1)*pi1(1 + i,-1)*Power(psi(1 + i,-1),4)*Power(r(1 + i,-1),2) - 3*alpha(-1 + i,-1)*Power(psi(-1 + i,-1),2)*Power(r(-1 + i,-1),2)*xi1(-1 + i,-1) +
					          3*alpha(1 + i,-1)*Power(psi(1 + i,-1),2)*Power(r(1 + i,-1),2)*xi1(1 + i,-1))/(dr*Power(psi(i,-1),4)*Power(r(i,-1),2))))/
					   (2.*((6*dr - dt*beta(-1 + i) + dt*beta(1 + i))*psi(i)*r(i) + 2*dt*beta(i)*(2*dr*psi(i) + 3*(-psi(-1 + i) + psi(1 + i))*r(i))));
				pi2_temp = (dr*dt*psi(i)*r(i)*((12*pi2(i,-1))/dt - 6*alpha(i)*DV(Power(phi1(i),2) + Power(phi2(i),2))*phi2(i)*Power(psi(i),2) -
					       (3*beta(-1 + i)*pi2(-1 + i)*Power(psi(-1 + i),4)*Power(r(-1 + i),2))/(dr*Power(psi(i),4)*Power(r(i),2)) + (3*beta(1 + i)*pi2(1 + i)*Power(psi(1 + i),4)*Power(r(1 + i),2))/(dr*Power(psi(i),4)*Power(r(i),2)) -
					       (3*alpha(-1 + i)*Power(psi(-1 + i),2)*Power(r(-1 + i),2)*xi2(-1 + i))/(dr*Power(psi(i),4)*Power(r(i),2)) + (3*alpha(1 + i)*Power(psi(1 + i),2)*Power(r(1 + i),2)*xi2(1 + i))/(dr*Power(psi(i),4)*Power(r(i),2)) +
					       (-2*beta(1 + i,-1)*pi2(i,-1)*Power(psi(i,-1),4)*Power(r(i,-1),2) - 6*dr*alpha(i,-1)*DV(Power(phi1(i),2) + Power(phi2(i),2),-1)*phi2(i,-1)*Power(psi(i,-1),6)*Power(r(i,-1),2) -
					          4*beta(i,-1)*pi2(i,-1)*Power(psi(i,-1),3)*r(i,-1)*(2*dr*psi(i,-1) + 3*(-psi(-1 + i,-1) + psi(1 + i,-1))*r(i,-1)) +
					          beta(-1 + i,-1)*(-3*pi2(-1 + i,-1)*Power(psi(-1 + i,-1),4)*Power(r(-1 + i,-1),2) + 2*pi2(i,-1)*Power(psi(i,-1),4)*Power(r(i,-1),2)) +
					          3*beta(1 + i,-1)*pi2(1 + i,-1)*Power(psi(1 + i,-1),4)*Power(r(1 + i,-1),2) - 3*alpha(-1 + i,-1)*Power(psi(-1 + i,-1),2)*Power(r(-1 + i,-1),2)*xi2(-1 + i,-1) +
					          3*alpha(1 + i,-1)*Power(psi(1 + i,-1),2)*Power(r(1 + i,-1),2)*xi2(1 + i,-1))/(dr*Power(psi(i,-1),4)*Power(r(i,-1),2))))/
					   (2.*((6*dr - dt*beta(-1 + i) + dt*beta(1 + i))*psi(i)*r(i) + 2*dt*beta(i)*(2*dr*psi(i) + 3*(-psi(-1 + i) + psi(1 + i))*r(i))));
				pi3_temp = (dr*dt*psi(i)*r(i)*((12*pi3(i,-1))/dt - (3*beta(-1 + i)*pi3(-1 + i)*Power(psi(-1 + i),4)*Power(r(-1 + i),2))/(dr*Power(psi(i),4)*Power(r(i),2)) +
					       (3*beta(1 + i)*pi3(1 + i)*Power(psi(1 + i),4)*Power(r(1 + i),2))/(dr*Power(psi(i),4)*Power(r(i),2)) - (3*alpha(-1 + i)*Power(psi(-1 + i),2)*Power(r(-1 + i),2)*xi3(-1 + i))/(dr*Power(psi(i),4)*Power(r(i),2)) +
					       (3*alpha(1 + i)*Power(psi(1 + i),2)*Power(r(1 + i),2)*xi3(1 + i))/(dr*Power(psi(i),4)*Power(r(i),2)) +
					       (-2*beta(1 + i,-1)*pi3(i,-1)*Power(psi(i,-1),4)*Power(r(i,-1),2) - 4*beta(i,-1)*pi3(i,-1)*Power(psi(i,-1),3)*r(i,-1)*(2*dr*psi(i,-1) + 3*(-psi(-1 + i,-1) + psi(1 + i,-1))*r(i,-1)) +
					          beta(-1 + i,-1)*(-3*pi3(-1 + i,-1)*Power(psi(-1 + i,-1),4)*Power(r(-1 + i,-1),2) + 2*pi3(i,-1)*Power(psi(i,-1),4)*Power(r(i,-1),2)) +
					          3*beta(1 + i,-1)*pi3(1 + i,-1)*Power(psi(1 + i,-1),4)*Power(r(1 + i,-1),2) - 3*alpha(-1 + i,-1)*Power(psi(-1 + i,-1),2)*Power(r(-1 + i,-1),2)*xi3(-1 + i,-1) +
					          3*alpha(1 + i,-1)*Power(psi(1 + i,-1),2)*Power(r(1 + i,-1),2)*xi3(1 + i,-1))/(dr*Power(psi(i,-1),4)*Power(r(i,-1),2))))/
					   (2.*((6*dr - dt*beta(-1 + i) + dt*beta(1 + i))*psi(i)*r(i) + 2*dt*beta(i)*(2*dr*psi(i) + 3*(-psi(-1 + i) + psi(1 + i))*r(i))));
				xi1_temp = (-((dt*alpha(-1 + i)*pi1(-1 + i))/(dr*Power(psi(-1 + i),2))) + (dt*alpha(1 + i)*pi1(1 + i))/(dr*Power(psi(1 + i),2)) - (dt*alpha(-1 + i,-1)*pi1(-1 + i,-1))/(dr*Power(psi(-1 + i,-1),2)) +
					     (dt*alpha(1 + i,-1)*pi1(1 + i,-1))/(dr*Power(psi(1 + i,-1),2)) - (dt*beta(-1 + i)*xi1(-1 + i))/dr + (dt*beta(1 + i)*xi1(1 + i))/dr - (dt*beta(-1 + i,-1)*xi1(-1 + i,-1))/dr + 4*xi1(i,-1) +
					     (dt*beta(1 + i,-1)*xi1(1 + i,-1))/dr)/4.;
				xi2_temp = (-((dt*alpha(-1 + i)*pi2(-1 + i))/(dr*Power(psi(-1 + i),2))) + (dt*alpha(1 + i)*pi2(1 + i))/(dr*Power(psi(1 + i),2)) - (dt*alpha(-1 + i,-1)*pi2(-1 + i,-1))/(dr*Power(psi(-1 + i,-1),2)) +
					     (dt*alpha(1 + i,-1)*pi2(1 + i,-1))/(dr*Power(psi(1 + i,-1),2)) - (dt*beta(-1 + i)*xi2(-1 + i))/dr + (dt*beta(1 + i)*xi2(1 + i))/dr - (dt*beta(-1 + i,-1)*xi2(-1 + i,-1))/dr + 4*xi2(i,-1) +
					     (dt*beta(1 + i,-1)*xi2(1 + i,-1))/dr)/4.;
				xi3_temp = (-((dt*alpha(-1 + i)*pi3(-1 + i))/(dr*Power(psi(-1 + i),2))) + (dt*alpha(1 + i)*pi3(1 + i))/(dr*Power(psi(1 + i),2)) - (dt*alpha(-1 + i,-1)*pi3(-1 + i,-1))/(dr*Power(psi(-1 + i,-1),2)) +
					     (dt*alpha(1 + i,-1)*pi3(1 + i,-1))/(dr*Power(psi(1 + i,-1),2)) - (dt*beta(-1 + i)*xi3(-1 + i))/dr + (dt*beta(1 + i)*xi3(1 + i))/dr - (dt*beta(-1 + i,-1)*xi3(-1 + i,-1))/dr + 4*xi3(i,-1) +
					     (dt*beta(1 + i,-1)*xi3(1 + i,-1))/dr)/4.;
				error = fabs(xi1_temp - xi1_ptr[i]) + fabs(pi1_temp - pi1_ptr[i]) + fabs(xi2_temp - xi2_ptr[i]) + fabs(pi2_temp - pi2_ptr[i])
									+ fabs(xi3_temp - xi3_ptr[i]) + fabs(pi3_temp - pi3_ptr[i]);
				bConvergence = bConvergence && (error < ERROR_TOLERANCE);
				xi1_ptr[i] = xi1_temp;
				pi1_ptr[i] = pi1_temp;
				xi2_ptr[i] = xi2_temp;
				pi2_ptr[i] = pi2_temp;
				xi3_ptr[i] = xi3_temp;
				pi3_ptr[i] = pi3_temp;
			}

			// Right boundary
			xi1_temp = (-(dt*xi1(-2 + i)) + 4*dt*xi1(-1 + i) - dt*xi1(-2 + i,-1) + 4*dt*xi1(-1 + i,-1) + 2*dr*xi1(i,-1) - 3*dt*xi1(i,-1))/(2*dr + 3*dt);
			pi1_temp = (-(dt*pi1(-2 + i)) + 4*dt*pi1(-1 + i) - dt*pi1(-2 + i,-1) + 4*dt*pi1(-1 + i,-1) + 2*dr*pi1(i,-1) - 3*dt*pi1(i,-1))/(2*dr + 3*dt);
			xi2_temp = (-(dt*xi2(-2 + i)) + 4*dt*xi2(-1 + i) - dt*xi2(-2 + i,-1) + 4*dt*xi2(-1 + i,-1) + 2*dr*xi2(i,-1) - 3*dt*xi2(i,-1))/(2*dr + 3*dt);
			pi2_temp = (-(dt*pi2(-2 + i)) + 4*dt*pi2(-1 + i) - dt*pi2(-2 + i,-1) + 4*dt*pi2(-1 + i,-1) + 2*dr*pi2(i,-1) - 3*dt*pi2(i,-1))/(2*dr + 3*dt);
			xi3_temp = (-(dt*xi3(-2 + i)) + 4*dt*xi3(-1 + i) - dt*xi3(-2 + i,-1) + 4*dt*xi3(-1 + i,-1) + 2*dr*xi3(i,-1) - 3*dt*xi3(i,-1))/(2*dr + 3*dt);
			pi3_temp = (-(dt*pi3(-2 + i)) + 4*dt*pi3(-1 + i) - dt*pi3(-2 + i,-1) + 4*dt*pi3(-1 + i,-1) + 2*dr*pi3(i,-1) - 3*dt*pi3(i,-1))/(2*dr + 3*dt);
			error = fabs(xi1_temp - xi1_ptr[i]) + fabs(pi1_temp - pi1_ptr[i]) + fabs(xi2_temp - xi2_ptr[i]) + fabs(pi2_temp - pi2_ptr[i])
					+ fabs(xi3_temp - xi3_ptr[i]) + fabs(pi3_temp - pi3_ptr[i]);
			bConvergence = bConvergence && (error < ERROR_TOLERANCE);
			xi1_ptr[i] = xi1_temp;
			pi1_ptr[i] = pi1_temp;
			xi2_ptr[i] = xi2_temp;
			pi2_ptr[i] = pi2_temp;
			xi3_ptr[i] = xi3_temp;
			pi3_ptr[i] = pi3_temp;

			// Check residuals and make sure that the solution stays within truncation error
			// Elliptic equations
			bool bResidualSmall = true;

			if (!b_use_minkowski)
			{
				IteratePsi(&bResidualSmall);
				bConvergence = bConvergence && bResidualSmall;

				IterateBeta(&bResidualSmall);
				bConvergence = bConvergence && bResidualSmall;

				IterateAlpha(&bResidualSmall);
				bConvergence = bConvergence && bResidualSmall;
			}

			UpdatePhi(&bResidualSmall);
			bConvergence = bConvergence && bResidualSmall;

			//OutputState_Debug(iterations);
		}

		// Check residual
		if (time % output_time_skip == 0 && var_output_mask[SPARSE_MASK_RESIDUAL_INDEX])
		{
			double* psi_residual = new double[gridsize];
			char name[NAME_MAX_LENGTH];
			strncpy(name, psi_residual_name.c_str(), NAME_MAX_LENGTH);
			ComputeResidualPsiEvolution(psi_residual);
			gft_out_bbox(name, time * dt, &gridsize, 1, bbox, psi_residual);
			delete [] psi_residual;
		}
	}

	void EvolveStep()
	{
		Evolve();
		memcpy(var_buffer, var_buffer_new, sizeof(double) * gridsize * var_count);
	}

	void ComputeMassAspect(double mass_aspect[])
	{
	    int i = 0;

	    mass_aspect[i] = (Power(-beta(i) + (r(i)*(-3*beta(i) + 4*beta(1 + i) - beta(2 + i)))/(2.*dr),2)*Power(psi(i),6)*r(i))/(18.*Power(alpha(i),2)) -
	    		   ((-3*psi(i) + 4*psi(1 + i) - psi(2 + i))*Power(r(i),2)*(-3*psi(i)*r(i) + 4*psi(1 + i)*r(1 + i) - psi(2 + i)*r(2 + i)))/(2.*Power(dr,2));

	    for (; i < gridsize - 1; ++i)
	    {
	        mass_aspect[i] = (Power(-beta(i) + (r(i)*(-beta(-1 + i) + beta(1 + i)))/(2.*dr),2)*Power(psi(i),6)*r(i))/(18.*Power(alpha(i),2)) -
	        		   ((-psi(-1 + i) + psi(1 + i))*Power(r(i),2)*(-(psi(-1 + i)*r(-1 + i)) + psi(1 + i)*r(1 + i)))/(2.*Power(dr,2));
	    }

	    mass_aspect[i] = (Power(-beta(i) + (r(i)*(beta(-2 + i) - 4*beta(-1 + i) + 3*beta(i)))/(2.*dr),2)*Power(psi(i),6)*r(i))/(18.*Power(alpha(i),2)) -
		   ((psi(-2 + i) - 4*psi(-1 + i) + 3*psi(i))*Power(r(i),2)*(psi(-2 + i)*r(-2 + i) - 4*psi(-1 + i)*r(-1 + i) + 3*psi(i)*r(i)))/(2.*Power(dr,2));
	}

	void Integrate(double output[], double input[], const size_t length)
	{
		size_t i = 0;

		output[i] = 0;

		for (; i < length; ++i)
			output[i] = output[i - 1] + (input[i - 1] + input[i]) * dr / 2.;
	}

	void OutputState()
	{
		char name[NAME_MAX_LENGTH];

		//ComputeMassAspect(mass_aspect);
		//gft_out_bbox(mass_aspect_name, time, &gridsize, 1, bbox, mass_aspect);

		//Integrate(field, &var_buffer[VAR_XI], gridsize);
		//gft_out_bbox(field_name, time, &gridsize, 1, bbox, field);

		for (int i = 0; i < var_count && i < OUTPUT_COUNT; ++i)
		{
			if (var_output_mask[i])
			{
				strncpy(name, var_names[i].c_str(), NAME_MAX_LENGTH);
				gft_out_bbox(name, time * dt, &gridsize, 1, bbox, &var_buffer[gridsize * i]);
			}
		}
		//gft_out_bbox(pi_name2, time, &gridsize, 1, bbox, &var_buffer[gridsize]);
		//gft_out_bbox(mass_aspect_name, time, &gridsize, 1, bbox, mass_aspect);
	}

	void OutputState_Debug(const double time_value)
	{
		char name[NAME_MAX_LENGTH];

		for (int i = 0; i < var_count && i < OUTPUT_COUNT; ++i)
		{
			strncpy(name, var_names_debug[i].c_str(), NAME_MAX_LENGTH);
			gft_out_bbox(name, time_value, &gridsize, 1, bbox, &var_buffer_new[gridsize * i]);
		}
	}

public:
	Application(char* pParamFile) : gridsize(129), Nt(129), var_count(VAR_COUNT),
		b_use_minkowski(0), init_data_type(DATA_TYPE_LEFT_MOVING_WAVE_PACKET1),
		output_time_skip(1), sparse_output_level(0), lambda(0.5), phi_amp(1.0), phi_delta(1.0), phi_r0(0.0),
		phi_perturb_amp(0.0), phi_perturb_delta(1.0), phi_perturb_r0(0.0), radius(100.0),
		event_horizon_radius(2.0), mass(1.0), init_data_omega(0.0), pFile(pParamFile), dr(0.0), dt(0.0),
		phi_init(nullptr), var_buffer(nullptr), var_buffer_new(nullptr), mass_aspect(nullptr),
		field(nullptr), time(0)
	{
		memset(output_tag, 0, sizeof(output_tag));
		memcpy(var_output_mask, var_output_masks[0], SPARSE_MASK_SIZE);
	}

	~Application()
	{
		delete [] phi_init;
		delete [] var_buffer;
		delete [] var_buffer_new;
		delete [] mass_aspect;
		delete [] field;
	}

	void Run()
	{
		LoadParameters();
		InitVariables();
		OutputState();

		cout << "Integrating equation " << Nt << " steps..." << endl;

		for (time = 1; time <= Nt; ++time)
		{
			cout << "time = " << time << endl;
			EvolveStep();
			if (time % output_time_skip == 0)
				OutputState();
		}

		cout << "Finished!" << endl;
	}
};

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        cout << "Incorrect number of parameters!" << endl;
        return EXIT_FAILURE;
    }

    Application app(argv[1]);
    app.Run();

    return EXIT_SUCCESS;
}
