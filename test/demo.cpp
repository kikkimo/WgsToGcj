#include "WGS2GCJ.h"
#include <iostream>
#include <vector>
#include <omp.h>
#include <assert.h>
#include <time.h>
#include <conio.h>

/**
 *  \brief Check if the point in china roughly
 *
 *  \param [in] lon: longitude of geodetic point, unit degree
 *  \param [in] lat: latitude of geodetic point,unit degree
 *  \return Returns true if point in china
 *  \time 15:22:01 2020/06/12
 */
bool OutOfChina(const double& lon, const double& lat) {
	return  !(72.004 <= lon && lon <= 137.8347 &&
		0.8293 <= lat && lat <= 55.8271);
}


void TestCase1()
{
	double lonWgs = 116.39123343289631;
	double latWgs = 39.9072885060602;
	std::cout.precision(16);
	std::cout << "WGS84 Point: (" << lonWgs << ", " << latWgs << ")\n";
	auto [lonGcj, latGcj] = Wgs2Gcj(lonWgs, latWgs);
	std::cout << "GCJ-02 Point [Wgs2Gcj]: (" << lonGcj << ", " << latGcj << ")\n";
	//simple linear iteration 
	auto [lonWgsNS, latWgsNS] = Gcj2Wgs_SimpleIteration(lonGcj, latGcj);
	std::cout << "WGS84 Point [simple linear iteration]: (" << lonWgsNS << ", " << latWgsNS << ")\n";

	//numerically differentiated cost function
	auto [lonWgsND, latWgsND] = Gcj2Wgs_NumbericDiff(lonGcj, latGcj);
	std::cout << "WGS84 Point [numerically derivation]: (" << lonWgsND << ", " << latWgsND << ")\n";

	//analytical differentiated cost function
	auto [lonWgsNAn, latWgsNAn] = Gcj2Wgs_AnalyticDiff(lonGcj, latGcj);
	std::cout << "WGS84 Point [analytical derivation]: (" << lonWgsNAn << ", " << latWgsNAn << ")\n";

#ifdef _USE_CERES
	//auto differentiable, use ceres
	auto [lonWgsNA, latWgsNA] = Gcj2Wgs_AutoDiff(lonGcj, latGcj);
	std::cout << "WGS84 Point [auto differentiable]: (" << lonWgsNA << ", " << latWgsNA << ")\n";
#endif
}

void TestCase2()
{
	//range of china
	std::vector<double> wgslonRange, wgslatRange;
	for (double wgslon = 72.004; wgslon <= 137.8347; wgslon += 0.01)
	{
		for (double wgslat = 0.8293; wgslat <= 55.8271; wgslat += 0.01)
		{
			if (OutOfChina(wgslon, wgslat))
				continue;
			wgslonRange.push_back(wgslon);
			wgslatRange.push_back(wgslat);
		}
	}
	size_t size = wgslonRange.size();
	assert(size == wgslatRange.size());
	double lfXMin = DBL_MAX;
	double lfXMax = 0.0l;
	double lfXAverage = 0.0l;
	double lfYMin = DBL_MAX;
	double lfYMax = 0.0l;
	double lfYAverage = 0.0l;
	int counts = 0;
	long t_begin = clock();
	#pragma omp parallel for
	for (long i = 0; i < size; ++i)
	{
		auto [gcj02lon, gcj02lat] = Wgs2Gcj(wgslonRange[i], wgslatRange[i]);
		auto [wgslonN, wgslatN] = Gcj2Wgs_SimpleIteration(gcj02lon, gcj02lat);
		double lfXError = wgslonN - wgslonRange[i];
		double lfYError = wgslatN - wgslatRange[i];
		double lfXE = fabs(lfXError);
		double lfYE = fabs(lfYError);
		#pragma omp critical
		{
			lfXMin = fabs(lfXMin) < lfXE ? lfXMin : lfXError;
			lfXMax = fabs(lfXMax) > lfXE ? lfXMax : lfXError;
			lfYMin = fabs(lfYMin) < lfYE ? lfYMin : lfYError;
			lfYMax = fabs(lfYMax) > lfYE ? lfYMax : lfYError;
			lfXAverage += lfXE * lfXE;
			lfYAverage += lfYE * lfYE;
			counts++;
		}
	}
	lfXAverage /= counts;
	lfYAverage /= counts;
	lfXAverage = sqrt(lfXAverage);
	lfYAverage = sqrt(lfYAverage);
	long t_end = clock();
	std::cout << "| SimpleIteration: " << counts << " case finished!" << std::endl;
	std::cout << "| minimum X error: " << lfXMin << std::endl;
	std::cout << "| maximum X error: " << lfXMax << std::endl;
	std::cout << "| average X error: " << lfXAverage << std::endl;
	std::cout << "| minimum Y error: " << lfYMin << std::endl;
	std::cout << "| maximum Y error: " << lfYMax << std::endl;
	std::cout << "| average Y error: " << lfYAverage << std::endl;
	std::cout << "| time use: " << t_end - t_begin << " ms" << std::endl;

	lfXMin = DBL_MAX;
	lfXMax = 0.0l;
	lfXAverage = 0.0l;
	lfYMin = DBL_MAX;
	lfYMax = 0.0l;
	lfYAverage = 0.0l;
	counts = 0;
	t_begin = clock();
	#pragma omp parallel for 
	for (long i = 0; i < size; ++i)
	{
		auto [gcj02lon, gcj02lat] = Wgs2Gcj(wgslonRange[i], wgslatRange[i]);
		auto [wgslonN, wgslatN] = Gcj2Wgs_AnalyticDiff(gcj02lon, gcj02lat);
		double lfXError = wgslonN - wgslonRange[i];
		double lfYError = wgslatN - wgslatRange[i];
		double lfXE = fabs(lfXError);
		double lfYE = fabs(lfYError);
		#pragma omp critical
		{
			lfXMin = fabs(lfXMin) < lfXE ? lfXMin : lfXError;
			lfXMax = fabs(lfXMax) > lfXE ? lfXMax : lfXError;
			lfYMin = fabs(lfYMin) < lfYE ? lfYMin : lfYError;
			lfYMax = fabs(lfYMax) > lfYE ? lfYMax : lfYError;
			lfXAverage += lfXE * lfXE;
			lfYAverage += lfYE * lfYE;
			counts++;
		}
	}
	lfXAverage /= counts;
	lfYAverage /= counts;
	lfXAverage = sqrt(lfXAverage);
	lfYAverage = sqrt(lfYAverage);
	t_end = clock();
	std::cout << "| analytical derivation: " << counts << " case finished!" << std::endl;
	std::cout << "| minimum X error: " << lfXMin << std::endl;
	std::cout << "| maximum X error: " << lfXMax << std::endl;
	std::cout << "| average X error: " << lfXAverage << std::endl;
	std::cout << "| minimum Y error: " << lfYMin << std::endl;
	std::cout << "| maximum Y error: " << lfYMax << std::endl;
	std::cout << "| average Y error: " << lfYAverage << std::endl;
	std::cout << "| time use: " << t_end - t_begin << " ms" << std::endl;

	lfXMin = DBL_MAX;
	lfXMax = 0.0l;
	lfXAverage = 0.0l;
	lfYMin = DBL_MAX;
	lfYMax = 0.0l;
	lfYAverage = 0.0l;
	counts = 0;
	t_begin = clock();
	#pragma omp parallel for
	for (long i = 0; i < size; ++i)
	{
		auto [gcj02lon, gcj02lat] = Wgs2Gcj(wgslonRange[i], wgslatRange[i]);
		auto [wgslonN, wgslatN] = Gcj2Wgs_NumbericDiff(gcj02lon, gcj02lat);
		double lfXError = wgslonN - wgslonRange[i];
		double lfYError = wgslatN - wgslatRange[i];
		double lfXE = fabs(lfXError);
		double lfYE = fabs(lfYError);
		#pragma omp critical
		{
			lfXMin = fabs(lfXMin) < lfXE ? lfXMin : lfXError;
			lfXMax = fabs(lfXMax) > lfXE ? lfXMax : lfXError;
			lfYMin = fabs(lfYMin) < lfYE ? lfYMin : lfYError;
			lfYMax = fabs(lfYMax) > lfYE ? lfYMax : lfYError;
			lfXAverage += lfXE * lfXE;
			lfYAverage += lfYE * lfYE;
			counts++;
		}
	}
	lfXAverage /= counts;
	lfYAverage /= counts;
	lfXAverage = sqrt(lfXAverage);
	lfYAverage = sqrt(lfYAverage);
	t_end = clock();
	std::cout << "| numerically differentiated: " << counts << " case finished!" << std::endl;
	std::cout << "| minimum X error: " << lfXMin << std::endl;
	std::cout << "| maximum X error: " << lfXMax << std::endl;
	std::cout << "| average X error: " << lfXAverage << std::endl;
	std::cout << "| minimum Y error: " << lfYMin << std::endl;
	std::cout << "| maximum Y error: " << lfYMax << std::endl;
	std::cout << "| average Y error: " << lfYAverage << std::endl;
	std::cout << "| time use: " << t_end - t_begin << " ms" << std::endl;
}

int main()
{
	TestCase1();
	TestCase2();
	getch();
	return 0;
}