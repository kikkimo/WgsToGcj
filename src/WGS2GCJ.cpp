#include <iostream>
#include <cmath>

#include "WGS2GCJ.h"

#ifdef _USE_CERES
#include <ceres/ceres.h>
#endif

//Beijing54 Geodetic coordinate system (Krasovsky reference ellipsoid)
constexpr double kKRASOVSKY_A = 6378245.0;				// equatorial radius [unit: meter]
constexpr double kKRASOVSKY_B = 6356863.0187730473;	// polar radius
constexpr double kKRASOVSKY_ECCSQ = 6.6934216229659332e-3; // first eccentricity squared
constexpr double kKRASOVSKY_ECC2SQ = 6.7385254146834811e-3; // second eccentricity squared
constexpr double PI = 3.14159265358979323846;   //π

constexpr double kDEG2RAD = PI / 180.0;
constexpr double kRAD2DEG = 180.0 / PI;

/**
 *  \brief Angle unit transform, degree to radian
 *
 *  \param [in] deg: angle in degrees
 *  \return Returns angle in radians
 *  \time 15:21:22 2020/06/12
 */
constexpr inline double Deg2Rad(const double deg) {
	return deg * kDEG2RAD;
}

/**
 *  \brief Angle unit transform, radian to degree
 *
 *  \param [in] rad: angle in radians
 *  \return Returns angle in degrees
 *  \time 15:21:01 2020/06/12
 */
constexpr inline double Rad2Deg(const double rad) {
	return rad * kRAD2DEG;
}
    
/**
 *  \brief Get geodetic offset used by GCJ-02 
 *
 *  \param [in] wgs84lon: longitude in WGS84 coordinate system [unit:degree] 
 *  \param [in] wgs84lat: latitude in WGS84 coordinate system [unit:degree] 
 *  \return Returns a pair of geodetic offset used by GCJ-02
 *  \time 15:28:33 2020/06/12
 */
std::pair<double,double> GetGeodeticOffset(const double& wgs84lon,const double& wgs84lat)
{
    //get geodetic offset relative to 'center china'
    double lon0 = wgs84lon - 105.0;
    double lat0 = wgs84lat - 35.0;

    //generate an pair offset roughly in meters
	double lon1 = 300.0 + lon0 + 2.0 * lat0 + 0.1 * lon0 * lon0 + 0.1 * lon0 * lat0 + 0.1 * sqrt(fabs(lon0));
	lon1 = lon1 + (20.0 * sin(6.0 * lon0 * PI) + 20.0 * sin(2.0 * lon0 * PI)) * 2.0 / 3.0;
	lon1 = lon1 + (20.0 * sin(lon0 * PI) + 40.0 * sin(lon0 / 3.0 * PI)) * 2.0 / 3.0;
	lon1 = lon1 + (150.0 * sin(lon0 / 12.0 * PI) + 300.0 * sin(lon0 * PI / 30.0)) * 2.0 / 3.0;
    double lat1 = -100.0 + 2.0 * lon0 + 3.0 * lat0 + 0.2 * lat0 * lat0 + 0.1 * lon0 * lat0 + 0.2 * sqrt(fabs(lon0));
    lat1 = lat1 + (20.0 * sin(6.0 * lon0 * PI) + 20.0 * sin(2.0 * lon0 * PI)) * 2.0 / 3.0;
    lat1 = lat1 + (20.0 * sin(lat0 * PI) + 40.0 * sin(lat0 / 3.0 * PI)) * 2.0 / 3.0;
    lat1 = lat1 + (160.0 * sin(lat0 / 12.0 * PI) + 320.0 * sin(lat0 * PI / 30.0)) * 2.0 / 3.0;

    //latitude in radian
    double B = Deg2Rad(wgs84lat);
    double sinB = sin(B), cosB = cos(B);
    double W = sqrt(1 - kKRASOVSKY_ECCSQ * sinB * sinB);
    double N = kKRASOVSKY_A / W;

    //geodetic offset used by GCJ-02
    double lon2 = Rad2Deg(lon1 / (N * cosB));
    double lat2 = Rad2Deg(lat1 * W * W / (N * (1 - kKRASOVSKY_ECCSQ)));
    return {lon2, lat2};
}

/**
 *  \brief Covert geodetic coordinate in WGS84 coordinate system to geodetic coordinate 
 *         in GCJ-02 coordinate system
 *
 *  \param [in] wgs84lon: longitude in WGS84 coordinate system [unit:degree]
 *  \param [in] wgs84lat: latitude in WGS84 coordinate system [unit:degree]
 *  \return Returns geodetic coordinate in GCJ-02 coordinate system
 *  \time 15:47:38 2020/06/12
 */
std::pair<double, double> Wgs2Gcj(const double& wgs84lon, const double& wgs84lat)
{
    auto [dlon, dlat] = GetGeodeticOffset(wgs84lon, wgs84lat);
    double gcj02lon = wgs84lon + dlon;
    double gcj02lat = wgs84lat + dlat;
    return { gcj02lon, gcj02lat };
}

/**
 *  \brief Covert geodetic coordinate in GCJ-02 coordinate system to geodetic coordinate
 *         in WGS84 coordinate system
 *
 *  \param [in] gcj02lon: longitude in GCJ-02 coordinate system [unit:degree]
 *  \param [in] gcj02lat: latitude in GCJ-02 coordinate system [unit:degree]
 *  \return Returns geodetic coordinate in WGS84 coordinate system
 *  \remark simple linear iteration
 *	\detail 
 *  \time 15:51:13 2020/06/12
 */
std::pair<double, double> Gcj2Wgs_SimpleIteration(const double& gcj02lon,
    const double& gcj02lat)
{
    auto [lon0, lat0] = Wgs2Gcj(gcj02lon, gcj02lat);
    int iterCounts = 0;
    while (++iterCounts < 1000)
    {
        auto[lon1, lat1] = Wgs2Gcj(lon0, lat0);
		double dlon = lon1 - gcj02lon;
		double dlat = lat1 - gcj02lat;
        lon1 = lon0 - dlon;
        lat1 = lat0 - dlat;
        //1.0e-9 degree corresponding to 0.1mm
        if (fabs(dlon) < 1.0e-9 && fabs(dlat) < 1.0e-9)
            break;
        lon0 = lon1;
        lat0 = lat1;
    }
    return {lon0 , lat0};
}

#ifdef _USE_CERES  //_USE_CERES

class AutoDiffCostFunc
{
public:
    AutoDiffCostFunc(const double lon, const double lat) :
        mlonGcj(lon), mlatGcj(lat) {}
	template <typename T>
    bool operator() (const T* const lonWgs, const T* const latWgs, T* residuals) const
    {
		//get geodetic offset relative to 'center china'
		T lon0 = lonWgs[0] - T(105.0);
		T lat0 = latWgs[0] - T(35.0);

		//generate an pair offset roughly in meters
		T lon1 = T(300.0) + lon0 + T(2.0) * lat0 + T(0.1) * lon0 * lon0 + T(0.1) * lon0 * lat0
            + T(0.1) * ceres::sqrt(ceres::abs(lon0));
		lon1 = lon1 + (T(20.0) * ceres::sin(T(6.0) * lon0 * T(PI)) + T(20.0) * ceres::sin(T(2.0) * lon0 * T(PI))) * T(2.0) / T(3.0);
		lon1 = lon1 + (T(20.0) * ceres::sin(lon0 * T(PI)) + T(40.0) * ceres::sin(lon0 / T(3.0) * T(PI))) * T(2.0) / T(3.0);
		lon1 = lon1 + (T(150.0) * ceres::sin(lon0 / T(12.0) * T(PI)) + T(300.0) * ceres::sin(lon0 * T(PI) / T(30.0))) * T(2.0) / T(3.0);
		T lat1 = T(-100.0) + T(2.0) * lon0 + T(3.0) * lat0 + T(0.2) * lat0 * lat0 + T(0.1) * lon0 * lat0
            + T(0.2) * ceres::sqrt(ceres::abs(lon0));
		lat1 = lat1 + (T(20.0) * ceres::sin(T(6.0) * lon0 * T(PI)) + T(20.0) * ceres::sin(T(2.0) * lon0 * T(PI))) * T(2.0) / T(3.0);
		lat1 = lat1 + (T(20.0) * ceres::sin(lat0 * T(PI)) + T(40.0) * ceres::sin(lat0 / T(3.0) * T(PI))) * T(2.0) / T(3.0);
		lat1 = lat1 + (T(160.0) * ceres::sin(lat0 / T(12.0) * T(PI)) + T(320.0) * ceres::sin(lat0 * T(PI) / T(30.0))) * T(2.0) / T(3.0);

		//latitude in radian
		T B = latWgs[0] * T(kDEG2RAD);
		T sinB = ceres::sin(B), cosB = ceres::cos(B);
        T W = ceres::sqrt(T(1) - T(kKRASOVSKY_ECCSQ) * sinB * sinB);
        T N = T(kKRASOVSKY_A) / W;

		//geodetic offset used by GCJ-02
		T lon2 = T(kRAD2DEG) * lon1 / (N * cosB);
		T lat2 = T(kRAD2DEG) * (lat1 * W * W / (N * (1 - kKRASOVSKY_ECCSQ)));

        //residuals
        residuals[0] = lonWgs[0] + lon2 - mlonGcj;
        residuals[1] = latWgs[0] + lat2 - mlatGcj;
        return true;
    }

private:
    double mlonGcj;
    double mlatGcj;
};

/**
 *  \brief Covert geodetic coordinate in GCJ-02 coordinate system to geodetic coordinate
 *         in WGS84 coordinate system
 *
 *  \param [in] gcj02lon: longitude in GCJ-02 coordinate system [unit:degree]
 *  \param [in] gcj02lat: latitude in GCJ-02 coordinate system [unit:degree]
 *  \return Returns geodetic coordinate in WGS84 coordinate system
 *  \remark the encryption formula is known and use an auto-differentiable cost function
 *  \time 15:51:13 2020/06/12
 */
std::pair<double, double> Gcj2Wgs_AutoDiff(const double& gcj02lon,
    const double& gcj02lat)
{
	ceres::Problem * poProblem = new ceres::Problem;
    AutoDiffCostFunc* pCostFunc = new AutoDiffCostFunc(gcj02lon,gcj02lat);

    double wgslon =  gcj02lon , wgslat =  gcj02lat;
	poProblem->AddResidualBlock(new ceres::AutoDiffCostFunction<AutoDiffCostFunc, 2, 1, 1>(pCostFunc),
        nullptr,
        &wgslon,
        &wgslat);

	ceres::Solver::Options options;
	options.max_num_iterations = 30;
	options.linear_solver_type = ceres::DENSE_QR;
	options.minimizer_progress_to_stdout = false;
	options.gradient_tolerance = 1e-16;
	options.function_tolerance = 1e-12;
	options.parameter_tolerance = 1e-14;
	ceres::Solver::Summary summary;
	ceres::Solve(options, poProblem, &summary);
	delete poProblem;		//auto free memory of cost function "pCostFunc"
    return { wgslon, wgslat };
}
#endif // _USE_CERES

/**
 *  \brief Calculate the partial derivatives with respect to estimated longitude in WGS84
 *
 *  \param [in] wgs84lon: estimated longitude in WGS84 coordinate system [unit:degree]
 *  \param [in] wgs84lat: estimated latitude in WGS84 coordinate system [unit:degree]
 *  \param [in] dlon: delta longitude (close to zero) in WGS84 coordinate system [unit:degree], 
 *  \return Returns partial derivatives with respect to estimated longitude in WGS84
 *  \time 20:26:16 2020/06/13
 */
std::pair<double, double> GetPartialDerivative_Lon(const double& wgs84lon, const double& wgs84lat, const double& dlon)
{
	double lonBk = wgs84lon + dlon;
	double lonFw = wgs84lon - dlon;
	auto [gcjlonBk, gcjlatBk] = Wgs2Gcj(lonBk, wgs84lat);
	auto [gcjlonFw, gcjlatFw] = Wgs2Gcj(lonFw, wgs84lat);
	double dlongcj_dlonwgs = (gcjlonBk - gcjlonFw) / (dlon * 2.0);
	double dlatgcj_dlonwgs = (gcjlatBk - gcjlatFw) / (dlon * 2.0);
	return { dlongcj_dlonwgs , dlatgcj_dlonwgs };
}

/**
 *  \brief Calculate the partial derivatives with respect to estimated latitude in WGS84
 *
 *  \param [in] wgs84lon: estimated longitude in WGS84 coordinate system [unit:degree]
 *  \param [in] wgs84lat: estimated latitude in WGS84 coordinate system [unit:degree]
 *  \param [in] dlat: delta latitude (close to zero) in WGS84 coordinate system [unit:degree],
 *  \return Returns partial derivatives with respect to estimated latitude in WGS84
 *  \time 20:26:25 2020/06/13
 */
std::pair<double, double> GetPartialDerivative_Lat(const double& wgs84lon, const double& wgs84lat, const double& dlat)
{
	double latBk = wgs84lat + dlat;
	double latFw = wgs84lat - dlat;
	auto [gcjlonBk, gcjlatBk] = Wgs2Gcj(wgs84lon, latBk);
	auto [gcjlonFw, gcjlatFw] = Wgs2Gcj(wgs84lon, latFw);
	double dlongcj_dlatwgs = (gcjlonBk - gcjlonFw) / (dlat * 2.0);
	double dlatgcj_dlatwgs = (gcjlatBk - gcjlatFw) / (dlat * 2.0);
	return { dlongcj_dlatwgs , dlatgcj_dlatwgs };
}

/**
 *  \brief Covert geodetic coordinate in GCJ-02 coordinate system to geodetic coordinate
 *         in WGS84 coordinate system
 *
 *  \param [in] gcj02lon: longitude in GCJ-02 coordinate system [unit:degree]
 *  \param [in] gcj02lat: latitude in GCJ-02 coordinate system [unit:degree]
 *  \return Returns geodetic coordinate in WGS84 coordinate system
 *  \remark the encryption formula is unknown but we can covert point in WGS84 to point
 *          in GCJ-02 with an API,then use the numerical derivation method to solve the 
 *          problem
 *	\detail Assuming the encryption formula is
 *
 *			gcj02lon = Wgs2Gcj(wgs84lon, wgs84lat)
 *			gcj02lat = Wgs2Gcj(wgs84lon, wgs84lat)
 *
 *	 In the rectification process, (wgs84lon, wgs84lat) are unknown items. Obviously,
 *   this is a system of nonlinear equations.
 *
 *   The linear formed error functions of forward intersection come from
 *   consideration of a Taylor series expansion.
 *           V = AX - b
 *    here:
 *    V: The residuals of the observed values
 *    A: The jacobian matrix
 *    X: The modification of the unknown items
 *    b: The constant terms of the error functions
 *
 *    Then the error functions written in vector form are:
 *    | V_lon | = | dlongcj_dlonwgs  dlongcj_dlatwgs |  | d_lonwgs | - | l_lon |
 *    | V_lat | = |         0        dlatgcj_dlatwgs |  | d_latwgs | - | l_lat |
 *    here:
 *    l_lon = longcj - longcj'                 // the modification of longcj
 *    l_lat = latgcj - latgcj'                 // the modification of latgcj
 *    longcj : the observed longitude in GCJ-02
 *    latgcj : the observed latitude in GCJ-02
 *    longcj' = Wgs2Gcj(wgs84lon',wgs84lat')    // estimated longitude in GCJ-02
 *    latgcj' = Wgs2Gcj(wgs84lon',wgs84lat')    // estimated latitude in GCJ-02
 *    wgs84lon' : estimated longitude in WGS84
 *    wgs84lat' : estimated latitude in WGS84
 *    d_lonwgs : unknown items
 *    d_latwgs : unknown items
 *    wgs84lon = wgs84lon' + d_lonwgs                           // update
 *    wgs84lat = wgs84lat' + d_latwgs
 *
 *	  let V = [V_lon V_lat]T = 0, then
 *	  d_latwgs = (l_lon * dlatgcj_dlonwgs - l_lat * dlongcj_dlonwgs) /
 *			(dlongcj_dlatwgs * dlatgcj_dlonwgs - dlatgcj_dlatwgs * dlongcj_dlonwgs)
 *	  d_lonwgs = (l_lon - dlongcj_dlatwgs * d_latwgs) / dlongcj_dlonwgs
 *
 *    This iterative procedure is repeated until X= [d_lonwgs d_latwgs]T are
 *    sufficiently small.
 *  \time 17:42:01 2020/06/12
 */
std::pair<double, double> Gcj2Wgs_NumbericDiff(const double& gcj02lon,
	const double& gcj02lat)
{
	double wgs84lon = gcj02lon, wgs84lat = gcj02lat;
	int nIterCount = 0;
	double tol = 1e-9;
	while (++nIterCount < 1000)
	{
		auto [dlongcj_dlonwgs, dlatgcj_dlonwgs] = GetPartialDerivative_Lon(wgs84lon, wgs84lat, tol);
		auto [dlongcj_dlatwgs, dlatgcj_dlatwgs] = GetPartialDerivative_Lat(wgs84lon, wgs84lat, tol);

		auto [gcj02lonEst, gcj02latEst] = Wgs2Gcj(wgs84lon, wgs84lat);
		double l_lon = gcj02lon - gcj02lonEst;
		double l_lat = gcj02lat - gcj02latEst;
		double d_latwgs = (l_lon * dlatgcj_dlonwgs - l_lat * dlongcj_dlonwgs) /
			(dlongcj_dlatwgs * dlatgcj_dlonwgs - dlatgcj_dlatwgs * dlongcj_dlonwgs);
		double d_lonwgs = (l_lon - dlongcj_dlatwgs * d_latwgs) / dlongcj_dlonwgs;

		if (fabs(d_latwgs) < tol && fabs(d_lonwgs) < tol)
			break;
		wgs84lon = wgs84lon + d_lonwgs;
		wgs84lat = wgs84lat + d_latwgs;
	}
	return { wgs84lon, wgs84lat };
}

/**
 *  \brief Covert geodetic coordinate in GCJ-02 coordinate system to geodetic coordinate
 *         in WGS84 coordinate system
 *
 *  \param [in] gcj02lon: longitude in GCJ-02 coordinate system [unit:degree]
 *  \param [in] gcj02lat: latitude in GCJ-02 coordinate system [unit:degree]
 *  \return Returns geodetic coordinate in WGS84 coordinate system
 *  \remark The encryption formula is known,and use the analytical derivation method to 
 *			solve the problem with high precision.
 *	\detail Assuming the encryption formula is 
 *
 *			gcj02lon = Wgs2Gcj(wgs84lon, wgs84lat)
 *			gcj02lat = Wgs2Gcj(wgs84lon, wgs84lat)
 *
 *	 In the rectification process, (wgs84lon, wgs84lat) are unknown items. Obviously,
 *   this is a system of nonlinear equations.
 *
 *   The linear formed error functions of forward intersection come from
 *   consideration of a Taylor series expansion.
 *           V = AX - b
 *    here:
 *    V: The residuals of the observed values
 *    A: The jacobian matrix
 *    X: The modification of the unknown items
 *    b: The constant terms of the error functions
 *
 *    Then the error functions written in vector form are:
 *    | V_lon | = | dlongcj_dlonwgs  dlongcj_dlatwgs |  | d_lonwgs | - | l_lon |
 *    | V_lat | = | dlatgcj_dlonwgs  dlatgcj_dlatwgs |  | d_latwgs | - | l_lat |
 *    here:
 *    l_lon = longcj - longcj'                 // the modification of longcj
 *    l_lat = latgcj - latgcj'                 // the modification of latgcj
 *    longcj : the observed longitude in GCJ-02
 *    latgcj : the observed latitude in GCJ-02
 *    longcj' = Wgs2Gcj(wgs84lon',wgs84lat')    // estimated longitude in GCJ-02
 *    latgcj' = Wgs2Gcj(wgs84lon',wgs84lat')    // estimated latitude in GCJ-02
 *    wgs84lon' : estimated longitude in WGS84
 *    wgs84lat' : estimated latitude in WGS84
 *    d_lonwgs : unknown items
 *    d_latwgs : unknown items
 *    wgs84lon = wgs84lon' + d_lonwgs                           // update
 *    wgs84lat = wgs84lat' + d_latwgs
 *
 *	  let V = [V_lon V_lat]T = 0, then
 *	  d_latwgs = (l_lon * dlatgcj_dlonwgs - l_lat * dlongcj_dlonwgs) /
 *			(dlongcj_dlatwgs * dlatgcj_dlonwgs - dlatgcj_dlatwgs * dlongcj_dlonwgs)
 *	  d_lonwgs = (l_lon - dlongcj_dlatwgs * d_latwgs) / dlongcj_dlonwgs
 *
 *    This iterative procedure is repeated until X= [d_lonwgs d_latwgs]T are
 *    sufficiently small.
 *  \time 01:54:46 2020/06/13
 */
std::pair<double, double> Gcj2Wgs_AnalyticDiff(const double& gcj02lon,
	const double& gcj02lat)
{
	double wgs84lon = gcj02lon, wgs84lat = gcj02lat;
	int nIterCount = 0;
	while (++nIterCount < 1000)
	{
		//get geodetic offset relative to 'center china'
		double lon0 = wgs84lon - 105.0;
		double lat0 = wgs84lat - 35.0;

		//generate an pair offset roughly in meters
		double lon1 = 300.0 + lon0 + 2.0 * lat0 + 0.1 * lon0 * lon0 + 0.1 * lon0 * lat0 + 0.1 * sqrt(fabs(lon0));
		lon1 = lon1 + (20.0 * sin(6.0 * lon0 * PI) + 20.0 * sin(2.0 * lon0 * PI)) * 2.0 / 3.0;
		lon1 = lon1 + (20.0 * sin(lon0 * PI) + 40.0 * sin(lon0 / 3.0 * PI)) * 2.0 / 3.0;
		lon1 = lon1 + (150.0 * sin(lon0 / 12.0 * PI) + 300.0 * sin(lon0 * PI / 30.0)) * 2.0 / 3.0;
		double lat1 = -100.0 + 2.0 * lon0 + 3.0 * lat0 + 0.2 * lat0 * lat0 + 0.1 * lon0 * lat0 + 0.2 * sqrt(fabs(lon0));
		lat1 = lat1 + (20.0 * sin(6.0 * lon0 * PI) + 20.0 * sin(2.0 * lon0 * PI)) * 2.0 / 3.0;
		lat1 = lat1 + (20.0 * sin(lat0 * PI) + 40.0 * sin(lat0 / 3.0 * PI)) * 2.0 / 3.0;
		lat1 = lat1 + (160.0 * sin(lat0 / 12.0 * PI) + 320.0 * sin(lat0 * PI / 30.0)) * 2.0 / 3.0;

		double g_lon0 = 0;
		if (lon0 > 0)
			g_lon0 = 0.05 / sqrt(lon0);
		else
			if (lon0 < 0)
				g_lon0 = -0.05 / sqrt(-lon0);
			else
				g_lon0 = 0;

		double PIlon0 = PI * lon0, PIlat0 = PI * lat0;
		double dlon1_dlonwgs = 1 + 0.2 * lon0 + 0.1 * lat0 + g_lon0
			+ ((120 * PI * cos(6 * PIlon0) + 40 * PI * cos(2 * PIlon0))
				+ (20 * PI * cos(PIlon0) + 40 * PI / 3.0 * cos(PIlon0 / 3.0))
				+ (12.5 * PI * cos(PIlon0 / 12.0) + 10 * PI * cos(PIlon0 / 30.0))) * 2.0 / 3.0;
		double dlon1_dlatwgs = 2 + 0.1 * lon0;

		double dlat1_dlonwgs = 2 + 0.1 * lat0 + 2 * g_lon0
			+ (120 * PI * cos(6 * PIlon0) + 40 * PI * cos(2 * PIlon0)) * 2.0 / 3.0;
		double dlat1_dlatwgs = 3 + 0.4 * lat0 + 0.1 * lon0
			+ ((20 * PI * cos(PIlat0) + 40.0 * PI / 3.0 * cos(PIlat0 / 3.0))
				+ (40 * PI / 3.0 * cos(PIlat0 / 12.0) + 32.0 * PI / 3.0 * cos(PIlat0 / 30.0))) * 2.0 / 3.0;

		//latitude in radian
		double B = Deg2Rad(wgs84lat);
		double sinB = sin(B), cosB = cos(B);
		double WSQ = 1 - kKRASOVSKY_ECCSQ * sinB * sinB;
		double W = sqrt(WSQ);
		double N = kKRASOVSKY_A / W;

		double dW_dlatwgs = -PI * kKRASOVSKY_ECCSQ * sinB * cosB / (180.0 * W);
		double dN_dlatwgs = -kKRASOVSKY_A * dW_dlatwgs / WSQ;

		double PIxNxCosB = PI * N * cosB;
		double dlongcj_dlonwgs = 1.0 + 180.0 * dlon1_dlonwgs / PIxNxCosB;
		double dlongcj_dlatwgs = 180 * dlon1_dlatwgs / PIxNxCosB -
			180 * lon1 * PI * (dN_dlatwgs * cosB - PI * N * sinB / 180.0) / (PIxNxCosB * PIxNxCosB);

		double PIxNxSubECCSQ = PI * N * (1 - kKRASOVSKY_ECCSQ);
		double dlatgcj_dlonwgs = 180 * WSQ * dlat1_dlonwgs / PIxNxSubECCSQ;
		double dlatgcj_dlatwgs = 1.0 + 180 * (N * (dlat1_dlatwgs * WSQ + 2.0 * lat1 * W * dW_dlatwgs) - lat1 * WSQ * dN_dlatwgs) /
			(N * PIxNxSubECCSQ);

		auto [gcj02lonEst, gcj02latEst] = Wgs2Gcj(wgs84lon, wgs84lat);
		double l_lon = gcj02lon - gcj02lonEst;
		double l_lat = gcj02lat - gcj02latEst;

		double d_latwgs = (l_lon * dlatgcj_dlonwgs - l_lat * dlongcj_dlonwgs) /
			(dlongcj_dlatwgs * dlatgcj_dlonwgs - dlatgcj_dlatwgs * dlongcj_dlonwgs);
		double d_lonwgs = (l_lon - dlongcj_dlatwgs * d_latwgs) / dlongcj_dlonwgs;

		if (fabs(d_latwgs) < 1.0e-9 && fabs(d_lonwgs) < 1.0e-9)
			break;
		wgs84lon = wgs84lon + d_lonwgs;
		wgs84lat = wgs84lat + d_latwgs;
	}
	return { wgs84lon, wgs84lat };
}



