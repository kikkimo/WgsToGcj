/***************************************************************************
 *
 * This class WGStoGCJ implements geodetic coordinate transform between geodetic 
 * coordinate in WGS84 coordinate system and geodetic coordinate in GCJ-02 
 * coordinate system.
 * 
 * details: https://blog.csdn.net/gudufuyun/article/details/106738942
 *
 * Copyright (c) 2020.  All rights reserved.
 * 
 *  kikkimo                                              
 *  School of Remote Sensing  and Information Engineering,
 *	WuHan University,
 *	Wuhan, Hubei, P.R. China. 430079
 *  fywhu@outlook.com
 *
 ***************************************************************************
 *																		   *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#pragma once
#include <utility>

#ifdef _WGSTOGCJ_API
#if (defined _WIN32 || defined WINCE || defined __CYGWIN__)
#define WGSTOGCJ_API_EXPORTS __declspec(dllexport)
#elif defined __GNUC__ && __GNUC__ >= 4
#define WGSTOGCJ_API_EXPORTS __attribute__((visibility("default")))
#endif
#endif	

#ifndef WGSTOGCJ_API_EXPORTS
#define WGSTOGCJ_API_EXPORTS
#endif

/**
 *  \brief Covert geodetic coordinate in WGS84 coordinate system to geodetic coordinate
 *         in GCJ-02 coordinate system
 *
 *  \param [in] wgs84lon: longitude in WGS84 coordinate system [unit:degree]
 *  \param [in] wgs84lat: latitude in WGS84 coordinate system [unit:degree]
 *  \return Returns geodetic coordinate in GCJ-02 coordinate system
 *  \time 15:47:38 2020/06/12
 */
WGSTOGCJ_API_EXPORTS std::pair<double, double>
Wgs2Gcj(const double& wgs84lon, const double& wgs84lat);

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
std::pair<double, double> __declspec(dllexport)
Gcj2Wgs_SimpleIteration(const double& gcj02lon, const double& gcj02lat);

#ifdef _USE_CERES
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
WGSTOGCJ_API_EXPORTS std::pair<double, double>
Gcj2Wgs_AutoDiff(const double& gcj02lon, const double& gcj02lat);
#endif

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
WGSTOGCJ_API_EXPORTS std::pair<double, double>
Gcj2Wgs_NumbericDiff(const double& gcj02lon, const double& gcj02lat);

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
WGSTOGCJ_API_EXPORTS std::pair<double, double>
Gcj2Wgs_AnalyticDiff(const double& gcj02lon, const double& gcj02lat);


