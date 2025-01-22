/* 
This collection of functions calculates the Critical Wind Speed for a stand using
the ForestGales roughness algorithm. The data needs to go through the fg_rou_dataprep_Norway_cpp
function first in order to gather all the required variables if they are not present.
*/


#include <Rcpp.h>

// TREE FUNCTIONS
/*
' Convert top to mean tree height (and viceversa), calculate crown dimensions (width and depth) from dbh and tree height, calculate diameter at base of the tree for Deflection Loading Factor, and calculate "equivalent mean height" from stand top height for the TMC method.
param dbh Diameter of the stem at breast height, i.e. 1.3m above the ground (cm). Depending on the method used ('roughness' or 'TMC') this can be either the arithmetic average of the dbh of all the trees in the stand, or the dbh of an individual tree.
param ht Tree height. Depending on the method used ('roughness' or 'TMC'), this can be either the mean tree in the stand, or each individual tree (m).
param cr_depth Length of the tree crown (m).
param top_ht Dominant tree height in the stand (m).
param mean_ht Mean tree height in the stand. For the TMC method, arithmetic mean height of the trees in the stand (m).
param equivalent_mean_ht Height of stand-level wind moment absorption hypothesised for the TMC method. Calculated from stand top height using a linear relationship derived from the mean height - top height regressions values across all tree species currently included in the model.
param param0_height Species-specific parameter for conversions between mean height and top height. Offset in linear regression equation.
param param1_height Species-specific parameter for conversions between mean height and top height. Multiplier in linear regression equation.
param param0_cr_width Species-specific parameter for calculation of crown width from dbh. Offset in linear regression equation.
param param1_cr_width Species-specific parameter for calculation of crown width from dbh. Multiplier in linear regression equation.
param param0_cr_depth Species-specific parameter for calculation of crown length from tree height. Offset in linear regression equation.
param param1_cr_depth Species-specific parameter for calculation of crown length from tree height. Multiplier in linear regression equation.
name Tree_Dimensions_Functions
title Tree Dimensions Functions
*/


// [[Rcpp::export]]
double diam_base_fun_cpp(double dbh, double ht, double cr_depth) {
	double wind_action_point = ht - (cr_depth/2);
	if(wind_action_point < 1.3 * 2.0) {
		wind_action_point = 1.3 * 2.0;
	}
	// dbase and dbh in cm in this formula
	double dbase = (dbh / pow((wind_action_point - 1.3), 0.333)) * pow(wind_action_point, 0.333);
	return dbase;
}

// DRAG
// [[Rcpp::export]]
double drag_fun_cpp(double uguess, double n_drag, double c_drag, double drag_upper_limit) {
	
	double drag;
	
	if(uguess < 10.0){
		drag = c_drag * pow(10.0, -n_drag);
	} else if (uguess > drag_upper_limit){
		drag = c_drag * pow(drag_upper_limit, -n_drag);
	} else {
		drag = c_drag * pow(uguess, -n_drag);
	}
	return drag;
}

// CANOPY BREADTH
// [[Rcpp::export]]
double canopy_breadth_fun_cpp(double cr_width, double uguess, double n_drag, double c_drag, double drag_upper_limit) {
	
	double canopy_breadth = cr_width * drag_fun_cpp(uguess, n_drag, c_drag, drag_upper_limit)/2.0;
	
	return canopy_breadth;
}


// MAX GAP FACTOR
/*
Effect of gap size on maximum bending moment
title Maximum Gap Factor
param gap_size Length of the upwind gap (m).
param ht Height of the tree. In the 'roughness' method, this is stand mean height. In the TMC method, this is 'equivalent mean height' (m).
return \code{max_gap_factor}, the effect of gap alone on the maximum bending moment exerted by the wind on a tree.
*/
// [[Rcpp::export]]
/* Common to Roughness and TMC Methods */

double max_gap_factor_fun_cpp(double gap_size, double ht) {
	
	double g_h;
	if(gap_size/ht > 10.0){
		g_h = 10.0;
	} else {
		g_h = gap_size/ht;
	}
	double max_gap_factor = pow(g_h/10.0, 0.3467);
	
	return max_gap_factor;
}

// EDGE GAP GUST FACTOR - ROUGHNESS VERSION
/*
Combined effect of gap size, edge effect, and gustiness on maximum bending moment for the roughness method.
title Edge Gap Gust Factor.
param spacing Mean distance between trees (m).
param mean_ht The mean tree height in the stand (m).
param dist_edge Distance of tree from the upwind edge (m).
param gap_size Length of the upwind gap (m).
return \code{edge_gap_gust_factor}, the combined effect of edge, gap, and gustiness on the maximum bending moment exerted by the wind on a tree. Used in the roughness method.
*/
// [[Rcpp::export]]

double edge_gap_gust_factor_fun_cpp(double spacing, double mean_ht, double dist_edge, double gap_size, Rcpp::DataFrame fgr_constants) {
		
	double s_h = spacing / mean_ht;
	if (s_h < 0.075) {
		s_h = 0.075;
	} else if(s_h > 0.45){
		s_h = 0.45;
	}
	double tree_heights_inside_forest = fgr_constants["tree_heights_inside_forest"];
	double dist_edge_lim = mean_ht * tree_heights_inside_forest;
	if(dist_edge > dist_edge_lim){
		dist_edge = dist_edge_lim;
	}
	double max_gap_factor = max_gap_factor_fun_cpp(gap_size, mean_ht);
	
	double edge_gap_gust_factor = ((2.7193 * s_h - 0.061) + (-1.273 * s_h + 0.9701) * pow((1.1127 * s_h + 0.0311), (tree_heights_inside_forest*mean_ht/mean_ht)) +
	((-1.273 * s_h + 0.9701) * (pow((1.1127 * s_h + 0.0311), (dist_edge/mean_ht)) - pow((1.1127 * s_h + 0.0311), (tree_heights_inside_forest*mean_ht/mean_ht)))) *
	max_gap_factor) /
    ((0.68 * s_h - 0.0385) + (-0.68 * s_h + 0.4785) * pow((1.7239 * s_h +0.0316), (tree_heights_inside_forest*mean_ht/mean_ht)));

	return edge_gap_gust_factor;
}

// EDGE GAP FACTOR - TMC VERSION
/*
Combined effect of gap size and edge effect on maximum bending moment for the single-tree method.
title Edge Gap Factor.
param spacing Mean distance between trees (m).
param equivalent_mean_ht Equivalent mean stand height: the level in the stand responsible for most of the momentum absorption (m).
param dist_edge Distance of tree from the upwind edge (m).
param gap_size Length of the upwind gap (m).
return \code{edge_gap_factor}, the combined effect of edge and gap on the maximum bending moment exerted by the wind on a tree. Used in the tmc method.
*/
// [[Rcpp::export]]
double edge_gap_factor_fun_cpp(double spacing, double equivalent_mean_ht, double dist_edge, double gap_size, Rcpp::DataFrame fgr_constants) {
	
	double s_h = spacing / equivalent_mean_ht;
	if (s_h < 0.075) {
		s_h = 0.075;
	} else if(s_h > 0.45){
		s_h = 0.45;
	}
	
	double tree_heights_inside_forest = fgr_constants["tree_heights_inside_forest"];
	double dist_edge_lim = equivalent_mean_ht * tree_heights_inside_forest;
	
	if (dist_edge > dist_edge_lim) {
		dist_edge = dist_edge_lim;
	}
	
	double max_gap_factor = max_gap_factor_fun_cpp(gap_size, equivalent_mean_ht);
	double edge_gap_factor = ((2.7193 * s_h - 0.061) + (-1.273 * s_h + 0.9701) * pow((1.1127 * s_h + 0.0311), (tree_heights_inside_forest*equivalent_mean_ht/equivalent_mean_ht)) +
	((-1.273 * s_h + 0.9701) * (pow((1.1127 * s_h + 0.0311), (dist_edge/equivalent_mean_ht)) - 
	pow((1.1127 * s_h + 0.0311), (tree_heights_inside_forest*equivalent_mean_ht/equivalent_mean_ht)))) *
	max_gap_factor) /
	((2.7193 * s_h - 0.061) + (-1.273 * s_h + 0.9701) * pow((1.1127 * s_h + 0.0311), (tree_heights_inside_forest * equivalent_mean_ht/equivalent_mean_ht)));
	
	return edge_gap_factor;
}

// TURNING COEFFICIENTS FUNCTIONS - TMC VERSION
/*
Calculate the turning moment coefficient.
param dbh Diameter of the stem at breast height, i.e. 1.3m above the ground (cm).
param ht Individual tree height (m).
param ci Competition Index (\code{bal}, \code{heg}, \code{none}) used.
param ci_value Value of \code{ci}.
title Turning Moment Coefficient Functions
*/
// [[Rcpp::export]]
double tc_intercept_fun_cpp(double dbh, double ht, std::string ci, double ci_value) {
	
	double tc_intercept = 0.0;
	
	/* Note that dbh is in meters ! */
	if (ci == "none") {
		tc_intercept = -9.64 + 113.78 * pow(dbh/100.0, 2) * ht;
	}
	/* Bal */
	else if (ci == "bal") {
		tc_intercept = 40.58 - 0.5337 * ci_value + 111.79 * pow(dbh/100.0, 2) * ht - 0.776 * pow(dbh/100.0, 2) * ht * ci_value;
	}
	/* Hegyi */
	else if (ci == "heg") {
		tc_intercept = 50.958 - 7.734 * ci_value + 123.819 * pow(dbh/100.0, 2) * ht - 26.535 * pow(dbh/100.0, 2) * ht * ci_value;
	}
	if (tc_intercept == 0) {
		tc_intercept = 9.64 + 113.78 * pow(dbh/100.0, 2) * ht;
	}
	
	return tc_intercept;
}

// [[Rcpp::export]]
double tc_zero_intercept_fun_cpp(double dbh, double ht, std::string ci, double ci_value) {
	
	double tc_zero_intercept = 0;
	
	/* Note that dbh is in meters ! */
	if(ci == "none") {
		tc_zero_intercept = 111.91 * pow(dbh/100.0, 2) * ht;
	} 
	/* Bal */
	else if(ci == "bal") {
		tc_zero_intercept = 0.130 * ci_value + 116.304 * pow(dbh/100.0, 2) * ht - 0.617 * pow(dbh/100.0, 2) * ht * ci_value;
	} 
	/* Hegyi */
	else if(ci == "heg") {
		tc_zero_intercept = 3.86 * ci_value + 124.252 * pow(dbh/100.0, 2) * ht - 17.285 * pow(dbh/100.0, 2) * ht * ci_value;
	}	
	if(tc_zero_intercept == 0) {
		tc_zero_intercept = 111.91 * pow(dbh/100.0, 2) * ht;
	}
	
	return tc_zero_intercept;
}
// [[Rcpp::export]]
double tc_zero_intercept_new_bal_fun_cpp(double dbh, double ht, std::string ci, double ci_value) {
	
	double tc_zero_intercept_new_bal = 0;
	
	/* Note that dbh is in meters ! */
	if(ci == "none") {
		tc_zero_intercept_new_bal = 111.91 * pow(dbh/100.0, 2) * ht;
	} 
	/* Bal */
	else if (ci == "bal") {
		tc_zero_intercept_new_bal = 0.274 * ci_value;
	} 
	/* Hegyi */
	else if (ci == "heg") {
		tc_zero_intercept_new_bal = 3.86 * ci_value + 124.252 * pow(dbh/100.0, 2) * ht - 17.285 * pow(dbh/100.0, 2) * ht * ci_value;
	}
	if (tc_zero_intercept_new_bal == 0) {
		tc_zero_intercept_new_bal = 111.91 * pow(dbh/100.0, 2) * ht;
	}
		
  return tc_zero_intercept_new_bal;
}

// [[Rcpp::export]]
double tc_zero_intercept_fun_balBA_cpp(double dbh, double ht, std::string ci, double ci_value) {

	double tc_zero_intercept = 0;
	
	/* Note that dbh is in meters ! */
	if(ci == "none") {
		tc_zero_intercept = 111.573 * pow(dbh/100.0, 2) * ht;
	} 
	/* Bal */
	else if (ci == "bal") {
		tc_zero_intercept = 113.54 * pow(dbh/100.0, 2) * ht - 20.494 * ci_value;
	} 
	/* Hegyi */
	else if (ci == "heg") {
		tc_zero_intercept = 3.86 * ci_value + 124.252 * pow(dbh/100.0, 2) * ht - 17.285 * pow(dbh/100.0, 2) * ht * ci_value;
	}
	if (tc_zero_intercept == 0) {
		tc_zero_intercept = 111.573 * pow(dbh/100.0, 2) * ht;
	}
		
  return tc_zero_intercept;
}



// CRITICAL MOMENTS

// Calculate the critical moments for breakage and overturning

// ROUGHNESS VERSION
/*
param dbh Diameter of the stem at breast height, i.e. 1.3m above the ground (cm). Depending on the method used ('roughness' or 'TMC') this can be either the arithmetic average of the dbh of all the trees in the stand, or the dbh of an individual tree.
param ht  Tree height (either mean tree height in the stand or individual tree height).
param mor Modulus of Rupture of green wood (MPa).
param fknot Knot factor (dimensionless).
param c_reg Regression coefficient of uprooting moment against stem weight (N m kg-1).
param stem_density Density of green wood of the stem (kg m-3).
param stem_vol Volume of the tree stem of the mean tree in the stand (m3). For the roughness method, this is stem volume of the mean tree. For the TMC method, this is individual tree stem voume.
name Critical_Moments_Functions
title Critical Moments Functions
*/
// [[Rcpp::export]]

double critical_moment_breakage_rou_cpp(double dbh, double mor, double fknot) {
	
	/* Note dbh in meters */
	double breaking_moment = mor * fknot * M_PI * pow((dbh/100.0),3) / 32.0;

	return breaking_moment;
}
// TMC VERSION
// [[Rcpp::export]]
double critical_moment_breakage_tmc_cpp(double dbh, double ht, double cr_depth, double mor, double fknot){
	
	/* Note dbh in meters */
	double dbase = diam_base_fun_cpp(dbh, ht, cr_depth);
	double breaking_moment = mor * fknot * M_PI * (pow(dbase/100.0, 3))/32.0;
	
	return breaking_moment;
}

// [[Rcpp::export]]
double critical_moment_overturning_cpp(double c_reg, double stem_density, double stem_vol) {
	double overturning_moment = c_reg * stem_density * stem_vol;
	
	return overturning_moment;
}


//LAMBDA CAPITAL
/*
Calculates the frontal area of the streamlined canopy per ground area.
title Lambda Capital Function.
param cr_width Width of the canopy in windless conditions (m).
param cr_depth Length of the canopy in windless conditions (m).
param spacing Mean distance between trees (m).
param uguess Critical wind speed at canopy top calculated with the roughness or single-tree method (m s-1).
param n_drag N parameter of the drag coefficient formula (dimensionless).
param c_drag C parameter of the drag coefficient formula (dimensionless).
param drag_upper_limit Maximum wind speed used during the experiments from which \code{n_drag} and \code{c_drag} were derived (m*s-1).
return \code{lambdacapital}, the frontal area of the streamlined crown under per ground area.
 The frontal area of the streamlined canopy per ground area (defined by spacing). Effectively, LambdaCapital is drag per unit ground area 
*/
// [[Rcpp::export]]
double lambdacapital_fun_cpp(double cr_width, double cr_depth, double spacing, double uguess, double n_drag, double c_drag, double drag_upper_limit) {
	
	double canopy_breadth = canopy_breadth_fun_cpp(cr_width, uguess, n_drag, c_drag, drag_upper_limit);
	double lambdacapital = 2.0 * (canopy_breadth * cr_depth / pow(spacing, 2));
	
	return lambdacapital;
}

// GAMMA SOLVED
/*
Calculates the ratio between the critical wind speed at canopy top and the friction velocity (uh / u*).
title Gamma Solved Function.
param cr_width Width of the canopy in windless conditions (m).
param cr_depth Length of the canopy in windless conditions (m).
param spacing Mean distance between trees (m).
param uguess Critical wind speed at canopy top calculated with the roughness or single-tree method (m s-1).
param n_drag N parameter of the drag coefficient formula (dimensionless).
param c_drag C parameter of the drag coefficient formula (dimensionless).
param drag_upper_limit Maximum wind speed used during the experiments from which \code{n_drag} and \code{c_drag} were derived (m*s-1).
return \code{gammasolved}, the ratio of critical wind speed at canopy top over friction velocity (uh / u*).
*/
// [[Rcpp::export]]

double gammasolved_fun_cpp(double cr_width, double cr_depth, double spacing, double uguess, double n_drag, double c_drag, double drag_upper_limit, Rcpp::DataFrame fgr_constants) {
	
	double gammasolved;
	double lambdacapital = lambdacapital_fun_cpp(cr_width, cr_depth, spacing, uguess, n_drag, c_drag, drag_upper_limit);
	double cs = fgr_constants["cs"];
	double cr = fgr_constants["cr"];
	
	if(lambdacapital > 0.6) {
		gammasolved = 1.0/sqrt(cs + cr * 0.3);
	} else {
		gammasolved = 1.0/sqrt(cs + cr * lambdacapital/2.0);
	}
	
	return gammasolved;
}


// ZPD
/*
Calculates the height of the zero plane displacement.
Zero Plane Displacement Function.
param cr_width Width of the canopy in windless conditions (m).
param cr_depth Length of the canopy in windless conditions (m).
param spacing Mean distance between trees (m).
param uguess Critical wind speed at canopy top calculated with the roughness or single-tree method (m s-1).
param n_drag N parameter of the drag coefficient formula (dimensionless).
param c_drag C parameter of the drag coefficient formula (dimensionless).
param drag_upper_limit Maximum wind speed used during the experiments from which \code{n_drag} and \code{c_drag} were derived (m*s-1).
param ht Height of the tree. In the 'roughness' method, this is stand mean height. In the TMC method, this is 'equivalent mean height' (m).
returns zpd, the height of the zero plane displacement (m).
-- Calculation of Zero-plane Displacement (d) and Aerodynamic Roughness (z0) from "Simplified Expressions for vegetation roughness and
zero-plane displacement as functions of canopy height and area index" by M.R. Raupach (1994)
*/
// [[Rcpp::export]]
double zpd_fun_cpp(double cr_width, double cr_depth, double spacing, double uguess, double n_drag, double c_drag, double drag_upper_limit, double ht, Rcpp::DataFrame fgr_constants) {
	double cd1 = fgr_constants["cd1"];
	
	double lambda1 = lambdacapital_fun_cpp(cr_width, cr_depth, spacing, uguess, n_drag, c_drag, drag_upper_limit);
	double zpd = (1.0 - ((1.0 - exp( - sqrt(cd1 * lambda1))) / sqrt(cd1 * lambda1))) * ht;
	return zpd;
}


// Z0
/*
' Calculates the length (m) of the canopy surface roughness.
Canopy Roughness Function.
param cr_width Width of the canopy in windless conditions (m).
param cr_depth Length of the canopy in windless conditions (m).
param spacing Mean distance between trees (m).
param uguess Critical wind speed at canopy top calculated with the roughness or single-tree method (m s-1).
param n_drag N parameter of the drag coefficient formula (dimensionless).
param c_drag C parameter of the drag coefficient formula (dimensionless).
param drag_upper_limit Maximum wind speed used during the experiments from which \code{n_drag} and \code{c_drag} were derived (m*s-1).
param ht Height of the tree. In the 'roughness' method, this is stand mean height. In the TMC method, this is 'equivalent mean height' (m).
returns z0, the length (m) of the surface roughness of the tree canopy.
-- Calculation of Zero-plane Displacement (d) and Aerodynamic Roughness (z0) from "Simplified Expressions for vegetation roughness and
zero-plane displacement as functions of canopy height and area index" by M.R. Raupach (1994)
*/
// [[Rcpp::export]]

double z0_fun_cpp(double cr_width, double cr_depth, double spacing, double uguess, double n_drag, double c_drag, double drag_upper_limit, double ht, Rcpp::DataFrame fgr_constants) {
	double k = fgr_constants["k"];
	double cw = fgr_constants["cw"];
	
	double zpd1 = zpd_fun_cpp(cr_width, cr_depth, spacing, uguess, n_drag, c_drag, drag_upper_limit, ht, fgr_constants);
	double gamma1 = gammasolved_fun_cpp(cr_width, cr_depth, spacing, uguess, n_drag, c_drag, drag_upper_limit, fgr_constants);
	double z0 = (ht - zpd1) * exp((-k * gamma1) + (log(cw) - 1.0 + pow(cw, -1)));
	return z0;
}


// TURNING MOMENT RATIOS FUNCTIONS - TMC VERSION
/*
Calculate the ratio between the turning moment coefficient before and after thinning in the stand.
param spacing_before Mean spacing of trees in the stand before any thinning (m).
param spacing_current Current mean spacing of trees in the stand (m).
param years_since_thin Number of years after the latest thinning.
param cr_width Width of the crown of the "equivalent mean tree" in the stand (m).
param cr_depth Length of the crown of the "equivalent mean tree" in the stand (m).
param uguess Critical wind speed at canopy top calculated with the roughness or single-tree method (m s-1).
param n_drag N parameter of the drag coefficient formula (dimensionless).
param c_drag C parameter of the drag coefficient formula (dimensionless).
param drag_upper_limit Maximum wind speed used during the experiments from which \code{n_drag} and \code{c_drag} were derived (m*s-1).
param ht Equivalent mean stand height: the level in the stand responsible for most of the momentum absorption (m).
param ci Competition Index (\code{bal}, \code{heg}, \code{none}) used.
name Turning_Moment_Ratios
title Turning Moment Ratio Functions
*/
// [[Rcpp::export]]
double tm_ratio_cpp(double spacing_before, double spacing_current, int years_since_thin, double cr_width, double cr_depth, double uguess, 
double n_drag, double c_drag, double drag_upper_limit, double ht, Rcpp::DataFrame fgr_constants) {
	
	double tmr_full;
	double d_before = zpd_fun_cpp(cr_width, cr_depth, spacing_before,
	uguess, n_drag, c_drag, drag_upper_limit, ht, fgr_constants);
	double d_current = zpd_fun_cpp(cr_width, cr_depth, spacing_current, 
	uguess, n_drag, c_drag, drag_upper_limit, ht, fgr_constants);
    double z0_before = z0_fun_cpp(cr_width, cr_depth, spacing_before, 
	uguess, n_drag, c_drag, drag_upper_limit, ht, fgr_constants);
	double z0_current = z0_fun_cpp(cr_width, cr_depth, spacing_current, 
	uguess, n_drag, c_drag, drag_upper_limit, ht, fgr_constants);
	
	if (years_since_thin > 5.0){tmr_full = 1.0;}
	else {
		tmr_full = (d_current/d_before) * pow(spacing_current/spacing_before, 2) * 
		pow(log((ht - d_before)/z0_before)/log((ht - d_current)/z0_current), 2) * 
		(1.0 - years_since_thin/5.0) + (years_since_thin/5.0);
	}
	
	return tmr_full;
}

// [[Rcpp::export]]
double tm_ratio_ci_cpp(double spacing_before, double spacing_current, int years_since_thin, double cr_width, double cr_depth, double uguess, 
double n_drag, double c_drag, double drag_upper_limit, double ht, std::string ci, Rcpp::DataFrame fgr_constants) {
	
	double tmrci_full;
	double d_before = zpd_fun_cpp(cr_width, cr_depth, spacing_before,
	uguess, n_drag, c_drag, drag_upper_limit, ht, fgr_constants);
	double d_current = zpd_fun_cpp(cr_width, cr_depth, spacing_current, 
	uguess, n_drag, c_drag, drag_upper_limit, ht, fgr_constants);
    double z0_before = z0_fun_cpp(cr_width, cr_depth, spacing_before, 
	uguess, n_drag, c_drag, drag_upper_limit, ht, fgr_constants);
	double z0_current = z0_fun_cpp(cr_width, cr_depth, spacing_current, 
	uguess, n_drag, c_drag, drag_upper_limit, ht, fgr_constants);
	
	if (years_since_thin > 5.0){tmrci_full = 1.0;}
	else if (ci == "none") {
		tmrci_full = (d_current/d_before) * pow(spacing_current/spacing_before, 2) * 
		pow(log((ht - d_before)/z0_before)/log((ht - d_current)/z0_current), 2) * 
		(1.0 - years_since_thin/5.0) + (years_since_thin/5.0);
	}
	else {tmrci_full = 1.0;}
	
	return tmrci_full;
}

// [[Rcpp::export]]
double tm_ratio_simple_cpp(double spacing_before, double spacing_current, int years_since_thin) {
	
	double tmr_simple;
	
	if(years_since_thin > 5.0) {tmr_simple = 1.0;} 
	else {
		tmr_simple = (0.99 * spacing_current/spacing_before) * (1 - years_since_thin/5.0) + (years_since_thin/5.0);
	}
	return tmr_simple;
}

// [[Rcpp::export]]
double tm_ratio_simple_ci_cpp(double spacing_before, double spacing_current, int years_since_thin, std::string ci) {
	
	double tmrci_simple;
	
	if(years_since_thin > 5.0) {tmrci_simple = 1.0;} 
	else if (ci == "none") {
		tmrci_simple = (0.99 * spacing_current/spacing_before) * (1 - years_since_thin/5.0) + (years_since_thin/5.0);
	} 
	else {tmrci_simple = 1.0;}
	
	return tmrci_simple;
}

// BENDING MOMENT
/*
Calculates the bending moment on a tree resulting from the load applied by the wind.
title Bending moment on a tree under wind loading
param dbh Arithmetic average of the diameter of the stem at breast height, i.e. 1.3m above the ground, for all trees in the stand (cm).
param ht Mean tree height in the stand (m).
param cr_width Width of the tree crown (m).
param cr_depth Depth of the tree crown (m).
param spacing Mean distance between trees (m).
param dist_edge Distance of the the tree from the upwind edge (m).
param gap_size Length of the upwind gap (m).
param uguess Critical wind speed at canopy top calculated with the roughness or single-tree method (m s-1).
param n_drag N parameter of the drag coefficient formula (dimensionless).
param c_drag C parameter of the drag coefficient formula (dimensionless).
param drag_upper_limit Maximum wind speed used during the experiments from which \code{n_drag} and \code{c_drag} were derived (m*s-1).
return Applied bending moment resulting from the wind load (N).
example bending_moment_rou(25, 20, 6.5, 12, 2.2, 0, 0, 15.6, 0.51, 2.35, 25)
*/
// [[Rcpp::export]]
double bending_moment_rou_cpp(double dbh, double ht, double cr_width, double cr_depth, double spacing, double dist_edge,
double gap_size, double uguess, double n_drag, double c_drag, double drag_upper_limit, double ro, Rcpp::DataFrame fgr_constants,
double aerodynamic_ht) {
	double k = fgr_constants["k"];
	
	double zpd_1 = zpd_fun_cpp(cr_width, cr_depth, spacing, uguess, n_drag, c_drag, drag_upper_limit, aerodynamic_ht, fgr_constants);
	double edge_1 = edge_gap_gust_factor_fun_cpp(spacing, ht, dist_edge, gap_size, fgr_constants);
	double z0_1 = z0_fun_cpp(cr_width, cr_depth, spacing, uguess, n_drag, c_drag, drag_upper_limit, aerodynamic_ht, fgr_constants);
		
	double bm_rou = zpd_1 * ro * edge_1 * pow((spacing * uguess * k)/log((aerodynamic_ht - zpd_1)/z0_1), 2);
	
	return bm_rou;
}

// DEFLECTION LOADING FACTOR
/*
Calculate the Deflection Loading Factor to account for the additional moment provided by the weight of the stem, crown (and snow, when present)
param bm Applied bending moment resulting from the wind load (N).
param ht Tree height. Depending on the method used ('roughness' or 'TMC'), this can be either the mean tree in the stand, or each individual tree (m).
param dbh Diameter of the stem at breast height, i.e. 1.3m above the ground (cm). Depending on the method used ('roughness' or 'TMC') this can be either the arithmetic average of the dbh of all the trees in the stand, or the dbh of an individual tree.
param cr_depth Length of the tree crown (m).
param cr_width Width of the tree crown (m).
param moe Modulus of Elasticity of green wood (MPa).
param stem_density Density of green wood of the stem (kg m-3).
param crown_density Density of of the tree crown (kg m-3).
param stem_density Density of green wood of the stem (kg m-3).
param stem_vol Volume of the tree stem of the mean tree in the stand (m3). For the roughness method, this is stem volume of the mean tree. For the TMC method, this is individual tree stem voume.
param snow_depth Depth of layer of snow on tree crown (cm).
param snow_density Density of snow (kg m-3).
param x Height along the tree stem where deflection is to be calculated (m).
param lever_arm Length of the lever arm (m). Typically, it is equal tree height.
param fow Applied wind loading in the \code{deflection_fun} function. It is a placeholder for the output of the \code{force_of_wind_fun} function.
param pull_height Height along the tree stem where the wind loading is applied (m).
name DFL_Functions
title Deflection Loading Factor Functions
*/
// [[Rcpp::export]]
double force_of_wind_fun_cpp(double bm, double ht, double cr_depth) {
	double force_of_wind = bm / (ht - cr_depth/2.0);
	return force_of_wind;
}

// [[Rcpp::export]]
double r_fun_cpp(double lever_arm, double pull_height) {
	//Ratio of distance from top of tree where force is applied to tree tree height (ht)
	double r_value = (lever_arm - pull_height) / lever_arm;
	return r_value;
}

// [[Rcpp::export]]
double i_fun_cpp(double dbh, double ht, double cr_depth) {
	//Second area moment of Inertia calculated at tree base:
	// dbase in m in this formula
	double i_value = M_PI * (pow(diam_base_fun_cpp(dbh/100.0, ht, cr_depth), 4) / 64.0);
	return i_value;
}

// [[Rcpp::export]]
double deflection_fun_cpp(double x, double lever_arm, double fow, double dbh, double ht, double cr_depth, double pull_height, double moe, Rcpp::DataFrame fgr_constants) {
	//Gardiner 1989 'Mechanical characteristics of Sitka Spruce' Eq. 5: deflection at distance x from top of tree. fow is ForceOfWind
	double an = fgr_constants["an"];
	double bn = fgr_constants["bn"];
	double cn = fgr_constants["cn"];

	double hh = an * bn * cn;
	double i_value = i_fun_cpp(dbh, ht, cr_depth);
	double deflection = ((fow * pow(lever_arm, 3)) / (moe* i_value * hh)) * (an * pow((x/lever_arm), cn) -
	r_fun_cpp(lever_arm, pull_height) * cn * pow((x/lever_arm), bn) +
	r_fun_cpp(lever_arm, pull_height) * bn * cn * (x/lever_arm) -
	an * cn * (x/lever_arm) + an * bn -
	r_fun_cpp(lever_arm, pull_height) * an * cn);
	
	return deflection;
}

// [[Rcpp::export]]
double dlf_fun_Norway_cpp(double bm, double ht, double cr_depth, double cr_width, double cr_prarea, double crown_off, double stem_vol, double dbh, double moe, double crown_density, double stem_density, double snow_load, Rcpp::DataFrame fgr_constants) {
	double grav = fgr_constants["grav"];
	
	double fow = force_of_wind_fun_cpp(bm, ht, cr_depth);
	/* Displacement at mid-crown height */
	double deflection_1 = deflection_fun_cpp(cr_depth/2.0, ht, fow, dbh, ht, cr_depth, (ht - cr_depth/2.0), moe, fgr_constants); 
	/* Deflection ¾ of the way down the stem.  Assumption that centre of stem mass is at this height */
	double deflection_2 = deflection_fun_cpp((3.0/4.0) * ht, ht, fow, dbh, ht, cr_depth, (ht - cr_depth/2.0), moe, fgr_constants);
	/* Weight of crown and snow */
	double weight_crown = ((((cr_prarea * cr_depth)/3.0) * crown_density) + (cr_prarea * snow_load)) * grav;
	/* Weight of the stem */
	double weight_stem = stem_vol * stem_density * grav;
	
	/* Note that BM is in both numerator and denominator, so the magnitude of the BM does not have any impact on DLF. It was decided to be kept in the equation for completeness */
	double ratioaddmomenttoBM  = ((crown_off + deflection_1) * weight_crown + deflection_2 * weight_stem)/bm;
	if (ratioaddmomenttoBM >= 1){ratioaddmomenttoBM = 0.999999;}
	double DLF_calc = 1.0/(1.0 - ratioaddmomenttoBM);
		
	return DLF_calc;	
}



// CWS BREAKAGE
/*
Calculates the critical wind speed for breakage using the roughness method.
title Critical wind speed for breakage - Roughness method.
param mean_ht Arithmetic mean height of the trees in the stand (m).
param mean_dbh Mean dbh of all trees in the stand (cm). Dbh is diameter at breast height, measured at 1.3m above the ground. Essential.
param spacing Mean distance between trees (m).
param dist_edge Distance from the upwind edge (m).
param gap_size Length of the upwind gap (m).
param mean_cr_width Width of the crown of the mean tree in the stand (m).
param mean_cr_depth Depth of the crown of the mean tree in the stand (m).
param moe Modulus of Elasticity of green wood (MPa).
param mor Modulus of Rupture of green wood (MPa).
param fknot Knot factor (dimensionless).
param n_drag N parameter of the drag coefficient formula (dimensionless).
param c_drag C parameter of the drag coefficient formula (dimensionless).
param drag_upper_limit Maximum wind speed used during the experiments from which \code{n_drag} and \code{c_drag} were derived (m*s-1).
param stem_vol Stem volume of the mean tree in the stand (m^3).
param stem_density Density of green wood of the stem (kg m-3).
param crown_density Crown density of the mean tree in the stand (kg m-3).
param snow_depth Depth of layer of snow on tree crown (cm).
param snow_density Density of snow (kg m-3).
param ro Air density (kg m-3).
return A list including: The critical wind speed (m s-1) for breakage at zero plane displacement height for the roughness method;
The applied bending moment; The zero plane displacement height; The ratio of critical wind speed at canopy top over friction velocity
(uh / u*); The deflection loading factor; The critical breaking moment; The combined effect of edge, gap, and gust on the applied
bending moment.
*/
// [[Rcpp::export]]
Rcpp::DataFrame uh_breakage_rou_Norway_cpp(Rcpp::NumericVector mean_ht, Rcpp::NumericVector mean_dbh, Rcpp::NumericVector spacing, Rcpp::NumericVector dist_edge, Rcpp::NumericVector gap_size, 
Rcpp::NumericVector mean_cr_width, Rcpp::NumericVector mean_cr_depth, Rcpp::NumericVector mean_cr_prarea, Rcpp::NumericVector crown_off, Rcpp::NumericVector moe, Rcpp::NumericVector mor, Rcpp::NumericVector fknot, Rcpp::NumericVector n_drag, Rcpp::NumericVector c_drag, Rcpp::NumericVector drag_upper_limit,
Rcpp::NumericVector stem_vol, Rcpp::NumericVector stem_density, Rcpp::NumericVector crown_density, Rcpp::NumericVector snow_load, double ro, Rcpp::DataFrame fgr_constants, Rcpp::NumericVector aerodynamic_ht) {
	
	int n = mean_ht.size();
	Rcpp::NumericVector edge_gap_gust_factor(n);
	Rcpp::NumericVector breaking_moment(n);
	Rcpp::NumericVector bm_rou(n);
	Rcpp::NumericVector dlf_calc(n);
	Rcpp::NumericVector dlf_used(n);
	Rcpp::NumericVector zpd(n);
	Rcpp::NumericVector gammasolved(n);
	Rcpp::NumericVector uh_b(n);

	double wind_precision = fgr_constants["wind_precision"];
	double dlf = fgr_constants["dlf"];
	double dlf_threshold = fgr_constants["dlf_threshold"];
	
		
	for (int i = 0; i < n; ++i){		
		// initialise the windspeed guesses
		double uguess = 25.0;
		double uguess1 = uguess;
		uh_b[i] = uguess/2.0;	
		
		edge_gap_gust_factor[i] = edge_gap_gust_factor_fun_cpp(spacing[i], mean_ht[i], dist_edge[i], gap_size[i], fgr_constants);
		breaking_moment[i] = critical_moment_breakage_rou_cpp(mean_dbh[i], mor[i], fknot[i]);	
		
		while(abs(uguess1 - uh_b[i]) > wind_precision){
			uguess1 = uguess;			
			bm_rou[i] = bending_moment_rou_cpp(mean_dbh[i], mean_ht[i], mean_cr_width[i], mean_cr_depth[i], spacing[i], dist_edge[i], gap_size[i], uguess, n_drag[i], c_drag[i], drag_upper_limit[i], ro, fgr_constants, aerodynamic_ht[i]);
			dlf_calc[i] = dlf_fun_Norway_cpp(bm_rou[i], mean_ht[i], mean_cr_depth[i], mean_cr_width[i], mean_cr_prarea[i], crown_off[i], stem_vol[i], mean_dbh[i], moe[i], crown_density[i], stem_density[i], snow_load[i], fgr_constants);
			dlf_used[i] = dlf_calc[i];
			if(dlf_used[i] < 0.0) {dlf_used[i] = dlf;} else if(dlf_used[i] > dlf_threshold) {dlf_used[i] = dlf_threshold;}
			zpd[i] = zpd_fun_cpp(mean_cr_width[i], mean_cr_depth[i], spacing[i], uguess, n_drag[i], c_drag[i], drag_upper_limit[i], aerodynamic_ht[i], fgr_constants);
			if(zpd[i] < 1.3){zpd[i] = 1.35;}
			gammasolved[i] = gammasolved_fun_cpp(mean_cr_width[i], mean_cr_depth[i], spacing[i], uguess, n_drag[i], c_drag[i], drag_upper_limit[i], fgr_constants);
			uh_b[i] = (1.0/spacing[i]) * sqrt(breaking_moment[i] /  (ro * edge_gap_gust_factor[i] * (zpd[i] - 1.3) * dlf_used[i])) * gammasolved[i];
			uguess = uh_b[i];
		}
	}
	
	Rcpp::DataFrame uh_b_results = Rcpp::DataFrame::create( Rcpp::_["uh_b"] = uh_b, Rcpp::_["bm_rou"] = bm_rou, 
	Rcpp::_["gammasolved"] = gammasolved, Rcpp::_["dlf_calc"] = dlf_calc, 
	Rcpp::_["breaking_moment"] = breaking_moment, Rcpp::_["edge_gap_gust_factor"] = edge_gap_gust_factor, 
	Rcpp::_["dlf_used"] = dlf_used, Rcpp::_["zpd_b"] = zpd) ;
	return uh_b_results;
	
}

// CWS BREAKAGE - TMC VERSION
// note I removed the dlf_threshold as an input variable here, instead modifications of this variable should be done directly in the fgr_constants file.

// [[Rcpp::export]]
Rcpp::DataFrame uh_breakage_tmc_Norway_cpp(Rcpp::NumericVector tree_ht, Rcpp::NumericVector dbh, Rcpp::NumericVector cr_depth, Rcpp::NumericVector cr_width, 
Rcpp::NumericVector cr_prarea, Rcpp::NumericVector crown_off, 
Rcpp::NumericVector spacing_current, Rcpp::NumericVector spacing_before, Rcpp::IntegerVector years_since_thin, Rcpp::NumericVector dist_edge, Rcpp::NumericVector gap_size, 
Rcpp::NumericVector moe, Rcpp::NumericVector mor, Rcpp::NumericVector fknot, 
Rcpp::NumericVector stem_vol, Rcpp::NumericVector stem_density, Rcpp::NumericVector crown_density, Rcpp::NumericVector snow_load, 
Rcpp::StringVector ci, Rcpp::NumericVector ci_value, Rcpp::NumericVector n_drag, Rcpp::NumericVector c_drag, Rcpp::NumericVector drag_upper_limit,
Rcpp::NumericVector equivalent_mean_ht, Rcpp::NumericVector stand_cr_depth, Rcpp::NumericVector stand_cr_width, Rcpp::DataFrame fgr_constants,
Rcpp::NumericVector aerodynamic_ht){
	
	int n = tree_ht.size();
	Rcpp::NumericVector edge_gap_factor(n);
	Rcpp::NumericVector breaking_moment(n);
	Rcpp::NumericVector tmc(n);
	Rcpp::NumericVector tmr_full(n);
	Rcpp::NumericVector dlf_calc(n);
	Rcpp::NumericVector dlf_used(n);
	Rcpp::NumericVector uguess(n);
	Rcpp::NumericVector uguess1(n);
	Rcpp::NumericVector uh_b(n);
	Rcpp::NumericVector bm_tmc(n);

	double wind_precision = fgr_constants["wind_precision"];
	double dlf = fgr_constants["dlf"];
	double dlf_threshold = fgr_constants["dlf_threshold"];
	
	for (int i = 0; i < n; ++i){
		std::string s_ci = Rcpp::as<std::string>(ci[i]);		
		// initialise the windspeed guesses
		breaking_moment[i] = critical_moment_breakage_tmc_cpp(dbh[i], tree_ht[i], cr_depth[i], mor[i], fknot[i]);
		edge_gap_factor[i] = edge_gap_factor_fun_cpp(spacing_current[i], equivalent_mean_ht[i], dist_edge[i], gap_size[i], fgr_constants);
		tmc[i] = tc_zero_intercept_fun_balBA_cpp(dbh[i], tree_ht[i], s_ci, ci_value[i]);
		if(tmc[i] < 1.0){tmc[i] = 1.0;}		
		uguess[i] = sqrt(breaking_moment[i]/tmc[i]);
		uguess1[i] = uguess[i];
		uh_b[i] = uguess[i]/2.0;	
		
		while(abs(uguess1[i] - uh_b[i]) > wind_precision){
			uguess1[i] = uguess[i];
			bm_tmc[i] = tmc[i] * pow(uguess[i], 2);
			dlf_calc[i] = dlf_fun_Norway_cpp(bm_tmc[i], tree_ht[i], cr_depth[i], cr_width[i], cr_prarea[i], crown_off[i], stem_vol[i], dbh[i], moe[i], crown_density[i], stem_density[i], snow_load[i], fgr_constants);
			dlf_used[i] = dlf_calc[i];
			if(dlf_used[i] < 1.0) {dlf_used[i] = dlf;} else if(dlf_used[i] > dlf_threshold) {dlf_used[i] = dlf_threshold;}
			tmr_full[i] = tm_ratio_cpp(spacing_before[i], spacing_current[i], years_since_thin[i], stand_cr_width[i], stand_cr_depth[i], uguess[i], n_drag[i], c_drag[i], drag_upper_limit[i], aerodynamic_ht[i], fgr_constants);
			uh_b[i] = sqrt(breaking_moment[i]/(tmc[i] * dlf_used[i] * tmr_full[i] * edge_gap_factor[i]));
			uguess[i] = uh_b[i];
		}
	}
	
	Rcpp::DataFrame uh_b_results = Rcpp::DataFrame::create( Rcpp::_["uh_b"] = uh_b, Rcpp::_["tmc"] = tmc, 
	Rcpp::_["tmr_full"] = tmr_full, Rcpp::_["dlf_calc"] = dlf_calc, 
	Rcpp::_["breaking_moment"] = breaking_moment, Rcpp::_["edge_gap_factor"] = edge_gap_factor, 
	Rcpp::_["dlf_used"] = dlf_used) ;
	return uh_b_results;
}
	

// [[Rcpp::export]]
Rcpp::DataFrame uh_breakage_tmc_tmr_simple_Norway_cpp(Rcpp::NumericVector tree_ht, Rcpp::NumericVector dbh, Rcpp::NumericVector cr_depth, Rcpp::NumericVector cr_width, 
Rcpp::NumericVector cr_prarea, Rcpp::NumericVector crown_off, 
Rcpp::NumericVector spacing_current, Rcpp::NumericVector spacing_before, Rcpp::IntegerVector years_since_thin, Rcpp::NumericVector dist_edge, Rcpp::NumericVector gap_size, 
Rcpp::NumericVector equivalent_mean_ht, Rcpp::NumericVector moe, Rcpp::NumericVector mor, Rcpp::NumericVector fknot, 
Rcpp::NumericVector stem_vol, Rcpp::NumericVector stem_density, Rcpp::NumericVector crown_density, Rcpp::NumericVector snow_load, 
Rcpp::StringVector ci, Rcpp::NumericVector ci_value, Rcpp::DataFrame fgr_constants){
	
	int n = tree_ht.size();
	Rcpp::NumericVector edge_gap_factor(n);
	Rcpp::NumericVector breaking_moment(n);
	Rcpp::NumericVector tmc(n);
	Rcpp::NumericVector tmr_simple(n);
	Rcpp::NumericVector dlf_calc(n);
	Rcpp::NumericVector dlf_used(n);
	Rcpp::NumericVector uguess(n);
	Rcpp::NumericVector uguess1(n);
	Rcpp::NumericVector uh_b(n);
	Rcpp::NumericVector bm_tmc(n);

	double wind_precision = fgr_constants["wind_precision"];
	double dlf = fgr_constants["dlf"];
	double dlf_threshold = fgr_constants["dlf_threshold"];
	
	for (int i = 0; i < n; ++i){	
		std::string s_ci = Rcpp::as<std::string>(ci[i]);		
		// initialise the windspeed guesses
		breaking_moment[i] = critical_moment_breakage_tmc_cpp(dbh[i], tree_ht[i], cr_depth[i], mor[i], fknot[i]);
		edge_gap_factor[i] = edge_gap_factor_fun_cpp(spacing_current[i], equivalent_mean_ht[i], dist_edge[i], gap_size[i], fgr_constants);
		tmc[i] = tc_zero_intercept_fun_balBA_cpp(dbh[i], tree_ht[i], s_ci, ci_value[i]);
		if(tmc[i] < 1.0){tmc[i] = 1.0;}		
		uguess[i] = sqrt(breaking_moment[i]/tmc[i]);
		uguess1[i] = uguess[i];
		uh_b[i] = uguess[i]/2.0;	
		
		while(abs(uguess1[i] - uh_b[i]) > wind_precision){
			uguess1[i] = uguess[i];
			bm_tmc[i] = tmc[i] * pow(uguess[i], 2);
			dlf_calc[i] = dlf_fun_Norway_cpp(bm_tmc[i], tree_ht[i], cr_depth[i], cr_width[i], cr_prarea[i], crown_off[i], stem_vol[i], dbh[i], moe[i], crown_density[i], stem_density[i], snow_load[i], fgr_constants);
			dlf_used[i] = dlf_calc[i];
			if(dlf_used[i] < 1.0) {dlf_used[i] = dlf;} else if(dlf_used[i] > dlf_threshold) {dlf_used[i] = dlf_threshold;}
			tmr_simple[i] = tm_ratio_simple_cpp(spacing_before[i], spacing_current[i], years_since_thin[i]);
			uh_b[i] = sqrt(breaking_moment[i]/(tmc[i] * dlf_used[i] * tmr_simple[i] * edge_gap_factor[i]));
			uguess[i] = uh_b[i];
		}
	}
	
	Rcpp::DataFrame uh_b_results = Rcpp::DataFrame::create( Rcpp::_["uh_b"] = uh_b, Rcpp::_["tmc"] = tmc, 
	Rcpp::_["tmr_simple"] = tmr_simple, Rcpp::_["dlf_calc"] = dlf_calc, 
	Rcpp::_["breaking_moment"] = breaking_moment, Rcpp::_["edge_gap_factor"] = edge_gap_factor, 
	Rcpp::_["dlf_used"] = dlf_used) ;
	return uh_b_results;
}


// CWS OVERTURNING
/*
Calculates the critical wind speed for overturning using the roughness method.
title Critical wind speed for overturning - Roughness method.
param mean_ht Arithmetic mean height of the trees in the stand (m).
param mean_dbh Mean dbh of all trees in the stand (cm). Dbh is diameter at breast height, measured at 1.3m above the ground. Essential.
param spacing Mean distance between trees (m).
param dist_edge Distance from the upwind edge (m).
param gap_size Length of the upwind gap (m).
param mean_cr_width Width of the crown of the mean tree in the stand (m).
param mean_cr_depth Depth of the crown of the mean tree in the stand (m).
param moe Modulus of Elasticity of green wood (MPa).
param c_reg Regression coefficients of uprooting moment against stem weight (N m kg-1).
param n_drag N parameter of the drag coefficient formula (dimensionless).
param c_drag C parameter of the drag coefficient formula (dimensionless).
param drag_upper_limit Maximum wind speed used during the experiments from which \code{n_drag} and \code{c_drag} were derived (m*s-1).
param stem_vol Stem volume of the mean tree in the stand (m^3).
param stem_density Density of green wood of the stem (kg m-3).
param crown_density Crown density of the mean tree in the stand (kg m-3).
param snow_depth Depth of layer of snow on tree crown (cm).
param snow_density Density of snow (kg m-3).
param ro Air density (kg m-3).
return A list including: The critical wind speed (m s-1) for overturning at zero plane displacement height for the roughness method;
The applied bending moment; The zero plane displacement height; The ratio of critical wind speed at canopy top over friction velocity
(uh / u*); The deflection loading factor; The critical overturning moment; The combined effect of edge, gap, and gust on the applied
bending moment.
*/

// [[Rcpp::export]]
Rcpp::DataFrame uh_overturning_rou_Norway_cpp(Rcpp::NumericVector mean_ht, Rcpp::NumericVector mean_dbh, Rcpp::NumericVector spacing, Rcpp::NumericVector dist_edge, Rcpp::NumericVector gap_size, 
Rcpp::NumericVector mean_cr_width, Rcpp::NumericVector mean_cr_depth, Rcpp::NumericVector mean_cr_prarea, Rcpp::NumericVector crown_off, Rcpp::NumericVector moe, Rcpp::NumericVector c_reg, Rcpp::NumericVector n_drag, Rcpp::NumericVector c_drag, Rcpp::NumericVector drag_upper_limit,
Rcpp::NumericVector stem_vol, Rcpp::NumericVector stem_density, Rcpp::NumericVector crown_density, Rcpp::NumericVector snow_load, double ro, Rcpp::DataFrame fgr_constants, Rcpp::NumericVector aerodynamic_ht) {
	
	int n = mean_ht.size();
	Rcpp::NumericVector bm_rou(n);
	Rcpp::NumericVector dlf_calc(n);
	Rcpp::NumericVector dlf_used(n);
	Rcpp::NumericVector zpd(n);
	Rcpp::NumericVector gammasolved(n);
	Rcpp::NumericVector uh_o(n);
	Rcpp::NumericVector edge_gap_gust_factor(n);
	Rcpp::NumericVector overturning_moment(n);

	double wind_precision = fgr_constants["wind_precision"];
	double dlf = fgr_constants["dlf"];
	double dlf_threshold = fgr_constants["dlf_threshold"];
		
	for (int i = 0; i < n; ++i){
		edge_gap_gust_factor[i] = edge_gap_gust_factor_fun_cpp(spacing[i], mean_ht[i], dist_edge[i], gap_size[i], fgr_constants);
		overturning_moment[i] = critical_moment_overturning_cpp(c_reg[i], stem_density[i], stem_vol[i]);	
		
		// initialise the windspeed guesses
		double uguess = 25.0;
		double uguess1 = uguess;
		uh_o[i] = uguess/2.0;		
		
		while(abs(uguess1 - uh_o[i]) > wind_precision){
			uguess1 = uguess;
			bm_rou[i] = bending_moment_rou_cpp(mean_dbh[i], mean_ht[i], mean_cr_width[i], mean_cr_depth[i], spacing[i], dist_edge[i], gap_size[i], uguess, n_drag[i], c_drag[i], drag_upper_limit[i], ro, fgr_constants, aerodynamic_ht[i]);
			dlf_calc[i] = dlf_fun_Norway_cpp(bm_rou[i], mean_ht[i], mean_cr_depth[i], mean_cr_width[i], mean_cr_prarea[i], crown_off[i], stem_vol[i], mean_dbh[i], moe[i], crown_density[i], stem_density[i], snow_load[i], fgr_constants);
			dlf_used[i] = dlf_calc[i];
			if(dlf_used[i] < 0.0) {dlf_used[i] = dlf;} else if(dlf_used[i] > dlf_threshold) {dlf_used[i] = dlf_threshold;}
			zpd[i] = zpd_fun_cpp(mean_cr_width[i], mean_cr_depth[i], spacing[i], uguess, n_drag[i], c_drag[i], drag_upper_limit[i], aerodynamic_ht[i], fgr_constants);
			if(zpd[i] < 1.3){zpd[i] = 1.35;}
			gammasolved[i] = gammasolved_fun_cpp(mean_cr_width[i], mean_cr_depth[i], spacing[i], uguess, n_drag[i], c_drag[i], drag_upper_limit[i], fgr_constants);
			uh_o[i] = (1.0/spacing[i]) * sqrt(overturning_moment[i] /  (ro * edge_gap_gust_factor[i] * zpd[i] * dlf_used[i])) * gammasolved[i];
			uguess = uh_o[i];
		}
	}
	
	Rcpp::DataFrame uh_o_results = Rcpp::DataFrame::create( Rcpp::_["uh_o"] = uh_o, Rcpp::_["bm_rou"] = bm_rou, 
	Rcpp::_["gammasolved"] = gammasolved, Rcpp::_["dlf_calc"] = dlf_calc, 
	Rcpp::_["overturning_moment"] = overturning_moment, Rcpp::_["edge_gap_gust_factor"] = edge_gap_gust_factor, 
	Rcpp::_["dlf_used"] = dlf_used, Rcpp::_["zpd_o"] = zpd) ;
	return uh_o_results;
	
}

// CWS OVERTURNING - TMC VERSION
// [[Rcpp::export]]
Rcpp::DataFrame uh_overturning_tmc_Norway_cpp(Rcpp::NumericVector tree_ht, Rcpp::NumericVector dbh, Rcpp::NumericVector cr_depth, Rcpp::NumericVector cr_width, 
Rcpp::NumericVector cr_prarea, Rcpp::NumericVector crown_off, 
Rcpp::NumericVector spacing_current, Rcpp::NumericVector spacing_before, Rcpp::IntegerVector years_since_thin, Rcpp::NumericVector dist_edge, Rcpp::NumericVector gap_size, 
Rcpp::NumericVector moe, Rcpp::NumericVector c_reg, Rcpp::NumericVector stem_vol, Rcpp::NumericVector stem_density, Rcpp::NumericVector crown_density, Rcpp::NumericVector snow_load, 
Rcpp::StringVector ci, Rcpp::NumericVector ci_value, Rcpp::NumericVector n_drag, Rcpp::NumericVector c_drag, Rcpp::NumericVector drag_upper_limit,
Rcpp::NumericVector equivalent_mean_ht, Rcpp::NumericVector stand_cr_depth, Rcpp::NumericVector stand_cr_width, Rcpp::DataFrame fgr_constants,
Rcpp::NumericVector aerodynamic_ht){
	
	int n = tree_ht.size();
	Rcpp::NumericVector edge_gap_factor(n);
	Rcpp::NumericVector overturning_moment(n);
	Rcpp::NumericVector tmc(n);
	Rcpp::NumericVector tmr_full(n);
	Rcpp::NumericVector dlf_calc(n);
	Rcpp::NumericVector dlf_used(n);
	Rcpp::NumericVector uguess(n);
	Rcpp::NumericVector uguess1(n);
	Rcpp::NumericVector uh_o(n);
	Rcpp::NumericVector bm_tmc(n);

	double wind_precision = fgr_constants["wind_precision"];
	double dlf = fgr_constants["dlf"];
	double dlf_threshold = fgr_constants["dlf_threshold"];
	
	for (int i = 0; i < n; ++i){		
		std::string s_ci = Rcpp::as<std::string>(ci[i]);
		// initialise the windspeed guesses
		overturning_moment[i] = critical_moment_overturning_cpp(c_reg[i], stem_density[i], stem_vol[i]);
		edge_gap_factor[i] = edge_gap_factor_fun_cpp(spacing_current[i], equivalent_mean_ht[i], dist_edge[i], gap_size[i], fgr_constants);
		tmc[i] = tc_zero_intercept_fun_balBA_cpp(dbh[i], tree_ht[i], s_ci, ci_value[i]);
		if(tmc[i] < 1.0){tmc[i] = 1.0;}		
		uguess[i] = sqrt(overturning_moment[i]/tmc[i]);
		uguess1[i] = uguess[i];
		uh_o[i] = uguess[i]/2.0;	
		
		while(abs(uguess1[i] - uh_o[i]) > wind_precision){
			uguess1[i] = uguess[i];
			bm_tmc[i] = tmc[i] * pow(uguess[i], 2);
			dlf_calc[i] = dlf_fun_Norway_cpp(bm_tmc[i], tree_ht[i], cr_depth[i], cr_width[i], cr_prarea[i], crown_off[i], stem_vol[i], dbh[i], moe[i], crown_density[i], stem_density[i], snow_load[i], fgr_constants);
			dlf_used[i] = dlf_calc[i];
			if(dlf_used[i] < 1.0) {dlf_used[i] = dlf;} else if(dlf_used[i] > dlf_threshold) {dlf_used[i] = dlf_threshold;}
			tmr_full[i] = tm_ratio_cpp(spacing_before[i], spacing_current[i], years_since_thin[i], stand_cr_width[i], stand_cr_depth[i], uguess[i], n_drag[i], c_drag[i], drag_upper_limit[i], aerodynamic_ht[i], fgr_constants);
			uh_o[i] = sqrt(overturning_moment[i]/(tmc[i] * dlf_used[i] * tmr_full[i] * edge_gap_factor[i]));
			uguess[i] = uh_o[i];
		}
	}
	
	Rcpp::DataFrame uh_o_results = Rcpp::DataFrame::create( Rcpp::_["uh_o"] = uh_o, Rcpp::_["tmc"] = tmc, 
	Rcpp::_["tmr_full"] = tmr_full, Rcpp::_["dlf_calc"] = dlf_calc, 
	Rcpp::_["overturning_moment"] = overturning_moment, Rcpp::_["edge_gap_factor"] = edge_gap_factor, 
	Rcpp::_["dlf_used"] = dlf_used) ;
	return uh_o_results;
}


// [[Rcpp::export]]
Rcpp::DataFrame uh_overturning_tmc_tmr_simple_cpp(Rcpp::NumericVector tree_ht, Rcpp::NumericVector dbh, Rcpp::NumericVector cr_depth, Rcpp::NumericVector cr_width, 
Rcpp::NumericVector cr_prarea, Rcpp::NumericVector crown_off, 
Rcpp::NumericVector spacing_current, Rcpp::NumericVector spacing_before, Rcpp::IntegerVector years_since_thin, Rcpp::NumericVector dist_edge, Rcpp::NumericVector gap_size, 
Rcpp::NumericVector equivalent_mean_ht, Rcpp::NumericVector moe, Rcpp::NumericVector c_reg, Rcpp::NumericVector stem_vol, Rcpp::NumericVector stem_density, Rcpp::NumericVector crown_density, 
Rcpp::NumericVector snow_load, Rcpp::StringVector ci, Rcpp::NumericVector ci_value, Rcpp::DataFrame fgr_constants){
	
	int n = tree_ht.size();
	Rcpp::NumericVector edge_gap_factor(n);
	Rcpp::NumericVector overturning_moment(n);
	Rcpp::NumericVector tmc(n);
	Rcpp::NumericVector tmr_simple(n);
	Rcpp::NumericVector dlf_calc(n);
	Rcpp::NumericVector dlf_used(n);
	Rcpp::NumericVector uguess(n);
	Rcpp::NumericVector uguess1(n);
	Rcpp::NumericVector uh_o(n);
	Rcpp::NumericVector bm_tmc(n);

	double wind_precision = fgr_constants["wind_precision"];
	double dlf = fgr_constants["dlf"];
	double dlf_threshold = fgr_constants["dlf_threshold"];
	
	for (int i = 0; i < n; ++i){	
		std::string s_ci = Rcpp::as<std::string>(ci[i]);
		// initialise the windspeed guesses
		overturning_moment[i] = critical_moment_overturning_cpp(c_reg[i], stem_density[i], stem_vol[i]);
		edge_gap_factor[i] = edge_gap_factor_fun_cpp(spacing_current[i], equivalent_mean_ht[i], dist_edge[i], gap_size[i], fgr_constants);
		tmc[i] = tc_zero_intercept_fun_balBA_cpp(dbh[i], tree_ht[i], s_ci, ci_value[i]);
		if(tmc[i] < 1.0){tmc[i] = 1.0;}		
		uguess[i] = sqrt(overturning_moment[i]/tmc[i]);
		uguess1[i] = uguess[i];
		uh_o[i] = uguess[i]/2.0;	
		
		while(abs(uguess1[i] - uh_o[i]) > wind_precision){
			uguess1[i] = uguess[i];
			bm_tmc[i] = tmc[i] * pow(uguess[i], 2);
			dlf_calc[i] = dlf_fun_Norway_cpp(bm_tmc[i], tree_ht[i], cr_depth[i], cr_width[i], cr_prarea[i], crown_off[i], stem_vol[i], dbh[i], moe[i], crown_density[i], stem_density[i], snow_load[i], fgr_constants);
			dlf_used[i] = dlf_calc[i];
			if(dlf_used[i] < 1.0) {dlf_used[i] = dlf;} else if(dlf_used[i] > dlf_threshold) {dlf_used[i] = dlf_threshold;}
			tmr_simple[i] = tm_ratio_simple_cpp(spacing_before[i], spacing_current[i], years_since_thin[i]);
			uh_o[i] = sqrt(overturning_moment[i]/(tmc[i] * dlf_used[i] * tmr_simple[i] * edge_gap_factor[i]));
			uguess[i] = uh_o[i];
		}
	}
	
	Rcpp::DataFrame uh_o_results = Rcpp::DataFrame::create( Rcpp::_["uh_o"] = uh_o, Rcpp::_["tmc"] = tmc, 
	Rcpp::_["tmr_simple"] = tmr_simple, Rcpp::_["dlf_calc"] = dlf_calc, 
	Rcpp::_["overturning_moment"] = overturning_moment, Rcpp::_["edge_gap_factor"] = edge_gap_factor, 
	Rcpp::_["dlf_used"] = dlf_used) ;
	return uh_o_results;
}

// CRITICAL SNOW LOAD
/* The critical snow load can be calculated after some simplifications made to the deflection loading factor. The simplifications calculates the additional moment created by a displaced crown.
In short, after considering the contribution from the stem weight to be negligible, the deflection function and the bending moment of the wind of the trees can be simplified assuming the wind acts
at the centre of the tree crown.
If the denominator of the simplified deflection loading factor is zero, then the tree will collapse under its own weight without any wind loading, and thus a critical snow load can be calculated.
See Critical Wind Speed and Critical Snow Load.pdf in T:\Aktive\DSU\3737\53280_Future_Forests\Akt2_Climate_Disturbances\T2_2_Snow and  https://doi.org/10.1080/02827581.2023.2289660 for an application
of this model. 
*/
// [[Rcpp::export]]
double csw_Norway_cpp(double tree_ht, double tree_dbh, double crown_weight, double cr_depth, double moe, double grav) {
	// Note that we consider the crown offset to be negligible for the calculation of the critical snow load due to the
	// simplifications made to the DLF to obtain the critical snow load.
	// Note that DBH must be specified in meters
	double param_a;
	double param_b;
	double csw;
	
	param_a = (64 * pow(tree_ht, 2) * grav)/(M_PI * moe * pow(tree_dbh/100, 4));
	param_b = 5.0 * (cr_depth/(2 * tree_ht)) - 5.96 * pow(cr_depth/(2 * tree_ht), 0.6) - 0.71 * pow(cr_depth/(2 * tree_ht), 2) + 1.67;
	
	csw = 1/(param_a * param_b) - crown_weight;
	return csw;	
}


// ELEVATE TO ANENOMETER HEIGHT
/*
Converts the critical wind speeds to the wind speed 10 meters above zero plane displacement (for anemometers' data). The wind is assumed to follow a logarithmic profile.
title Elevate
param uh Critical wind speed at canopy top calculated with the roughness or single-tree method (m s-1).
param z0 Aerodynamic roughness length of the canopy (m).
param d Zero plane displacement height of the canopy (m).
param ht Height at which the critical wind speed is calculated (m). Note that this differs between the 'roughness' (stand mean height) and 'TMC' (1.05 * stand top height) methods.
return A critical wind speed that can be compared with met data.
*/
// [[Rcpp::export]]
double elevate_cpp(double uh, double z0, double d, double ht, double ht_above_d) {
	
	double u_elev = uh * log(ht_above_d/z0)/log((ht - d)/z0);
	return u_elev;
}

// Calculate breakage at different height along the trunk
// [[Rcpp::export]]
Rcpp::DataFrame uh_breakage_heightz_cpp(Rcpp::NumericVector uh_b, Rcpp::NumericVector tree_ht, Rcpp::NumericVector cr_depth, Rcpp::NumericVector dbh, double m = 2.0/3.0){
	
	int nsize = tree_ht.size();
	Rcpp::NumericVector hwind(nsize);
	Rcpp::NumericVector hmin_cws(nsize);
	Rcpp::NumericVector d0(nsize);
	Rcpp::NumericVector diam_htz(nsize);
	Rcpp::NumericVector uh_b_min(nsize);
	
	for (int i = 0; i < nsize; ++i){
		hwind[i] = tree_ht[i] - (cr_depth[i]/2.0);
		if (hwind[i] < 1.3 * 2.0) {hwind[i] = 1.3 * 2;}
		hmin_cws[i] = tree_ht[i] - (cr_depth[i]/2) * (1 + (1/(3 * m - 1)));
		d0[i] = diam_base_fun_cpp(dbh[i], tree_ht[i], cr_depth[i]);
		diam_htz[i] = d0[i] * pow((tree_ht[i] - hmin_cws[i])/tree_ht[i], m);
		uh_b_min[i] = uh_b[i] * sqrt(pow(diam_htz[i], 3)/pow(dbh[i], 3)) * sqrt((hwind[i] - 1.3)/(hwind[i] - hmin_cws[i]));
	}
	
	Rcpp::DataFrame uh_b_min_results = Rcpp::DataFrame::create( Rcpp::_["uh_b_min"] = uh_b_min, Rcpp::_["hmin_cws"] = hmin_cws, 
	Rcpp::_["diam_htz"] = diam_htz) ;
	return uh_b_min_results;
}


// ANNUAL EXCEEDANCE PROBABILITY - TO DO


// MAIN FUNCTION - ROUGHNESS VERSION

// [[Rcpp::export]]    
Rcpp::List fg_rou_Norway_cpp(Rcpp::List inputdata_full, Rcpp::DataFrame fgr_constants, Rcpp::DataFrame species_parameters,
Rcpp::Nullable<Rcpp::NumericVector> full_output_ = R_NilValue, std::string breakage_basecanopy = "no"){
	// Check that the data has been properly prepared and that all the variables are present
	std::cout << "Make sure to have run your data through the data preparation procedure first."<<std::endl;
	
	// We extract all the required variables to run through the ForestGales function
	Rcpp::DataFrame required(inputdata_full["required"]);
	Rcpp::DataFrame desirable(inputdata_full["desirable"]);
	Rcpp::DataFrame advanced(inputdata_full["advanced"]);
	Rcpp::DataFrame others(inputdata_full["others"]);
	
	// Identifier variables
	Rcpp::StringVector stand_id = required["stand_id"];
	Rcpp::StringVector species = required["species"];
	
	
	// Extract the required variables
	Rcpp::NumericVector mean_ht = desirable["mean_ht"];
	Rcpp::NumericVector aerodynamic_ht = desirable["aerodynamic_ht"];
	Rcpp::NumericVector mean_dbh = required["mean_dbh"];
	Rcpp::NumericVector spacing = required["spacing"];
	Rcpp::NumericVector dist_edge = desirable["dist_edge"];
	Rcpp::NumericVector gap_size = desirable["gap_size"];
	Rcpp::NumericVector mean_cr_width = desirable["mean_cr_width"];
	Rcpp::NumericVector mean_cr_depth = desirable["mean_cr_depth"];
	Rcpp::NumericVector mean_cr_prarea = desirable["mean_cr_prarea"];
	Rcpp::NumericVector crown_off = desirable["crown_off"];
	Rcpp::NumericVector elev_ht = advanced["elev_ht"];
	Rcpp::NumericVector ht_above_d = advanced["ht_above_d"];
	Rcpp::NumericVector moe = advanced["moe"];
	Rcpp::NumericVector mor = advanced["mor"];
	Rcpp::NumericVector c_reg = advanced["c_reg"];
	Rcpp::NumericVector fknot = advanced["fknot"];
	Rcpp::NumericVector n_drag = advanced["n_drag"];
	Rcpp::NumericVector c_drag = advanced["c_drag"];
	Rcpp::NumericVector drag_upper_limit = advanced["drag_upper_limit"];
	Rcpp::NumericVector stem_vol = advanced["stem_vol"];
	Rcpp::NumericVector stem_density = advanced["stem_density"];
	Rcpp::NumericVector crown_density = advanced["crown_density"];
	Rcpp::NumericVector snow_load = others["snow_load"];
	
	
	double ro = fgr_constants["ro_default"];
	
	// Get the length of the dataframe
	//int n = stand_id.size();
	
	// Calculate the critical wind speeds for breakage and for overturning
	Rcpp::DataFrame uh_b_results = uh_breakage_rou_Norway_cpp(mean_ht, mean_dbh, spacing, 
	dist_edge, gap_size, mean_cr_width, mean_cr_depth, mean_cr_prarea, crown_off, moe, 
	mor, fknot, n_drag, c_drag, drag_upper_limit, stem_vol, 
	stem_density, crown_density, snow_load, ro,
	fgr_constants, aerodynamic_ht);
	
	
	Rcpp::DataFrame uh_o_results = uh_overturning_rou_Norway_cpp(mean_ht, mean_dbh, spacing, 
	dist_edge, gap_size, mean_cr_width, mean_cr_depth, mean_cr_prarea, crown_off, moe, 
	c_reg, n_drag, c_drag, drag_upper_limit, stem_vol, 
	stem_density, crown_density, snow_load, ro,
	fgr_constants, aerodynamic_ht);
	
	// Extract the variables from these results
	Rcpp::NumericVector uh_b = uh_b_results["uh_b"];
	Rcpp::NumericVector uh_o = uh_o_results["uh_o"];
	Rcpp::NumericVector zpd_b = uh_b_results["zpd_b"];
	Rcpp::NumericVector zpd_o = uh_o_results["zpd_o"];
	
	// Elevate the value
	int n = uh_b.size();
	Rcpp::NumericVector z0_b(n);
	Rcpp::NumericVector z0_o(n);
	Rcpp::NumericVector u_elev_b(n);
	Rcpp::NumericVector u_elev_o(n);
	Rcpp::NumericVector u_damage(n);
	Rcpp::NumericVector u_elev_damage(n);
	Rcpp::StringVector mode_of_damage(n);
	Rcpp::NumericVector lambdacapital_b(n);
	Rcpp::NumericVector drag_b(n);
	Rcpp::NumericVector lambdacapital_o(n);
	Rcpp::NumericVector drag_o(n);
	
	for (int i = 0; i < n; ++i){
		z0_b[i] = z0_fun_cpp(mean_cr_width[i], mean_cr_depth[i], spacing[i], uh_b[i], 
		n_drag[i], c_drag[i], drag_upper_limit[i], aerodynamic_ht[i], fgr_constants);
		z0_o[i] = z0_fun_cpp(mean_cr_width[i], mean_cr_depth[i], spacing[i], uh_o[i],
		n_drag[i], c_drag[i], drag_upper_limit[i], aerodynamic_ht[i], fgr_constants);

		u_elev_b[i] = elevate_cpp(uh_b[i], z0_b[i], zpd_b[i], elev_ht[i], ht_above_d[i]);
		u_elev_o[i] = elevate_cpp(uh_o[i], z0_o[i], zpd_o[i], elev_ht[i], ht_above_d[i]);
		
		if(uh_b[i] < uh_o[i]){u_damage[i] = uh_b[i];} else {u_damage[i] = uh_o[i];} 
		if(u_elev_b[i] < u_elev_o[i]){u_elev_damage[i] = u_elev_b[i];} else {u_elev_damage[i] = u_elev_o[i];} 
		
		if(u_elev_b[i] < u_elev_o[i]){mode_of_damage[i] = "Breakage";} else {mode_of_damage[i] = "Overturning";}
		
		
		lambdacapital_b[i] = lambdacapital_fun_cpp(mean_cr_width[i], mean_cr_depth[i], 
		spacing[i], uh_b[i], n_drag[i], c_drag[i], drag_upper_limit[i]);
		drag_b[i] = drag_fun_cpp(uh_b[i], n_drag[i], c_drag[i], drag_upper_limit[i]);
		
		lambdacapital_o[i] = lambdacapital_fun_cpp(mean_cr_width[i], mean_cr_depth[i], 
		spacing[i], uh_o[i], n_drag[i], c_drag[i], drag_upper_limit[i]);
		drag_o[i] = drag_fun_cpp(uh_o[i], n_drag[i], c_drag[i], drag_upper_limit[i]);
	}
	
	// If breakage at the base of the stem
	Rcpp::DataFrame output_data_Bmin;
	if(breakage_basecanopy == "yes"){
		
		Rcpp::DataFrame uh_b_min_results = uh_breakage_heightz_cpp(uh_b, mean_ht, mean_cr_depth, mean_dbh);
		Rcpp::NumericVector uh_b_min = uh_b_min_results["uh_b_min"];
		Rcpp::NumericVector hmin_cws = uh_b_min_results["hmin_cws"];
		Rcpp::NumericVector diam_htz = uh_b_min_results["diam_htz"];	

		Rcpp::NumericVector zpd_b_min(n);
		Rcpp::NumericVector z0_b_min(n);
		Rcpp::NumericVector u_elev_b_min(n);
		Rcpp::NumericVector lambdacapital_b_min(n);
		Rcpp::NumericVector gammasolved_b_min(n);
		Rcpp::NumericVector drag_b_min(n);
		
		for (int i = 0; i < n; ++i){
			zpd_b_min[i] = zpd_fun_cpp(mean_cr_width[i], mean_cr_depth[i], spacing[i], uh_b_min[i], 
			n_drag[i], c_drag[i], drag_upper_limit[i], aerodynamic_ht[i], fgr_constants);
			z0_b_min[i] = z0_fun_cpp(mean_cr_width[i], mean_cr_depth[i], spacing[i], uh_b_min[i], 
			n_drag[i], c_drag[i], drag_upper_limit[i], aerodynamic_ht[i], fgr_constants);

			u_elev_b_min[i] = elevate_cpp(uh_b_min[i], z0_b_min[i], zpd_b_min[i], elev_ht[i], ht_above_d[i]);
			lambdacapital_b_min[i] = lambdacapital_fun_cpp(mean_cr_width[i], mean_cr_depth[i], 
			spacing[i], uh_b_min[i], n_drag[i], c_drag[i], drag_upper_limit[i]);
			drag_b_min[i] = drag_fun_cpp(uh_b_min[i], n_drag[i], c_drag[i], drag_upper_limit[i]);
			gammasolved_b_min = gammasolved_fun_cpp(mean_cr_width[i], mean_cr_depth[i],spacing[i], uh_b_min[i], 
			n_drag[i], c_drag[i], drag_upper_limit[i], fgr_constants);
			if(u_elev_damage[i] > u_elev_b_min[i]){
				u_elev_damage[i] = u_elev_b_min[i];
				u_damage[i] = uh_b_min[i];
				mode_of_damage[i] = "Breakage_base_canopy";
			}
		}
		output_data_Bmin = Rcpp::DataFrame::create(Rcpp::_["u_elev_b_min"] = Rcpp::clone(u_elev_b_min), Rcpp::_["uh_b_min"] = uh_b_min_results["uh_b_min"],
		Rcpp::_["zpd_b_min"] = Rcpp::clone(zpd_b_min), Rcpp::_["z0_b_min"] = Rcpp::clone(z0_b_min), Rcpp::_["drag_b_min"] = Rcpp::clone(drag_b_min),
		Rcpp::_["gammasolved_b_min"] = Rcpp::clone(gammasolved_b_min), Rcpp::_["lambdacapital_b_min"] = Rcpp::clone(lambdacapital_b_min),
		Rcpp::_["u_elev_damage"] = Rcpp::clone(u_elev_damage), Rcpp::_["mode_of_damage"] = Rcpp::clone(mode_of_damage));		
	}
	
	// Calculate probability of breakage, overturning, and damage - to do
	
	// Add stems per hectare
	others["sph"] = round(10000.0/(pow(spacing,2)), 0);
	
	
	Rcpp::List output_data;	
	
	// Choose whether to return full dataframe - including the input data - or just the results dataframe
	if (full_output_.isNotNull()) {
        Rcpp::NumericVector full_output(full_output_);        // casting to underlying type NumericVector
		output_data["required"] = Rcpp::clone(required);
		output_data["desirable"] = Rcpp::clone(desirable);
		output_data["advanced"] = Rcpp::clone(advanced);
		output_data["others"] = Rcpp::clone(others);
		output_data["results_BO"] = Rcpp::DataFrame::create(Rcpp::_["u_elev_b"] = Rcpp::clone(u_elev_b), Rcpp::_["uh_b"] = uh_b_results["uh_b"],
		//Rcpp::_["prob_b"] = Rcpp::clone(prob_b), 
		Rcpp::_["zpd_b"] = Rcpp::clone(zpd_b),
		Rcpp::_["z0_b"] = Rcpp::clone(z0_b), Rcpp::_["drag_b"] = Rcpp::clone(drag_b),
		Rcpp::_["gammasolved_b"] = uh_b_results["gammasolved"], Rcpp::_["lambdacapital_b"] = Rcpp::clone(lambdacapital_b),
		Rcpp::_["breaking_moment"] = uh_b_results["breaking_moment"],
		Rcpp::_["u_elev_o"] = Rcpp::clone(u_elev_o), Rcpp::_["uh_o"] = uh_o_results["uh_o"],
		//Rcpp::_["prob_o"] = Rcpp::clone(prob_o), 
		Rcpp::_["zpd_o"] = Rcpp::clone(zpd_o),
		Rcpp::_["z0_o"] = Rcpp::clone(z0_o), Rcpp::_["drag_o"] = Rcpp::clone(drag_o),
		Rcpp::_["gammasolved_o"] = uh_o_results["gammasolved"], Rcpp::_["lambdacapital_o"] = Rcpp::clone(lambdacapital_o),
		Rcpp::_["overturning_moment"] = uh_o_results["overturning_moment"]);
		
		if(breakage_basecanopy == "yes"){
			output_data["results_Bmin"] = Rcpp::clone(output_data_Bmin);
			output_data["results_damage"] =  Rcpp::DataFrame::create(Rcpp::_["mode_of_damage"] = output_data_Bmin["mode_of_damage"], 
			Rcpp::_["u_elev_damage"] = output_data_Bmin["u_elev_damage"], Rcpp::_["u_damage"] = output_data_Bmin["u_damage"], 
			Rcpp::_["bm_rou"] = uh_b_results["bm_rou"], Rcpp::_["dlf_calc"] = uh_b_results["dlf_calc"],
			Rcpp::_["dlf_used"] = uh_b_results["dlf_used"], Rcpp::_["edge_gap_gust_factor"] = uh_b_results["edge_gap_gust_factor"]);	
		} else {
			output_data["results_damage"] =  Rcpp::DataFrame::create(
			Rcpp::_["mode_of_damage"] = Rcpp::clone(mode_of_damage), Rcpp::_["u_elev_damage"] = Rcpp::clone(u_elev_damage),
			Rcpp::_["u_damage"] = Rcpp::clone(u_damage), 
			//Rcpp::_["prob_damage"] = Rcpp::clone(prob_damage),
			Rcpp::_["bm_rou"] = uh_b_results["bm_rou"], Rcpp::_["dlf_calc"] = uh_b_results["dlf_calc"],
			Rcpp::_["dlf_used"] = uh_b_results["dlf_used"], Rcpp::_["edge_gap_gust_factor"] = uh_b_results["edge_gap_gust_factor"]);
		}		
	} else {
		if(breakage_basecanopy == "yes"){
			output_data["results"] = Rcpp::DataFrame::create(Rcpp::_["stand_id"] = Rcpp::clone(stand_id), Rcpp::_["species"] = Rcpp::clone(species),
			Rcpp::_["u_elev_damage"] = output_data_Bmin["u_elev_damage"], Rcpp::_["mode_of_damage"] = output_data_Bmin["mode_of_damage"],
			Rcpp::_["u_elev_b"] = Rcpp::clone(u_elev_b), Rcpp::_["u_elev_b_min"] = output_data_Bmin["u_elev_b_min"], Rcpp::_["u_elev_o"] = Rcpp::clone(u_elev_o),
			Rcpp::_["max_stem_weight_warning"] = others["max_stem_weight_warning"]);			
		} else {
			output_data["results"] = Rcpp::DataFrame::create(Rcpp::_["stand_id"] = Rcpp::clone(stand_id), Rcpp::_["species"] = Rcpp::clone(species),
			Rcpp::_["u_elev_damage"] = Rcpp::clone(u_elev_damage), Rcpp::_["mode_of_damage"] = Rcpp::clone(mode_of_damage),
			Rcpp::_["u_elev_b"] = Rcpp::clone(u_elev_b), Rcpp::_["u_elev_o"] = Rcpp::clone(u_elev_o),
			//Rcpp::_["prob_b"] = Rcpp::clone(prob_b), Rcpp::_["prob_o"] = Rcpp::clone(prob_o), Rcpp::_["prob_damage"] = Rcpp::clone(prob_damage),
			Rcpp::_["max_stem_weight_warning"] = others["max_stem_weight_warning"]);
		}
	}
	
	return output_data;		
}




// MAIN FUNCTION - TMC VERSION

// [[Rcpp::export]]    
Rcpp::List fg_tmc_Norway_cpp(Rcpp::List inputdata_full, Rcpp::DataFrame fgr_constants, Rcpp::DataFrame species_parameters,
Rcpp::Nullable<Rcpp::NumericVector> full_output_ = R_NilValue, std::string breakage_basecanopy = "no"){
	// Check that the data has been properly prepared and that all the variables are present
	std::cout << "Make sure to have run your data through the data preparation procedure first."<<std::endl;
	
	// We extract all the required variables to run through the ForestGales function
	Rcpp::DataFrame required(inputdata_full["required"]);
	Rcpp::DataFrame desirable(inputdata_full["desirable"]);
	Rcpp::DataFrame advanced(inputdata_full["advanced"]);
	Rcpp::DataFrame others(inputdata_full["others"]);
	
	// Identifier variables
	Rcpp::StringVector stand_id = required["stand_id"];
	Rcpp::StringVector tree_id = required["tree_id"];
	Rcpp::StringVector species = required["species"];
	
	
	// Extract the required variables
	// Required
	Rcpp::NumericVector tree_ht = required["tree_ht"];
	Rcpp::NumericVector dbh = required["dbh"];
	Rcpp::NumericVector spacing_current = required["spacing_current"];
	Rcpp::NumericVector stand_mean_dbh = required["stand_mean_dbh"];
	Rcpp::NumericVector stand_top_ht = required["stand_top_ht"];
	Rcpp::NumericVector equivalent_mean_ht = required["equivalent_mean_ht"];
	// Desirable
	Rcpp::NumericVector cr_width = desirable["cr_width"];
	Rcpp::NumericVector cr_depth = desirable["cr_depth"];
	Rcpp::NumericVector cr_prarea = desirable["cr_prarea"];
	Rcpp::NumericVector crown_off = desirable["crown_off"];
	Rcpp::NumericVector stand_cr_width = desirable["stand_cr_width"];
	Rcpp::NumericVector stand_cr_depth = desirable["stand_cr_depth"];
	Rcpp::NumericVector dist_edge = desirable["dist_edge"];
	Rcpp::NumericVector gap_size = desirable["gap_size"];
	// Advanced
	Rcpp::StringVector ci = advanced["ci"];
	Rcpp::NumericVector ci_value = advanced["ci_value"];
	Rcpp::NumericVector spacing_before = advanced["spacing_before"];
	Rcpp::IntegerVector years_since_thin = advanced["years_since_thin"];
	Rcpp::NumericVector elev_ht = advanced["elev_ht"];
	Rcpp::NumericVector aerodynamic_ht = advanced["aerodynamic_ht"];
	Rcpp::NumericVector ht_above_d = advanced["ht_above_d"];
	Rcpp::NumericVector moe = advanced["moe"];
	Rcpp::NumericVector mor = advanced["mor"];
	Rcpp::NumericVector c_reg = advanced["c_reg"];
	Rcpp::NumericVector fknot = advanced["fknot"];
	Rcpp::NumericVector n_drag = advanced["n_drag"];
	Rcpp::NumericVector c_drag = advanced["c_drag"];
	Rcpp::NumericVector drag_upper_limit = advanced["drag_upper_limit"];
	Rcpp::NumericVector stem_vol = advanced["stem_vol"];
	Rcpp::NumericVector stem_density = advanced["stem_density"];
	Rcpp::NumericVector crown_density = advanced["crown_density"];
	Rcpp::NumericVector snow_load = others["snow_load"];
		
	// Calculate the critical wind speeds for breakage and for overturning	
	Rcpp::DataFrame uh_b_results = uh_breakage_tmc_Norway_cpp(tree_ht, dbh, cr_depth, 
	cr_width, cr_prarea,  crown_off, spacing_current, spacing_before, years_since_thin, 
	dist_edge, gap_size, moe, mor, fknot, stem_vol, stem_density, 
	crown_density, snow_load, ci, ci_value, 
	n_drag, c_drag, drag_upper_limit, equivalent_mean_ht, 
	stand_cr_depth, stand_cr_width, fgr_constants, aerodynamic_ht);
	
	Rcpp::DataFrame uh_o_results = uh_overturning_tmc_Norway_cpp(tree_ht, dbh, cr_depth, 
	cr_width, cr_prarea,  crown_off, spacing_current, spacing_before, years_since_thin, 
	dist_edge, gap_size, moe, c_reg, stem_vol, stem_density, 
	crown_density, snow_load, ci, ci_value, 
	n_drag, c_drag, drag_upper_limit, equivalent_mean_ht, 
	stand_cr_depth, stand_cr_width, fgr_constants, aerodynamic_ht);
	
	// Extract the variables from these results
	Rcpp::NumericVector uh_b = uh_b_results["uh_b"];
	Rcpp::NumericVector uh_o = uh_o_results["uh_o"];	
	
	// Elevate the value
	int n = uh_b.size();
	Rcpp::NumericVector zpd_b(n);
	Rcpp::NumericVector zpd_o(n);
	Rcpp::NumericVector z0_b(n);
	Rcpp::NumericVector z0_o(n);
	Rcpp::NumericVector lambdacapital_b(n);
	Rcpp::NumericVector lambdacapital_o(n);
	Rcpp::NumericVector gammasolved_b(n);
	Rcpp::NumericVector gammasolved_o(n);
	Rcpp::NumericVector canopybreadth_b(n);
	Rcpp::NumericVector canopybreadth_o(n);
	Rcpp::NumericVector drag_b(n);
	Rcpp::NumericVector drag_o(n);
	
	Rcpp::NumericVector u_elev_b(n);
	Rcpp::NumericVector u_elev_o(n);
	Rcpp::NumericVector u_damage(n);
	Rcpp::NumericVector u_elev_damage(n);
	Rcpp::StringVector mode_of_damage(n);
	
	for (int i = 0; i < n; ++i){
		zpd_b[i] = zpd_fun_cpp(stand_cr_width[i], stand_cr_depth[i], spacing_current[i], 
		uh_b[i], n_drag[i], c_drag[i], drag_upper_limit[i], aerodynamic_ht[i], 
		fgr_constants);		
		z0_b[i] = z0_fun_cpp(stand_cr_width[i], stand_cr_depth[i], spacing_current[i], uh_b[i], 
		n_drag[i], c_drag[i], drag_upper_limit[i], aerodynamic_ht[i], fgr_constants);
		lambdacapital_b[i] = lambdacapital_fun_cpp(stand_cr_width[i], stand_cr_depth[i], 
		spacing_current[i], uh_b[i], n_drag[i], c_drag[i], drag_upper_limit[i]);
		gammasolved_b[i] = gammasolved_fun_cpp(stand_cr_width[i], stand_cr_depth[i], 
		spacing_current[i], uh_b[i], n_drag[i], c_drag[i], drag_upper_limit[i], 
		fgr_constants);
		canopybreadth_b[i] = canopy_breadth_fun_cpp(stand_cr_width[i], uh_b[i], 
		n_drag[i], c_drag[i], drag_upper_limit[i]);
		drag_b[i] = drag_fun_cpp(uh_b[i], n_drag[i], c_drag[i], drag_upper_limit[i]);
		
		zpd_o[i] = zpd_fun_cpp(stand_cr_width[i], stand_cr_depth[i], spacing_current[i], 
		uh_o[i], n_drag[i], c_drag[i], drag_upper_limit[i], aerodynamic_ht[i], 
		fgr_constants);		
		z0_o[i] = z0_fun_cpp(stand_cr_width[i], stand_cr_depth[i], spacing_current[i], uh_o[i], 
		n_drag[i], c_drag[i], drag_upper_limit[i], aerodynamic_ht[i], fgr_constants);
		lambdacapital_o[i] = lambdacapital_fun_cpp(stand_cr_width[i], stand_cr_depth[i], 
		spacing_current[i], uh_o[i], n_drag[i], c_drag[i], drag_upper_limit[i]);
		gammasolved_o[i] = gammasolved_fun_cpp(stand_cr_width[i], stand_cr_depth[i], 
		spacing_current[i], uh_o[i], n_drag[i], c_drag[i], drag_upper_limit[i], 
		fgr_constants);
		canopybreadth_o[i] = canopy_breadth_fun_cpp(stand_cr_width[i], uh_o[i], 
		n_drag[i], c_drag[i], drag_upper_limit[i]);		
		drag_o[i] = drag_fun_cpp(uh_o[i], n_drag[i], c_drag[i], drag_upper_limit[i]);

		u_elev_b[i] = elevate_cpp(uh_b[i], z0_b[i], zpd_b[i], elev_ht[i], ht_above_d[i]);
		u_elev_o[i] = elevate_cpp(uh_o[i], z0_o[i], zpd_o[i], elev_ht[i], ht_above_d[i]);
		
		if(uh_b[i] < uh_o[i]){u_damage[i] = uh_b[i];} else {u_damage[i] = uh_o[i];} 
		if(u_elev_b[i] < u_elev_o[i]){u_elev_damage[i] = u_elev_b[i];} else {u_elev_damage[i] = u_elev_o[i];} 
		
		if(u_elev_b[i] < u_elev_o[i]){mode_of_damage[i] = "Breakage";} else {mode_of_damage[i] = "Overturning";}		
		
	}
	
	// Calculate probability of breakage, overturning, and damage - to do
	
	// Add stems per hectare
	others["sph_before"] = round(10000.0/(pow(spacing_before, 2)), 0);
	others["sph_current"] = round(10000.0/(pow(spacing_current, 2)), 0);
	
	// If breakage at the base of the stem
	Rcpp::DataFrame output_data_Bmin;
	if(breakage_basecanopy == "yes"){
		
		Rcpp::DataFrame uh_b_min_results = uh_breakage_heightz_cpp(uh_b, tree_ht, cr_depth, dbh);
		Rcpp::NumericVector uh_b_min = uh_b_min_results["uh_b_min"];
		Rcpp::NumericVector hmin_cws = uh_b_min_results["hmin_cws"];
		Rcpp::NumericVector diam_htz = uh_b_min_results["diam_htz"];	

		Rcpp::NumericVector zpd_b_min(n);
		Rcpp::NumericVector z0_b_min(n);
		Rcpp::NumericVector u_elev_b_min(n);
		Rcpp::NumericVector lambdacapital_b_min(n);
		Rcpp::NumericVector gammasolved_b_min(n);
		Rcpp::NumericVector drag_b_min(n);
		
		for (int i = 0; i < n; ++i){
			zpd_b_min[i] = zpd_fun_cpp(stand_cr_width[i], stand_cr_depth[i], spacing_current[i], uh_b_min[i], 
			n_drag[i], c_drag[i], drag_upper_limit[i], aerodynamic_ht[i], fgr_constants);
			z0_b_min[i] = z0_fun_cpp(stand_cr_width[i], stand_cr_depth[i], spacing_current[i], uh_b_min[i], 
			n_drag[i], c_drag[i], drag_upper_limit[i], aerodynamic_ht[i], fgr_constants);

			u_elev_b_min[i] = elevate_cpp(uh_b_min[i], z0_b_min[i], zpd_b_min[i], elev_ht[i], ht_above_d[i]);
			lambdacapital_b_min[i] = lambdacapital_fun_cpp(stand_cr_width[i], stand_cr_depth[i], 
			spacing_current[i], uh_b_min[i], n_drag[i], c_drag[i], drag_upper_limit[i]);
			drag_b_min[i] = drag_fun_cpp(uh_b_min[i], n_drag[i], c_drag[i], drag_upper_limit[i]);
			gammasolved_b_min = gammasolved_fun_cpp(stand_cr_width[i], stand_cr_depth[i],spacing_current[i], uh_b_min[i], 
			n_drag[i], c_drag[i], drag_upper_limit[i], fgr_constants);
			if(u_elev_damage[i] > u_elev_b_min[i]){
				u_elev_damage[i] = u_elev_b_min[i];
				u_damage[i] = uh_b_min[i];
				mode_of_damage[i] = "Breakage_base_canopy";
			}
		}
		
		output_data_Bmin = Rcpp::DataFrame::create(Rcpp::_["u_elev_b_min"] = Rcpp::clone(u_elev_b_min), Rcpp::_["uh_b_min"] = uh_b_min_results["uh_b_min"],
		Rcpp::_["zpd_b_min"] = Rcpp::clone(zpd_b_min), Rcpp::_["z0_b_min"] = Rcpp::clone(z0_b_min), Rcpp::_["drag_b_min"] = Rcpp::clone(drag_b_min),
		Rcpp::_["gammasolved_b_min"] = Rcpp::clone(gammasolved_b_min), Rcpp::_["lambdacapital_b_min"] = Rcpp::clone(lambdacapital_b_min),
		Rcpp::_["u_elev_damage"] = Rcpp::clone(u_elev_damage), Rcpp::_["mode_of_damage"] = Rcpp::clone(mode_of_damage));		
	}
	
	
	// Choose whether to return full dataframe - including the input data - or just the results dataframe
	Rcpp::List output_data;	
	
	if (full_output_.isNotNull()) {
        Rcpp::NumericVector full_output(full_output_);        // casting to underlying type NumericVector
		output_data["required"] = Rcpp::clone(required);
		output_data["desirable"] = Rcpp::clone(desirable);
		output_data["advanced"] = Rcpp::clone(advanced);
		output_data["others"] = Rcpp::clone(others);
		output_data["results_BO"] =  Rcpp::DataFrame::create(Rcpp::_["u_elev_b"] = Rcpp::clone(u_elev_b), Rcpp::_["uh_b"] = uh_b_results["uh_b"],
		//Rcpp::_["prob_b"] = Rcpp::clone(prob_b), 
		Rcpp::_["zpd_b"] = Rcpp::clone(zpd_b), Rcpp::_["z0_b"] = Rcpp::clone(z0_b), Rcpp::_["drag_b"] = Rcpp::clone(drag_b),
		Rcpp::_["gammasolved_b"] = Rcpp::clone(gammasolved_b), Rcpp::_["lambdacapital_b"] = Rcpp::clone(lambdacapital_b), 
		Rcpp::_["canopybreadth_b"] = Rcpp::clone(canopybreadth_b), Rcpp::_["breaking_moment"] = uh_b_results["breaking_moment"],
		Rcpp::_["u_elev_o"] = Rcpp::clone(u_elev_o), Rcpp::_["uh_o"] = uh_o_results["uh_o"],
		//Rcpp::_["prob_o"] = Rcpp::clone(prob_o), 
		Rcpp::_["zpd_o"] = Rcpp::clone(zpd_o), Rcpp::_["z0_o"] = Rcpp::clone(z0_o), Rcpp::_["drag_o"] = Rcpp::clone(drag_o),
		Rcpp::_["gammasolved_o"] = Rcpp::clone(gammasolved_o), Rcpp::_["lambdacapital_o"] = Rcpp::clone(lambdacapital_o),
		Rcpp::_["canopybreadth_o"] = Rcpp::clone(canopybreadth_o), Rcpp::_["overturning_moment"] = uh_o_results["overturning_moment"]);
		
		if(breakage_basecanopy == "yes"){
			output_data["results_Bmin"] = Rcpp::clone(output_data_Bmin);
			output_data["results_damage"] =  Rcpp::DataFrame::create(Rcpp::_["mode_of_damage"] = output_data_Bmin["mode_of_damage"], 
			Rcpp::_["u_elev_damage"] = output_data_Bmin["u_elev_damage"], Rcpp::_["u_damage"] = output_data_Bmin["u_damage"], 
			Rcpp::_["tmc"] = uh_b_results["tmc"], Rcpp::_["tmr_full"] = uh_b_results["tmr_full"],
			Rcpp::_["dlf_calc"] = uh_b_results["dlf_calc"], Rcpp::_["dlf_used"] = uh_b_results["dlf_used"]);	
		} else {
			output_data["results_damage"] =  Rcpp::DataFrame::create(
			Rcpp::_["mode_of_damage"] = Rcpp::clone(mode_of_damage), Rcpp::_["u_elev_damage"] = Rcpp::clone(u_elev_damage),
			Rcpp::_["u_damage"] = Rcpp::clone(u_damage), 
			//Rcpp::_["prob_damage"] = Rcpp::clone(prob_damage),
			Rcpp::_["tmc"] = uh_b_results["tmc"], Rcpp::_["tmr_full"] = uh_b_results["tmr_full"],
			Rcpp::_["dlf_calc"] = uh_b_results["dlf_calc"], Rcpp::_["dlf_used"] = uh_b_results["dlf_used"]);
		}		
	} else {
		if(breakage_basecanopy == "yes"){
			output_data["results"] = Rcpp::DataFrame::create(Rcpp::_["stand_id"] = Rcpp::clone(stand_id), Rcpp::_["tree_id"] = Rcpp::clone(tree_id), 
			Rcpp::_["species"] = Rcpp::clone(species),
			Rcpp::_["u_elev_damage"] = output_data_Bmin["u_elev_damage"], Rcpp::_["mode_of_damage"] = output_data_Bmin["mode_of_damage"],
			Rcpp::_["u_elev_b"] = Rcpp::clone(u_elev_b), Rcpp::_["u_elev_b_min"] = output_data_Bmin["u_elev_b_min"], Rcpp::_["u_elev_o"] = Rcpp::clone(u_elev_o),
			Rcpp::_["max_stem_weight_warning"] = others["max_stem_weight_warning"]);			
		} else {
			output_data["results"] = Rcpp::DataFrame::create(Rcpp::_["stand_id"] = Rcpp::clone(stand_id), Rcpp::_["tree_id"] = Rcpp::clone(tree_id), 
			Rcpp::_["species"] = Rcpp::clone(species),
			Rcpp::_["u_elev_damage"] = Rcpp::clone(u_elev_damage), Rcpp::_["mode_of_damage"] = Rcpp::clone(mode_of_damage),
			Rcpp::_["u_elev_b"] = Rcpp::clone(u_elev_b), Rcpp::_["u_elev_o"] = Rcpp::clone(u_elev_o),
			//Rcpp::_["prob_b"] = Rcpp::clone(prob_b), Rcpp::_["prob_o"] = Rcpp::clone(prob_o), Rcpp::_["prob_damage"] = Rcpp::clone(prob_damage),
			Rcpp::_["max_stem_weight_warning"] = others["max_stem_weight_warning"]);
		}
	}
	
	return output_data;	
}