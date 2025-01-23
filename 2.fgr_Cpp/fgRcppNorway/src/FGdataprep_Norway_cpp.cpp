
#include <Rcpp.h>

// Note we focus on the three main tree - NS (Norway spruce), SP (Scots pine) and BI (birch) - so some of the functions are only defined for those (namely stem volume)
// [[Rcpp::export]]                                                                                                                                           
bool contains_cpp(std::string s, Rcpp::DataFrame L) {                                                                                                                  
    Rcpp::CharacterVector nv = L.names();                                                                                                                     
    for (int i=0; i<nv.size(); i++) {                                                                                                                         
        if (std::string(nv[i]) == s) {                                                                                                                        
            return true;                                                                                                                                      
        }                                                                                                                                                     
    }                                                                                                                                                         
    return false;                                                                                                                                             
}

// [[Rcpp::export]]
Rcpp::NumericVector extract_create_numcolumns_cpp(std::string s, Rcpp::DataFrame DF, int n){
	
	Rcpp::NumericVector out(n);
	if(contains_cpp(s, DF)){
		out = DF[s];
	} else {
		for (int i = 0; i < n; ++i){
			out[i] = NA_REAL;
		}
	}
	
	return out;
}
// [[Rcpp::export]]
Rcpp::IntegerVector extract_create_intcolumns_cpp(std::string s, Rcpp::DataFrame DF, int n){

	Rcpp::IntegerVector out(n);
	if(contains_cpp(s, DF)){
		out = DF[s];
	} else {
		for (int i = 0; i < n; ++i){
			out[i] = NA_INTEGER;
		}
	}
	
	return out;
}


// [[Rcpp::export]]
Rcpp::StringVector extract_create_strcolumns_cpp(std::string s, Rcpp::DataFrame DF, int n){

	Rcpp::StringVector out(n);
	if(contains_cpp(s, DF)){
		out = DF[s];
	} else {
		for (int i = 0; i < n; ++i){
			out[i] = NA_STRING;
		}
	}
	
	return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector sp_params_Norway_fun_cpp(Rcpp::DataFrame species_parameters, Rcpp::String species, 
Rcpp::String location, Rcpp::String dbh_class, Rcpp::String leafstatus){
	
	// Subset the DataFrame species_parameter to the row corresponding to the species (location and dbh_class)  
	
	std::string rnames_id;	
	std::string species_u = std::string(species);
	std::string location_u = std::string(location);
	std::string dbh_class_u = std::string(dbh_class);
	std::string leafstatus_u = std::string(leafstatus);
	
	
	if(location_u.empty() && dbh_class_u.empty() && leafstatus_u.empty()){rnames_id = species_u;}
	else{		
		// create the identifier from the row.names
		char location_u1 = location_u.at(0);
		char dbh_class_u1 = dbh_class_u.at(0);
		if(species_u == "BI"){
			rnames_id = species_u + "." + location_u1 + "." + dbh_class_u1 + "." + leafstatus_u;
		} else {rnames_id = species_u + "." + location_u1 + "." + dbh_class_u1;}
	}
	// First get the size in number of columns of the dataframe
	int nCols = species_parameters.size();
	// Then create a NumericVector as output with nCols columns
	Rcpp::NumericVector sp_params(nCols);
	// Find which row index is the ID (species, location, dbh_class) of interest
	Rcpp::CharacterVector sp_rnames = species_parameters.attr("row.names");  
	// Set up storage for found IDs
	int id_loc = NA_INTEGER;  
	int sp_rnames_size = sp_rnames.size();
	// Loop through input and check if the input matches the target	
	for(int i = 0; i < sp_rnames_size; i++){
		if(sp_rnames[i] == rnames_id){
			id_loc = i;
		}
	}
	// Then fill the vector with the column values for that row
	// Except for the first columns (columns 0 to 4) which are strings
	for (int j = 5; j < nCols; j++) {
		Rcpp::NumericVector column = species_parameters[j] ;
		sp_params[j] = column[id_loc] ;
	}
	// add names
	sp_params.attr("names") = species_parameters.attr("names");
	return sp_params;  
}


// [[Rcpp::export]]
double top_ht_to_mean_ht_cpp(double param0_height, double param1_height, double top_ht){
	double mean_ht = param0_height + param1_height *top_ht;
	return mean_ht;
}

// [[Rcpp::export]]
double mean_ht_to_top_ht_cpp(double param0_height, double param1_height, double mean_ht){
	double top_ht = (mean_ht - param0_height) / param1_height;
	return top_ht;
}

// [[Rcpp::export]]
double eq_mean_ht_fun_cpp(double top_ht) {
	double equivalent_mean_ht = 1.03 * top_ht - 1.8;
	return equivalent_mean_ht;
}

// [[Rcpp::export]]
double canopy_width_fun_cpp(double param0_cr_width, double param1_cr_width, double dbh){
	double crown_width = param0_cr_width + param1_cr_width * dbh;
	return crown_width;
}

// [[Rcpp::export]]
double canopy_depth_Norway_fun_cpp(std::string species, Rcpp::NumericVector species_parameters, double dbh, double ht, std::string focus = "tree"){
	
	//Note that dbh and ht need to be in mm and cm respectively for these formulas
	
	double crown_depth;
	
	
	if(focus == "stand"){
		double param0_cr_depth = species_parameters["p_sacd0"];
		double param1_cr_depth = species_parameters["p_sacd1"];	
		double param2_cr_depth = species_parameters["p_sacd2"];
		double param3_cr_depth = species_parameters["p_sacd3"];
		double param4_cr_depth = species_parameters["p_sacd4"];
		crown_depth = param0_cr_depth + param1_cr_depth * dbh + param2_cr_depth * ht + param3_cr_depth * pow(dbh, 2) + param4_cr_depth * pow(ht, 2);
		crown_depth = pow(crown_depth, 2)/10;
	} else {
		double param0_cr_depth = species_parameters["p_cd0"];
		double param1_cr_depth = species_parameters["p_cd1"];	
		double param2_cr_depth = species_parameters["p_cd2"];
		double param3_cr_depth = species_parameters["p_cd3"];
		double param4_cr_depth = species_parameters["p_cd4"];
		if(species == "SP"){
			crown_depth = ht - exp(param0_cr_depth + param1_cr_depth * dbh + param2_cr_depth * ht + param4_cr_depth * pow(ht, 2));
		} else {
			crown_depth = ht - pow(param0_cr_depth + param1_cr_depth * dbh + param2_cr_depth * ht + param3_cr_depth * pow(dbh, 2) + param4_cr_depth * pow(ht, 2), 2);
		}			
		crown_depth = crown_depth/10;
	}	
	return crown_depth;
}


// DEFAULT VOLUME FUNCTIONS
/*
Calculates the stem volume of a tree.
title Stem Volume Function.
param species Tree species.
param dbh Diameter of the stem at breast height, i.e. 1.3m above the ground (cm). Depending on the method used ('roughness' or 'TMC') this can be either the arithmetic average of the dbh of all the trees in the stand, or the dbh of an individual tree.
param ht Tree height. Depending on the method used ('roughness' or 'TMC'), this can be either the mean tree in the stand, or each individual tree (m).
return The volume of the stem (m3). Applies the method in Fonweban et al. (2012) to calculate stem volume.
*/


// [[Rcpp::export]]
double stem_volume_fonweban_cpp(double dbh, double ht, Rcpp::NumericVector species_parameters){
	/* From Fonweban et al. 2012. Eq 4 */	
	double param0_tvf = species_parameters["param0_vol"];
	double param1_tvf = species_parameters["param1_vol"];
	double param2_tvf = species_parameters["param2_vol"];
	double stem_vol = param0_tvf * (pow(dbh, param1_tvf) * pow(ht, param2_tvf));	
	return stem_vol;
}

// [[Rcpp::export]]
/* applies the method in Honer 1967 to calculate stem volume.*/
double stem_volume_quebec_cpp(double dbh, double ht, Rcpp::NumericVector species_parameters) {
	//From JC Ruel
	double param0_tvq = species_parameters["param0_vol"];
	double param1_tvq = species_parameters["param1_vol"];
	double param2_tvq = species_parameters["param2_vol"];
	double stem_vol = 0.004389 * pow(dbh, 2) * pow((1- 0.04365 * param0_tvq), 2) / (param1_tvq + 0.3048 * param2_tvq/ht);
	return stem_vol;
}

// [[Rcpp::export]]
/* applies the method in Laasasenaho 1982 to calculate stem volume.*/
double stem_volume_laasasenaho_cpp(double dbh, double ht, Rcpp::NumericVector species_parameters) {
	double param0_tvl = species_parameters["param0_vol"];
	double param1_tvl = species_parameters["param1_vol"];
	double param2_tvl = species_parameters["param2_vol"];
	double param3_tvl = species_parameters["param3_vol"];
	double param4_tvl = species_parameters["param4_vol"];
	double stem_vol = (pow(dbh, param1_tvl) * pow(ht, param2_tvl) * pow((ht - 1.3), param3_tvl) * exp(param0_tvl + param4_tvl * dbh)) / 1000.0;
	return stem_vol;
}

// [[Rcpp::export]]
/* calculates stem volume for Japanese Larch grown in Japan.*/
double stem_volume_japanese_larch_japan_cpp(double dbh, double ht) {
	/* parameters for this volume function need to be kept here rather than in the species_parameters df because
	they depend on the value of dbh, and as such they need to be evaluated. */
	double param1_tvjl;
	double param2_tvjl;
	double param3_tvjl;
	double param4_tvjl;
	if (dbh < 11.0) {
		param1_tvjl = 0.77430;
		param2_tvjl = 1.87385;
		param3_tvjl = 0.94852;
		param4_tvjl = 1.0026;
	} else if (dbh < 21.0) {
		param1_tvjl = 0.58495;
		param2_tvjl = 1.96416;
		param3_tvjl = 1.04523;
		param4_tvjl = 1.0027;
	} else if (dbh < 31.0) {
		param1_tvjl = 0.67205;
		param2_tvjl = 1.84173;
		param3_tvjl = 1.11080;
		param4_tvjl = 1.0012;
	} else {
		param1_tvjl = 0.79071;
		param2_tvjl = 1.74034;
		param3_tvjl = 1.13316;
		param4_tvjl = 1.0016;
	}  
	double stem_vol = param4_tvjl * pow(10,(-5 + param1_tvjl + param2_tvjl * log10(dbh) + param3_tvjl * log10(ht)));
	return stem_vol;
}


// [[Rcpp::export]]
double stem_vol_andretreslag_fun_cpp(std::string species, double dbh, double ht, Rcpp::NumericVector species_parameters) {
	
	double stem_vol_f = 0;	
	std::string species_u = species.substr(0);
	
	if(species == "SS" || species == "CP" || species == "LP" || species == "EL" || species == "HL" || species == "DF" || species == "JL" ||
	species == "NF" || species == "GF" || species == "WH" || species == "BE" || species == "OK" || species == "MP" || species == "RP" || species == "EG" || species_u == "U") {
		// MP=U1; RP=U2; EG=U3
		stem_vol_f = stem_volume_fonweban_cpp(dbh, ht, species_parameters);
	}
	if(species == "JLJ") {
		stem_vol_f = stem_volume_japanese_larch_japan_cpp(dbh, ht);
	}
	if(species == "WS" || species == "BS" || species == "BF" || species == "JP") {
		// WS=U4; BS=U5; BF=U6; JP=U7
		stem_vol_f = stem_volume_quebec_cpp(dbh, ht, species_parameters);
	}
	return stem_vol_f;
}
// [[Rcpp::export]]
double stem_vol_NSNorway_fun_cpp(std::string location, std::string dbh_class, double dbh, double ht, Rcpp::NumericVector species_parameters){
	double stem_vol_f = 0;	
	
	double param0_tvl = species_parameters["param0_vol"];
	double param1_tvl = species_parameters["param1_vol"];
	double param2_tvl = species_parameters["param2_vol"];
	double param3_tvl = species_parameters["param3_vol"];
	double param4_tvl = species_parameters["param4_vol"];
	double param5_tvl = species_parameters["param5_vol"];
	
	if (location == "annen"){
		if(dbh_class == "medium"){
			stem_vol_f = (param0_tvl + param1_tvl * dbh *  pow(ht, 2) + param2_tvl * pow(ht, 2) + param3_tvl * dbh * ht + param4_tvl * ht + param5_tvl * dbh)/1000;
		} else {
			stem_vol_f = (param0_tvl + param1_tvl * pow(dbh, 2) * ht + param2_tvl * dbh * pow(ht, 2) + param3_tvl * pow(ht, 2) + param4_tvl * dbh * ht)/1000;
		}
	} else {
		stem_vol_f = (param0_tvl * pow(ht, param1_tvl) * pow(dbh, param2_tvl) * pow(ht - 1.3, param3_tvl) * pow(dbh + 40, param4_tvl))/1000;
	}
	
	return stem_vol_f;
}
// [[Rcpp::export]]
double stem_vol_SPNorway_fun_cpp(std::string location, std::string dbh_class, double dbh, double ht, Rcpp::NumericVector species_parameters){
	double stem_vol_f = 0;	
	
	double param0_tvl = species_parameters["param0_vol"];
	double param1_tvl = species_parameters["param1_vol"];
	double param2_tvl = species_parameters["param2_vol"];
	double param3_tvl = species_parameters["param3_vol"];
	double param4_tvl = species_parameters["param4_vol"];
	
	if (location == "annen"){
		if(dbh_class == "liten"){
			stem_vol_f = (param0_tvl + param1_tvl * pow(dbh, 2) + param2_tvl * pow(dbh, 2) * ht + param3_tvl * dbh * pow(ht, 2))/1000;
		} else {
			stem_vol_f = (param0_tvl + param1_tvl * pow(dbh, 2) + param2_tvl * pow(dbh, 2) * ht + param3_tvl * pow(dbh, 2) * (3.17935 + (1.02890 * dbh) - (0.27023 * (dbh / ht))))/1000;
		}
	} else {
		stem_vol_f = param0_tvl * pow(ht, param1_tvl) * pow(dbh, param2_tvl) * pow(ht - 1.3, param3_tvl) * pow(dbh + 100, param4_tvl)/1000;
	}
	
	return stem_vol_f;
}
// [[Rcpp::export]]
double stem_vol_BINorway_fun_cpp(std::string location, std::string dbh_class, double dbh, double ht, Rcpp::NumericVector species_parameters){
	double stem_vol_f = 0;	
	
	double param0_tvl = species_parameters["param0_vol"];
	double param1_tvl = species_parameters["param1_vol"];
	double param2_tvl = species_parameters["param2_vol"];
	double param3_tvl = species_parameters["param3_vol"];
	double param4_tvl = species_parameters["param4_vol"];
	double param5_tvl = species_parameters["param5_vol"];
	
	stem_vol_f = (param0_tvl * (param1_tvl + (param2_tvl * pow(dbh, 2)) + (param3_tvl * pow(dbh, 2) * ht) + (param4_tvl * dbh * pow(ht, 2)) + (param5_tvl * pow(ht, 2))))/1000;
	
	return stem_vol_f;
}
// [[Rcpp::export]]
double stem_vol_Norway_fun_cpp(std::string species, std::string location, std::string dbh_class, double dbh, double ht, Rcpp::NumericVector species_parameters) {
	
	double stem_vol_f = 0;	
	
	if(species == "NS"){
		stem_vol_f = stem_vol_NSNorway_fun_cpp(location, dbh_class, dbh, ht, species_parameters);
	} else if(species == "SP"){
		stem_vol_f = stem_vol_SPNorway_fun_cpp(location, dbh_class, dbh, ht, species_parameters);
	} else if(species == "BI"){
		stem_vol_f = stem_vol_BINorway_fun_cpp(location, dbh_class, dbh, ht, species_parameters);
	}
	
	return stem_vol_f;
}

// MAIN FUNCTION FOR DATA PREPARATION - ROUGHNESS METHOD
// [[Rcpp::export]]
Rcpp::List fg_rou_dataprep_Norway_cpp(Rcpp::DataFrame inputdata, Rcpp::DataFrame fgr_constants, Rcpp::DataFrame species_parameters, std::string season = "winter",
std::string country = "fgr"){
	// 1. Check essentials (species, mean_ht&top_ht, mean_dbh, spacing):
	//first step check that all the essential variables are present in the inputdata DataFrame
	//essential variables are the following: stand_id, species, mean_dbh, spacing
	try
	{
		Rcpp::StringVector stand_id = inputdata["stand_id"];
		Rcpp::StringVector species = inputdata["species"];
		Rcpp::NumericVector mean_dbh = inputdata["mean_dbh"]; 
		Rcpp::NumericVector spacing = inputdata["spacing"]; 
	}
	catch (...)
	{
		std::cout << "Missing some key variables! Check that stand_id, species, mean_dbh and spacing are in your dataframe"<<std::endl;
		return EXIT_FAILURE;
	}	
	
	// If all variables are present, load them
	Rcpp::StringVector stand_id = inputdata["stand_id"];
	Rcpp::StringVector species = inputdata["species"];
	Rcpp::NumericVector mean_dbh = inputdata["mean_dbh"]; 
	Rcpp::NumericVector spacing = inputdata["spacing"];
	spacing = round(spacing, 1);
	// Check if any of the variables contain NA	
	Rcpp::LogicalVector stand_id_na = is_na(stand_id);
	Rcpp::LogicalVector species_na = is_na(species);
	Rcpp::LogicalVector mean_dbh_na = is_na(mean_dbh);
	Rcpp::LogicalVector spacing_na = is_na(spacing);
	if(is_true(any(stand_id_na)) || is_true(any(species_na)) || 
	is_true(any(mean_dbh_na)) || is_true(any(spacing_na))){
		std::cout<<"One or more of the key variables contains NA! Check your input data"<<std::endl;
		return EXIT_FAILURE;
	}
	// Get the length of the dataframe
	int n = stand_id.size();
	
	// do the same for the height variables and check that they don't both contain NAs (i.e. they are both absent from the initial dataframe or they are present but missing
	// data)
	Rcpp::NumericVector top_ht = extract_create_numcolumns_cpp("top_ht", inputdata, n);
	Rcpp::NumericVector mean_ht = extract_create_numcolumns_cpp("mean_ht", inputdata, n);
	for (int i = 0; i < n; ++i){
		bool top_ht_i_na = Rcpp::NumericVector::is_na(top_ht[i]);
		bool mean_ht_i_na = Rcpp::NumericVector::is_na(mean_ht[i]);
		if(top_ht_i_na && mean_ht_i_na){
			std::cout<<"Both top_ht and mean_ht are missing in some/all stands, please check input data"<<std::endl;
			return EXIT_FAILURE;
		}
	}
	
	// 2. Extract the desirable variables if they are present
	// If a variable is not present, create it with default calculations, if present and contains NAs, we will replace the NAs with the default values
	// fylke - the fylke ( = region) is important for calculating appropriate stem volume as the equations may differ. If not present it will be defaulted to alle
	Rcpp::IntegerVector fylke;
	fylke = extract_create_intcolumns_cpp("fylke", inputdata, n);
	// create the location vector for the stem volume and crown dimension functions (only for NS, SP and BI)
	Rcpp::StringVector location(n);
	for (int i = 0; i < n; ++i){
		if(species[i] == "NS" || species[i] == "SP"){
			if(fylke[i] != 11 && fylke[i] != 12 && fylke[i] != 14 && fylke[i] != 15 ){location[i] = "annen";}
			else {location[i] = "vest";}
		} else if(species[i] == "BI"){location[i] = "alle";}
		else {location[i] = NA_STRING;}
	}
	// dbh_class - required for appropriate stem volume functions for NS, SP and BI
	Rcpp::StringVector dbh_class(n);
	for (int i = 0; i < n; ++i){
		if((location[i] == "vest") || (species[i] == "BI")){dbh_class[i] = "alle";}
		else if(species[i] == "SP"){
			if(mean_dbh[i] < 11.1){dbh_class[i] = "liten";}
			else {dbh_class[i] = "stor";}
		}
		else if(species[i] == "NS"){
			if(mean_dbh[i] < 10.1){dbh_class[i] = "liten";}
			else if(mean_dbh[i] >= 10.1 && mean_dbh[i] < 12.9){dbh_class[i] = "medium";}
			else{dbh_class[i] = "stor";}
		}
		else {dbh_class[i] = NA_STRING;}
	}	
	// leafstatus - for birch in Norway only
	Rcpp::StringVector leafstatus(n);
	for (int i = 0; i < n; ++i){
		if(species[i] == "BI"){
			if(season == "summer"){leafstatus[i] = "lo";}
			else {leafstatus[i] = "ll";}
		}
		else {leafstatus[i] = NA_STRING;}
	}		
	// aerodynamic_ht
	Rcpp::NumericVector	aerodynamic_ht;
	aerodynamic_ht = extract_create_numcolumns_cpp("aerodynamic_ht", inputdata, n);
	// mean_cr_width
	Rcpp::NumericVector	mean_cr_width;
	mean_cr_width = extract_create_numcolumns_cpp("mean_cr_width", inputdata, n);
	// mean_cr_depth
	Rcpp::NumericVector	mean_cr_depth;
	mean_cr_depth = extract_create_numcolumns_cpp("mean_cr_depth", inputdata, n);
	
	
	// soil_group
	Rcpp::IntegerVector	soil_group;
	soil_group = extract_create_intcolumns_cpp("soil_group", inputdata, n);
	// rooting
	Rcpp::IntegerVector	rooting;
	rooting = extract_create_intcolumns_cpp("rooting", inputdata, n);
	// Gap_size
	Rcpp::NumericVector	gap_size;
	gap_size = extract_create_numcolumns_cpp("gap_size", inputdata, n);
	// dist_edge
	Rcpp::NumericVector dist_edge;
	dist_edge = extract_create_numcolumns_cpp("dist_edge", inputdata, n);
	
	// 3. Advanced variables
	// ht_above_d
	Rcpp::NumericVector	ht_above_d;
	ht_above_d = extract_create_numcolumns_cpp("ht_above_d", inputdata, n);
	// Elev_ht
	Rcpp::NumericVector	elev_ht;
	elev_ht = extract_create_numcolumns_cpp("elev_ht", inputdata, n);
	// MOE
	Rcpp::NumericVector	moe;
	moe = extract_create_numcolumns_cpp("moe", inputdata, n);
	// MOR
	Rcpp::NumericVector	mor;
	mor = extract_create_numcolumns_cpp("mor", inputdata, n);	
	// stem_density
	Rcpp::NumericVector	stem_density;
	stem_density = extract_create_numcolumns_cpp("stem_density", inputdata, n);
	// crown_density
	Rcpp::NumericVector	crown_density;
	crown_density = extract_create_numcolumns_cpp("crown_density", inputdata, n);
	// fknot
	Rcpp::NumericVector	fknot;
	fknot = extract_create_numcolumns_cpp("fknot", inputdata, n);
	// c_drag
	Rcpp::NumericVector	c_drag;
	c_drag = extract_create_numcolumns_cpp("c_drag", inputdata, n);
	// n_drag
	Rcpp::NumericVector	n_drag;
	n_drag = extract_create_numcolumns_cpp("n_drag", inputdata, n);
	// drag_upper_limit
	Rcpp::NumericVector	drag_upper_limit;
	drag_upper_limit = extract_create_numcolumns_cpp("drag_upper_limit", inputdata, n);
	// c_reg
	Rcpp::NumericVector	c_reg;
	c_reg = extract_create_numcolumns_cpp("c_reg", inputdata, n);	
	// crown_offset
	Rcpp::NumericVector	crown_off;
	crown_off = extract_create_numcolumns_cpp("crown_off", inputdata, n);	
	// stem_vol
	Rcpp::NumericVector	stem_vol;
	stem_vol = extract_create_numcolumns_cpp("stem_vol", inputdata, n);
	// crown_vol
	Rcpp::NumericVector	crown_vol;
	crown_vol = extract_create_numcolumns_cpp("crown_vol", inputdata, n);
	// stem_weight
	Rcpp::NumericVector stem_weight(n);
	stem_weight = extract_create_numcolumns_cpp("stem_weight", inputdata, n);
	// max_stem_weight_warning
	Rcpp::StringVector max_stem_weight_warning(n);
	for (int i = 0; i < n; ++i){max_stem_weight_warning[i] = NA_STRING;}
	//Initiate snow_load
	Rcpp::NumericVector snow_load(n);
	
	//4. Wind variables - to do 
	
	//5. Replace the NA values with default calculations
	
	double tree_heights_inside_forest = fgr_constants["tree_heights_inside_forest"];
	for (int i = 0; i < n; ++i){
		std::string s_species = Rcpp::as<std::string>(species[i]);
		std::string s_location = Rcpp::as<std::string>(location[i]);
		std::string s_dbhclass = Rcpp::as<std::string>(dbh_class[i]);
		std::string s_leaf = Rcpp::as<std::string>(leafstatus[i]);
		// select the appropriate species_parameter row
		// if the species is not present or the species 2-letter identifier not correct, display a message
		try{Rcpp::NumericVector sp_params = sp_params_Norway_fun_cpp(species_parameters, s_species, s_location, s_dbhclass, s_leaf);} catch (...) {
			std::cout << "This species is not present in the species_parameters dataframe or its 2-letter identifier is not correct." << std::endl;
			return EXIT_FAILURE;
		}
		// if the species is other than NS, BI or SP, use the default fgr species parameters, otherwise use the country-specific ones
		Rcpp::NumericVector sp_params = sp_params_Norway_fun_cpp(species_parameters, s_species, s_location, s_dbhclass, s_leaf);
		
		// Fill in the variables
		
		// Desirable variables
		if (Rcpp::NumericVector::is_na(top_ht[i])) {top_ht[i] = mean_ht_to_top_ht_cpp(sp_params["param0_height"], sp_params["param1_height"], mean_ht[i]);}
		if (Rcpp::NumericVector::is_na(mean_ht[i])) {mean_ht[i] = top_ht_to_mean_ht_cpp(sp_params["param0_height"], sp_params["param1_height"], top_ht[i]);}
		if (Rcpp::NumericVector::is_na(aerodynamic_ht[i])) {aerodynamic_ht[i] = top_ht[i];}
		if (Rcpp::NumericVector::is_na(mean_cr_width[i])) {
			mean_cr_width[i] = canopy_width_fun_cpp(sp_params["param0_cr_width"], sp_params["param1_cr_width"], mean_dbh[i]);
			mean_cr_width[i] = std::round(mean_cr_width[i]/0.01)*0.01;}
		if (Rcpp::NumericVector::is_na(mean_cr_depth[i])) {
			mean_cr_depth[i] = canopy_depth_Norway_fun_cpp(s_species, sp_params, 10 * mean_dbh[i], 10 * top_ht[i], "stand");
			if(mean_cr_depth[i] > mean_ht[i]){mean_cr_depth[i] = mean_ht[i];}
			mean_cr_depth[i] = std::round(mean_cr_depth[i]/0.01)*0.01;}
		if (Rcpp::IntegerVector::is_na(soil_group[i])) {soil_group[i] = 1;}
		if (Rcpp::IntegerVector::is_na(rooting[i])) {rooting[i] = 3;}
		if (Rcpp::NumericVector::is_na(gap_size[i])) {gap_size[i] = 0;}
		if (Rcpp::NumericVector::is_na(dist_edge[i])) {dist_edge[i] = mean_ht[i] * tree_heights_inside_forest;}
		
		
				
		// Advanced variables		
		std::string soil_rd = "soil_" + std::to_string(int(soil_group[i])) + "_rd_" + std::to_string(int(rooting[i]));
		if (Rcpp::NumericVector::is_na(elev_ht[i])) {elev_ht[i] = aerodynamic_ht[i];} //Default height to calculate elevated_crit is mean_ht
		if (Rcpp::NumericVector::is_na(ht_above_d[i])) {ht_above_d[i] = 10;}
		if (Rcpp::NumericVector::is_na(moe[i])) {moe[i] = sp_params["moe"];}
		if (Rcpp::NumericVector::is_na(mor[i])) {mor[i] = sp_params["mor"];}
		if (Rcpp::NumericVector::is_na(stem_density[i])) {stem_density[i] = sp_params["stem_density"];}
		if (Rcpp::NumericVector::is_na(crown_density[i])) {crown_density[i] = sp_params["canopy_density"];}
		if (Rcpp::NumericVector::is_na(fknot[i])) {fknot[i] = sp_params["fknot"];}
		if (Rcpp::NumericVector::is_na(c_drag[i])) {c_drag[i] = sp_params["c_drag"];}
		if (Rcpp::NumericVector::is_na(n_drag[i])) {n_drag[i] = sp_params["n_drag"];}
		if (Rcpp::NumericVector::is_na(drag_upper_limit[i])) {drag_upper_limit[i] = sp_params["drag_upper_limit"];}
		if (Rcpp::NumericVector::is_na(c_reg[i])) {c_reg[i] = sp_params["c_reg_" + soil_rd];}
		if (Rcpp::NumericVector::is_na(crown_off[i])) {crown_off[i] = 0;}
		
		
		if (Rcpp::NumericVector::is_na(stem_vol[i])) {
			if(s_species == "BI" || s_species == "NS" || s_species == "SP"){
				stem_vol[i] = stem_vol_Norway_fun_cpp(s_species, s_location, s_dbhclass, mean_dbh[i], mean_ht[i], sp_params);
			} else {
				stem_vol[i] = stem_vol_andretreslag_fun_cpp(s_species, mean_dbh[i], mean_ht[i], sp_params);
			}
		}
		
		if (Rcpp::NumericVector::is_na(crown_vol[i])) {crown_vol[i] = 1.0/3.0 * M_PI * mean_cr_depth[i] * pow(mean_cr_width[i]/2, 2);}
		
		stem_weight[i] = stem_vol[i] * stem_density[i];
		
		if(stem_weight[i] > sp_params["msw_" + soil_rd]) {max_stem_weight_warning[i] = "Warning: msw";} 
		
		// Wind variables - TO ADD		
		
	}
	
	// crown projected area:
	Rcpp::NumericVector	mean_cr_prarea;
	mean_cr_prarea = M_PI * pow((mean_cr_width/2.0), 2);
	// SNOW VARIABLES
	// Choose whether using snow_depth and density or snow_load - snow load need to be present in the input dataframe to be able to use it
	if(contains_cpp("snow_load", inputdata)){
		snow_load = inputdata["snow_load"];
	} else {
		// snow_depth
		Rcpp::NumericVector	snow_depth(n);
		snow_depth = extract_create_numcolumns_cpp("snow_depth", inputdata, n);
		for (int i = 0; i < n; ++i){
			if (Rcpp::NumericVector::is_na(snow_depth[i])) {snow_depth[i] = 0;}
		}
		// snow_density
		Rcpp::NumericVector	snow_density(n);
		snow_density = extract_create_numcolumns_cpp("snow_density", inputdata, n);
		for (int i = 0; i < n; ++i){
			if (Rcpp::NumericVector::is_na(snow_density[i])) {snow_density[i] = fgr_constants["snow_density_default"];}
		}
		// snow_load
		snow_load = snow_depth * snow_density;
	}
	
	// crown weight
	Rcpp::NumericVector crown_weight(n); 
	crown_weight = crown_vol * crown_density;
	// snow_weight
	Rcpp::NumericVector snow_weight;
	snow_weight = snow_load * mean_cr_prarea;
	

	Rcpp::List output_data;
	output_data["required"] = Rcpp::DataFrame::create(Rcpp::_["stand_id"] = Rcpp::clone(stand_id), Rcpp::_["species"] = Rcpp::clone(species),
	Rcpp::_["fylke"] = Rcpp::clone(fylke), Rcpp::_["mean_dbh"] = Rcpp::clone(mean_dbh), 
	Rcpp::_["top_ht"] = Rcpp::clone(top_ht), Rcpp::_["spacing"] = Rcpp::clone(spacing));	
	output_data["desirable"] = Rcpp::DataFrame::create(Rcpp::_["mean_ht"] = Rcpp::clone(mean_ht), Rcpp::_["aerodynamic_ht"] = Rcpp::clone(aerodynamic_ht),
	Rcpp::_["mean_cr_width"] = Rcpp::clone(mean_cr_width), Rcpp::_["mean_cr_depth"] = Rcpp::clone(mean_cr_depth),
	Rcpp::_["mean_cr_prarea"] = Rcpp::clone(mean_cr_prarea), Rcpp::_["crown_off"] = Rcpp::clone(crown_off),
	Rcpp::_["soil_group"] = Rcpp::clone(soil_group), Rcpp::_["rooting"] = Rcpp::clone(rooting),
	Rcpp::_["gap_size"] = Rcpp::clone(gap_size),
	Rcpp::_["dist_edge"] = Rcpp::clone(dist_edge));	
	output_data["advanced"] = Rcpp::DataFrame::create(Rcpp::_["elev_ht"] = Rcpp::clone(elev_ht), Rcpp::_["ht_above_d"] = Rcpp::clone(ht_above_d),
	Rcpp::_["moe"] = Rcpp::clone(moe), Rcpp::_["mor"] = Rcpp::clone(mor), 
	Rcpp::_["stem_density"] = Rcpp::clone(stem_density), Rcpp::_["crown_density"] = Rcpp::clone(crown_density), 
	Rcpp::_["fknot"] = Rcpp::clone(fknot), Rcpp::_["c_drag"] = Rcpp::clone(c_drag), 
	Rcpp::_["n_drag"] = Rcpp::clone(n_drag), Rcpp::_["drag_upper_limit"] = Rcpp::clone(drag_upper_limit),
	Rcpp::_["c_reg"] = Rcpp::clone(c_reg), 	
	Rcpp::_["stem_vol"] = Rcpp::clone(stem_vol), Rcpp::_["crown_vol"] = Rcpp::clone(crown_vol), 
	Rcpp::_["stem_weight"] = Rcpp::clone(stem_weight), Rcpp::_["crown_weight"] = Rcpp::clone(crown_weight));
	output_data["others"] = Rcpp::DataFrame::create(Rcpp::_["season"] = season, Rcpp::_["country"] = country, 
	Rcpp::_["snow_load"] = Rcpp::clone(snow_load), Rcpp::_["snow_weight"] = Rcpp::clone(snow_weight), 
	Rcpp::_["max_stem_weight_warning"] = Rcpp::clone(max_stem_weight_warning));	
	
	return output_data;
}


// MAIN FUNCTION FOR DATA PREPARATION - TMC METHOD
// [[Rcpp::export]]
Rcpp::List fg_tmc_dataprep_Norway_cpp(Rcpp::DataFrame inputdata, Rcpp::DataFrame fgr_constants, Rcpp::DataFrame species_parameters, std::string season = "winter",
                                     std::string country = "fgr"){
	// 1. Check essentials (stand_id, tree_id, species, tree_ht, dbh, spacing_current, stand_mean_dbh, stand_top_ht):
	//first step check that all the essential variables are present in the inputdata DataFrame
	//essential variables are the following: stand_id, tree_id, species, tree_ht, dbh, spacing_current, stand_mean_dbh, stand_top_ht
	try
	{
		Rcpp::StringVector stand_id = inputdata["stand_id"];
		Rcpp::StringVector tree_id = inputdata["tree_id"];
		Rcpp::StringVector species = inputdata["species"];
		Rcpp::NumericVector tree_ht = inputdata["tree_ht"]; 
		Rcpp::NumericVector dbh = inputdata["dbh"]; 
		Rcpp::NumericVector spacing_current = inputdata["spacing_current"];
		Rcpp::NumericVector stand_mean_dbh = inputdata["stand_mean_dbh"];
		Rcpp::NumericVector stand_top_ht = inputdata["stand_top_ht"];
	}
	catch (...)
	{
		std::cout << "Missing some key variables! Check that stand_id, tree_id, species, tree_ht, dbh, spacing_current, stand_mean_dbh and stand_top_ht are in your dataframe"<<std::endl;
		return EXIT_FAILURE;
	}	
	
	// If all variables are present, load them
	Rcpp::StringVector stand_id = inputdata["stand_id"];
	Rcpp::StringVector tree_id = inputdata["tree_id"];
	Rcpp::StringVector species = inputdata["species"];
	Rcpp::NumericVector tree_ht = inputdata["tree_ht"]; 
	Rcpp::NumericVector dbh = inputdata["dbh"]; 
	Rcpp::NumericVector spacing_current = inputdata["spacing_current"];
	Rcpp::NumericVector stand_mean_dbh = inputdata["stand_mean_dbh"];
	Rcpp::NumericVector stand_top_ht = inputdata["stand_top_ht"];
	spacing_current = round(spacing_current, 1);
	// Check if any of the variables contain NA	
	Rcpp::LogicalVector stand_id_na = is_na(stand_id);
	Rcpp::LogicalVector tree_id_na = is_na(tree_id);
	Rcpp::LogicalVector species_na = is_na(species);
	Rcpp::LogicalVector tree_ht_na = is_na(tree_ht);
	Rcpp::LogicalVector dbh_na = is_na(dbh);
	Rcpp::LogicalVector spacing_current_na = is_na(spacing_current);
	Rcpp::LogicalVector stand_mean_dbh_na = is_na(stand_mean_dbh);
	Rcpp::LogicalVector stand_top_ht_na = is_na(stand_top_ht);
	if(is_true(any(stand_id_na)) || is_true(any(tree_id_na)) || is_true(any(species_na)) || is_true(any(tree_ht_na)) || 
	is_true(any(dbh_na)) || is_true(any(spacing_current_na)) || is_true(any(stand_mean_dbh_na)) || is_true(any(stand_top_ht_na))){
		std::cout<<"One or more of the key variables contains NA! Check your input data"<<std::endl;
		return EXIT_FAILURE;
	}
	// Get the length of the dataframe
	int n = stand_id.size();
	
	// equivalent_mean_ht
	Rcpp::NumericVector	equivalent_mean_ht;
	equivalent_mean_ht = extract_create_numcolumns_cpp("equivalent_mean_ht", inputdata, n);
	// predominant species
	Rcpp::StringVector	predominant_species;
	predominant_species = extract_create_strcolumns_cpp("predominant_species", inputdata, n);
	for (int i = 0; i < n; ++i){
		if (Rcpp::StringVector::is_na(predominant_species[i])) {predominant_species[i] = species[i];}
	}
	
	// 1.2 Get the location, dbh class and leaf status
	// fylke - the fylke is important for calculating appropriate stem volume. If not present it will be defaulted to alle
	Rcpp::IntegerVector fylke;
	fylke = extract_create_intcolumns_cpp("fylke", inputdata, n);
	// create the location vector for the stem volume and crown dimension functions (only for NS, SP and BI)
	Rcpp::StringVector location(n);
	for (int i = 0; i < n; ++i){
		if(species[i] == "NS" || species[i] == "SP"){
			if(fylke[i] != 11 && fylke[i] != 12 && fylke[i] != 14 && fylke[i] != 15 ){location[i] = "annen";}
			else {location[i] = "vest";}
		} else if(species[i] == "BI"){location[i] = "alle";}
		else {location[i] = NA_STRING;}
	}
	// dbh_class - required for appropriate stem volume functions for NS, SP and BI
	Rcpp::StringVector dbh_class(n);
	for (int i = 0; i < n; ++i){
		if((location[i] == "vest") || (species[i] == "BI")){dbh_class[i] = "alle";}
		else if(species[i] == "SP"){
			if(dbh[i] < 11.1){dbh_class[i] = "liten";}
			else {dbh_class[i] = "stor";}
		}
		else if(species[i] == "NS"){
			if(dbh[i] < 10.1){dbh_class[i] = "liten";}
			else if(dbh[i] >= 10.1 && dbh[i] < 12.9){dbh_class[i] = "medium";}
			else{dbh_class[i] = "stor";}
		}
		else {dbh_class[i] = NA_STRING;}
	}		
	// stand_dbh_class - required for appropriate stem volume functions for NS, SP and BI for stand level
	Rcpp::StringVector stand_dbh_class(n);
	for (int i = 0; i < n; ++i){
		if((location[i] == "vest") || (predominant_species[i] == "BI")){stand_dbh_class[i] = "alle";}
		else if(predominant_species[i] == "SP"){
			if(stand_mean_dbh[i] < 11.1){stand_dbh_class[i] = "liten";}
			else {stand_dbh_class[i] = "stor";}
		}
		else if(predominant_species[i] == "NS"){
			if(stand_mean_dbh[i] < 10.1){stand_dbh_class[i] = "liten";}
			else if(stand_mean_dbh[i] >= 10.1 && stand_mean_dbh[i] < 12.9){stand_dbh_class[i] = "medium";}
			else{stand_dbh_class[i] = "stor";}
		}
		else {stand_dbh_class[i] = NA_STRING;}
	}		
	// leafstatus - for birch in Norway only
	Rcpp::StringVector leafstatus(n);
	for (int i = 0; i < n; ++i){
		if(species[i] == "BI"){
			if(season == "summer"){leafstatus[i] = "lo";}
			else {leafstatus[i] = "ll";}
		}
		else {leafstatus[i] = NA_STRING;}
	}	
	// stand_leafstatus - for birch in Norway only
	Rcpp::StringVector stand_leafstatus(n);
	for (int i = 0; i < n; ++i){
		if(predominant_species[i] == "BI"){
			if(season == "summer"){stand_leafstatus[i] = "lo";}
			else {stand_leafstatus[i] = "ll";}
		}
		else {stand_leafstatus[i] = NA_STRING;}
	}
			
	// 2. Extract the desirable variables if they are present
	// If a variable is not present, create it with default calculations, if present and contains NAs, we will replace the NAs with the default values
	// cr_width
	Rcpp::NumericVector	cr_width;
	cr_width = extract_create_numcolumns_cpp("cr_width", inputdata, n);
	// cr_depth
	Rcpp::NumericVector	cr_depth;
	cr_depth = extract_create_numcolumns_cpp("cr_depth", inputdata, n);
	// crown_off
	Rcpp::NumericVector	crown_off;
	crown_off = extract_create_numcolumns_cpp("crown_off", inputdata, n);
	// stand_cr_width
	Rcpp::NumericVector	stand_cr_width;
	stand_cr_width = extract_create_numcolumns_cpp("stand_cr_width", inputdata, n);
	// stand_cr_depth
	Rcpp::NumericVector	stand_cr_depth;
	stand_cr_depth = extract_create_numcolumns_cpp("stand_cr_depth", inputdata, n);
	// soil_group
	Rcpp::IntegerVector	soil_group;
	soil_group = extract_create_intcolumns_cpp("soil_group", inputdata, n);
	// rooting
	Rcpp::IntegerVector	rooting;
	rooting = extract_create_intcolumns_cpp("rooting", inputdata, n);
	// gap_size
	Rcpp::NumericVector	gap_size;
	gap_size = extract_create_numcolumns_cpp("gap_size", inputdata, n);
	// dist_edge
	Rcpp::NumericVector	dist_edge;
	dist_edge = extract_create_numcolumns_cpp("dist_edge", inputdata, n);	
	
	// 3. Advanced variables
	// ci
	Rcpp::StringVector	ci;
	ci = extract_create_strcolumns_cpp("ci", inputdata, n);
	// ci_value
	Rcpp::NumericVector	ci_value;
	ci_value = extract_create_numcolumns_cpp("ci_value", inputdata, n);
	// spacing_before
	Rcpp::NumericVector	spacing_before;
	spacing_before = extract_create_numcolumns_cpp("spacing_before", inputdata, n);
	// years_since_thin
	Rcpp::IntegerVector	years_since_thin;
	years_since_thin = extract_create_intcolumns_cpp("years_since_thin", inputdata, n);
	// elev_ht
	Rcpp::NumericVector	elev_ht;
	elev_ht = extract_create_numcolumns_cpp("elev_ht", inputdata, n);
	// aerodynamic_ht
	Rcpp::NumericVector	aerodynamic_ht;
	aerodynamic_ht = extract_create_numcolumns_cpp("aerodynamic_ht", inputdata, n);
	// ht_above_d
	Rcpp::NumericVector	ht_above_d;
	ht_above_d = extract_create_numcolumns_cpp("ht_above_d", inputdata, n);
	// MOE
	Rcpp::NumericVector	moe;
	moe = extract_create_numcolumns_cpp("moe", inputdata, n);
	// MOR
	Rcpp::NumericVector	mor;
	mor = extract_create_numcolumns_cpp("mor", inputdata, n);
	// stem_density
	Rcpp::NumericVector	stem_density;
	stem_density = extract_create_numcolumns_cpp("stem_density", inputdata, n);
	// crown_density
	Rcpp::NumericVector	crown_density;
	crown_density = extract_create_numcolumns_cpp("crown_density", inputdata, n);
	// fknot
	Rcpp::NumericVector	fknot;
	fknot = extract_create_numcolumns_cpp("fknot", inputdata, n);
	// c_drag
	Rcpp::NumericVector	c_drag;
	c_drag = extract_create_numcolumns_cpp("c_drag", inputdata, n);
	// n_drag
	Rcpp::NumericVector	n_drag;
	n_drag = extract_create_numcolumns_cpp("n_drag", inputdata, n);
	// drag_upper_limit
	Rcpp::NumericVector	drag_upper_limit;
	drag_upper_limit = extract_create_numcolumns_cpp("drag_upper_limit", inputdata, n);
	// c_reg
	Rcpp::NumericVector	c_reg;
	c_reg = extract_create_numcolumns_cpp("c_reg", inputdata, n);
	// stem_vol
	Rcpp::NumericVector	stem_vol;
	stem_vol = extract_create_numcolumns_cpp("stem_vol", inputdata, n);
	// crown_vol
	Rcpp::NumericVector	crown_vol;
	crown_vol = extract_create_numcolumns_cpp("crown_vol", inputdata, n);
	// stem_weight
	Rcpp::NumericVector stem_weight(n);
	stem_weight = extract_create_numcolumns_cpp("stem_weight", inputdata, n);
	// max_stem_weight_warning
	Rcpp::CharacterVector max_stem_weight_warning(n);
	for (int i = 0; i < n; ++i){max_stem_weight_warning[i] = NA_STRING;}
	//Initiate snow_load
	Rcpp::NumericVector snow_load(n);
	
	
	// 4. Wind climate variables - to do
	
	
	//5. Replace the NA values with default calculations
	// Desirable variables	
	double tree_heights_inside_forest = fgr_constants["tree_heights_inside_forest"];
	for (int i = 0; i < n; ++i){
		std::string s_species = Rcpp::as<std::string>(species[i]);
		std::string s_predomsp = Rcpp::as<std::string>(predominant_species[i]);
		std::string s_location = Rcpp::as<std::string>(location[i]);
		std::string s_dbhclass = Rcpp::as<std::string>(dbh_class[i]);
		std::string s_sdbhclass = Rcpp::as<std::string>(stand_dbh_class[i]);
		std::string s_leaf = Rcpp::as<std::string>(leafstatus[i]);
		std::string s_sleaf = Rcpp::as<std::string>(stand_leafstatus[i]);
		// select the appropriate species_parameter row
		// if the species is not present or the species 2-letter identifier not correct, display a message
		try{Rcpp::NumericVector sp_params = sp_params_Norway_fun_cpp(species_parameters, s_species, s_location, s_dbhclass, s_leaf);} catch (...) {
			std::cout << "This species is not present in the species_parameters dataframe or its 2-letter identifier is not correct." << std::endl;
			return EXIT_FAILURE;
		}
		// if the species is other than NS, BI or SP, use the default fgr species parameters, otherwise use the NIBIO-Norway ones
		Rcpp::NumericVector sp_params = sp_params_Norway_fun_cpp(species_parameters, s_species, s_location, s_dbhclass, s_leaf);
		
		
		// select the appropriate species_parameter row for the surrounding trees - predominant species
		// if the species is not present or the species 2-letter identifier not correct, display a message
		try{Rcpp::NumericVector predomsp_params = sp_params_Norway_fun_cpp(species_parameters, s_predomsp, s_location, s_sdbhclass, s_sleaf);} catch (...) {
			std::cout << "This species is not present in the species_parameters dataframe or its 2-letter identifier is not correct." << std::endl;
			return EXIT_FAILURE;
		}
		Rcpp::NumericVector predomsp_params = sp_params_Norway_fun_cpp(species_parameters, s_predomsp, s_location, s_sdbhclass, s_sleaf);
		
		// Calculate equivalent_mean_ht
		equivalent_mean_ht[i] = eq_mean_ht_fun_cpp(stand_top_ht[i]);	
		
		if (Rcpp::NumericVector::is_na(cr_width[i])) {
			cr_width[i] = canopy_width_fun_cpp(sp_params["param0_cr_width"], sp_params["param1_cr_width"], dbh[i]);
			cr_width[i] = std::round(cr_width[i]/0.01)*0.01;}
		if (Rcpp::NumericVector::is_na(cr_depth[i])) {
			cr_depth[i] = canopy_depth_Norway_fun_cpp(s_species, sp_params, 10 * dbh[i], 10 * tree_ht[i], "tree");
			if(cr_depth[i] > tree_ht[i]){cr_depth[i] = tree_ht[i];}
			cr_depth[i] = std::round(cr_depth[i]/0.01)*0.01;}
		if (Rcpp::NumericVector::is_na(stand_cr_width[i])) {
			stand_cr_width[i] = canopy_width_fun_cpp(predomsp_params["param0_cr_width"], predomsp_params["param1_cr_width"], stand_mean_dbh[i]);
			stand_cr_width[i] = std::round(stand_cr_width[i]/0.01)*0.01;}
		if (Rcpp::NumericVector::is_na(stand_cr_depth[i])) {
			stand_cr_depth[i] = canopy_depth_Norway_fun_cpp(s_predomsp, predomsp_params, 10 * stand_mean_dbh[i], 10 * equivalent_mean_ht[i], "stand");
			if(stand_cr_depth[i] > equivalent_mean_ht[i]){stand_cr_depth[i] = equivalent_mean_ht[i];}
			stand_cr_depth[i] = std::round(stand_cr_depth[i]/0.01)*0.01;}
		if (Rcpp::IntegerVector::is_na(soil_group[i])) {soil_group[i] = 1;}
		if (Rcpp::IntegerVector::is_na(rooting[i])) {rooting[i] = 3;}
		if (Rcpp::NumericVector::is_na(gap_size[i])) {gap_size[i] = 0;}
		if (Rcpp::NumericVector::is_na(dist_edge[i])) {dist_edge[i] = equivalent_mean_ht[i] * tree_heights_inside_forest;}
		
		// Advanced variables
		std::string soil_rd = "soil_" + std::to_string(int(soil_group[i])) + "_rd_" + std::to_string(int(rooting[i]));
		if (Rcpp::StringVector::is_na(ci[i])) {ci[i] = "none";}	
		// we don't replace the NA values for ci_value since the default with ci = "none" are NA values		
		if (Rcpp::NumericVector::is_na(spacing_before[i])) {spacing_before[i] = spacing_current[i];}		
		if (Rcpp::IntegerVector::is_na(years_since_thin[i])) {years_since_thin[i] = 5;}	
		if (years_since_thin[i] > 5){years_since_thin[i] = 5;} else if(years_since_thin[i] < 0){years_since_thin[i] = 0;}
		if (Rcpp::NumericVector::is_na(elev_ht[i])) {elev_ht[i] = 1.05 * stand_top_ht[i];} 		
		if (Rcpp::NumericVector::is_na(aerodynamic_ht[i])) {aerodynamic_ht[i] = stand_top_ht[i];}	
		if (Rcpp::NumericVector::is_na(ht_above_d[i])) {ht_above_d[i] = 10;}
		if (Rcpp::NumericVector::is_na(moe[i])) {moe[i] = sp_params["moe"];}
		if (Rcpp::NumericVector::is_na(mor[i])) {mor[i] = sp_params["mor"];}
		if (Rcpp::NumericVector::is_na(stem_density[i])) {stem_density[i] = sp_params["stem_density"];}
		if (Rcpp::NumericVector::is_na(crown_density[i])) {crown_density[i] = sp_params["canopy_density"];}
		if (Rcpp::NumericVector::is_na(crown_off[i])) {crown_off[i] = 0;}
		if (Rcpp::NumericVector::is_na(fknot[i])) {fknot[i] = sp_params["fknot"];}
		if (Rcpp::NumericVector::is_na(c_drag[i])) {c_drag[i] = predomsp_params["c_drag"];}
		if (Rcpp::NumericVector::is_na(n_drag[i])) {n_drag[i] = predomsp_params["n_drag"];}
		if (Rcpp::NumericVector::is_na(drag_upper_limit[i])) {drag_upper_limit[i] = predomsp_params["drag_upper_limit"];}	
		if (Rcpp::NumericVector::is_na(c_reg[i])) {c_reg[i] = sp_params["c_reg_" + soil_rd];}		
		if (Rcpp::NumericVector::is_na(stem_vol[i])) {
			if(s_species == "BI" || s_species == "NS" || s_species == "SP"){
				stem_vol[i] = stem_vol_Norway_fun_cpp(s_species, s_location, s_dbhclass, dbh[i], tree_ht[i], sp_params);
			} else {
				stem_vol[i] = stem_vol_andretreslag_fun_cpp(s_species, dbh[i], tree_ht[i], sp_params);
			}
		}
		if (Rcpp::NumericVector::is_na(crown_vol[i])) {crown_vol[i] = 1.0/3.0 * M_PI * cr_depth[i] * pow(cr_width[i]/2, 2);}
		
		stem_weight[i] = stem_vol[i] * stem_density[i];		
		//msw[i] = sp_params["msw_soil_" + std::to_string(int(soil_group[i])) + "_rd_" + std::to_string(int(rooting[i]))];
		if(stem_weight[i] > sp_params["msw_" + soil_rd]) {max_stem_weight_warning[i] = "Warning: msw";} 
		
		// Wind variables - to do
		
	}
	
	// crown projected area:
	Rcpp::NumericVector	cr_prarea;
	cr_prarea = M_PI * pow((cr_width/2.0), 2);
	
	// crown weight
	Rcpp::NumericVector crown_weight(n); 
	crown_weight = crown_vol * crown_density;
	
	// SNOW VARIABLES
	// Choose whether using snow_depth and density or snow_load - snow load need to be present in the input dataframe to be able to use it
	if(contains_cpp("snow_load", inputdata)){
		snow_load = inputdata["snow_load"];
	} else {
		// snow_depth
		Rcpp::NumericVector	snow_depth(n);
		snow_depth = extract_create_numcolumns_cpp("snow_depth", inputdata, n);
		for (int i = 0; i < n; ++i){
			if (Rcpp::NumericVector::is_na(snow_depth[i])) {snow_depth[i] = 0;}
		}
		// snow_density
		Rcpp::NumericVector	snow_density(n);
		snow_density = extract_create_numcolumns_cpp("snow_density", inputdata, n);
		for (int i = 0; i < n; ++i){
			if (Rcpp::NumericVector::is_na(snow_density[i])) {snow_density[i] = fgr_constants["snow_density_default"];}
		}
		// snow_load
		snow_load = snow_depth * snow_density;
	}
	
	// snow_weight
	Rcpp::NumericVector snow_weight;
	snow_weight = snow_load * cr_prarea;
	
	spacing_before = round(spacing_before, 1);
	Rcpp::List output_data;
	output_data["required"] = Rcpp::DataFrame::create(Rcpp::_["stand_id"] = Rcpp::clone(stand_id), Rcpp::_["tree_id"] = Rcpp::clone(tree_id), 
	Rcpp::_["species"] = Rcpp::clone(species), Rcpp::_["predominant_species"] = Rcpp::clone(predominant_species),
	Rcpp::_["tree_ht"] = Rcpp::clone(tree_ht), Rcpp::_["dbh"] = Rcpp::clone(dbh), Rcpp::_["spacing_current"] = Rcpp::clone(spacing_current),
	Rcpp::_["stand_mean_dbh"] = Rcpp::clone(stand_mean_dbh), Rcpp::_["stand_top_ht"] = Rcpp::clone(stand_top_ht), 
	Rcpp::_["equivalent_mean_ht"] = Rcpp::clone(equivalent_mean_ht));
	output_data["desirable"] = Rcpp::DataFrame::create(Rcpp::_["cr_width"] = Rcpp::clone(cr_width), Rcpp::_["cr_depth"] = Rcpp::clone(cr_depth),
	Rcpp::_["cr_prarea"] = Rcpp::clone(cr_prarea), Rcpp::_["crown_off"] = Rcpp::clone(crown_off),
	Rcpp::_["stand_cr_width"] = Rcpp::clone(stand_cr_width), Rcpp::_["stand_cr_depth"] = Rcpp::clone(stand_cr_depth),
	Rcpp::_["soil_group"] = Rcpp::clone(soil_group), Rcpp::_["rooting"] = Rcpp::clone(rooting),
	Rcpp::_["gap_size"] = Rcpp::clone(gap_size), Rcpp::_["dist_edge"] = Rcpp::clone(dist_edge));
	output_data["advanced"] = Rcpp::DataFrame::create(Rcpp::_["ci"] = Rcpp::clone(ci), Rcpp::_["ci_value"] = Rcpp::clone(ci_value),
	Rcpp::_["spacing_before"] = Rcpp::clone(spacing_before), Rcpp::_["years_since_thin"] = Rcpp::clone(years_since_thin),	
	Rcpp::_["elev_ht"] = Rcpp::clone(elev_ht), Rcpp::_["aerodynamic_ht"] = Rcpp::clone(aerodynamic_ht), Rcpp::_["ht_above_d"] = Rcpp::clone(ht_above_d),
	Rcpp::_["moe"] = Rcpp::clone(moe), Rcpp::_["mor"] = Rcpp::clone(mor), 
	Rcpp::_["stem_density"] = Rcpp::clone(stem_density), Rcpp::_["crown_density"] = Rcpp::clone(crown_density), 
	Rcpp::_["fknot"] = Rcpp::clone(fknot), Rcpp::_["c_drag"] = Rcpp::clone(c_drag), 
	Rcpp::_["n_drag"] = Rcpp::clone(n_drag), Rcpp::_["drag_upper_limit"] = Rcpp::clone(drag_upper_limit),
	Rcpp::_["c_reg"] = Rcpp::clone(c_reg), 	
	Rcpp::_["stem_vol"] = Rcpp::clone(stem_vol), Rcpp::_["crown_vol"] = Rcpp::clone(crown_vol), 
	Rcpp::_["stem_weight"] = Rcpp::clone(stem_weight), Rcpp::_["crown_weight"] = Rcpp::clone(crown_weight));
	output_data["others"] = Rcpp::DataFrame::create(Rcpp::_["season"] = season, Rcpp::_["country"] = country, 
	Rcpp::_["snow_load"] = Rcpp::clone(snow_load), Rcpp::_["snow_weight"] = Rcpp::clone(snow_weight), 
	Rcpp::_["max_stem_weight_warning"] = Rcpp::clone(max_stem_weight_warning));
	
	return output_data;
}