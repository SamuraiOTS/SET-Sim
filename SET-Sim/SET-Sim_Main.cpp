

// ============================================
// Defining the Inclusions
// ============================================
#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>
#include <vector>
using namespace std;


// ============================================
// PHYSICAL CONSTANTS
// ============================================
const double G = 6.674e-11;              // Gravitational constant (m^3/kg·s^2)
const double c = 299792458;              // Speed of light (m/s)
const double pi = 3.14159265359;         // Pi

// ============================================
// EARTH PARAMETERS
// ============================================
double M = 5.972e24;                     // Earth's mass (kg)
double omega = 7.292e-5;                 // Earth's angular velocity (rad/s)
double g_0 = 9.807;                      // Earth's surface gravity (m/s^2)
double r_0 = 6.371e6;                    // Earth's radius (m)
double r_GEO = 4.2164e7;                 // Radius from Earth's center to GEO (m)
double T = 86400;                        // Orbital period for GEO (seconds)

// ============================================
// MATERIAL PROPERTIES
// ============================================
double material_density;                 // Material density (kg/m^3)
double tensile_strength;                 // Tensile strength of material (Pa)
double yield_strength;                   // Yield strength (Pa)
double sigma_max;                        // Maximum tensile stress (Pa)
double safety_factor = 2.0;              // Safety factor (typically 2-3)
double FOM;                              // Figure of Merit (km)

// ============================================
// CABLE GEOMETRY
// ============================================
double r;                                // Distance from Earth's center (m)
double r_top;                            // Top of cable/counterweight position (m)
double cross_area;                       // Cross-sectional area at point r (m^2)
double A_0;                              // Cross-sectional area at Earth's surface (m^2)
double A_GEO;                            // Cross-sectional area at GEO (m^2)
double A_top = 1e-4;                     // Start with 1 cm² at the top where stress is minimal
double L_c = 4.96e6;                     // Characteristic length for Earth (m)
double total_cable_length;               // Total cable length (m)
double taper_ratio;                      // Taper ratio (A_GEO/A_surface)
double cumulative_tension = 0;           // For the calculation of weight -> cross section

// ============================================
// FORCES AND STRESSES
// ============================================
double g_eff;                            // Effective gravity at radius r (m/s^2)
double tension;                          // Cable tension at point r (N)
double T_GEO;                            // Tension at geostationary orbit (N)
double stress;                           // Stress at point of interest (Pa)
int failure_height;                      // Failure height if stress exceeds sigma_max
double max_stress = 0;

// ============================================
// MASS CALCULATIONS
// ============================================
double M_cable;                          // Total cable mass (kg)
double M_uniform;                        // Mass with uniform cross-section (kg)
double MR;                               // Mass ratio (tapered to uniform)
double M_cw;                             // Counterweight mass (kg)
double r_cw;                             // Counterweight radius from Earth's center (m)

// ============================================
// CLIMBER PARAMETERS
// ============================================
double climber_mass;                     // Climber mass (kg)
double v_climber;                        // Velocity of climber along cable (m/s)
double climb_velo;                       // Climbing velocity (m/s)
double a_Coriolis;                       // Coriolis acceleration for climbers (m/s^2)
double E_GEO;                            // Gravitational potential energy to GEO (J)
double P;                                // Power requirement for climber (W)
double P_losses;                         // Friction and efficiency losses (W)

// ============================================
// DYNAMIC PROPERTIES
// ============================================
double f_1;                              // Fundamental frequency (Hz)
double longitudinal_velo;                // Longitudinal wave velocity (m/s)

// ============================================
// ENVIRONMENTAL FORCES
// ============================================
double F_solar;                          // Solar radiation pressure force (N)
double solar_const = 1361;               // Solar constant (W/m^2)
double albedo_coeff;                     // Albedo coefficient (0-1)
double angle_of_incidence;               // Angle of incidence for solar radiation (rad)
double F_drag;                           // Atmospheric drag force (N)
double atmos_density;                    // Atmospheric density (kg/m^3)
double drag_coef;                        // Drag coefficient
double relative_velo;                    // Relative velocity for drag calc (m/s)

// ============================================
// Defining the Structs
// ============================================
struct MaterialProperties {
    string name;                        // kg/m^3
    double material_density;            // Pa (ultimate tensile strength)
    double tensile_strength;            // Pa (yield strength for sigma_max)
    double yield_strength;              // Pa (maximum allowable stress with safety factor)
    double sigma_max;                   // Safety factor applied
    double safety_factor;               // Figure of Merit (tensile_strength/density) in km
    double FOM;
	double price_per_kg;                // $/kg
	double totalCost;                   // total cost of the cable
}matProp;

struct locationProperties {
    string name;
    double latitude;                    // degrees
    double longitude;                   // degrees
    double elevation;                   // meters above sea level
    double atmos_density;               // kg/m^3
    double wind_speed;                  // m/s
    double seismic_activity;            // Richter scale
}locProp;

struct locationEffects {
	int    sectnum;					   // Section number/index
	double coriolis_parameter;          // Coriolis parameter at the location (1/s)
	double centrifugal_effect;          // Centrifugal effect at the location (m/s^2)
	double gravity_effect;              // Local gravity effect at the location (m/s^2)
    
}locEff;

//store as many properties of each section as needed 
struct sectionProperties {

	int sectNum;                        // Section number/index
	double radius;                      // Distance from Earth's center (m)  ====== use this for calculations :(  ======
    double length;                      // Length of the section (m)
    double cross_area;                  // Cross-sectional area (m^2)
    double mass;                        // Mass of the section (kg)
    double tension;                     // Tension at the section (N)
	double stress;                      // Stress at the section (Pa)
}sectProp;

//this is to hold the effects applied on each section
struct sectionEffects {

	int sectNum;                        // Section number/index
    double wind_load;                   // Wind load on the section (N)
    double seismic_load;                // Seismic load on the section (N)
    double thermal_expansion;           // Thermal expansion effect (m)
	double climber_load;                // Load from climbers on the section (N)
	double counterweight_load;          // Load from counterweight on the section (N)
	double solar_radiation_load;        // Load from solar radiation pressure on the section (N)
	double gravity_effect;              // Effect of local gravity variations (m/s^2)

}sectEff;

// ============================================
// Defining the arrays
// ============================================
vector<double> cableLength;
vector<double> cableWidth;
vector<double> cableWeight;
vector<double> cableSegmentWeight;
vector<double> latitudeList;
vector<double> longitudeList;
vector<double> elevationList;
vector<double> atmos_densityList;
vector<double> wind_speedList;
vector<double> seismic_activityList;
vector<locationProperties> locationList;     // to hold multiple location properties (i don't know wether to use a struct or just multiple arrays here.)
vector<sectionProperties> sectionList;       // to hold multiple section properties (dude, structs are life atp)
vector<sectionEffects> sectionEffectsList;   // to hold multiple section effects 
vector<locationEffects> locationEffectsList; // to hold multiple location effects

// ============================================
// Settings
// ============================================
// This is to define any settings like DEBUG that helps me understand where i messed up
// Made it verbose for no actual reason. check debug function to see details
int debug = 0;


/*
A = all debug level (TODO)
M = Main function debug level
L = Location debug level
C = Cable debug level

*/


// ============================================
// TODO LIST
// ============================================
// Add in wind loading, seismic activity, and other environmental factors
// Add in climber effects
// Add in location effects
// Add in dynamic properties
// Add in counterweight effects
// Add in solar and atmospheric drag effects
// Make sure that each section is stored in the sectionProperties struct for future use
// Make elevatorGeometry function that calls this one and adds in the other effects I.E. extra cables, climbers, counterweight, etc.
// Flesh out any other functions that are not yet
// Make a GUI (maybe?) -> Python? 
	//Full graphics and visualization of the cable and climber system need stress graphs and other data visualizations



// ============================================
// Actual gosh darn code finally
// ============================================


//calculate all of the base orbital mechanics
double orbitalMechanics(double r = r_GEO) {
    r_GEO = pow((G * M * pow(T, 2) / (4 * pow(pi, 2))), (1.0 / 3.0)); // calculate the radius from center of earth to Geostationary orbit

    //to calculate valid effective gravity at point r
    g_eff = (G * (M / pow(r, 2))) - (pow(omega, 2) * r); // calculating gravity is (earths gravity - centrifugal force) should be close to 0 at r_GEO and 9.8 at r_0

    cout << "r_GEO: " << r_GEO << endl;
    cout << "g_eff: " << g_eff;
    return g_eff;
}


//this is to add the variables to the construct
void materialProperties(string n, double density, double tensile, double yield_strength, double sf = 2.0) {
    matProp.name = n;
    matProp.material_density = density;
    matProp.tensile_strength = tensile;
    matProp.yield_strength = yield_strength;
    matProp.safety_factor = sf;
    //calculations for finding the max strength before fail accounting for the safetyfactor and FOM
    matProp.sigma_max = (yield_strength / sf);
    matProp.FOM = (tensile / (density * g_0) / 1000.0);
    
}

//section properties function to make custom setions if needed
void makeSectionProperties(int sectNum, double density, double tensile, double yield, double sf) {
    sectProp.sectNum = sectNum;
    //this is to define custom sections for different possible compositions
	//TODO: make this work with the rest of the program
	//TODO: add in the ability to define segment connections at specific heights i.e. GEO, 100km, etc. OR if a user wants to use multiple materials, at material connection, average the two materials' properties for a more generalized transition
    sectionList.push_back(sectProp);

}

//this is to define the material properties
void materialDefinition() {
    cout << "what material are you looking to use? S - steel || C - copper || N - Carbon Nanotubes || T - Titanium || O - other: " << endl;
    char choice;
    cin >> choice;
    switch (toupper(choice)) {
        case 'S':
            cout << "You have chosen Steel." << endl;
            material_density = 7850.0;                 // Material density (kg/m³)
            tensile_strength = 1.52e9;                 // Tensile strength of material (Pa)
            yield_strength = 1.24e9;                   // yield strength (Pa)
            safety_factor = 2.0;                       // Safety factor (typically 2-3)
            materialProperties("Steel", material_density, tensile_strength, yield_strength, safety_factor);
            cout << "Global variables: " << endl;
            cout << "Mat_density: " << material_density << endl << "Tensile_strength: " << tensile_strength << endl << "Yield_strength: " << yield_strength << endl << "Safety Factor: " << safety_factor << endl;
            cout << "Calculated: " << endl << "FOM: " << matProp.FOM << endl << "Sigma Max: " << matProp.sigma_max << endl;
            break;
        case 'C': {
            string copperOption;
            cout << "You have chosen Copper." << endl;
            cout << "Is the copper cold-worked? (y/N): ";
            cin.ignore(); // Clear leftover newline from previous input
            getline(cin, copperOption);
            material_density = 8960.0;                 // Material density (kg/m³)
            safety_factor = 2.0;                       // Safety factor (typically 2-3)
            if (!copperOption.empty() && (copperOption[0] == 'y' || copperOption[0] == 'Y')) {
                tensile_strength = 3.5e8;              // Tensile strength of material (Pa)
                yield_strength = 4.0e8;                // yield strength (Pa)
            }
            else {
                tensile_strength = 2.0e9;             // Tensile strength of material (Pa)
                yield_strength = 8.0e7;               // yield strength (Pa)
            }
            materialProperties("Copper", material_density, tensile_strength, yield_strength, safety_factor);
            cout << "Global variables: " << endl;
            cout << "Mat_density: " << material_density << endl << "Tensile_strength: " << tensile_strength << endl << "Yield_strength: " << yield_strength << endl << "Safety Factor: " << safety_factor << endl;
            cout << "Calculated: " << endl << "FOM: " << matProp.FOM << endl << "Sigma Max: " << matProp.sigma_max;
            break;
        }
        case 'N':
            cout << "You have chosen Carbon Nanotubes." << endl;
            material_density = 1400;                 // Material density (kg/m³)
            tensile_strength = 1.5e11;               // Tensile strength of material (Pa)
            yield_strength = 6.0e10;                 // yield strength (Pa)
            safety_factor = 2.0;                     // Safety factor (typically 2-3)
            materialProperties("Carbon Nanotubes", material_density, tensile_strength, yield_strength, safety_factor);
            cout << "Global variables: " << endl;
            cout << "Mat_density: " << material_density << endl << "Tensile_strength: " << tensile_strength << endl << "Yield_strength: " << yield_strength << endl << "Safety Factor: " << safety_factor << endl;
            cout << "Calculated: " << endl << "FOM: " << matProp.FOM << endl << "Sigma Max: " << matProp.sigma_max << endl;
            break;
        case 'T':
            cout << "You have chosen Titanium." << endl;
            material_density = 4506.0;                 // Material density (kg/m³)
            tensile_strength = 1e8;                    // Tensile strength of material (Pa)
            yield_strength = 1.4e8;                    // yield strength (Pa)
            safety_factor = 2.0;                       // Safety factor (typically 2-3)
            materialProperties("Titanium", material_density, tensile_strength, yield_strength, safety_factor);
            cout << "Global variables: " << endl;
            cout << "Mat_density: " << material_density << endl << "Tensile_strength: " << tensile_strength << endl << "Yield_strength: " << yield_strength << endl << "Safety Factor: " << safety_factor << endl;
            cout << "Calculated: " << endl << "FOM: " << matProp.FOM << endl << "Sigma Max: " << matProp.sigma_max << endl;
            break;
        case 'O':
            cout << "You have chosen to define these yourself." << endl;
            cout << "Material Density (kg/m^3): ";
            cin >> material_density;          // Material density (kg/m³)
            cout << endl << "Tensile Strength (Pa): ";
            cin >> tensile_strength;          // Tensile strength of material (Pa)
            cout << endl << "yield_strengths (Pa): ";
            cin >> yield_strength;            // Maximum tensile stress (Pa)
            cout << endl << "Safety Factor (typically 2-3) (2): ";
            cin >> safety_factor;             // Safety factor (typically 2-3)
            break;
        }
}

//This is for cableGeometry to check whether or not the cable has snapped at any point.
bool isBroken(int num_segments) {

    for (int i = 0; i < num_segments; i++) {
        double stress_at_point = cableWeight[i] / cableWidth[i];

        //track max stress
        if (stress_at_point > max_stress) {
            max_stress = stress_at_point;
        }

        //check if stress exceeds material limit
        if (stress_at_point > matProp.sigma_max) {
            failure_height = i;
            return true;
        }

    }
    return false;
}

//This is to calculate location effects
void locationEffect(int i) {
    double local_g;
    double coriolis_p;
    double local_omega;
    double velo_0;
    double toRAD = pi / 180.0;
    double gravity_variation;

    int i = sectProp.sectNum; //section number to calculate at - change as needed
    //calculate sin of latitude in radians
    double sinLatitude = sin(locProp.latitude * toRAD);
    double cosLatitude = cos(locProp.latitude * toRAD);

    //calculate local g based on latitude
    local_g = orbitalMechanics(sectionList[i].radius);
    locationEffectsList[i].gravity_effect += local_g;

    //local coriolis parameter
    coriolis_p = 2 * omega * sinLatitude;
	locationEffectsList[i].coriolis_parameter = coriolis_p;

    //ground velocity at latitude
    local_omega = pow(omega, 2) * sectionList[i].radius * cosLatitude;

	locationEffectsList[i].centrifugal_effect = local_omega;

	gravity_variation = g_0 * (1 + 0.0053024 * pow(sinLatitude, 2) - 0.0000058 * pow(sin(2 * locProp.latitude * toRAD), 2));
   
}

//This is to get the location
void locationChoice() {
    cout << "Specific location (Y) or compare for best results (N)? (y/N): ";
    string locChoice;
    cin.ignore();
    getline(cin, locChoice);
    if (locChoice == "") locChoice = "N"; //default to N if empty
    if (locChoice == "y" || locChoice == "Y") {
        cout << "Enter latitude (-90 to 90): ";
        double latitude;
        double longitude;
        cin >> latitude;
        cout << "Enter longitude (-180 to 180) (0 by default): ";
        cin >> longitude;

        //put it in the locProp struct
        locProp.latitude = latitude;
        locProp.longitude = longitude;
    }
    else { // the illusion of choice lmao
        cout << "Did you want to compare specific locations or use the defualt of: " << endl;
        cout << "1. Equator (0deg)" << endl;
        cout << "2. 28.5N (Cape Canaveral)" << endl;
        cout << "3. 28.5N (Kennedy Space Center)" << endl;
        cout << "4. 34N (Vandenberg Space Force Base)" << endl;
        cout << "5. 45.6N (Baikonur Cosmodrome)" << endl;
        cout << "6. 5.2N (Guiana Space Centre)" << endl;
        cout << "Or you can choose to enter your own latitudes as well." << endl;
        cout << "Either enter 'default' or 'custom': ";
        string locType;
        getline(cin, locType);

        if (locType == "default" || locType == "D" || locType == "d") {
            //location list has name, latitude, longitude, elevation, atmos_density, wind_speed, seismic activity
            locationList.push_back({ "Equator", 0.0, 0.0, 0.0, 1.225, 5.0, 0.0 });
            locationList.push_back({ "Cape Canaveral", 28.5, -80.6, 3.0, 1.225, 10.0, 0.0 });
            locationList.push_back({ "Kennedy Space Center", 28.5, -80.6, 3.0, 1.225, 10.0, 0.0 });
            locationList.push_back({ "Vandenberg Space Force Base", 34.7, -120.6, 112.0, 1.225, 8.0, 1.5 });
            locationList.push_back({ "Baikonur Cosmodrome", 45.6, 63.3, 90.0, 1.225, 7.0, 0.5 });
            locationList.push_back({ "Guiana Space Centre", 5.2, -52.8, 10.0, 1.225, 6.0, 0.0 });

            cout << "Using default locations." << endl;
        }
        else { //choose your own adventure style
            cout << "Enter latitudes separated by spaces (end with -999): ";
            double lat;
            while (true) {
                cout << "Latitude: ";
                cin >> lat;
                if (lat == -999) break;
                latitudeList.push_back(lat);
            }
            cout << "Would you like to assign Longitudes as well? (0 by default) (y/N): ";
            string lonChoice;
            getline(cin, lonChoice);
            if (lonChoice == "y" || lonChoice == "Y") {
                double longitude;
                for (size_t i = 0; i < latitudeList.size(); i++) {
                    cout << "Enter longitude for latitude " << latitudeList[i] << ": ";
                    cin >> longitude;
                    longitudeList.push_back(longitude);
                }
            }
            else { // defaults to 0 if not specified
                for (size_t i = 0; i < latitudeList.size(); i++) {
                    longitudeList.push_back(0.0);
                }

            }
            cout << "Would you like to assign any other variables as well? elevation, atmos_density, wind_speed, seismic activity (0 by default) (y/N): ";
            string otherChoice;
            getline(cin, otherChoice);
            if (otherChoice == "y" || otherChoice == "Y") {
                double elevation;
                double atmos_density;
                double wind_speed;
                double seismic_activity;
                cout << "Would you like to assign elevation? (y/N): ";
                string elevChoice;
                getline(cin, elevChoice);
                if (elevChoice == "y" || elevChoice == "Y") {
                    for (size_t i = 0; i < latitudeList.size(); i++) {
                        cout << "Enter elevation for latitude " << latitudeList[i] << ": ";
                        cin >> elevation;
                        elevationList.push_back(elevation);
                    }
                }
                cout << "Would you like to assign atmospheric density? (y/N): ";
                string atmosChoice;
                getline(cin, atmosChoice);
                if (atmosChoice == "y" || atmosChoice == "Y") {
                    for (size_t i = 0; i < latitudeList.size(); i++) {
                        cout << "Enter atmospheric density for latitude " << latitudeList[i] << ": ";
                        cin >> atmos_density;
                        atmos_densityList.push_back(atmos_density);
                    }
                }
                cout << "Would you like to assign wind speed? (y/N): ";
                string windChoice;
                getline(cin, windChoice);
                if (windChoice == "y" || windChoice == "Y") {
                    for (size_t i = 0; i < latitudeList.size(); i++) {
                        cout << "Enter wind speed for latitude " << latitudeList[i] << ": ";
                        cin >> wind_speed;
                        wind_speedList.push_back(wind_speed);
                    }
                }
                cout << "Would you like to assign seismic activity? (y/N): ";
                string seismicChoice;
                getline(cin, seismicChoice);
                if (seismicChoice == "y" || seismicChoice == "Y") {
                    for (size_t i = 0; i < latitudeList.size(); i++) {
                        cout << "Enter seismic activity for latitude " << latitudeList[i] << ": ";
                        cin >> seismic_activity;
                        seismic_activityList.push_back(seismic_activity);
                    }
                }
            }
            else { // defaults to 0 (atmos is 1.225) ifnot specified
                for (size_t i = 0; i < latitudeList.size(); i++) {
                    longitudeList.push_back(0.0);
                    elevationList.push_back(0.0);
                    atmos_densityList.push_back(1.225);
                    wind_speedList.push_back(0.0);
                    seismic_activityList.push_back(0.0);
                }
                for (size_t i = 0; i < latitudeList.size(); i++) {
                    locationList.push_back({ "Custom Location", latitudeList[i], longitudeList[i], elevationList[i], atmos_densityList[i], wind_speedList[i], seismic_activityList[i] });
                }
            }
        }
    }
    cout << "====== Location Summary ======" << endl;
    if (debug == 1) {
        for (size_t i = 0; i < locationList.size(); i++) {
            cout << "Location: " << i + 1 << endl << "Name: " << locationList[i].name << endl << "Latitude: " << locationList[i].latitude << endl << "Longitude: " << locationList[i].longitude << endl << "Atmospheric_density: " << locationList[i].atmos_density << endl << "Elevation: " << locationList[i].elevation << endl << "Wind Speed: " << locationList[i].wind_speed << endl << endl;
        }
    }
}

//climberGeometry is used for calculating and storing the data of the cross-section, weight, length of the climber system
void climberGeometry() {
}

//climberDynamics is used for calculating and storing the data of the forces, power, velocity, acceleration of the climber system
void climberDynamics() {
}

//climberEffects is used for calculating and storing the data of the effects of the climber system on the cable
void climberEffects() {
}

//counterweightGeometry is used for calculating and storing the data of the cross-section, weight, length of the counterweight system
void counterweightGeometry() {
}

//counterweightDynamics is used for calculating and storing the data of the forces, power, velocity, acceleration of the counterweight system
void counterweightDynamics() {
}
//counterweightEffects is used for calculating and storing the data of the effects of the counterweight system on the cable
void counterweightEffects() {
}
//cableGeometry is used for calculating and storing the data of the cross-section, weight, length of the cable system
void cableGeometry() {
    //Init variables

    bool cable_failed = false;
    double g_eff_local;



    //cable extends 50% beyond GEO for counterweight
    r_top = r_GEO * 1.5;

    total_cable_length = r_top - r_0;

    cout << endl << "Total cable length: " << total_cable_length / 1000 << " km" << endl;

    double dr = 1000; // 1 km segments
    int num_segments = total_cable_length / dr;

    cout << "Number of segments " << num_segments << endl;

    //resize the vectors to hold all results
    cableLength.resize(num_segments);
    cableWidth.resize(num_segments);
    cableWeight.resize(num_segments);
    cableSegmentWeight.resize(num_segments);
    locationEffectsList.resize(num_segments);

    //set A_TOP - cross sectional area of the top of the cable
    double A_top = 1e-4; // 1cm^2
    cout << "Initial area at top (A_top: " << A_top * 10000 << " cm2" << endl;

    //pre-calculate cumulative g_eff for optimization
    //creates an array where cumulative_g_eff[i] = integram from 0 to i

    vector<double> cumulative_g_eff(num_segments + 1);
    cumulative_g_eff[0] = 0; //base case

    for (int i = 0; i < num_segments; i++) {
        double r_i = r_0 + (i * dr);
        double g_eff_i = (G * M / pow(r_i, 2)) - pow(omega, 2) * r_i;
        cumulative_g_eff[i + 1] = cumulative_g_eff[i] + g_eff_i * dr;
    }

    cout << "Starting the heavy lifting..." << endl;

    //main calculation loop (top to bottom)
    //reset cumulative tension - starting at 0 at the top
    cumulative_tension = 0;
    bool use_location_effects = !locationList.empty() || (locProp.latitude != 0);
    //loop from top to bottom
    for (int i = num_segments - 1; i >= 0; i--) {
        //effective radius from earth's center
        double current_r = r_0 + (i * dr);
		locationEffect(int i); // calculate location effects for this section
        //effective gravity at this point (negative beyond GEO)
        if (use_location_effects) {
            locationEffect(); // calculate location effects for this section
            g_eff_local = (G * M / pow(current_r, 2)) - pow(omega, 2) * current_r + locationEffectsList[i].gravity_effect;
        }
        else {
            g_eff_local = (G * M / pow(current_r, 2)) - pow(omega, 2) * current_r;
        }
        //simple area from tension
        double required_area_simple;
        if (cumulative_tension > 0) {
            required_area_simple = cumulative_tension / matProp.sigma_max;
        }
        else {
            required_area_simple = A_top;
        }

        //exponential taper formula
        //get integral from current point i to TOP
        //uses pre-calculated cumulative array for 0(1) lookup
        double g_eff_integral = cumulative_g_eff[num_segments] - cumulative_g_eff[i];

        //apply exponential taper formula
        //accounts for the accumulated gravitational stress from here to the top
        double taper_factor = exp((material_density * abs(g_eff_integral)) / matProp.sigma_max);
        double tapered_area = A_top * taper_factor;

        //choose the safer area
        double final_area = max(required_area_simple, tapered_area);

        //ensure minimum area constraint
        if (final_area < A_top) {
            final_area = A_top;
        }
        //calculate the segment weight
        double segment_weight = material_density * final_area * dr * abs(g_eff_local);

        //update cumulative tension
        cumulative_tension += segment_weight;

        //store results
        cableLength[i] = current_r - r_0;
        cableWidth[i] = final_area;
        cableWeight[i] = cumulative_tension;
        cableSegmentWeight[i] = segment_weight;

        //store section properties
        sectionList.resize(num_segments); // ensure the sectionList is the correct size
        sectionList[i].sectNum = i;
        sectionList[i].radius = current_r;
        sectionList[i].length = dr;
        sectionList[i].cross_area = final_area;

        if (debug == 1) {
            //progress indicator
            if (i % 1000 == 0) {
                cout << "Height: " << (current_r - r_0) / 1000 << " km, " << endl
                    << "Area: " << final_area * 10000 << " cm^2, " << endl
                    << "Tension: " << cumulative_tension / 1e6 << " MN" << endl << endl;
            }
        }
    }
    //check for failure
    cout << "Checking for cable integrity..." << endl;
    cable_failed = isBroken(num_segments);
    //report results
    if (cable_failed) {
        cout << endl << "====== Cable Integrity Check ======" << endl;
        cout << endl << "Cable Failure" << endl;
        cout << "Failure at height: " << cableLength[failure_height] / 1000 << " km" << endl;
        cout << "Stress at failure: " << (cableWeight[failure_height] / cableWidth[failure_height]) / 1e9 << " GPa" << endl;
        cout << "Material Limit: " << matProp.sigma_max / 1e9 << " GPa" << endl;
    }
    else {
        cout << endl << "====== Cable Integrity Check ======" << endl;
        cout << endl << "SUCCESS: The cable STANDS STRONG" << endl;
        cout << "Maximum stress: " << max_stress / 1e9 << " GPa" << endl;
        cout << "Material limit: " << matProp.sigma_max / 1e9 << " GPa" << endl;
        cout << "Safety margin: " << ((matProp.sigma_max - max_stress) / matProp.sigma_max * 100) << "%" << endl;

    }

    //additional stats
    double total_mass = 0;
    double max_area = 0;
    double min_area = 1e10; // large to start

    for (int i = 0; i < num_segments; i++) {
        //mass of segment
        total_mass += material_density * cableWidth[i] * dr;

        //track min/max area
        if (cableWidth[i] > max_area) max_area = cableWidth[i];
        if (cableWidth[i] < min_area) min_area = cableWidth[i];
    }



    cout << endl << "====== Cable Stats ======" << endl;
    cout << "Total cable mass: " << total_mass / 1e6 << " metric tons" << endl;
    cout << "Minimum area (at top): " << min_area * 10000 << " cm^2" << endl;
    cout << "Maximum area (near Earth): " << max_area * 10000 << " cm^2" << endl;
    cout << "Taper ratio: " << max_area / min_area << endl;
    cout << "Cable extends from " << r_0 / 1000 << " km to " << r_top / 1000 << " km from Earth's center" << endl;
    cout << "Done" << endl;

}

//elevatorDynamics is used for calculating and storing the data of the forces, power, velocity, acceleration of the full elevator system
void elevatorDynamics() {
}

//elevatorEffects is used for calculating and storing the data of the effects of the full elevator system on the cable
void elevatorEffects() {
}
//elevatorGeometry is used for calculating and storing the data of the cross-section, weight, length of the full elevator system
void elevatorGeometry() {
}

//This is to calculate price per kg
void pricePerKg() {
    if (matProp.name == "Steel") {
        matProp.price_per_kg = 0.5; // $0.5 per kg
    }
    else if (matProp.name == "Copper") {
        matProp.price_per_kg = 10.0; // $10.0 per kg
    }
    else if (matProp.name == "Carbon Nanotubes") {
        matProp.price_per_kg = 400.0; // $400 per kg (current estimate, can vary widely)
    }
    else if (matProp.name == "Titanium") {
        matProp.price_per_kg = 15.0; // $15.0 per kg
    }
    else {
        matProp.price_per_kg = 20.0; // default for other materials
    }
	cout << endl << "====== Material Cost Analysis ======" << endl;
	cout << endl << "Price per kg of " << matProp.name << ": $" << matProp.price_per_kg << " per kg" << endl;
    //calculate total cable cost
    double total_mass = 0;
    for (double segment_weight : cableSegmentWeight) {
        total_mass += segment_weight / g_0; // convert weight (N) to mass (kg)
    }
    double total_cost = total_mass * matProp.price_per_kg;
	cout << "Estimated total cable cost: $" << total_cost / 1e6 << " million" << endl;
    matProp.totalCost = total_cost;
}




//This is for the main function
int main(){
    if (debug ==2) cout << endl << "calling materialDefinition()" << endl;
    materialDefinition();
    if (debug == 2) cout << endl << "calling locationChoice()" << endl;
    locationChoice();
    if (debug == 2) cout << endl << "calling locationEffects()" << endl;
    locationEffect();
    if (debug == 2) cout << endl << "calling orbitalMechanics()" << endl;
    orbitalMechanics();
    if (debug == 2) cout << endl << "calling cableGeometry()" << endl;
    cableGeometry();
    if (debug == 2) cout << endl << "calling pricePerKg()" << endl;
    pricePerKg();



    
    return 0;
}

//the data flow should be as follows:
//user defines materials OR chooses to define sections -> user chooses location or compares multiple
//all of the geometry and dynamics functions are called to calculated the effects of each part of the system
//all the effects are calculated and applied to the sectionEffects struct
//elevator properties are calculated and effects are applied to simulate stress, tension, and other factors on the cable

