



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
}matProp;
// ============================================
// Defining the arrays
// ============================================
vector<double> cableLength;
vector<double> cableWidth;
vector<double> cableWeight;
vector<double> cableSegmentWeight;


// ============================================
// Settings
// ============================================
// This is to define any settings like DEBUG that helps me understand where i messed up
bool debug = 0;

// ============================================
// Actual mfn code
// ============================================



//calculate all of the base orbital mechanics
void orbitalMechanics() {
    r_GEO = pow((G * M * pow(T, 2) / (4 * pow(pi, 2))), (1.0 / 3.0)); // calculate the radius from center of earth to Geostationary orbit

    r = r_GEO; //to calculate valid effective gravity at point r
    g_eff = (G * (M / pow(r, 2))) - (pow(omega, 2) * r); // calculating gravity is (earths gravity - centrifugal force) should be close to 0 at r_GEO and 9.8 at r_0

    cout << "r_GEO: " << r_GEO << endl;
    cout << "g_eff: " << g_eff;
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
    case 'C':
        char copperOption;
        cout << "You have chosen Copper." << endl;
        cout << "Is the copper cold-worked? (y/N): ";
        cin >> copperOption;
        material_density = 8960.0;                 // Material density (kg/m³)
        safety_factor = 2.0;                       // Safety factor (typically 2-3)
        if (toupper(copperOption) == 'Y') {
            tensile_strength = 3.5e8;              // Tensile strength of material (Pa)
            yield_strength = 3.2e9;                // yield strength (Pa)
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
    case 'N':
        cout << "You have chosen Carbon Nanotubes." << endl;
        material_density = 1400;                 // Material density (kg/m³)
        tensile_strength = 1.5e11;               // Tensile strength of material (Pa)
        yield_strength = 1.0e11;                 // yield strength (Pa)
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

//cableGeometry is used for calculating and storing the data of the cross-section, weight, length of the elevator
void cableGeometry() {
    //Init variables

    bool cable_failed = false;
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

    //loop from top to bottom
    for (int i = num_segments - 1; i >= 0; i--) {
        //effective radius from earth's center
        double current_r = r_0 + (i * dr);

        //effective gravity at this point (negative beyond GEO)
        double g_eff_local = (G * M / pow(current_r, 2)) - pow(omega, 2) * current_r;

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
        cout << endl << "Cable Failure" << endl;
        cout << "Failure at height: " << cableLength[failure_height] / 1000 << " km" << endl;
        cout << "Stress at failure: " << (cableWeight[failure_height] / cableWidth[failure_height]) / 1e9 << " GPa" << endl;
        cout << "Material Limit: " << matProp.sigma_max / 1e9 << " GPa" << endl;
    }
    else {
        cout << endl << "SUCCESS: The cable STANDS STRONG" << endl;
        cout << "Maximum stress: " << max_stress / 1e9 << " GPa" << endl;
        cout << " Material limit: " << matProp.sigma_max / 1e9 << " GPa" << endl;
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



int main()
{
    if (debug == 1) cout << endl << "calling materialDefinition()" << endl;
    materialDefinition();
    if (debug == 1) cout << endl << "calling orbitalMechanics()" << endl;
    orbitalMechanics();
    if (debug == 1) cout << endl << "calling cableGeometry()" << endl;
    cableGeometry();

    return 0;
}

