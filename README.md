# SET-Sim
## Space Elevator Tether Simulator

Space elevator cable stress calculator. 

## What This Is

C++ program that calculates whether a cable extending from Earth's surface to beyond geostationary orbit would break under its own weight. It usually does. Unless you use carbon nanotubes. Then it might survive.

## The Physics

The cable experiences varying forces along its length:
- Gravity decreases with altitude (inverse square law)
- Centrifugal force increases with altitude (Earth's rotation)
- At GEO (42,164 km): forces balance
- Below GEO: net downward force
- Above GEO: net upward force

Each segment must support all weight above it. The program calculates this cumulative tension and determines required cross-sectional area at each point.

## Current Implementation

- Numerical integration using 1km segments from surface to 1.5×GEO radius
- Exponential taper calculation based on material properties
- Prefix sum optimization for O(n) complexity instead of O(n²)
- Material database: steel, titanium, copper, carbon nanotubes, or custom input
- Failure analysis with stress/strain calculations

## Build

```bash
g++ -o orbital-bushido main.cpp -std=c++11
./orbital-bushido
```

## Output

The program will tell you:
- Maximum stress location and magnitude
- Required taper ratio (base area / top area)
- Total cable mass
- Failure point if applicable

## Known Limitations

- Assumes uniform material properties
- Ignores wind loading and orbital debris
- No thermal expansion calculations
- Assumes perfect manufacturing (no defects)
- Single cable design (no redundancy)

## Future Work

Planned features that haven't been implemented yet:

- **Latitude/longitude placement** - Ecuador vs other locations, accounting for local gravity variations
- **CSV import/export** - For parameter sweeps and result analysis
- **Materials database integration** - search and pull from material science databases like The Materials Project for data integrity
- **Climber dynamics** - Power requirements, Coriolis effects, transit time calculations
- **Weather modeling** - Wind loads at various altitudes, particularly jet stream
- **Oscillation modes** - Natural frequencies, damping requirements, resonance avoidance
- **Redundant cable designs** - Multiple smaller cables vs single large cable
- **Cost modeling** - Material costs, launch costs, construction logistics
- **Anchor station specs** - Foundation requirements, ocean platform vs land-based
- **Power beaming** - Laser power delivery to climbers, efficiency calculations

## Why This Exists

My dream is that a space elevator will one day be a reality and wanted to dive into current hurdles stopping it. The math speaks "not likely". This program shows why. Carbon nanotubes with 150 GPa tensile strength might work if we can manufacture them in 60,000 km lengths without defects.

## License

MIT. Do whatever you want with it.
