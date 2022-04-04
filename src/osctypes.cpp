struct position {
    public:
    double x = 0;
    double y = 0;
    double z = 0;
};

struct orbparam {
    public:
    long int sma; // semi major axis (m)
    double ecc; // eccentricity
    double mna; // mean anomaly
    double aop; // argument of periapsis (rad)
    double inc; // inclination (rad)
    double asc; // longitude of the ascending node (rad)
};

struct orbrot {
    // Orbital reference rotation is prograde positive x and anti-normal positive z, this is the offset from that rotation
    public:
    double q;
    double x;
    double y;
    double z;
};

struct powermodel {
    // Using vec indices to indicate state, 0 for low power/off, 1 for idle, 2 for in use/max load, >2 custom. Usage in W
    std::vector<double> pstates; 
};