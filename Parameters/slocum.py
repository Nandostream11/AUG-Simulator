class SLOCUM_PARAMS:
    class GLIDER_CONFIG:
        HULL_MASS = 40.0
        FIXED_POINT_MASS = 0.0
        BALLAST_MASS = 1.0
        INT_MOVABLE_MASS = 9.0

        FLUID_DISP_MASS = 50.0

        VEH_DENSITY = 0.0  # change

        BODY_LEN = 1.5  # meters
        RADIUS = 0.11  # meters

        # Added mass terms
        MF1 = 5  # kg
        MF2 = 60  # kg
        MF3 = 70  # kg

        # Added inertia terms
        J1 = 4  # kg.m2
        J2 = 12  # kg.m2
        J3 = 11  # kg.m2

    class HYDRODYNAMICS:
        KL = 132.5
        KL0 = 0.0

        KD = 25
        KD0 = 2.15

        KM = -100
        KM0 = 0.0

        KOmega1 = 50
        KOmega2 = 50

    class VARIABLES:
        GLIDE_ANGLE = 25  # degrees
        SPEED = 0.3  # speed of glider in vertical plane in m/s
        BALLAST_RATE = 0.025  # Ballast rate input in kg/s

        # Position of internal mass mm
        rp1 = None  # This will vary
        rp2 = 0.0
        rp3 = 0.05  # fixed, in meters

        # Position of ballast mass mb
        # Fixed, assumes ballast mass does not move
        rb1 = 0.0
        rb2 = 0.0
        rb3 = 0.0

        # Position of fixed point mass mw
        rw1 = 0.0
        rw2 = 0.0
        rw3 = 0.0

        PHI = 0.0  # degrees
        THETA = 45.0  # degrees
        PSI = 0.0  # degrees

    class CONTROLS:
        wp1 = 0.01  # 0.005
        wp2 = 0.0
        wp3 = 0.0

        wb1 = 0.0
        wb2 = 0.0
        wb3 = 0.0

        ww1 = 0.0
        ww2 = 0.0
        ww3 = 0.0
