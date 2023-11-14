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
        KL = 132.5  # 135
        KL0 = 0.0

        KD = 25  # 45
        KD0 = 2.15  # 2

        K_beta = 20

        KM = -100  # -50
        KM0 = 0.0

        K_MY = 100
        K_MR = -60

        KOmega11 = -20  # Kq1
        KOmega12 = -60  # Kq2
        KOmega13 = -20  # Kq3
        KOmega21 = 0
        KOmega22 = 0
        KOmega23 = 0

    class VARIABLES:
        GLIDE_ANGLE = 25  # degrees
        SPEED = 0.8  # speed of glider in vertical plane in m/s
        BALLAST_RATE = 0.001  # Ballast rate input in kg/s

        # Position of internal mass mm
        rp1 = None  # This will vary
        rp2 = 0.02  # fixed, in meters. Change this to change roll (banked turns)
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

        PHI = 45  # roll angle in degrees
        THETA = 25  # pitch angle in degrees
        PSI = 0.0  # yaw angle in degrees

        BETA = 1.0  # sideslip angle in degrees
        
        RUDDER = 10.0 # rudder angle

    class CONTROLS:
        wp1 = 0.001  # m/s2
        wp2 = 0.0
        wp3 = 0.0

        wb1 = 0.0
        wb2 = 0.0
        wb3 = 0.0

        ww1 = 0.0
        ww2 = 0.0
        ww3 = 0.0
