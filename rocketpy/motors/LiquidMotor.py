# -*- coding: utf-8 -*-

__author__ = "Giovani Hidalgo Ceotto, Oscar Mauricio Prada Ramirez, João Lemes Gribel Soares, Mateus Stano Junqueira, Pedro Henrique Marinho Bressan, Patrick Bales, Lakshman Peri, Gautam Yarramreddy, Curtis Hu, and William Bradford"
__copyright__ = "Copyright 20XX, RocketPy Team"
__license__ = "MIT"

try:
    from functools import cached_property
except ImportError:
    from rocketpy.tools import cached_property

from rocketpy.Function import Function, funcify_method
from .Motor import Motor


class LiquidMotor(Motor):
    """Class to specify characteristics and useful operations for Liquid
    motors. This class inherits from the Motor class.
    """

    def __init__(
        self,
        thrustSource,
        dry_mass,
        center_of_dry_mass,
        dry_inertia,
        nozzleRadius,
        burn_time=None,
        nozzlePosition=0,
        reshapeThrustCurve=False,
        interpolationMethod="linear",
        coordinateSystemOrientation="nozzleToCombustionChamber",
    ):
        """Initialize LiquidMotor class, process thrust curve and geometrical
        parameters and store results.

        Parameters
        ----------
        thrustSource : int, float, callable, string, array
            Motor's thrust curve. Can be given as an int or float, in which
            case the thrust will be considered constant in time. It can
            also be given as a callable function, whose argument is time in
            seconds and returns the thrust supplied by the motor in the
            instant. If a string is given, it must point to a .csv or .eng file.
            The .csv file shall contain no headers and the first column must
            specify time in seconds, while the second column specifies thrust.
            Arrays may also be specified, following rules set by the class
            Function. See help(Function). Thrust units are Newtons.
        burn_time: float, tuple of float, optional
            Motor's burn time.
            If a float is given, the burn time is assumed to be between 0 and the
            given float, in seconds.
            If a tuple of float is given, the burn time is assumed to be between
            the first and second elements of the tuple, in seconds.
            If not specified, automatically sourced as the range between the first- and
            last-time step of the motor's thrust curve. This can only be used if the
            motor's thrust is defined by a list of points, such as a .csv file, a .eng
            file or a Function instance whose source is a list.
        dry_mass : int, float
            The total mass of the motor structure, including chambers
            and tanks, when it is empty and does not contain any propellant.
        center_of_dry_mass : int, float
            The position, in meters, of the motor's center of mass with respect
            to the motor's coordinate system when it is devoid of propellant.
            See `Motor.coordinateSystemOrientation`.
        dry_inertia : tuple, list
            Tuple or list containing the motor's dry mass inertia tensor
            components, in kg*m^2. This inertia is defined with respect to the
            the `center_of_dry_mass` position.
            Assuming e_3 is the rocket's axis of symmetry, e_1 and e_2 are
            orthogonal and form a plane perpendicular to e_3, the dry mass
            inertia tensor components must be given in the following order:
            (I_11, I_22, I_33, I_12, I_13, I_23), where I_ij is the
            component of the inertia tensor in the direction of e_i x e_j.
            Alternatively, the inertia tensor can be given as (I_11, I_22, I_33),
            where I_12 = I_13 = I_23 = 0.
        nozzleRadius : int, float
            Motor's nozzle outlet radius in meters.
        nozzlePosition : float
            Motor's nozzle outlet position in meters, specified in the motor's
            coordinate system. See `Motor.coordinateSystemOrientation` for
            more information.
        reshapeThrustCurve : boolean, tuple, optional
            If False, the original thrust curve supplied is not altered. If a
            tuple is given, whose first parameter is a new burn out time and
            whose second parameter is a new total impulse in Ns, the thrust
            curve is reshaped to match the new specifications. May be useful
            for motors whose thrust curve shape is expected to remain similar
            in case the impulse and burn time varies slightly. Default is
            False.
        interpolationMethod : string, optional
            Method of interpolation to be used in case thrust curve is given
            by data set in .csv or .eng, or as an array. Options are 'spline'
            'akima' and 'linear'. Default is "linear".
        coordinateSystemOrientation : string, optional
            Orientation of the motor's coordinate system. The coordinate system
            is defined by the motor's axis of symmetry. The origin of the
            coordinate system  may be placed anywhere along such axis, such as at the
            nozzle area, and must be kept the same for all other positions specified.
            Options are "nozzleToCombustionChamber" and "combustionChamberToNozzle".
            Default is "nozzleToCombustionChamber".
        """
        super().__init__(
            thrustSource,
            dry_mass,
            center_of_dry_mass,
            dry_inertia,
            nozzleRadius,
            burn_time,
            nozzlePosition,
            reshapeThrustCurve,
            interpolationMethod,
            coordinateSystemOrientation,
        )

        self.positioned_tanks = []

    @funcify_method("Time (s)", "Exhaust Velocity (m/s)")
    def exhaustVelocity(self):
        """Computes the exhaust velocity of the motor from its mass flow
        rate and thrust.

        Parameters
        ----------
        t : float
            Time in seconds.

        Returns
        -------
        self.exhaustVelocity : Function
            Gas exhaust velocity of the motor.
        """
        return self.thrust / (-1 * self.massFlowRate)

    @funcify_method("Time (s)", "Propellant Mass (kg)")
    def propellantMass(self):
        """Evaluates the mass of the motor as the sum of each tank mass.

        Parameters
        ----------
        t : float
            Time in seconds.

        Returns
        -------
        Function
            Mass of the motor, in kg.
        """
        propellantMass = 0

        for positioned_tank in self.positioned_tanks:
            propellantMass += positioned_tank.get("tank").mass

        return propellantMass

    @cached_property
    def propellantInitialMass(self):
        """Property to store the initial mass of the propellant.

        Returns
        -------
        float
            Initial mass of the propellant, in kg.
        """
        return self.propellantMass(0)

    @funcify_method("Time (s)", "mass flow rate (kg/s)", extrapolation="zero")
    def massFlowRate(self):
        """Evaluates the mass flow rate of the motor as the sum of each tank
        mass flow rate.

        Parameters
        ----------
        t : float
            Time in seconds.

        Returns
        -------
        Function
            Mass flow rate of the motor, in kg/s.

        See Also
        --------
        `Motor.totalMassFlowRate` :
            Calculates the total mass flow rate of the motor assuming
            constant exhaust velocity.
        """
        massFlowRate = Function(0)

        for positioned_tank in self.positioned_tanks:
            massFlowRate += positioned_tank.get("tank").netMassFlowRate

        return massFlowRate

    @funcify_method("Time (s)", "center of mass (m)")
    def centerOfPropellantMass(self):
        """Evaluates the center of mass of the motor from each tank center of
        mass and positioning. The center of mass height is measured relative to
        the motor nozzle.

        Parameters
        ----------
        t : float
            Time in seconds.

        Returns
        -------
        Function
            Center of mass of the motor, in meters.
        """
        totalMass = 0
        massBalance = 0

        for positioned_tank in self.positioned_tanks:
            tank = positioned_tank.get("tank")
            tankPosition = positioned_tank.get("position")
            totalMass += tank.mass
            massBalance += tank.mass * (tankPosition + tank.centerOfMass)

        return massBalance / totalMass

    @funcify_method("Time (s)", "Inertia I_11 (kg m²)")
    def propellant_I_11(self):
        """Inertia tensor 11 component of the propellant, the inertia is
        relative to the e_1 axis, centered at the instantaneous propellant
        center of mass.

        Parameters
        ----------
        t : float
            Time in seconds.

        Returns
        -------
        float
            Propellant inertia tensor 11 component at time t.

        Notes
        -----
        The e_1 direction is assumed to be the direction perpendicular to the
        motor body axis.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Moment_of_inertia#Inertia_tensor
        """
        I_11 = Function(0)
        centerOfMass = self.centerOfPropellantMass

        for positioned_tank in self.positioned_tanks:
            tank = positioned_tank.get("tank")
            tankPosition = positioned_tank.get("position")
            I_11 += (
                tank.inertiaTensor
                + tank.mass * (tankPosition + tank.centerOfMass - centerOfMass) ** 2
            )

        return I_11

    @funcify_method("Time (s)", "Inertia I_22 (kg m²)")
    def propellant_I_22(self):
        """Inertia tensor 22 component of the propellant, the inertia is
        relative to the e_2 axis, centered at the instantaneous propellant
        center of mass.

        Parameters
        ----------
        t : float
            Time in seconds.

        Returns
        -------
        float
            Propellant inertia tensor 22 component at time t.

        Notes
        -----
        The e_2 direction is assumed to be the direction perpendicular to the
        motor body axis, and perpendicular to e_1.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Moment_of_inertia#Inertia_tensor
        """
        return self.I_11

    @funcify_method("Time (s)", "Inertia I_33 (kg m²)")
    def propellant_I_33(self):
        """Inertia tensor 33 component of the propellant, the inertia is
        relative to the e_3 axis, centered at the instantaneous propellant
        center of mass.

        Parameters
        ----------
        t : float
            Time in seconds.

        Returns
        -------
        float
            Propellant inertia tensor 33 component at time t.

        Notes
        -----
        The e_3 direction is assumed to be the axial direction of the rocket
        motor.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Moment_of_inertia#Inertia_tensor
        """
        return 0

    @funcify_method("Time (s)", "Inertia I_12 (kg m²)")
    def propellant_I_12(self):
        return 0

    @funcify_method("Time (s)", "Inertia I_13 (kg m²)")
    def propellant_I_13(self):
        return 0

    @funcify_method("Time (s)", "Inertia I_23 (kg m²)")
    def propellant_I_23(self):
        return 0

    def addTank(self, tank, position):
        """Adds a tank to the rocket motor.

        Parameters
        ----------
        tank : Tank
            Tank object to be added to the rocket motor.
        position : float
            Position of the tank relative to the motor nozzle, in meters.
            The position is measured from the nozzle tip to the tank
            geometry reference zero point.
        """
        self.positioned_tanks.append({"tank": tank, "position": position})

    def allInfo(self):
        """Prints out all data and graphs available about the Motor.

        Parameters
        ----------
        None

        Return
        ------
        None
        """
        # Print nozzle details
        print("Nozzle Details")
        print("Nozzle Radius: " + str(self.nozzleRadius) + " m")

        # Print motor details
        print("\nMotor Details")
        print("Total Burning Time: " + str(self.burnDuration) + " s")
        print(
            "Total Propellant Mass: "
            + "{:.3f}".format(self.propellantInitialMass)
            + " kg"
        )
        print(
            "Average Propellant Exhaust Velocity: "
            + "{:.3f}".format(self.exhaustVelocity.average(*self.burn_time))
            + " m/s"
        )
        print("Average Thrust: " + "{:.3f}".format(self.averageThrust) + " N")
        print(
            "Maximum Thrust: "
            + str(self.maxThrust)
            + " N at "
            + str(self.maxThrustTime)
            + " s after ignition."
        )
        print("Total Impulse: " + "{:.3f}".format(self.totalImpulse) + " Ns")

        # Show plots
        print("\nPlots")
        self.thrust.plot(*self.burn_time)
        self.totalMass.plot(*self.burn_time)
        self.massFlowRate.plot(*self.burn_time)
        self.exhaustVelocity.plot(*self.burn_time)
        self.centerOfMass.plot(*self.burn_time)
        self.I_11.plot(*self.burn_time)
        self.I_22.plot(*self.burn_time)
        self.I_33.plot(*self.burn_time)
        self.I_12.plot(*self.burn_time)
        self.I_13.plot(*self.burn_time)
        self.I_23.plot(*self.burn_time)
