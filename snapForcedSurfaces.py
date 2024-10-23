from numpy import (cos, sin, zeros, zeros_like,
                   pi, vstack, linspace, array,
                   arctan, sqrt, nan)
from scipy.integrate import solve_bvp, quad


class snapFixedRadius:
    """This is a wrapper class for the solver of the system of differential
    equations that describe the shape of the interface of a 2D droplet sitting
    on a smooth sinusoidal surface of amplitude 'a'.

    The system of differential equations is
    \begin{align}
        \frac{\mathrm{d}\phi}{\mathrm{d}s} & = k0 + k1 x \\
        \frac{\mathrm{d}x}{\mathrm{d}s} & = \cos \phi \\
        \frac{\mathrm{d}y}{\mathrm{d}s} & = \sin \phi
    \end{align}
    The system is solved using the solve_bvp from the scipy.integrate package.
    As the variable $s$ ranges from 0 to $S$, and $S$ is uknown, a change of 
    variables is necessary, thus using $\zeta = s/S$ and now $\zeta \in [0, 1]$
    keeping $S$ as one of the parameters to be found by solve_bvp.


    Class members
    -------------

    a : float
        Amplitude of the sinusoidal surface 

    Rb : float
        Base radius of the droplet 

    X : float
        Centre of the droplet, i.e., mid-point between contact points.

    qe : float
        Microscopic contact angle in radians

    guess_vars : list
        guess values for the [arc length, base curvature, forcing constant]

    guess_state : 3 by N array of floats
        
    sol : bunch scipy object 
        Object containing the solution of the boundary value problem as given
        by scipy.integrate.solve_bvp 

    msg : str
        Message to print if the solution to the system is yet to be calculated 


    Class methods
    -------------
    
    update_values 
    update_guess_state
    update_guess_vars
    calculate_shape
    get_interface_arclen
    get_force
    evaluate

    """

    def __init__(
        self, 
        contact_angle=None,
        droplet_centre=None,
        droplet_radius=None,
        surface_amplitude=None
    ):
        """Quick way to set up some of the configuration parameters of a droplet 
        on a sinusoidal surface.


        Parameters 
        ----------        

        contact_angle : float, optional
            The microscopic contact angle in radians of the liquid-gas interface
            on both contact points.

        droplet_centre : float, optional
            The centre position of the droplet in dimensionless units, i.e., pro-
            portional to the wavelenth.

        droplet_radius : float, optional
            The horizontal distance of the contact points from the centre in
            dimensionless units.

        surface_amplitude : float, optional
            The amplitude of the sinusoidal surface in dimensionless units.


        Returns
        -------

        None
        """
        update = (not contact_angle is None and
                  not droplet_centre is None and
                  not droplet_radius is None and
                  not surface_amplitude is None)
        if (update):
            self.update_values(contact_angle,
                               droplet_centre,
                               droplet_radius,
                               surface_amplitude)
        
        self.msg = "Solution to the system differential equations not obtained."

        self.update_guess_state()

    def _h(x):
        return self.a*(1.0 - cos(2.0*pi*x))

    def _dh(x):
        return 2.0*pi*self.a*sin(2.0*pi*x)

    def update_values(self, contact_angle, droplet_centre, droplet_radius, surface_amplitude):
        """An alternative way to set up some of the configuration parameters.
        This can be executed once the object is already initialised.


        Parameters
        ----------

        contact_angle : float, optional
            The microscopic contact angle in radians of the liquid-gas interface
            on both contact points.

        droplet_centre : float, optional
            The centre position of the droplet in dimensionless units, i.e., pro-
            portional to the wavelenth.

        droplet_radius : float, optional
            The horizontal distance of the contact points from the centre in
            dimensionless units.
 
        surface_amplitude : float, optional
            The amplitude of the sinusoidal surface in dimensionless units.
           

        Returns
        -------

        None
        """

        self.qe = contact_angle
        self.X = droplet_centre
        self.Rb = droplet_radius
        self.a = surface_amplitude
        self.guess_vars = [pi*self.Rb, 1.0/self.Rb, 0.0]
        return None

    def update_guess_state(self, some_state=None):
        """Used to set up a guess for the solution to the system of 
        differential equations for the shape of the droplet.


        Parameters
        ----------
        
        some_state : list, optional
            ...


        Returns
        -------

        None

        """

        self.guess_state = some_state
        return None

    def update_guess_vars(self, arc_len, curvature, force):
        """The unknown parameters guessed for an initial value to solve 
        the system of equations by the solve_bvp algorithm.


        Parameters
        ----------

        arc_len : float 
            The 'S' parameter, i.e., the total arch length.

        curvature : float 
            The 'k0' parameter, which stands as a the curvature if no
            external force is applied.

        force : float
            The 'k1' parameter, which stands for a body force, similar
            to gravity, that pushes sideways the droplet.


        Returns
        -------

        None

        """
        self.guess_vars = [arc_len, curvature, force]
        return None

    def calculate_shape(self):
        """Routine to calculate the shape of the interface and the unknown 
        constants by solving the non-linear system of differential equations
        and finding the constants such that the boundary conditions are 
        satisfied. This is based on the package scipy.integrate.solve_bvp.


        Parameters
        ----------

        None


        Returns
        -------

        None

        """

        def fun(x, y, p):
            """
            variables
            x -> \zeta : s / S  (arc length from the right contact line
            divided by the total arc length)

            y[0] -> \phi : the angle of the unit tangent vector.
            y[1] -> x : x coordinate
            y[2] -> y : y coordinate

            p[0] -> S : total arc length
            p[1] -> k0 : base curvature
            p[2] -> k1 : forcing coefficient
            """

            x0 = self.X + self.Rb
            xf = self.X - self.Rb

            return p[0] * vstack([
                p[1] + p[2] * (y[1] - self.X),
                cos(y[0]),
                sin(y[0])
            ])

        def bc(ya, yb, p):
            """
            The boundary conditions for the system.

            Right contact point:
            ya[0] -> phi(0)
            ya[1] -> x(0)
            ya[2] -> y(0)

            Left contact line:
            yb[0] -> phi(1)
            yb[1] -> x(1)
            yb[2] -> y(1)

            The cross-sectional area is to be found once the solution is obtained.

            Parameters
            ----------

            ya : 3-element list of floats
                The boundary conditions for the right contact line 

            yb : 3-element list of floats 
                The boundary conditions for the left contact line 

            p : 3-element list of floats 
                The value of the a priori unknown parameters 


            Returns
            -------

            eqn : 6-element np array
                the reminder of the 3-left and 3-right boundary conditions

            """

            x0 = self.X + self.Rb
            xf = self.X - self.Rb

            return array([
                ya[0] - pi + self.qe - arctan(self._dh(x0)),
                ya[1] - x0,
                ya[2] - self._h(x0),

                yb[0] - pi - self.qe - arctan(self._dh(xf)),
                yb[1] - xf,
                yb[2] - self._h(xf)
            ])

        if not self.guess_state is None:
            x = self.guess_state[0]
            y = self.guess_state[1]
        else:
            x = linspace(0, 1, 32)
            y = zeros((4, x.size))

            y[0] = pi - self.qe + 2 * self.qe * x
            y[1] = self.X + self.Rb * x
            y[2] = self.h(self.X) + self.Rb * x * (1 - x)

        self.sol = solve_bvp(fun,
                             bc,
                             x,
                             y,
                             p=self.guess_vars,
                             max_nodes=2**12,
                             verbose=0,
                             tol=1e-6)
        if not self.sol.success:
            print("Unable to find the shape of the interface.")

        return self.sol.success

    def get_interface_arclen(self):
        """Gives the total liquid-gas interface arc length once the solution
        to the differential equations is found.

        Parameters
        ----------

        None 

        
        Returns
        -------

        S : float 
            The arc length.
        """

        if self.sol.success:
            return self.sol.p[0]
        else:
            raise ValueError(
                msg 
              + "\nUnable to calculate the interface arc length."
            )

        return nan

    def get_force(self):
        """The force calculated from the solution to the boundary
        value problem.


        Parameters
        ----------

        None


        Returns
        -------

        k1 : float 
            The coefficient of the forcing term.

        """

        if self.sol.success:
            return self.sol.p[2]
        else:
            raise ValueError(
                msg
              + "\nThe force obained is as good as the initial guess,"
                "giving NaN."
            )

        return nan

    def evaluate(self, N=128):
        """Evaluates the functions that are solution to the system of
        differential equations at N points.

        Parameters
        ----------

        N : int, optional
            The number of points to evaluate the functions. By default
            N = 128.


        Returns
        -------
           
        solution : tuple of (1darray, ndarray of 3 by N array of floats)
            The evaluated functions:
            (
                \zeta,
                solution[:, 0] -> phi(\zeta)
                solution[:, 1] -> x(\zeta)
                solution[:, 2] -> y(\zeta)
            )
        """

        if self.sol.success:
            x_plot = linspace(0, 1, N)
            y_plot = self.sol.sol(x_plot)
            return (x_plot, y_plot)
        else:
            raise ValueError(
                msg
              + "\nUnable to compute the functions solution to the system"
                "of differential equations."
            )

        return None


