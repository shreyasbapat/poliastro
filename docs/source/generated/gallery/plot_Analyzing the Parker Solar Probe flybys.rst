.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_generated_gallery_plot_Analyzing the Parker Solar Probe flybys.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_generated_gallery_plot_Analyzing the Parker Solar Probe flybys.py:


=======================================
Analyzing the Parker Solar Probe flybys
=======================================

 1. Modulus of the exit velocity, some features of Orbit #2
 ----------------------------------------------------------

First, using the data available in the reports, we try to compute some of the properties of orbit #2.
This is not enough to completely define the trajectory, but will give us information later on in the process.



.. code-block:: python


    from astropy import units as u

    T_ref = 150 * u.day
    print(T_ref)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    150.0 d



.. code-block:: python


    from poliastro.bodies import Earth, Sun, Venus

    k = Sun.k
    print(k)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    Name   = Heliocentric gravitational constant
      Value  = 1.32712442099e+20
      Uncertainty  = 10000000000.0
      Unit  = m3 / s2
      Reference = IAU 2009 system of astronomical constants



.. code-block:: python


    import numpy as np







.. math::
    T = 2 \pi \sqrt{\frac{a^3}{\mu}} \Rightarrow a = \sqrt[3]{\frac{\mu T^2}{4 \pi^2}}



.. code-block:: python


    a_ref = np.cbrt(k * T_ref**2 / (4 * np.pi**2)).to(u.km)
    print(a_ref.to(u.au))





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    0.5524952607456924 AU


.. math::
   \varepsilon = -\frac{\mu}{r} + \frac{v^2}{2} = -\frac{\mu}{2a} \Rightarrow v = +\sqrt{\frac{2\mu}{r} - \frac{\mu}{a}}



.. code-block:: python


    energy_ref = (-k / (2 * a_ref)).to(u.J / u.kg)
    print(energy_ref)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    -802837548.527052 J / kg



.. code-block:: python


    from poliastro.twobody import Orbit
    from poliastro.util import norm
    from astropy.time import Time

    flyby_1_time = Time("2018-09-28", scale="tdb")
    print(flyby_1_time)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    2018-09-28 00:00:00.000



.. code-block:: python


    r_mag_ref = norm(Orbit.from_body_ephem(Venus, epoch=flyby_1_time).r)
    print(r_mag_ref.to(u.au))





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    0.7257313162319988 AU



.. code-block:: python


    v_mag_ref = np.sqrt(2 * k / r_mag_ref - k / a_ref)
    print(v_mag_ref.to(u.km / u.s))





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    28.967363509562784 km / s


2. Lambert arc between #0 and #1
--------------------------------

To compute the arrival velocity to Venus at flyby #1, we have the necessary data to solve the boundary value problem.



.. code-block:: python


    d_launch = Time("2018-08-11", scale="tdb")
    print(d_launch)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    2018-08-11 00:00:00.000



.. code-block:: python


    ss0 = Orbit.from_body_ephem(Earth, d_launch)
    ss1 = Orbit.from_body_ephem(Venus, epoch=flyby_1_time)

    tof = flyby_1_time - d_launch

    from poliastro import iod

    (v0, v1_pre), = iod.lambert(Sun.k, ss0.r, ss1.r, tof.to(u.s))
    print(v0)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    [ 9.59933726 11.29855172  2.92449333] km / s



.. code-block:: python

    print(v1_pre)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    [-16.98082099  23.30752839   9.13129077] km / s



.. code-block:: python


    print(norm(v1_pre))





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    30.24846495185759 km / s


3. Flyby #1 around Venus
------------------------

We compute a flyby using poliastro with the default value of the entry angle,
just to discover that the results do not match what we expected.



.. code-block:: python


    from poliastro.threebody.flybys import compute_flyby

    V = Orbit.from_body_ephem(Venus, epoch=flyby_1_time).v
    print(V)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    [ 648499.73735241 2695078.44750227 1171563.7170508 ] km / d



.. code-block:: python


    h = 2548 * u.km

    d_flyby_1 = Venus.R + h
    print(d_flyby_1.to(u.km))





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    8599.8 km



.. code-block:: python


    V_2_v_, delta_ = compute_flyby(v1_pre, V, Venus.k, d_flyby_1)

    print(norm(V_2_v_))





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    27.75533877003213 km / s


4. Optimization
---------------

Now we will try to find the value of $\theta$ that satisfies our requirements.



.. code-block:: python


    def func(theta):
        V_2_v, _ = compute_flyby(v1_pre, V, Venus.k, d_flyby_1, theta * u.rad)
        ss_1 = Orbit.from_vectors(Sun, ss1.r, V_2_v, epoch=flyby_1_time)
        return (ss_1.period - T_ref).to(u.day).value







There are two solutions:



.. code-block:: python

    import matplotlib.pyplot as plt
    # %matplotlib inline

    theta_range = np.linspace(0, 2 * np.pi)
    plt.figure()
    plt.plot(theta_range, [func(theta) for theta in theta_range])
    plt.axhline(0, color='k', linestyle="dashed")
    plt.show()




.. image:: /generated/gallery/images/sphx_glr_plot_Analyzing the Parker Solar Probe flybys_001.png
    :class: sphx-glr-single-img





.. code-block:: python

    print(func(0))





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    -9.142672330001131



.. code-block:: python


    print(func(1))





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    7.09811543934556



.. code-block:: python


    from scipy.optimize import brentq

    theta_opt_a = brentq(func, 0, 1) * u.rad
    print(theta_opt_a.to(u.deg))





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    38.59870925415651 deg



.. code-block:: python


    theta_opt_b = brentq(func, 4, 5) * u.rad
    print(theta_opt_b.to(u.deg))





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    279.34770004025205 deg



.. code-block:: python


    V_2_v_a, delta_a = compute_flyby(v1_pre, V, Venus.k, d_flyby_1, theta_opt_a)
    V_2_v_b, delta_b = compute_flyby(v1_pre, V, Venus.k, d_flyby_1, theta_opt_b)

    print(norm(V_2_v_a))





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    28.96736350956281 km / s



.. code-block:: python


    print(norm(V_2_v_b))





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    28.96736350956282 km / s


5. Exit Orbit
---------------

And finally, we compute orbit #2 and check that the period is the expected one.



.. code-block:: python


    ss01 = Orbit.from_vectors(Sun, ss1.r, v1_pre, epoch=flyby_1_time)
    print(ss01)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    0 x 1 AU x 18.8 deg (HCRS) orbit around Sun (☉) at epoch 2018-09-28 00:00:00.000 (TDB)


The two solutions have different inclinations, so we still have to find out which is the good one.
We can do this by computing the inclination over the ecliptic - however, as the original data was
in the International Celestial Reference Frame (ICRF), whose fundamental plane is parallel to the
Earth equator of a reference epoch, we have change the plane to the Earth ecliptic,
which is what the original reports use.



.. code-block:: python


    ss_1_a = Orbit.from_vectors(Sun, ss1.r, V_2_v_a, epoch=flyby_1_time)
    print(ss_1_a)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    0 x 1 AU x 25.0 deg (HCRS) orbit around Sun (☉) at epoch 2018-09-28 00:00:00.000 (TDB)



.. code-block:: python


    ss_1_b = Orbit.from_vectors(Sun, ss1.r, V_2_v_b, epoch=flyby_1_time)
    print(ss_1_b)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    0 x 1 AU x 13.1 deg (HCRS) orbit around Sun (☉) at epoch 2018-09-28 00:00:00.000 (TDB)


Let's define a function to do that quickly for us, using the get_frame function from poliastro.frames:



.. code-block:: python


    from astropy.coordinates import CartesianRepresentation
    from poliastro.frames import Planes, get_frame

    def change_plane(ss_orig, plane):
        """Changes the plane of the Orbit.

        """
        ss_orig_rv = ss_orig.frame.realize_frame(
            ss_orig.represent_as(CartesianRepresentation)
        )

        dest_frame = get_frame(ss_orig.attractor, plane, obstime=ss_orig.epoch)

        ss_dest_rv = ss_orig_rv.transform_to(dest_frame)
        ss_dest_rv.representation_type = CartesianRepresentation

        ss_dest = Orbit.from_vectors(
            ss_orig.attractor,
            r=ss_dest_rv.data.xyz,
            v=ss_dest_rv.data.differentials['s'].d_xyz,
            epoch=ss_orig.epoch,
            plane=plane,
        )
        return ss_dest

    print(change_plane(ss_1_a, Planes.EARTH_ECLIPTIC))





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    0 x 1 AU x 3.5 deg (HeliocentricEclipticJ2000) orbit around Sun (☉) at epoch 2018-09-28 00:00:00.000 (TDB)



.. code-block:: python


    print(change_plane(ss_1_b, Planes.EARTH_ECLIPTIC))





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    0 x 1 AU x 13.1 deg (HeliocentricEclipticJ2000) orbit around Sun (☉) at epoch 2018-09-28 00:00:00.000 (TDB)


Therefore, the correct option is the first one.



.. code-block:: python


    print(ss_1_a.period.to(u.day))





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    149.99999999999991 d



.. code-block:: python


    print(ss_1_a.a)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    82652114.57939689 km


And, finally, we plot the solution:



.. code-block:: python


    from poliastro.plotting import OrbitPlotter
    plt.figure()
    frame = OrbitPlotter()

    frame.plot(ss0, label=Earth)
    frame.plot(ss1, label=Venus)
    frame.plot(ss01, label="#0 to #1")
    frame.plot(ss_1_a, label="#1 to #2");



.. rst-class:: sphx-glr-horizontal


    *

      .. image:: /generated/gallery/images/sphx_glr_plot_Analyzing the Parker Solar Probe flybys_002.png
            :class: sphx-glr-multi-img

    *

      .. image:: /generated/gallery/images/sphx_glr_plot_Analyzing the Parker Solar Probe flybys_003.png
            :class: sphx-glr-multi-img




**Total running time of the script:** ( 0 minutes  7.887 seconds)


.. _sphx_glr_download_generated_gallery_plot_Analyzing the Parker Solar Probe flybys.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: plot_Analyzing the Parker Solar Probe flybys.py <plot_Analyzing the Parker Solar Probe flybys.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: plot_Analyzing the Parker Solar Probe flybys.ipynb <plot_Analyzing the Parker Solar Probe flybys.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.readthedocs.io>`_
