poliastro - Astrodynamics in Python
===================================

.. image:: _static/logo_text.png
   :width: 80%
   :align: center

**poliastro** is an open source (MIT) collection of
Python functions useful in Astrodynamics and Orbital
Mechanics, focusing on interplanetary applications.
It provides a simple and intuitive API and handles
physical quantities with units. Some of its awesome
features are:

* Analytical and numerical orbit propagation
* Conversion between position and velocity vectors and
  classical orbital elements
* Coordinate frame transformations
* Hohmann and bielliptic maneuvers computation
* Trajectory plotting
* Initial orbit determination (Lambert problem)
* Planetary ephemerides (using SPICE kernels via Astropy)
* Computation of Near-Earth Objects (NEOs)

And more to come!

poliastro is developed by an open, international community.
Release announcements and general discussion take place
on our `mailing list`_ and `chat`_.

.. _`mailing list`: https://groups.io/g/poliastro-dev
.. _`chat`: https://riot.im/app/#/room/#poliastro:matrix.org

.. include:: form.rst

.. figure:: _static/molniya.png
   :align: center
   :figwidth: 60%
   :alt: Molniya orbit

   Plot of a `Molniya orbit`_ around the Earth
   (\\(a = 26600\\,\\mathrm{km}, e = 0.75,
   i = 63.4 \\mathrm{{}^{\\circ}} \\)).

The `source code`_, `issue tracker`_ and `wiki`_ are hosted on
GitHub, and all contributions and feedback are more than welcome.
You can test poliastro in your browser using binder,
a cloud Jupyter notebook server:

.. image:: https://img.shields.io/badge/launch-binder-e66581.svg?style=flat-square
   :target: https://beta.mybinder.org/v2/gh/poliastro/poliastro/master?filepath=index.ipynb

.. _`source code`: https://github.com/poliastro/poliastro
.. _`issue tracker`: https://github.com/poliastro/poliastro/issues
.. _`wiki`: https://github.com/poliastro/poliastro/wiki/

poliastro works on recent versions of Python and
is released under the MIT license, hence allowing
commercial use of the library.

.. code-block:: python

    from poliastro.examples import molniya ......................
    from poliastro.plotting import plot

    plot(molniya)

.. include:: success.rst

----

.. _`Molniya orbit`: http://en.wikipedia.org/wiki/Molniya_orbit

Contents
--------

.. toctree::
   :maxdepth: 2

   about
   getting_started
   user_guide
   jupyter
   references
   api
   changelog
