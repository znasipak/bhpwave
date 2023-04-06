References
==========

This code is built upon the self-force formalism, which aims to describe
the dynamics and gravitational wave radiation of black hole binaries using perturbation theory.
Here we give a brief list of relevant references and citations in the self-force literature
that were useful and/or necessary for building :code:`bhpwave`. A comprehensive
review of black hole perturbation theory and the self-force formalism is provided in
`(Pound and Wardell; 2021)`_ along with explicit derivations for adiabatic inspirals
and waveforms. For some of the first numerical studies of
quasi-circular inspirals in Kerr spacetime, see `(Detweiler; 1978)`_, `(Kennefick; 1998)`_, and `(Hughes; 2000)`_.
For recent studies, see `(Taracchini, et al.; 2014)`_ and `(Gralla, Hughes, and Warburton; 2016)`_.
For numerical calculations of both multi-periodic waveform harmonics for both
snapshot and adiabatic inspirals, see `(Drasco and Hughes; 2006)`_ and `(Hughes, et al.; 2021)`_.
Recent open-source codes that implement these methods are discussed further in
the :ref:`external-codes` section below.

Citations
---------

* `(Teukolsky; 1973)`_
* `(Teukolsky and Press; 1974)`_
* `(Detweiler; 1978)`_
* `(Kennefick; 1998)`_
* `(Hughes; 2000)`_
* `(Mino; 2003)`_
* `(Drasco and Hughes; 2006)`_
* `(Hinderer and Flanagan; 2008)`_
* `(Taracchini, et al.; 2014)`_
* `(Gralla, Hughes, and Warburton; 2016)`_
* `(Gourgoulhon, et al.; 2019)`_
* `(Chua, et al.; 2020)`_
* `(Hughes, et al.; 2021)`_
* `(Katz, et al.; 2021)`_
* `(Pound and Wardell; 2021)`_

.. _external-codes:

External Waveform Codes
-----------------------

For general reference, we include other open-source tools that rely on black hole perturbation theory to generate
gravitational waveforms. These tools serve as great sources for comparison with :code:`bhpwave`.

*   `kerrgeodesic_gw`_: A SageMath package, called :code:`kerrgeodesic_gw`, that creates "snapshot" waveforms (see `(Drasco and Hughes; 2006)`_ for more details) for small bodies on circular
    geodesics around Kerr black holes. The waveforms are similar to those produced in :code:`bhpwave`, except
    they do not include backreaction so the inspiral of the small body is not modeled. This code is based
    on the work `(Gourgoulhon, et al.; 2019)`_ and the references therein.

*   `FastEMRIWaveforms`_: A Python package, called :code:`few`, that creates adiabatic waveforms for small bodies
    undergoing eccentric, equatorial inspirals around Schwarzschild massive black holes (MBHs).
    Therefore, :code:`bhpwave` and :code:`few` overlap when the Kerr spin paramter :math:`a=0` in :code:`bhpwave` and when
    the orbital eccentricity :math:`e_0 = 0.` in :code:`few`. This code is based on the work 
    `(Chua, et al.; 2020)`_ and `(Katz, et al.; 2021)`_ and the references therein.

.. _FastEMRIWaveforms: https://bhptoolkit.org/FastEMRIWaveforms/
.. _kerrgeodesic_gw: https://sagemanifolds.obspm.fr/kerrgeodesic_gw/reference/


.. _(Teukolsky; 1973): https://ui.adsabs.harvard.edu/abs/1973ApJ...185..635T/abstract
.. _(Teukolsky and Press; 1974): https://ui.adsabs.harvard.edu/abs/1974ApJ...193..443T/abstract
.. _(Detweiler; 1978): https://ui.adsabs.harvard.edu/abs/1978ApJ...225..687D/abstract
.. _(Kennefick; 1998): https://journals.aps.org/prd/abstract/10.1103/PhysRevD.58.064012
.. _(Hughes; 2000): https://arxiv.org/abs/gr-qc/9910091
.. _(Mino; 2003): https://arxiv.org/abs/gr-qc/0302075
.. _(Drasco and Hughes; 2006): https://arxiv.org/abs/gr-qc/9910091
.. _(Hinderer and Flanagan; 2008): https://arxiv.org/abs/0805.3337
.. _(Taracchini, et al.; 2014): https://arxiv.org/abs/1404.1819
.. _(Gralla, Hughes, and Warburton; 2016): https://arxiv.org/abs/1603.01221
.. _(Gourgoulhon, et al.; 2019): https://www.aanda.org/articles/aa/abs/2019/07/aa35406-19/aa35406-19.html
.. _(Chua, et al.; 2020): https://arxiv.org/abs/2008.06071
.. _(Hughes, et al.; 2021): https://arxiv.org/abs/2102.02713
.. _(Katz, et al.; 2021): https://arxiv.org/abs/2104.04582
.. _(Pound and Wardell; 2021): https://arxiv.org/abs/2101.04592