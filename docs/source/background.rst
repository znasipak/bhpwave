Background
==========

For those new to adiabatic waveforms, we give a brief overview of adiabatic inspirals within the
context of black hole perturbation theory and the self-force formalism, with a focus on 
quasi-circular inspirals. This also serves as a reference for the conventions and notation used
throughout the documentation.

Adiabatic inspirals
-------------------

We begin by focusing on the motion of a small body with mass :math:`mu` on a bound orbit about
a massive black hole (MBH) with mass :math:`M` and Kerr spin parameter :math:`a`. Assuming that
:math:`\mu \ll M`, then the leading-order motion is given by a geodesic.

In the restricted case that this leading-order motion is circular, then we can describe the geodesic 
in terms of its (Boyer-Lindquist) radius :math:`r_\mathrm{orb}`, 
its orbital energy :math:`E_\mathrm{orb}`, and its orbital frequency :math:`\Omega`. 
These can all be expressed in terms of the velocity parameter :math:`v^2 = M/r_\mathrm{orb}` and dimensionless
Kerr spin parameter :math:`q = a/M`,

.. math::
    E_\mathrm{orb} &= \frac{1 - 2 v^2 + q v^3}{\sqrt{1-3v^2+ 2qv^3}},
    \\
    \Omega &= \frac{v^3}{M(1 + q v^3)}.

Due to small body's it will emit gravitational waves, which radiate away 
energy. This is captured by the time-averaged GW flux :math:`\langle \dot{E} \rangle_\mathrm{GW}`.
Thus the small body cannot remain on a geodesic. In the adiabatic approximation,
we assume that the time-rate-of-change of orbital energy must be balanced by the GW flux. Therefore,
the orbit evolves according to the adiabatic equations of motion,

.. math::
    \frac{dE_\mathrm{orb}}{dt} &= - \langle \dot{E}\rangle_\mathrm{GW},
    \\
    \frac{d\Phi}{dt} &= \Omega,

where :math:`\Phi` is the orbital phase of :math:`\mu`. Consequently, the trajectory of the small 
body (in Boyer-Lindquist coordinates) is given by

.. math::
    x^\mu_p = (t_p,r_p,\theta_p,\phi_p) = (t, r_0(E_\mathrm{orb}(t)), \pi/2, \Phi(t)).

We consider the small body to rapidly plunge and merge with the MBH when :math:`r_0(t) = r_\mathrm{ISCO}`, where 
:math:`r_\mathrm{ISCO}` is the radius of the innermost stable circular orbit (ISCO), given by

.. math::
    r_\mathrm{ISCO} &= 3 + z_2 - \mathrm{sgn}(q)
    \sqrt{(3 - z_1)(3 + z_1 + 2z_2)},
    \\
    z_1 & = 1 + (1 - q^2)^{1/3}
    \left((1-q)^{1/3} + (1 + q)^{1/3} \right),
    \\
    z_2 &= \sqrt{3 q^2 + z_1^2},

where :math:`q<0` represents a retrograde orbit: an orbit where the orbital angular momentum of 
the small body and spin angular momentum of the MBH are anti-aligned.

Adiabatic waveforms
-------------------