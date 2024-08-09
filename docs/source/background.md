# Background

We give a brief overview of adiabatic inspirals and waveforms within the
context of black hole perturbation theory and the self-force formalism, with a focus on
quasi-circular inspirals. This also serves as a reference for the conventions and notation used
throughout the documentation in the {doc}`api/api` section. More in depth reviews of
adiabatic inspirals and waveforms can be found within the {doc}`ref` section.

## Adiabatic inspirals

```{image} images/quasicircular_cartoon.png
:align: left
:width: 250
```

First we focus on the motion of a small body with mass $\mu$ on a bound orbit about
a massive black hole (MBH) with mass $M$ and Kerr spin parameter $a$. Assuming that
$\mu \ll M$, then the leading-order motion is given by a geodesic.
In the restricted case that this leading-order motion is circular, then we can describe the geodesic
in terms of its (Boyer-Lindquist) radius $r_\mathrm{orb}$,
its orbital energy $E_\mathrm{orb}$, and its orbital frequency $\Omega.$
These can all be expressed in terms of the velocity parameter $v^2 = M/r_\mathrm{orb}$ and dimensionless
Kerr spin parameter $q = a/M$ [(Detweiler; 1978)],

$$
E_\mathrm{orb} &= \frac{1 - 2 v^2 + q v^3}{\sqrt{1-3v^2+ 2qv^3}}, \\ \Omega &= \frac{v^3}{M(1 + q v^3)}.
$$ (geoparams)

Due to small body's mass $\mu$, the binary emits gravitational waves, which radiate away
energy. This is captured by the time-averaged GW flux $\langle \dot{E} \rangle_\mathrm{GW}$.
(See [(Teukolsky and Press; 1974)], [(Detweiler; 1978)], [(Hughes; 2000)], [(Hughes, et al.; 2021)] for more details.)
Thus the small body cannot remain on a geodesic. In the adiabatic approximation,
we assume that the time-rate-of-change of orbital energy must be balanced by the GW flux [(Hughes; 2000)], [(Mino; 2003)]. Therefore,
the orbit evolves according to the adiabatic equations of motion,

$$
\frac{dE_\mathrm{orb}}{dt} &= - \langle \dot{E}\rangle_\mathrm{GW}, \\ \frac{d\Phi}{dt} &= \Omega,
$$ (fluxbalance)

where $\Phi$ is the orbital phase of $\mu$. The trajectory of the small
body (in Boyer-Lindquist coordinates) is then given by

$$
x^\mu_p = (t_p,r_p,\theta_p,\phi_p) = (t, r_\mathrm{orb}(E_\mathrm{orb}(t)), \pi/2, \Phi(t)).
$$ (traj)

We consider the small body to rapidly plunge and merge with the MBH when $r_\mathrm{orb}(t=t_\mathrm{merge}) = r_\mathrm{ISCO}$, where
$r_\mathrm{ISCO}$ is the radius of the innermost stable circular orbit (ISCO), given by

$$
r_\mathrm{ISCO} &= 3 + z_2 - \mathrm{sgn}(q) \sqrt{(3 - z_1)(3 + z_1 + 2z_2)}, \\ z_1 & = 1 + (1 - q^2)^{1/3} \left((1-q)^{1/3} + (1 + q)^{1/3} \right), \\ z_2 &= \sqrt{3 q^2 + z_1^2}.
$$ (isco)

Formally, the adiabatic approximation breaks down as the small body approaches the ISCO. However,
for small mass-ratios ($\mu/M \ll 1$) this transition from adiabatic inspiral to plunge is quite short,
and thus can be neglected when generating leading-order waveforms.
Note that $q<0$ represents a retrograde orbit: an orbit where the orbital angular momentum of
the small body and spin angular momentum of the MBH are anti-aligned.

## Adiabatic waveforms

From this adiabatic inspiral, we construct
the resulting adiabatic gravitational waveform. It can be expressed as a sum over many
harmonics with complex amplitudes $H_{lm}(t) = H_{lm}(r_\mathrm{orb}(t))$ and
phases $\Phi_m(t) = m\Phi(t)$ [(Hughes, et al.; 2021)] [(Pound and Wardell; 2021)],

$$
h(t, r_\mathrm{obs}, \theta_\mathrm{obs}, \phi_\mathrm{obs}) &= \frac{\mu}{r_\mathrm{obs}} \sum_{l=2}^\infty \sum_{|m|=1}^{l} H_{lm}(t) e^{i\Phi_m(t)} {}_{-2} Y_{lm}(\theta_\mathrm{obs}, \phi_\mathrm{obs}), \\ &= \sum_{l=2}^\infty \sum_{|m|=1}^{l} h_{lm}(t,r_\mathrm{obs}, \theta_\mathrm{obs}, \phi_\mathrm{obs}), \\ &= h_+(t) - i h_\times(t),
$$ (hsum_full)

where $(r_\mathrm{obs}, \theta_\mathrm{obs}, \phi_\mathrm{obs})$ is the position of the
observor in Boyer-Lindquist coordinates, ${}_{-2} Y_{lm}$ is the spin-weighted spherical
harmonic of spin-weight $-2$, and $h_+$ and $h_\times$ are the plus
and cross polarizations of the gravitational wave strain in the source frame (of the binary).
The slowly-varying amplitudes experience a change $\Delta H_{lm} \sim O(1)$ over
the inspiral, while the rapidly-varying phases experience a change $\Delta H_{lm} \sim O(M/\mu)$
[(Hinderer and Flanagan; 2008)].

## Simplified waveform summation

For circular orbits, the harmonic amplitudes possess the symmetry $H_{l,-m} = (-1)^{l}\bar{H}_{lm}$,
where the overbar denotes complex conjugation. Similarly, the spin-weighted spherical harmonics
possess the relation ${}_{s}Y_{l, -m} = (-1)^{s+m} {}_{-s} \bar{Y}_{lm}$. Combining these
and separating the complex amplitudes and spin-weighted spherical harmonics in terms of their magnitudes and phases,
i.e., $H_{lm}(t) = A_{lm}(t) e^{i\phi_{lm}(t)}$ and ${}_{s}Y_{lm}(\theta,\phi) = {}_{s}y_{lm}(\theta) e^{im\phi}$, we find

$$
h_{lm}^+ + h_{l, -m}^+ &= \frac{\mu}{r_\mathrm{obs}} A_{lm}({}_{-2} y_{lm} + (-1)^{l+m} {}_{+2}y_{lm})\cos(m\Phi - m\phi_\mathrm{obs} - \phi_{lm}), \\ & = H^+_{lm}, \\ h_{lm}^\times + h_{l, -m}^\times &= \frac{\mu}{r_\mathrm{obs}} A_{lm}({}_{-2} y_{lm} - (-1)^{l+m} {}_{+2}y_{lm})\sin(m\Phi - m\phi_\mathrm{obs} - \phi_{lm}), \\ & = H^\times_{lm}.
$$ (HlmPlusCross)

This reduces the mode-sum over the waveform harmonics and allows one to easily separate the calculation
of the plus and cross polarizations,

$$
h_{+,\times} = \sum_{l=2}^\infty \sum_{m = 1}^{l} H^{+,\times}_{lm}.
$$ (hsum)

In {code}`bhpwave`, waveforms are generated using the mode-sum in Eq. {eq}`hsum`. Therefore,
requesting a specific $(l, m)$ mode returns the combination $H^+_{lm} - i H^\times_{lm}$.
In other words, for a given value of $m$, {code}`bhpwave` performs an internal sum over both
the postive and negative values of $m$ by default.

## Transforming frames of reference

While Eq. {eq}`hsum_full` provides a full description of the adiabatic waveform in the
source frame, often we are interested in modeling the waveform measured by an observor in
their own frame. In {code}`bhpwave` we take the observor frame to be the solar system
barycenter (SSB) frame. Following the conventions of [(Katz, et al.; 2021)], an observor
in the SSB frame parametrizes a generic source in terms of the parameters

- $M$: the (redshifted) mass of the massive black hole
- $\mu$: the (redshifted) mass of the smaller compact object
- $q$: the dimensionless black hole spin
- $p_0$: the initial semi-latus rectum of the binary (a measure of radial separation)
- $e_0$: the initial orbital eccentricity of the binary
- $x_0$: cosine of the intial orbital inclination of the binary
- $D_L$: the luminosity distance to the source (which is equivalent to $r_\mathrm{obs}$)
- $q_{S}$: the polar angle of the source's sky location
- $\phi_{S}$: the azimuthal angle of the source's sky location
- $q_{K}$: the polar angle of the Kerr spin vector
- $\phi_{K}$: the azimuthal angle of the Kerr spin vector
- $\Phi_{\phi 0}$: the initial azimuthal position of the small compact object
- $\Phi_{r 0}$: the phase describing the initial radial position and velocity of the small compact object
- $\Phi_{\theta 0}$: Phase describing the initial polar position and velocity of the small compact object

For quasi-circular inspirals, we have the simplifications $p_0 = r_0$,
$e_0 = 0$, $x_0 = 1$, $\Phi_{r 0} = 0$,  $\Phi_{\theta 0} = 0$.

The SSB angle variables $(q_{S}, \phi_{S}, q_{K}, \phi_{K})$ are related
to the source frame angles $(\theta, \phi)$ via the transformations [(Katz, et al.; 2021)],

$$
\phi = -\frac{\pi}{2}, \qquad \cos\theta = -(\sin q_S\sin q_K \cos(\phi_S - \phi_K) + \cos q_S \cos q_K).
$$ (angleTransforms)

Additionally, the $z$-axes of the SSB and source frames will generally be misaligne. Consequently,
the plus and cross gravitational wave polarizations in the source frame will not be equivalent to the
plus and cross polarizations observed in the SSB frame. Instead the observed waveform $h_\mathrm{SSB}$
will experience a phase shift $h_\mathrm{SSB} = e^{2i\psi} h$, where the change in the polarization
angle $\psi$ is defined as,

$$
e^{i\psi} = \frac{\cos q_S\sin q_K \cos(\phi_S - \phi_K) - \sin q_S \cos q_K - i \sin q_K \sin(\phi_S - \phi_K)}{|\cos q_S\sin q_K \cos(\phi_S - \phi_K) - \sin q_S \cos q_K - i \sin q_K \sin(\phi_S - \phi_K)|}.
$$ (polarization)

## Frequency domain waveforms

In progress...

[(chua, et al.; 2020)]: https://arxiv.org/abs/2008.06071
[(detweiler; 1978)]: https://ui.adsabs.harvard.edu/abs/1978ApJ...225..687D/abstract
[(drasco and hughes; 2006)]: https://arxiv.org/abs/gr-qc/9910091
[(gourgoulhon, et al.; 2019)]: https://www.aanda.org/articles/aa/abs/2019/07/aa35406-19/aa35406-19.html
[(hinderer and flanagan; 2008)]: https://arxiv.org/abs/0805.3337
[(hughes, et al.; 2021)]: https://arxiv.org/abs/2102.02713
[(hughes; 2000)]: https://arxiv.org/abs/gr-qc/9910091
[(katz, et al.; 2021)]: https://arxiv.org/abs/2104.04582
[(kennefick; 1998)]: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.58.064012
[(mino; 2003)]: https://arxiv.org/abs/gr-qc/0302075
[(pound and wardell; 2021)]: https://arxiv.org/abs/2101.04592
[(teukolsky and press; 1974)]: https://ui.adsabs.harvard.edu/abs/1974ApJ...193..443T/abstract
[(teukolsky; 1973)]: https://ui.adsabs.harvard.edu/abs/1973ApJ...185..635T/abstract
