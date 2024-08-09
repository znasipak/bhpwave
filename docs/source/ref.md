# References

This code is built upon the self-force formalism, which aims to describe
the dynamics and gravitational wave radiation of black hole binaries using perturbation theory.
Here we give a brief list of relevant references and citations in the self-force literature
that were useful and/or necessary for building {code}`bhpwave`. A comprehensive
review of black hole perturbation theory and the self-force formalism is provided in
[(Pound and Wardell; 2021)] along with explicit derivations for adiabatic inspirals
and waveforms. For some of the first numerical studies of
quasi-circular inspirals in Kerr spacetime, see [(Detweiler; 1978)], [(Kennefick; 1998)], and [(Hughes; 2000)].
For recent studies, see [(Taracchini, et al.; 2014)] and [(Gralla, Hughes, and Warburton; 2016)].
For numerical calculations of both multi-periodic waveform harmonics for both
snapshot and adiabatic inspirals, see [(Drasco and Hughes; 2006)] and [(Hughes, et al.; 2021)].
Recent open-source codes that implement these methods are discussed further in
the {ref}`external-codes` section below.

## Citations

- [(Teukolsky; 1973)]
- [(Teukolsky and Press; 1974)]
- [(Detweiler; 1978)]
- [(Kennefick; 1998)]
- [(Hughes; 2000)]
- [(Mino; 2003)]
- [(Drasco and Hughes; 2006)]
- [(Hinderer and Flanagan; 2008)]
- [(Taracchini, et al.; 2014)]
- [(Gralla, Hughes, and Warburton; 2016)]
- [(Gourgoulhon, et al.; 2019)]
- [(Chua, et al.; 2020)]
- [(Hughes, et al.; 2021)]
- [(Katz, et al.; 2021)]
- [(Pound and Wardell; 2021)]

(external-codes)=

## External Waveform Codes

For general reference, we include other open-source tools that rely on black hole perturbation theory to generate
gravitational waveforms. These tools serve as great sources for comparison with {code}`bhpwave`.

- [kerrgeodesic_gw]: A SageMath package, called {code}`kerrgeodesic_gw`, that creates "snapshot" waveforms (see [(Drasco and Hughes; 2006)] for more details) for small bodies on circular
  geodesics around Kerr black holes. The waveforms are similar to those produced in {code}`bhpwave`, except
  they do not include backreaction so the inspiral of the small body is not modeled. This code is based
  on the work [(Gourgoulhon, et al.; 2019)] and the references therein.
- [FastEMRIWaveforms]: A Python package, called {code}`few`, that creates adiabatic waveforms for small bodies
  undergoing eccentric, equatorial inspirals around Schwarzschild massive black holes (MBHs).
  Therefore, {code}`bhpwave` and {code}`few` overlap when the Kerr spin paramter $a=0$ in {code}`bhpwave` and when
  the orbital eccentricity $e_0 = 0.$ in {code}`few`. This code is based on the work
  [(Chua, et al.; 2020)] and [(Katz, et al.; 2021)] and the references therein.

[(chua, et al.; 2020)]: https://arxiv.org/abs/2008.06071
[(detweiler; 1978)]: https://ui.adsabs.harvard.edu/abs/1978ApJ...225..687D/abstract
[(drasco and hughes; 2006)]: https://arxiv.org/abs/gr-qc/9910091
[(gourgoulhon, et al.; 2019)]: https://www.aanda.org/articles/aa/abs/2019/07/aa35406-19/aa35406-19.html
[(gralla, hughes, and warburton; 2016)]: https://arxiv.org/abs/1603.01221
[(hinderer and flanagan; 2008)]: https://arxiv.org/abs/0805.3337
[(hughes, et al.; 2021)]: https://arxiv.org/abs/2102.02713
[(hughes; 2000)]: https://arxiv.org/abs/gr-qc/9910091
[(katz, et al.; 2021)]: https://arxiv.org/abs/2104.04582
[(kennefick; 1998)]: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.58.064012
[(mino; 2003)]: https://arxiv.org/abs/gr-qc/0302075
[(pound and wardell; 2021)]: https://arxiv.org/abs/2101.04592
[(taracchini, et al.; 2014)]: https://arxiv.org/abs/1404.1819
[(teukolsky and press; 1974)]: https://ui.adsabs.harvard.edu/abs/1974ApJ...193..443T/abstract
[(teukolsky; 1973)]: https://ui.adsabs.harvard.edu/abs/1973ApJ...185..635T/abstract
[fastemriwaveforms]: https://bhptoolkit.org/FastEMRIWaveforms/
[kerrgeodesic_gw]: https://sagemanifolds.obspm.fr/kerrgeodesic_gw/reference/
