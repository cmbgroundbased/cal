@page church_model The atmospheric model

The atmospheric model implemented in CAL comes down to the famous S. Church paper, <em>Church, S.E., 1995. Predicting residual levels of atmospheric sky noise in ground-based observations of the Cosmic Background Radiation. Monthly Notices of the Royal Astronomical Society, 272(3), pp.551-569.</em>. In this page the main results are summed up, and discussed in reference with the code.

# The Atmosphere as a dispersive medium
The atmosphere is a dispersive medium, this means that its electrical permittivity \f$\varepsilon(\omega) \in \mathbb{C}\f$. So we can write

\f[\varepsilon(\omega) = \varepsilon_R(\omega) + i\,\varepsilon_I(\omega)\f]

where the real is proportional on the refractive index, while the image part is proportional to the atmospheric absorption coefficient

\f[ \varepsilon(\omega)_R = \sqrt{n}; \,\, \varepsilon(\omega)_I = \frac{\lambda \alpha}{4\pi} \f]

The contribute to complex permittivity is due mostly by the water vapor, and Oxygen molecules that are dissolved in the atmosphere. The Oxygen molecules are well mixed in the atmosphere and they affect the CMB observations, basically, with a constant load for constant elevation scans. The same is not true for what concerns the water vapor distributions. The water vapor molecules are distributed in the atmosphere with a spatial spectra density that follow a power-low distributions. If \f$ k \f$ represents the vectors of the spatial Fourier transform, the water vapor is distributed in a volume of atmosphere with this spatial power spectral density

\f[ P(k) = k ^ {-b}. \f]

The water vapor affects both the refractive index and the absorption coefficient because they are connected one to each other by the Kramers-Kroning relations.

The refractive index fluctuations affect, in case of single dish telescope, the telescope pointing, while the fluctuations in the absorption coefficient, since it's related to the atmospheric emission, are seen by the telescope with antenna temperature fluctuations.

\image{center} html "./img/emission_atm.jpg" "Absorption Coefficient"
