# JWST offsets for MIRI coronagraphic observations of R Aqr
Author: Massimo Marengo (2024)
Last update: 09/05/2024
Note: this software is in development and not yet completed

## Requirements ##
- numpy
- matplotlib
- astropy
- pySIAF
- coronographic_offset_ta

This last module has been developed by Jonathan Aguilar and is available at
https://github.com/STScI-MIRI/coronagraphy_offset_ta

## Description ##

This package provides a set of Jupyter notebooks to calculate the dx, dy
offsets for source acquisition from an offset source. The current notebook 
provides the calculation for R Aqr in two of the MIRI 4QPM coronagraphs.
A similar notebook for the PSF reference star is being prepared.

The calculations are done in the following steps:

- Define observing parameters
- Define science target and acquisition target astrometry
- Calculate updated astrometry for science and acquisition target
  at the beginning and end of observing window. This also provides the
  offset between the science target and the acquisition star, to check
  that this distance is within the Visit Splitting Distance.
- Calculate dx, dy offsets based on updated coordinates. This is to be done
  separately for each mask used, and for the beginning and end time of the
  observation window.

A separate notebook will be used to compare the position of the acquisition
stars with respect to the main PSF features (e.g. spikes) of the nearby bright
stars (R Aqr and the PSF star respectively). This will provide appropriate
constraints for the V3PA angle (and the observing time window) that will 
allow a safe target acquisition.

### Observing Parameters ###

This scripts are designed to allow experimenting with different obserning
windows (must be within the visibility window), to then use as constraints
for the observations. The following parameters are required:

- JWST_TIME1: beginning of JWST observing window
- JWST_TIME2: beginning of JWST observing window
- V3PA_TIME1: chosen V3PA at the beginning of JWST observing window
- V3PA_TIME2: chosen V3PA at the end of JWST observing window
- CORON_ID_A: # ID for first coronagraphic observation
- CORON_ID_B: # ID for second coronagraphic observation

For the R Aqr program we use `CORON_ID_A = '1140'` and `CORON_ID_B = '1550'`.

The V3PA angles need to be set within the available angles for a specified
oberving date. The allowed range can be obtained from APT: Visit Planner ->
Reports -> Visit -> Total Roll Analysis for Visit. We need to check with
technical support that this is really as it works.

### Science and acquisition target astronomy ###

For both science target we need to provide a SkyCoord object with RA, DEC,
PMRA, PMDEC, PARALLAX, date when the coordinates had been derived, and
a coordinate frame of reference (e.g. 'icrs') with reference epoch.

For R Aqr we rely on our VLBI astrometry. For the PSF reference and the 
acquisition stars we rely on Gaia DR3 astronomy.

### Products ###

The scripts are designed to provide the following parameters, that will be
used to update the observation planning in APT:

- dx, dy offsets for target acquisition (directly calculated by the notebooks).
  These go in APT in `Special Requirements` -> `Offset`.
- Optimal observation window for which the offsets have been calculated,
  and that satisfy the visibility and asteroid avoidance requirements, and for
  the change in offsets is small enough to not affect significantly the
  positioning of the coronagraph on the science source. There may also be
  some interplay between the obserbing window and the choice of V3PA. The
  observing window is set in APT in `Special Requirements` -> `Timing` ->
  `Between Dates`.
- V3PA range, to be choses such that the target acquisition star is clear from
  the main features (e.g. the diffraction spikes) of the science targets PSF.
  This is set in APT in `Special Requirements` -> `Position Angle` ->
  `PA Range`. Set checkbox on `V3PA`.

## Status ##

The calculation described above are available for R Aqr, and can be accessed
in the Jupyter notebook `r_aqr_miri_offset_ta.ipnb`.

We need to create a new notebook to check the position of the acquisition
stars with respect to the R Aqr (and, separately, the PSF star) PSF, using
the `webpsf` tools. This will allow to determine the relative position of
the acquisition star with respect to the main PSF features (e.g. spikes) of
the nearby bright source.

The whole procedure need then to be repeated for the PSF star.
