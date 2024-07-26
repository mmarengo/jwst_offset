"""
Calculate x and y offset for target acquisition with JWST.

This module provides functions to calculate the x and y offset for target
acquisition with JWST. The offset is calculated based on the position of the
science target (RA, DEC and proper motion), and the position of a Gaia DR3
acquisition star (also RA, DEC and proper motion). This routine is used when
direct acquisition of the science target is not possible because of saturation.

The offset is calculated with the following steps:
* Calculate the position of the science target at the time of the observation
* Calculate the position of the acquisition star at the time of the observation
* Calculate the offset between the two positions

The module is implemented using the astropy.coordinates module. It contains
the following routines
"""