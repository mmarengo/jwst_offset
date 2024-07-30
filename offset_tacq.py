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
the following routines:

>>> update_coords: Update the coordinates of a target given its RA, DEC and
    proper motion.
>>> calculate_offset: Calculate the offset between two targets given their
    coordinates.
"""

import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord

def update_coords(ref_coords, new_time):
    """
    Updates the coordinates of a source given its RA, DEC and proper motion.

    Parameters
    ----------

    ref_coords: SkyCoord
        The reference coordinates of the source. It should be an
        astropy.coordinates.SkyCoord object with the following attributes:
        - ra: The right ascension of the source in degrees.
        - dec: The declination of the source in degrees.
        - pm_ra_cosdec: The proper motion in right ascension in mas/yr.
        - pm_dec: The proper motion in declination in mas/yr.
        - ref_epoch: The epoch of the coordinates.
        - obstime: The time of the observation when ref_coords were measured.

    new_time: Time
        The time of the observation when the new coordinates are needed.
    
    Returns
    -------

    new_coords: SkyCoord
        The updated coordinates of the source at the new_time.
    """