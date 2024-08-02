# jwst_offset
Calculate offset for target acquisition with JWST/MIRI. This is the offset
calculated as radial offset and PA angle (East to North) from the Gaia
acquisition star to the target. At the moment this offset is not yet converter
into the X,Y offset required by APT (see note at the end of Jupiter notebook).

For now refer only to the Jupyter notebooks; once they are validated the
software will be converted into a standalone Pythin module.

* r_aqr_offsets.ipynb: calculate offset for R Aqr (science target)
* t_cet_offsets.ipynb: calculate offset for T Cet (PSF reference)

The second notebook is not yet available at the current time.
