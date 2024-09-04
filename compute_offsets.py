"""
Offset TA Tools
Author: Jonathan Aguilar
Latest update: Feb 8, 2024

USAGE
-----
This tool is meant to help users compute offsets in the case that they need to
perform TA on a target other than their science target. The data required are
the sky coordinates of the acquisition and science targets, and precise
telescope ROLL_REF angles, at the epoch of observation.

There are two ways to use this tool: running it as a script, or importing it as
a module. To run it as a script, users should edit the relevant values at the
bottom of this file, after `if __name__ == "__main__"`. The instructions for
editing each section are provided in comments next to the relevant lines of
code. Once the values have been entered, this can be run as a script (e.g.
`python compute_offsets.py`). It is suggested that users make a new copy of the
file for every unique set of parameters.

An example of using this file as an imported module is provided in
`example.py`. In this case, `compute_offsets.py` should be in the same
directory as the code performing the import.

The main function of `compute_offsets.py` is also called `compute_offsets()`.
It will print out the appropriate X and Y offsets to enter into APT, and
display some diagnostic plots if requested.


Requirements:
- numpy
- matplotlib
- astropy
- pySIAF
"""

import numpy as np

import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Polygon

from astropy.coordinates import SkyCoord
from astropy import units

import pysiaf
from pysiaf import Siaf


#------------------------------------------------------#
#----------------- Offset Computation -----------------#
#------------------------------------------------------#
miri = Siaf("MIRI")

def create_attmat(
        position : SkyCoord,
        aper : pysiaf.aperture.JwstAperture,
        pa : float,
        idl_offset : tuple[float, float] = (0., 0.),
        set_matrix : bool = False
) -> np.ndarray:
    """
    Create an attitude matrix for JWST when the reference point of a particular
    aperture is pointed at a given position for a specified PA

    Parameters
    ----------
    position : SkyCoord
      skycoord position on the sky
    aper : pysiaf.aperture.JwstAperture
      pySIAF-defined aperture object
    pa : float
      PA angle with respect to the V3 axis measured at the aperture reference
      point. This corresponds to the V3PA field in the APT PA range special
      requirement, and the ROLL_REF keyword in the data. This is *not* the
      PA_V3 keyword value, which is the PA angle of the V3 axis measured at the
      telescope boresight.
    idl_offset : tuple[float, float] = (0.0, 0.0)
      allows you to specify an arbitrary position in IDL coordinates that
      corresponds to the position
    set_matrix : bool = True
      if True, also set the matrix on the current aperture in addition to returning it
    Output
    ------
    attmat : np.ndarray
      matrix that pySIAF can use to specify the attitude of the telescope
    """
    v2, v3 = aper.idl_to_tel(idl_offset[0], idl_offset[1])
    # v2, v3 = aper.reference_point('tel')
    # compute the attitude matrix when we're pointing directly at the TA target
    attmat = pysiaf.utils.rotations.attitude_matrix(v2, v3,
                                                    ra=position.ra.deg,
                                                    dec=position.dec.deg,
                                                    pa=pa)
    if set_matrix == True:
        aper.set_attitude_matrix(attmat)
    return attmat


def sky_to_idl(
        stars : list[dict],
        aper : pysiaf.aperture.JwstAperture,
        pa : float,
        idl_offset : tuple[float, float] = (0., 0.)
) -> list[dict]:
    """
    Convert RA and Dec positions of a TA star and its target with an offset
    into a detector position (measured from the reference point, in arcsec)
    for a given PA of the V3 axis w.r.t. North.
    Assume the TA star is centered on the aperture (idl = (0, 0))

    Parameters
    ----------
    stars : list of {"label": label, "position": pos} elements to plot
      the first element is the ACQ target. The aperture is centerd on this target.
    aper : SIAF object for the aperture being used (e.g. MIRIM_MASK1550)
    pa : the PA in degrees of the V3 axis of the telescope (measured eastward of North) at the observation epoch

    Output
    ------
    idl_coords : dict {label: idl_tuple} of IDL coordinates for each target
      a dictionary of IDL x, y coordinates for the TA star and the science target
      the TA star coordinates should be very close to 0
    """
    acq_pos = stars[0]['position']
    attmat = create_attmat(acq_pos, aper, pa, idl_offset)
    aper.set_attitude_matrix(attmat)
    idl_coords = []
    # ta star - should be close to 0
    for star in stars:
        label = star['label']
        pos = star['position']
        idl_coords.append({'label': label, 'position': np.array(aper.sky_to_idl(pos.ra.deg, pos.dec.deg))})
    return idl_coords

def print_offset_information(
        slew_from : dict,
        slew_to: dict,
        idl_coords : list,
        sep : units.Quantity,
        pa : units.Quantity,
        offset : np.ndarray,
) -> None:
    """
    Print the information about offsets to the user.

    Parameters
    ----------
    define your parameters

    Output
    ------
    None - prints to screen

    """
    print_output = []
    print_output.append("Computing offset command values to slew from:")
    print_output.append(f"\t{slew_from['label']}")
    print_output.append(f"\t\t RA: \t {slew_from['position'].ra.degree}")
    print_output.append(f"\t\t Dec: \t {slew_from['position'].dec.degree}")
    print_output.append(f"to:")
    print_output.append(f"\t{slew_to['label']}")
    print_output.append(f"\t\t RA: \t {slew_to['position'].ra.degree}")
    print_output.append(f"\t\t Dec: \t {slew_to['position'].dec.degree}")
    print_output.append(f"\n")
    print_output.append(f"Separation and PA: {sep.mas:0.2f} mas, {pa.degree:0.2f} deg\n")

    print_output.append(f"\n")
    print_output.append(f"After TA but before slewing, the position of the ACQ star should be close to (0, 0):")
    print_output.append(f"\t" + ', '.join(f"{i:+0.3e}" for i in idl_coords[0]['position']) + " arcsec")
    print_output.append(f"... and the position of the SCI star is:")
    print_output.append(f"\t" + ', '.join(f"{i:+0.3e}" for i in idl_coords[1]['position']) + " arcsec")

    print_output.append(f"\n")
    print_output.append(f"When the ACQ star is centered, the SCI star is at:")
    print_output.append(f"\tdX: {idl_coords[1]['position'][0]:+2.6f} arcsec")
    print_output.append(f"\tdY: {idl_coords[1]['position'][1]:+2.6f} arcsec")

    if len(idl_coords) > 2:
        print_output.append(f"\n")
        print_output.append(f"Before the slew, here are the positions of the other targets provided:")
        for ic in idl_coords[2:]:
            print_output.append(f"{ic['label']}")
            print_output.append(f"\tdX: {ic['position'][0]:+2.6f} arcsec")
            print_output.append(f"\tdY: {ic['position'][1]:+2.6f} arcsec")

    print_output.append(f"\n")
    print_output.append(f"Sanity check: on-sky angular separation should be the same distance as that of the slew.")
    # print in nice columns
    sep_as_str = f"{sep:0.6f}"
    slew_mag_str = f"{np.linalg.norm(offset) * units.arcsec :0.6f}"
    check_str = [['Separation', sep_as_str], ['Slew magnitude', slew_mag_str]]
    for row in check_str:
        print_output.append("{: >20} {: >20}".format(*row))

    if len(idl_coords) > 2:
        print_output.append(f"\n")
        print_output.append(f"After the slew, here are the positions of other targets provided:")
        for ic in idl_coords[2:]:
            print_output.append(f"{ic['label']}")
            print_output.append(f"\tdX: {ic['position'][0]+offset[0]:+2.6f} arcsec")
            print_output.append(f"\tdY: {ic['position'][1]+offset[1]:+2.6f} arcsec")

    print_output.append("\n")
    print_output.append("Therefore, the commanded offsets that will move the coronagraph from the ACQ star to the SCI are:")
    print_output.append(f"\tdX: {offset[0]:+4.6f} arcsec")
    print_output.append(f"\tdY: {offset[1]:+4.6f} arcsec")
    print_output.append("\n")

    for line in print_output:
        print(line)


def compute_offsets(
        slew_from: dict,
        slew_to: dict,
        v3pa: float,
        coron_id : str,
        other_stars : list = [],
        verbose : int = 1,
        show_plots : bool = True,
        plot_full : bool = False,
        return_offsets : bool = False,
) -> np.ndarray :
    """
    Compute the slews for the TA sequences, print the offsets, and show the plots if requested.
    How it works:
    - Point the coronagraphic aperture (MASK or CORON) at the TA star by
      setting an attitude matrix for the V3PA value
    - The offset is the negative of the IDL coordinates of the SCI star
    - The rest of the machinery is basically just making verification plots

    Parameters
    ----------
    slew_from: dict
      A dictionary containing the label and position of the TA target, set by
      the user in compute_offsets.py
    slew_to: dict
      A dictionary containing the label and position of the science target, set by
      the user in compute_offsets.py
    v3pa: float
      The PA_V3 angle of the telescope for this observation
    coron_id : str
      Identifier for the desired coronagraphic subarray. Must be one of '1065',
      '1140', '1550', and 'LYOT'
    other_stars : list
      A list of dicts of other stars in the field, in the same format as slew_from/slew_to
    verbose : int = 1
      print diagnostics and offsets to screen.
      0 : nothing is printed
      1 : all output is printed
      2 : only the final IDl coordinates of the stars are printed
    show_plots : bool = True
      If True, display the diagnostic plots. If False, only print the offsets.
    plot_full : bool = True
      If True, plot the MIRI Imager footprint in addition to the coronagraphic
      subarray. This is useful for FULL array readout mode.
    return_offsets : bool = False
      If True, return an array of dx and dy offsets

    Output
    ------
    Prints offsets and shows plots. Returns a dict of floats

    """
    coron_ids = ['1065','1140','1550', 'LYOT']
    # make sure coron_id is valid
    try:
        assert(coron_id.upper() in coron_ids)
    except AssertionError:
        print(f"Error: bad value for `coron_id` ({coron_id})")
        print(f"Must be one of [{', '.join(coron_ids)}].")
        return np.array([np.nan, np.nan])
    coron_id = coron_id.upper()
    star_positions = [slew_from, slew_to] + other_stars
    labels = {'ACQ': slew_from['label'],
              'SCI': slew_to['label']}


    # Offsets to science target
    sep = star_positions[0]['position'].separation(star_positions[1]['position']).to(units.arcsec)
    pa = star_positions[0]['position'].position_angle(star_positions[1]['position']).to(units.deg)


    # There are two relevant apertures: MIRIM_MASK[XXXX], which is the entire subarray, and
    # MIRIM_CORON[XXXX], which is just the portion that gets illuminated
    # let's combine all the SIAF objects in a dict for convenience

    all_apers = {}
    all_apers['UR'] = miri[f'MIRIM_TA{coron_id}_UR']
    all_apers['CUR'] = miri[f'MIRIM_TA{coron_id}_CUR']
    all_apers['coro'] = miri[f'MIRIM_CORON{coron_id}']
    all_apers['mask'] = miri[f'MIRIM_MASK{coron_id}']
    if plot_full == True:
        # plot the MIRIM detector and the imager illuminated footprint
        # all_apers['full'] = miri['MIRIM_FULL']
        all_apers['imager'] = miri['MIRIM_ILLUM']


    idl_coords = sky_to_idl(star_positions,
                            all_apers['mask'],
                            v3pa)
    # The offset you apply is as if you were moving the science target - i.e.
    # the negative of its position
    offset = -1*np.array(idl_coords[1]['position'])

    if verbose == 1:
        print_offset_information(
            slew_from=slew_from,
            slew_to=slew_to,
            idl_coords=idl_coords,
            sep=sep,
            pa=pa,
            offset=offset
        )
    if verbose == 2:
        len_label = max(len(star['label']) for star in idl_coords)
        print("IDL coordinates of all stars after slew:")
        for star in idl_coords:
            label = star['label']
            idl = star['position']
            print(f"{label:{len_label}s}:\t{idl[0]+offset[0]:+0.10f}, {idl[1]+offset[1]:+0.10f}")
        print("")

    if show_plots == True:
        make_plots(
            all_apers,
            star_positions,
            v3pa,
            idl_coords,
            offset,
        )

    if return_offsets == True:
        return offset


def work_backwards(
        slew_from : dict,
        slew_to : dict,
        coron_id : str,
        v3pa : float,
        offset : tuple[float, float] = (0., 0.),
        other_stars : list = [],
) -> list[dict] :
    """
    Work backwards to find out where the slew_to star ended up, given sky
    coordinates for the targets, the commanded offset, and the v3pa value of
    the observation.

    Parameters
    ----------
    slew_from : dict
      label and coordinate of the ACQ target
    slew_to : dict
      label and coordinate of the SCI target
    coron_id : str
      Identifier for the desired coronagraphic subarray. Must be one of '1065',
      '1140', '1550', and 'LYOT'
    v3pa : float
      the v3pa angle of the telescope at the aperture reference position, in
      degrees
    offset : np.ndarray[float]
      the x and y offset commanded
    other_stars : list
      A list of dicts of other stars in the field that you might want to keep
      track of, in the same format as slew_from/slew_to

    Output
    ------
    idl_positions : list[dict]
      A list of positions in subarray IDL coordinates provided targets. Each
      list entry has format {'label': label, 'position': position}.

    """
    coron_id = coron_id.upper()
    coro = miri[f'MIRIM_CORON{coron_id}']
    mask = miri[f'MIRIM_MASK{coron_id}']

    star_positions = [slew_from, slew_to] + other_stars

    idl_coords = sky_to_idl(star_positions,
                            coro,
                            v3pa,
                            idl_offset=offset)

    return idl_coords

#--------------------------------------------#
#----------------- Plotting -----------------#
#--------------------------------------------#

def plot_apers(
        ax,
        attmat,
        aper_dict,
        format_dict={}
):
    """
    Helper function to plot the apertures for a given part of the TA sequence
    ax : axis to plot on
    attmat : attitude matrix
    aper_dict : dictionary of apertures to plot
    format_dict : aperture plot formatting parameters
    """
    for k, aper in aper_dict.items():
        aper.set_attitude_matrix(attmat)
        format_thisaper = format_dict.copy()
        if k == 'mask':
            format_thisaper['mark_ref'] = True
        if k == 'coro':
            # skip the illuminated region aperture, it's too crowded
            pass
        else:
            aper.plot(ax=ax, label=False, frame='sky', fill=False, **format_thisaper)

def plot_before_offset_slew(
        aper_dict : dict,
        idl_coords : list,
        star_positions : list =[],
        ax = None,
):
    """Plot the scene on the detector when you're pointed at the acquisition target"""
    # plot 1 : POV of the detector
    if ax is None:
        fig, ax = plt.subplots(1, 1, layout='constrained')
    else:
        fig = ax.get_figure()
    if star_positions != []:
        title = f"{star_positions[0]['label']} --> {star_positions[1]['label']}"
        fig.suptitle(title)
    ax.set_title("""Positions *before* offset slew""")
    frame = 'idl' # options are: tel (telescope), det (detector), sci (aperture)
    aper_dict['mask'].plot(ax=ax, label=False, frame=frame, c='C0')
    aper_dict['coro'].plot(ax=ax, label=False, frame=frame, c='C1', mark_ref=True)

    ax.scatter(0, 0,
               c='k',
               label=f"ACQ/{idl_coords[0]['label']}",
               marker='x',
               s=100)
    ax.scatter(*idl_coords[1]['position'],
               label=f"SCI/{idl_coords[1]['label']}",
               marker="*",
               c='k')
    for star in idl_coords[2:]:
        ax.scatter(*star['position'],
                   # c='k',
                   label=star['label'],
                   marker='.',
                   s=50)
    ax.add_artist(quad_boundaries(aper_dict['coro'], kwargs={'fc': 'grey'}))
    ax.legend()
    ax.set_aspect("equal")
    ax.grid(True, ls='--', c='grey', alpha=0.5)
    return fig


def plot_after_offset_slew(
        aper_dict : dict,
        idl_coords : list,
        offset : np.ndarray,
        star_positions : list =[],
        ax = None,
):
    """Plot the scene on the detector when you're pointed at the science target"""
    # plot 1 : POV of the detector
    if ax is None:
        fig, ax = plt.subplots(1, 1, layout='constrained')
    else:
        fig = ax.get_figure()
    if star_positions != []:
        title = f"{star_positions[0]['label']} --> {star_positions[1]['label']}"
        fig.suptitle(title)
    ax.set_title("""Positions *after* offset slew""")
    frame = 'idl' # options are: tel (telescope), det (detector), sci (aperture)
    aper_dict['mask'].plot(ax=ax, label=False, frame=frame, c='C0')
    aper_dict['coro'].plot(ax=ax, label=False, frame=frame, c='C1', mark_ref=True)

    ax.scatter(0 + offset[0], 0 + offset[1],
               c='k',
               label=f"ACQ/{idl_coords[0]['label']}",
               marker='x',
               s=100)
    ax.scatter(*(idl_coords[1]['position'] + offset),
               label=f"SCI/{idl_coords[1]['label']}",
               marker="*",
               c='k')
    for star in idl_coords[2:]:
        ax.scatter(*(star['position'] + offset),
                   # c='k',
                   label=star['label'],
                   marker='.',
                   s=50)
    ax.add_artist(quad_boundaries(aper_dict['coro'], kwargs={'fc': 'grey'}))
    ax.legend()
    ax.set_aspect("equal")
    ax.grid(True, ls='--', c='grey', alpha=0.5)
    return fig


def plot_detector_ta_sequence(
        aper_dict,
        ta_sequence,
        idl_coords,
        star_positions={},
        axes=None
) -> mpl.figure.Figure:
    """Plot the TA sequence as seen by the detector"""
    if axes is None:
        nrows = 1
        ncols = 4
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(2*ncols, 4*nrows),
                                 sharex=True, sharey=True, layout='constrained')
    else:
        fig = axes[0].get_figure()

    if star_positions != {}:
        title= f"{star_positions[0]['label']} --> {star_positions[1]['label']}\n"
    else:
        title=''
    fig.suptitle(title + f"TA sequence, as seen by the detector")

    # plot the SIAF apertures on every plot
    for ax in axes:
        aper_dict['mask'].plot(ax=ax, label=False, frame='det', fill=False, c='C0')
        ax.plot([], [], c='C0', label='Readout')
        aper_dict['coro'].plot(ax=ax, label=False, frame='det', mark_ref=True, fill=False, c='C1')
        ax.plot([], [], c='C1', label='Illuminated')
        # TA aperturesL
        for aper_id in ['UR', 'CUR']:
            aper_dict[aper_id].plot(ax=ax, label=False, frame='det', mark_ref=True, fill=False, c='C2')
        ax.plot([], [], c='C2', label='TA regions')

    # plot the positions of the stars at each step in the TA sequence
    # Outer TA
    ax = axes[0]
    ta_aper_id = "UR"
    ax.set_title("Step 1\n" + f"{ta_aper_id} TA region")

    # use the TA aperture object to convert coordinates
    ta_aper = aper_dict[ta_aper_id]
    acq_pos = ta_aper.idl_to_det(*ta_sequence[ta_aper_id][0]['position'])
    sci_pos = ta_aper.idl_to_det(*ta_sequence[ta_aper_id][1]['position'])
    acq_label = ta_sequence[ta_aper_id][0]['label']
    sci_label = ta_sequence[ta_aper_id][1]['label']

    ax.scatter(*acq_pos, 
               c='k', label=acq_label, marker='x', s=100)
    ax.scatter(*sci_pos,
               c='k', label=sci_label, marker='*', s=100)
    other_pos = {i['label']: ta_aper.idl_to_det(*i['position']) for i in ta_sequence[ta_aper_id][2:]}
    for label, pos in other_pos.items():
        ax.scatter(*pos, marker='.', s=50, label=label)

    # put the legend on this plot
    ax.legend(loc='best', ncol=1, fontsize='small', markerscale=0.7)

    # Inner TA
    ax = axes[1]
    ta_aper_id = 'CUR'
    ax.set_title("Step 2\n" + f"{ta_aper_id} TA region")
    # use the TA aperture object to convert coordinates
    ta_aper = aper_dict[ta_aper_id]
    ta_aper.plot(ax=ax, label=False, frame='det', mark_ref=True, fill=False, c='C2')

    acq_pos = ta_aper.idl_to_det(*ta_sequence[ta_aper_id][0]['position'])
    sci_pos = ta_aper.idl_to_det(*ta_sequence[ta_aper_id][1]['position'])
    ax.scatter(*acq_pos,
               c='k', label='ACQ', marker='x', s=100)
    ax.scatter(*sci_pos,
               c='k', label='SCI', marker='*', s=100)
    other_pos = [ta_aper.idl_to_det(*i['position']) for i in ta_sequence[ta_aper_id][2:]]
    for pos in other_pos:
        ax.scatter(*pos, marker='.', s=50)

#     # TA star centered
    ax = axes[2]
    ax.set_title("Step 3\n" + "TA star centered")
    # plot the final TA before the offset is applied
    aper = aper_dict['coro']
    ax.scatter(*aper.idl_to_det(*idl_coords[0]['position']),
               c='k', label='ACQ', marker='x', s=100)
    ax.scatter(*aper.idl_to_det(*idl_coords[1]['position']),
               c='k', label='SCI', marker='*', s=100)
    other_pos = [ta_aper.idl_to_det(*i['position']) for i in ta_sequence[ta_aper_id][2:]]
    for pos in other_pos:
        ax.scatter(*pos, marker='.', s=50)

    # Offset applied
    ax = axes[3]
    ax.set_title("Step 4\n" + "Offset applied")
    # apply the offset to the position
    aper  = aper_dict['coro']
    acq_label = ta_sequence['slew'][0]['label']
    acq_pos = aper.idl_to_det(*ta_sequence['slew'][0]['position'])
    sci_label = ta_sequence['slew'][1]['label']
    sci_pos = aper.idl_to_det(*ta_sequence['slew'][1]['position'])
    ax.scatter(*acq_pos, 
               c='k', label=f'ACQ/{acq_label}', marker='x', s=100)
    ax.scatter(*sci_pos, 
               c='k', label=f'SCI/{sci_label}', marker='*', s=100)
    for pos in ta_sequence['slew'][2:]:
        ax.scatter(*aper.idl_to_det(*pos['position']),
                   marker='.', s=50)

    for ax in axes:
        # plot customizations
        # ax.label_outer()
        ax.set_aspect('equal')
        ax.grid(True, ls='--', c='grey', alpha=0.5)

    return fig


def plot_sky_ta_sequence(aper_dict, star_positions, v3pa, offset, colors, axes=None):

    if axes is None:
        nrows = 1
        ncols = 4
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                                 layout='constrained',
                                 sharex=True, sharey=True)
    else:
        fig = axes[0].get_figure()

    targ_label= f"{star_positions[0]['label']} --> {star_positions[1]['label']}"
    fig.suptitle(targ_label)

    for ax in axes.ravel():
        acq_pos = (star_positions[0]['position'].ra.deg, star_positions[0]['position'].dec.deg)
        sci_pos = (star_positions[1]['position'].ra.deg, star_positions[1]['position'].dec.deg)
        other_pos = [(star['position'].ra.deg, star['position'].dec.deg) for star in star_positions[2:]]
        ax.scatter(*acq_pos,
                   c='k', label=f"ACQ/{star_positions[0]['label']}", marker='x', s=100)
        ax.scatter(*sci_pos,
                   c='k', label=f"SCI/{star_positions[1]['label']}", marker='*', s=100)
        for pos in other_pos:
            ax.scatter(*pos, marker='.', s=50)


    # We start TA in the outer TA region
    ax = axes[0]
    ax.set_title(f"Step 1\nUR TA region")

    # center the attitude matrix at the Outer TA ROI
    attmat = create_attmat(star_positions[0]['position'], aper_dict['UR'], v3pa)
    formatting = dict(c=colors[0], alpha=1, ls='-')
    plot_apers(ax, attmat, aper_dict, formatting)


    # Continue to step 2 of TA, in the inner TA region
    ax = axes[1]
    ax.set_title(f"Step 2\nCUR TA region")

    # center the attitude matrix at the Inner TA ROI
    attmat = create_attmat(star_positions[0]['position'], aper_dict['CUR'], v3pa)
    formatting = dict(c=colors[1], alpha=1, ls='-')
    plot_apers(ax, attmat, aper_dict, formatting)



    # plot the final TA before the offset is applied
    ax = axes[2]
    ax.set_title("Step 3\nCentered")

    # center the attitude matrix on the coronagraph reference position
    attmat = create_attmat(star_positions[0]['position'], aper_dict['coro'], v3pa)
    formatting = dict(c=colors[2], alpha=1, ls='-')
    plot_apers(ax, attmat, aper_dict, formatting)



    # Plot the apertures and sources after the offset slew
    ax = axes[3]
    # compute the ra, dec of the offset from the TA position
    attmat = create_attmat(star_positions[0]['position'], aper_dict['coro'], v3pa)
    aper_dict['coro'].set_attitude_matrix(attmat)
    ra, dec = aper_dict['coro'].idl_to_sky(*(-offset))
    tel_sky = SkyCoord(ra=ra, dec=dec, unit='deg', frame='icrs')
    # compute the new attitude matrix at the slew position
    attmat = create_attmat(tel_sky, aper_dict['coro'], v3pa)
    ax.set_title(f"Step 4\nOffset applied")#\nTel-Targ sep: {tel_sky.separation(star_positions['SCI']).to('mas'):0.2e}")
    formatting = dict(c=colors[3], alpha=1, ls='-')
    plot_apers(ax, attmat, aper_dict, formatting)



    for ax in axes:
        # plot customizations
        ax.set_ylabel("Dec [deg]")
        ax.set_xlabel("RA [deg]")
        # ax.label_outer()
        ax.set_aspect("equal") 
        ax.grid(True, ls='--', c='grey', alpha=0.5)    
        # fix x-axis labels
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
    # RA increases to the left
    # they share axes, so you just need to invert one of them
    axes.flat[0].invert_xaxis()

    return fig


def plot_observing_sequence(
        aper_dict : dict,
        ta_sequence : dict,
        v3pa : float,
        idl_coords : list,
        star_positions : list,
        offset : np.ndarray,
) -> mpl.figure.Figure :
    """
    Plot the observing sequence both from the detector POV and the sky POV.

    Parameters
    ----------
    aper_dict : dict
      a dictionary of the various Siaf apertures to plot
    ta_sequence : dict
      a dictionary containing the IDL coordinates of the targets in the science
      aperture at each stage of the TA sequence
    v3pa : float
      the v3pa angle of the aperture at the time of observation
    idl_coords : ??
      ??
    star_positions : list
      a list of sky positions for each star
    offset : np.ndarray[float]
      The x, y IDL offset that will be applied after TA

    Output
    ------
    fig : mpl.figure.Figure
      A figure containing the plot of each step of the TA sequence.
      The top row is in the detector frame, the bottom row is in the sky frame

    """
    nrows = 2 # top row: detector, bottom row: sky
    ncols = 4 # outer, inner, centered, slew
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                             layout='constrained',
                             sharex='row', sharey='row')

    # plot the sequence from the POV of the detector
    fig = plot_detector_ta_sequence(aper_dict, ta_sequence, idl_coords, star_positions, axes[0])
    # plot the sequence from the POV of the sky
    colors = mpl.cm.plasma(np.linspace(0.2, 0.9, 4))
    fig = plot_sky_ta_sequence(aper_dict, star_positions, v3pa, offset, colors, axes[1])

    return fig

def plot_sky_ta_sequence_one_axis(aper_dict, star_positions, v3pa, offset, colors):
    """Plot the TA sequence on the sky, all on one axis"""
    nrows = 1
    ncols = 1
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 5), layout='constrained')

    targ_label = f"{star_positions[0]['label']} --> {star_positions[1]['label']}"
    fig.suptitle(f"{targ_label}\nTA sequence, in RA/Dec")

    colors = mpl.cm.plasma_r(np.linspace(0.2, 0.9, 4))

    acq_pos = (star_positions[0]['position'].ra.deg, star_positions[0]['position'].dec.deg)
    sci_pos = (star_positions[1]['position'].ra.deg, star_positions[1]['position'].dec.deg)
    other_pos = [(star['position'].ra.deg, star['position'].dec.deg) for star in star_positions[2:]]
    ax.scatter(*acq_pos,
               c='k',
               label=f"ACQ/{star_positions[0]['label']}",
               marker='x',
               s=100)
    ax.scatter(*sci_pos,
               c='k',
               label=f"SCI/{star_positions[1]['label']}",
               marker='*',
               s=100)
    for pos in other_pos:
        ax.scatter(*pos,
                   marker='.',
                   s=50)


    # We start TA in the outer TA region
    attmat = create_attmat(star_positions[0]['position'], aper_dict['UR'], v3pa)
    formatting = dict(c=colors[0], alpha=1, ls='dotted')
    plot_apers(ax, attmat, aper_dict, formatting)
    ax.plot([], [], 
            **formatting,
            label='Step 1: Outer TA step')


    # Continue to step 2 of TA, in the inner TA region
    attmat = create_attmat(star_positions[0]['position'], aper_dict['CUR'], v3pa)
    formatting = dict(c=colors[1], alpha=1, ls='dashdot')
    plot_apers(ax, attmat, aper_dict, formatting)
    ax.plot([], [], 
            **formatting,
            label='Step 2: Inner TA step')    


    # plot the final TA before the offset is applied
    # the telescope is now pointing the center of the coronagraph at the TA star
    attmat = create_attmat(star_positions[0]['position'], aper_dict['coro'], v3pa)
    formatting = dict(c=colors[2], alpha=1, ls='dashed')
    plot_apers(ax, attmat, aper_dict, formatting)
    ax.plot([], [],
            **formatting,
            label='Step 3: Before offset')


    # the telescope now places the TA star at the commanded offset
    # note that you must CHANGE THE SIGN OF THE OFFSET to get the position of the reference point
    ra, dec = aper_dict['coro'].idl_to_sky(*(-offset))
    tel_sky = SkyCoord(ra=ra, dec=dec, unit='deg', frame='icrs')
    attmat = create_attmat(tel_sky, aper_dict['coro'], v3pa)

    formatting = dict(c=colors[3], alpha=1, ls='solid')
    plot_apers(ax, attmat, aper_dict, formatting)
    ax.plot([], [],
            **formatting,
            label='Step 4: After offset')


    # plot formatting
    ax.set_ylabel("Dec [deg]")
    ax.set_xlabel("RA [deg]")
    ax.set_aspect("equal")
    ax.grid(True, ls='--', c='grey', alpha=0.5)
    # fix x-axis labels
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
    ax.legend()
    # RA increases to the left
    ax.invert_xaxis()

    return fig


def make_plots(
        aper_dict : dict,
        star_positions : list,
        v3pa : float,
        idl_coords : list,
        offset : np.ndarray,
) -> None:
    """
    Generate the plots

    Parameters
    ----------
    aper_dict : dict,
      a dict of the SIAF apertures to plot
    star_positions : dict,
      dictionary of the sky coordinates of the ACQ and SCI targets
    idl_coords : list,
      list of dicts of the idl coordinates of each star in the list, at each phase of the TA sequence
    offset : np.ndarray,
      the amount of the slew if you want to draw it I guess
    Output
    ------
    Displays 4 plots

    """

    figures = []

    fig1, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6), layout='constrained')
    fig1 = plot_before_offset_slew(aper_dict, idl_coords, star_positions, ax=axes[0])
    fig1 = plot_after_offset_slew(aper_dict, idl_coords, offset, star_positions, ax=axes[1])
    figures.append(fig1)

    # Plot 2: The TA sequence on the detector
    # compute the positions in IDL coordinates at each step of the sequence
    ta_sequence = {}
    for aper_id in ['UR', 'CUR', 'coro']:
        aper = aper_dict[aper_id]
        ta_sequence[aper_id] = sky_to_idl(star_positions, 
                                          aper, 
                                          v3pa)
    ta_sequence['slew'] = [{'label': i['label'], 'position': np.array(i['position']) + offset} for i in ta_sequence['coro']]

    # Plot 2: The TA sequence in RA and Dec on a single plot
    colors = mpl.cm.plasma(np.linspace(0.2, 0.9, 4))
    fig2 = plot_sky_ta_sequence_one_axis(aper_dict, star_positions, v3pa, offset, colors)
    figures.append(fig2)

    # Plot 3: plot detector and sky POV on same figure
    fig3 = plot_observing_sequence(
        aper_dict, ta_sequence, v3pa, idl_coords,
        star_positions, offset
    )
    figures.append(fig3)

    # # now, show the plots in correct order
    for fig in figures[::-1]:
        fig.show()
    plt.show()



def quad_boundaries(aperture, kwargs={}):
    """
    Generate a polygon to plot the 4QPM quadrant boundaries. Stolen from the JWST
    coronagraphic visibility tool

    Parameters
    ----------
    aperture: a pysiaf.Siaf aperture for the 1065, 1140, or 1550 coronagraph
    kwargs : {} arguments to pass to Polygon

    Output
    ------
    mask : matplotlib.patches.Polygon object
    """

    y_angle = np.deg2rad(aperture.V3IdlYAngle)
    corners_x, corners_y = aperture.corners(to_frame='idl')
    min_x, min_y = np.min(corners_x), np.min(corners_y)
    max_x, max_y = np.max(corners_x), np.max(corners_y)

    width_arcsec = 0.33
    x_verts0 = np.array([
        min_x,
        -width_arcsec,
        -width_arcsec,
        width_arcsec,
        width_arcsec,
        max_x,
        max_x,
        width_arcsec,
        width_arcsec,
        -width_arcsec,
        -width_arcsec,
        min_x
    ])
    y_verts0 = np.array([
        width_arcsec,
        width_arcsec,
        max_y,
        max_y,
        width_arcsec,
        width_arcsec,
        -width_arcsec,
        -width_arcsec,
        min_y,
        min_y,
        -width_arcsec,
        -width_arcsec
    ])
    x_verts = np.cos(y_angle) * x_verts0 + np.sin(y_angle) * y_verts0
    y_verts = -np.sin(y_angle) * x_verts0 + np.cos(y_angle) * y_verts0

    verts = np.concatenate([x_verts[:, np.newaxis], y_verts[:, np.newaxis]], axis=1)
    mask = Polygon(verts, alpha=0.5, **kwargs)
    return mask


if __name__ == "__main__":

    ###############################
    ########## USER INPUT #########
    ###############################

    # Star positions - make sure to enter all values.
    # The coordinates should be given *at the time of observation*
    # By default, we expect the coordinates to be given in RA, Dec system,
    # in units of degrees, and in the ICRS frame.
    # Users who are comfortable with astropy.coordiantes.SkyCoord may set them
    # as they like.
    # The `label` is used for identifying each component in plots and print output

    # The "slew_to" variable stores the position of the final target of the observations.
    # The "slew_from" variable stores the position of the star that is used for TA.
    slew_from = { # HD 167855 A
        'label': 'A',
        'position': SkyCoord( 
            272.81369648, 69.25014279,
            unit='deg',
            frame='icrs',
        )
    }
    slew_to = { # HD 167855 B
        'label': 'B',
        'position': SkyCoord( 
            272.81221024, 69.24905315,
            unit='deg',
            frame='icrs',
        )
    }

    # Telescope V3PA
    # enter the PA angle of the *telescope* V3 axis, at the time of the observation
    v3pa = 320.074

    # Choose a coronagraph by assigning one of the following to `coron_id`:
    # 1065, 1140, 1550, LYOT
    coron_id = '1550'

    # Print output
    verbose = True

    # Plotting - set to False if you don't want to show plots
    show_plots = False

    ###############################
    ####### END USER INPUT ########
    ###############################

    offsets = compute_offsets(
        slew_from, slew_to, v3pa, coron_id,
        verbose=verbose,
        show_plots=show_plots,
        plot_full = True,
        return_offsets= True
    )

