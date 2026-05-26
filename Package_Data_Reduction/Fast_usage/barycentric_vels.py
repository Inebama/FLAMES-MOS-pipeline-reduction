import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.coordinates import AltAz
from astropy.coordinates import solar_system_ephemeris
from astropy.coordinates import get_body_barycentric_posvel
from astropy.coordinates import EarthLocation, SkyCoord

from astropy import units as u
from astroquery.simbad import Simbad

# --------------------------------------------------------------
# Function to query object coordinates from SIMBAD
# --------------------------------------------------------------
def get_target_coords(target_name):
    """
    Query SIMBAD for a target and return its ICRS coordinates as an Astropy SkyCoord object.
    No aliases are used; fails cleanly if the target is not found.
    """

    # Ensure RA/DEC in degrees (new standard field names in astroquery)
    custom_simbad = Simbad()
    custom_simbad.add_votable_fields('ra', 'dec')

    try:
        result = custom_simbad.query_object(target_name)
        if result is None:
            raise ValueError(f"Target '{target_name}' not found in SIMBAD.")

        # Handle both possible field name conventions
        for ra_key, dec_key in [('RA_d', 'DEC_d'), ('RA', 'DEC'), ('ra', 'dec')]:
            if ra_key in result.colnames and dec_key in result.colnames:
                ra = float(result[ra_key][0]) * u.deg
                dec = float(result[dec_key][0]) * u.deg
                return SkyCoord(ra=ra, dec=dec, frame='icrs')

        # If the RA/DEC fields are not present
        raise ValueError(f"RA/DEC fields not found in SIMBAD result for '{target_name}'.")

    except Exception as e:
        raise ValueError(f"Could not resolve target '{target_name}' in SIMBAD. Error: {e}")


# --------------------------------------------------------------
# Function to compute barycentric velocity correction
# --------------------------------------------------------------
# def barycentric_velocity(target_coord, times, location):

#     with solar_system_ephemeris.set('builtin'):
#         # Get observer's position and velocity relative to the Solar System barycenter
#         pos, vel = get_body_barycentric_posvel('earth', times)
#         obs_pos, obs_vel = location.get_gcrs_posvel(times)
#         total_vel = vel + obs_vel[1]  # Earth's + observatory's velocity

#     # Unit vector in target direction
#     direction = target_coord.icrs.represent_as('cartesian').get_xyz().value
#     direction /= np.linalg.norm(direction)

#     # Dot product gives projected velocity (positive = moving away)
#     barycorr = np.array([np.dot(total_vel[i].xyz.to(u.m/u.s).value, direction) for i in range(len(times))]) * u.m/u.s

#     return barycorr.to(u.km/u.s)

#full barycentric correction
def barycentric_velocity(target_coord, times, location):
    return target_coord.radial_velocity_correction(obstime=times, location=location)

# --------------------------------------------------------------
# Main function display and get several points at once
def plot_barycentric_velocity(
    target_name,
    mjd_range=None,
    n_points=100,
    mjd_array=None,
    observatory="Paranal",
    savefile=None,
    return_values=False,
):
    # Observatory list (manually define coordinates if not available)
    observatories = {
        "paranal": EarthLocation.of_site("Paranal Observatory"),
        "lasilla": EarthLocation.of_site("La Silla Observatory"),
        "ctio": EarthLocation.of_site("cerro tololo interamerican observatory"),  # correct lowercase name
        "campana": EarthLocation.of_site("Las Campanas Observatory"),
    }

    if observatory not in observatories:
        raise ValueError(f"Unknown observatory '{observatory}'. Choose one of: {list(observatories.keys())}")

    location = observatories[observatory]

    # Get coordinates of target
    target_coord = get_target_coords(target_name)

    # Define time sampling
    if mjd_array is not None:
        times = Time(mjd_array, format='mjd')
    elif mjd_range is not None:
        t_start, t_end = mjd_range
        times = Time(np.linspace(t_start, t_end, n_points), format='mjd')
    else:
        raise ValueError("You must provide either mjd_range=(start, end) or mjd_array=[...]")

    # Compute barycentric velocities
    velocities = barycentric_velocity(target_coord, times, location)

    # Plot
    plt.figure(figsize=(8, 5))
    plt.plot(times.mjd, velocities.value, 'o-', ms=3)
    plt.xlabel("MJD")
    plt.ylabel("Barycentric velocity correction (km/s)")
    plt.title(f"Barycentric velocity variation for {target_name} ({observatory})")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Save file
    if savefile:
        with open(savefile, 'w') as f:
            f.write(f"# Target: {target_name}\n")
            f.write(f"# Coordinates (ICRS): RA={target_coord.ra.deg:.6f} deg, DEC={target_coord.dec.deg:.6f} deg\n")
            f.write(f"# Observatory: {observatory}\n")
            f.write("# Columns: MJD   BarycentricVelocity_km/s\n")
            for t, v in zip(times.mjd, velocities.value):
                f.write(f"{t:.6f} {v:.6f}\n")
        print(f"Saved data to: {savefile}")
    if return_values:
        return times.mjd,velocities.value




# Compute orbital and rotational contributions separately
def compute_components(target_coord, times, location):
    """
    Compute the orbital (Earth barycentric), rotational (observatory), and total
    barycentric velocity components projected along the line of sight.
    Returns three arrays (v_orb, v_rot, v_total) in km/s.
    """
    from astropy.coordinates import solar_system_ephemeris, get_body_barycentric_posvel
    from astropy import units as u
    import numpy as np

    # Earth's barycentric velocity (orbital)
    with solar_system_ephemeris.set('builtin'):
        earth_pos, earth_vel = get_body_barycentric_posvel('earth', times)
    # Convert from AU/day to m/s and ensure shape (N,3)
    earth_vel_si = (earth_vel.xyz.to(u.m/u.s)).T.value  # shape (N,3)

    # Observatory velocity (rotation)
    gcrs_pos, gcrs_vel = location.get_gcrs_posvel(times)
    obs_vel_si = (gcrs_vel.xyz.to(u.m/u.s)).T.value  # shape (N,3)

    # Unit vector to target in ICRS
    uv = target_coord.icrs.represent_as('cartesian').get_xyz().value
    uv = uv / np.linalg.norm(uv)
    
    # Project velocities onto line of sight
    v_orb = np.dot(earth_vel_si, uv)
    v_rot = np.dot(obs_vel_si, uv)
    v_total = v_orb + v_rot

    v_orb = (v_orb * u.m / u.s).to(u.km / u.s)
    v_rot = (v_rot * u.m / u.s).to(u.km / u.s)
    v_total = (v_total * u.m / u.s).to(u.km / u.s)

    return v_orb, v_rot, v_total



import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.coordinates import EarthLocation, GCRS, SkyCoord
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric_posvel
from astropy import units as u

def components_plot(
    target_name="Alnitak",
    observatory="Paranal",
    mjd_start=60000.0,
    mjd_end=60001.0,
    n_points=400
):
    """
    Extended barycentric-component plot:
    - classical components projected on fixed ICRS direction
    - same components but projecting onto time-dependent GCRS/apparent direction
    - astropy full radial_velocity_correction
    - residuals showing what is missing (precession/nutation/aberration/relativity)
    """

    # --- Location ---
    if observatory.lower() == "paranal":
        loc = EarthLocation.of_site("Paranal Observatory")
    else:
        loc = EarthLocation.of_site(observatory)

    # --- Times and target ---
    times = Time(np.linspace(mjd_start, mjd_end, n_points), format='mjd')
    target = get_target_coords(target_name)   # assumes your existing simbad resolver

    # --- Earth barycentric and observatory GCRS velocities ---
    with solar_system_ephemeris.set('builtin'):
        _, earth_vel = get_body_barycentric_posvel('earth', times)

    # get_gcrs_posvel returns (pos, vel). vel.xyz is a (3, N) quantity (units m/s)
    gcrs_pos, gcrs_vel = loc.get_gcrs_posvel(times)

    # Normalize shapes to arrays of shape (N,3) in SI (m/s)
    earth_vel_si = earth_vel.xyz.to(u.m / u.s).T.value      # shape (N,3)
    obs_vel_si = gcrs_vel.xyz.to(u.m / u.s).T.value         # shape (N,3)

    # --- Unit direction vectors ---
    # ICRS (fixed direction, single vector)
    uv_icrs = target.icrs.represent_as('cartesian').get_xyz().value
    uv_icrs = np.array(uv_icrs).astype(float).ravel()
    uv_icrs = uv_icrs / np.linalg.norm(uv_icrs)

    # GCRS (time dependent apparent direction) -> transform target to GCRS at each obstime
    # This yields a SkyCoord with xyz arrays; we'll extract per-time unit vectors
    target_gcrs = target.transform_to(GCRS(obstime=times))
    # get_xyz gives shape (3, N)
    uv_gcrs_3xn = target_gcrs.cartesian.get_xyz().value
    # transpose to (N,3) and normalize each vector
    uv_gcrs = uv_gcrs_3xn.T.astype(float)
    norms = np.linalg.norm(uv_gcrs, axis=1)[:, None]
    uv_gcrs = uv_gcrs / norms

    # --- Projections (dot products) ---
    # Each projection is dot(velocity_vector, unit_vector) for each epoch (arrayized)
    v_earth_icrs = np.einsum('ij,j->i', earth_vel_si, uv_icrs)   # shape (N,)
    v_obs_icrs   = np.einsum('ij,j->i', obs_vel_si, uv_icrs)

    # time-dependent projections onto GCRS direction
    v_earth_gcrs = np.einsum('ij,ij->i', earth_vel_si, uv_gcrs)  # row-wise dot
    v_obs_gcrs   = np.einsum('ij,ij->i', obs_vel_si, uv_gcrs)

    # Classical sums
    v_classic_icrs = v_earth_icrs + v_obs_icrs
    v_classic_gcrs = v_earth_gcrs + v_obs_gcrs

    # --- Astropy full correction (most complete) ---
    # radial_velocity_correction returns Quantity (units: km/s)
    v_astropy = target.radial_velocity_correction(obstime=times, location=loc).to(u.km / u.s).value

    # Convert classical arrays (m/s) to km/s for plotting
    CONV = 1.0 / 1000.0
    v_earth_icrs_k = v_earth_icrs * CONV
    v_obs_icrs_k   = v_obs_icrs * CONV
    v_classic_icrs_k = v_classic_icrs * CONV

    v_earth_gcrs_k = v_earth_gcrs * CONV
    v_obs_gcrs_k   = v_obs_gcrs * CONV
    v_classic_gcrs_k = v_classic_gcrs * CONV

    # --- Residuals (in m/s for clarity) ---
    # convert astropy (km/s) to m/s for residuals
    v_astropy_ms = v_astropy * 1000.0

    # residuals between full astropy and the two classical approaches
    resid_icrs_ms = v_astropy_ms - (v_classic_icrs * 1.0)    # both in m/s
    resid_gcrs_ms = v_astropy_ms - (v_classic_gcrs * 1.0)

    # difference introduced by using time-dependent GCRS direction vs fixed ICRS direction
    delta_dir_ms = (v_classic_gcrs - v_classic_icrs) * 1.0   # in m/s

    # --- Plotting ---
    fig, axes = plt.subplots(4, 1, figsize=(11, 11), sharex=True,
                             gridspec_kw={'height_ratios':[2,1,1,1]})

    ax0 = axes[0]
    ax0.plot(times.mjd, v_astropy, label='Astropy full correction', lw=1.5)
    ax0.plot(times.mjd, v_classic_icrs_k, '--', label='Classical (ICRS direction)', lw=1)
    ax0.plot(times.mjd, v_classic_gcrs_k, ':', label='Classical (GCRS/apparent dir)', lw=1)
    ax0.set_ylabel('Velocity (km/s)')
    ax0.set_title(f"Barycentric components for {target_name} at {observatory}")
    ax0.legend()
    ax0.grid(True)

    # Second panel: split components (GCRS projection)
    ax1 = axes[1]
    ax1.plot(times.mjd, v_earth_gcrs_k, label='Earth barycentric (proj on GCRS)', lw=1)
    ax1.plot(times.mjd, v_obs_gcrs_k, label='Observer geocentric (proj on GCRS)', lw=1)
    ax1.set_ylabel('km/s')
    ax1.legend()
    ax1.grid(True)

    # Third panel: residuals (astropy - classical_gcrs) in m/s
    ax2 = axes[2]
    ax2.plot(times.mjd, resid_gcrs_ms, 'o-', ms=3)
    ax2.axhline(0, color='k', lw=0.6)
    ax2.set_ylabel('Astropy - classical (m/s)')
    ax2.grid(True)

    # Fourth panel: direction-change effect (GCRS vs ICRS) in m/s
    ax3 = axes[3]
    ax3.plot(times.mjd, delta_dir_ms, 'o-', ms=3, label='Projection change (GCRS - ICRS)')
    ax3.plot(times.mjd, (v_obs_gcrs - v_obs_icrs) * 1.0, 'x-', ms=3, label='Obs proj change (m/s)')
    ax3.axhline(0, color='k', lw=0.6)
    ax3.set_ylabel('m/s')
    ax3.set_xlabel('MJD')
    ax3.legend()
    ax3.grid(True)

    plt.tight_layout()
    plt.show()

    # --- Print summary amplitudes ---
    def amp(arr_ms):
        return np.nanmax(arr_ms) - np.nanmin(arr_ms)

    print("Amplitude summaries (peak-to-peak):")
    print(f"  Earth barycentric (GCRS, km/s): {np.ptp(v_earth_gcrs_k):.3f} km/s")
    print(f"  Observer geocentric (GCRS, km/s): {np.ptp(v_obs_gcrs_k):.3f} km/s")
    print(f"  Classical total (GCRS, km/s): {np.ptp(v_classic_gcrs_k):.3f} km/s")
    print(f"  Astropy full (km/s): {np.ptp(v_astropy):.3f} km/s")
    print(f"  Residual (astropy - classical_gcrs) peak-to-peak: {amp(resid_gcrs_ms):.1f} m/s")
    print(f"  Direction-change effect (GCRS-ICRS) peak-to-peak: {amp(delta_dir_ms):.1f} m/s")

    # Return arrays in case you want to save or inspect them programmatically
    out = dict(
        times=times,
        v_astropy=(v_astropy * u.km/u.s),
        v_classic_icrs=(v_classic_icrs_k * u.km/u.s),
        v_classic_gcrs=(v_classic_gcrs_k * u.km/u.s),
        v_earth_gcrs=(v_earth_gcrs_k * u.km/u.s),
        v_obs_gcrs=(v_obs_gcrs_k * u.km/u.s),
        resid_gcrs=(resid_gcrs_ms * u.m/u.s),
        delta_dir=(delta_dir_ms * u.m/u.s),
    )
    return out



# --------------------------------------------------------------
# Example usage
# --------------------------------------------------------------
if __name__ == "__main__":
    # Example getting Baricentric just
    # plot_barycentric_velocity(
    #     target_name="47 Tuc",
    #     mjd_range=(60000, 60365),
    #     n_points=1200,
    #     observatory="paranal",
    #     savefile=None,
    #     return_values=True
    # )

    components_plot(
        target_name="47 Tuc",
        observatory="Paranal",
        mjd_start=60000.0,
        mjd_end=60365,
        n_points=4000) # decompose rotation and translational contribution