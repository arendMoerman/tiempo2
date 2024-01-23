"""!
@file
File containing templates for input dictionaries.
"""

instrument = {
        "freqs_filt"    : "Range of frequencies of filterbank in GHz.",
        "R"             : "Resolving power f / df.",
        "eta_inst"      : "Efficiency of entire chip.",
        "freq_sample"   : "Readout frequency in Hertz.",
        }

telescope = {
        "Dtel"          : "Diameter of telescope in meters.",
        "Ttel"          : "Temperature of telescope in Kelvin.",
        "Tgnd"          : "Temperature of ground around telescope in Kelvin.",
        "eta_ap"        : "Aperture efficiency of telescope, as function of instrument frequencies. If a single number is given, assume same aperture efficiency across entire frequency range.",
        "eta_mir"       : "Mirror efficiency.",
        "eta_fwd"       : "Front-to-back efficiency.",
        "chop_mode"     : "How to chop. Can choose 'none', 'direct', 'abba'."
        "freq_chop"     : "Chopping frequency in Hertz. If None, no chopping.",
        "dAz_chop"      : "Angular separation between chopping paths.",
        }

atmosphere = {
        "Tatm"          : "Temperature of atmosphere in Kelvin.",
        "filename"      : "Name of file containing ARIS screen.",
        "path"          : "Path to ARIS file",
        "dx"            : "Gridsize of ARIS screen along x-axis in meter.",
        "dy"            : "Gridsize of ARIS screen along y-axis in meter.",
        "h_column"      : "Reference height of atmospheric column.",
        "v_wind"        : "Windspeed in meters per second.",
        "PWV0"          : "Mean PWV value in millimeters.",
        }

SZsource = {
        "type"          : "Type of source (SZ).",
        "Az"            : "Azimuthal lower and upper limits of source map in degrees.",
        "El"            : "Elevation lower and upper limits of source map in degrees.",
        "nAz"           : "Number of Azimuth points.",
        "nEl"           : "Number of Elevation points.",
        # MockSZ specific
        "Te"            : "Electron temperature of cluster gas in Kev.",
        "ne0"           : "Central electron density in # per square centimeter.",
        "beta"          : "Isothermal-beta structure coefficient.",
        "v_pec"         : "Peculiar cluster velocity, relative to CMB, in kilometers per second.",
        "rc"            : "Cluster core radius in kiloparsec.",
        "Da"            : "Angular diameter distance in megaparsec.",
        "freqs_src"     : "Range of frequencies over which to simulate source signal, in GHz.",
        }

Galsource = {
        "type"          : "Type of source (GalSpec).",
        "Az"            : "Azimuthal lower and upper limits of source map in degrees.",
        "El"            : "Elevation lower and upper limits of source map in degrees.",
        "nAz"           : "Number of Azimuth points.",
        "nEl"           : "Number of Elevation points.",
        # GalSpec specific
        "lum"           : "Luminosity in log(L_fir/L_sol).",
        "z"             : "redshift of galaxy.",
        "f_lo"          : "Lower frequency limit.",
        "f_hi"          : "Upper frequency limit.",
        "nfreqs"        : "Number of frequency bins, linearly spaced.",
        "lwidth"        : "Linewidth of spectral lines in km/s.",
        "COlines"       : "Kamenetzky or Rosenberg.",
        "lines"         : "Bonato or Spinoglio.",
        "mollines"      : "T/F, add molecular lines."
        }

load_source = {
        "path"          : "Path to saved source datacube.",
        "filename"      : "Name of saved source.",
        }

simparams = {
        "name_sim"      : "Name of simulation.",
        "t_obs"         : "Total observation time in seconds.",
        "nThreads"      : "Number of CPU threads to use.",
        "outDir"        : "Output directory of simulation. Automatically uses name_sim as name of file.",
        }
