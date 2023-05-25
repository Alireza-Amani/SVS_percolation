# <<< doc >>> -----------------------------------------------------------------
"""
This is where we provide:
    - required data or path to data
    - values for parameters that describe the system
    - values for parameters that are related to methods/techniques used in the
        article

All of the variables will be encapsulated in a dataclass.


"""
# _________________________________________________________________ <<< doc >>>

# <<< imports >>> -------------------------------------------------------------
from pathlib import Path
from dataclasses import dataclass, field
from collections import OrderedDict
from copy import deepcopy
import pandas as pd
import numpy as np


from svspaper_package.svs.svs_model import SVSModel
from svspaper_package.svs.prep_svs import ModelInputData
from svspaper_package.svs.perturb_run import PerturbAndRun
from svspaper_package.helpers.help_functions import (
    fill_hyprop_dict, fit_swrc, find_best_sample, prec_cr, calc_bias_std,
    get_ar1_param, create_perturbed_data, save_csv_as_met, read_field_meteo,
    read_main_met_source, create_enclosure_dataframe, calc_specific_humidity,
    create_ensemble_scenarios
)
from svspaper_package.svs.helper_functions import assert_file_exist



# _____________________________________________________________ <<< imports >>>

# <<< variable definition >>> -----------------------------------------------

# number of ensemble members
N_ENS = 5

# number of CPU cores to use
N_CORES = 5

# random seed
RANDOM_SEED = 1915


paths_dict = {
    # path to the HYPROP test data
    "hyprop_data": Path("./data/hyprop_data.csv"),

    # path to field weather station data
    "field_weather": Path("./data/average_field_stations.csv"),

    # path to the main meteo data that SVS uses
    "main_meteo": Path("./data/saint_g_ERA5.csv"),

    # folder to save the perturbed dataframes
    "met_ens": Path("./output/perturbed_met_files"),

    # a folder to host the SVS runs
    "svs_runs": Path("./output/svs_runs"),

    # path to the SVS executable
    'svs_executable': Path("./SVS_source_code/svs_exec")

}


# soil types information
soil_types_dict = {
    "CM": {  # CM: cover material
        # template name for the soil volumetric WC columns in HYPROP data
        "vol_swc_col": "CM_smpl#_VWC_percent",
        "suction_col": "CM_smpl#_suction_kPa",
        "parameters": {
            'sand': 66.2,
            'clay': 7.7,
            'ksat': 2.9e-5,
            'rhosoil': 1575.0,  # kg/m3
            # other params are inferred from the HYPROP data
        },
        "limits": OrderedDict(

            sand=[55.2, 78.5],
            clay=[5.2, 11.7],
            wsat=[0.292, 0.327],
            wfc=[0.047, 0.095],
            wwilt=[0.02, 0.063],
            user_wfcdp=[0.5, 0.99],  # of Wsat
            ksat=[7.00E-06, 5.40E-05],
            rhosoil=[1464, 1686],
        )

    },

    "AB": {
        "vol_swc_col": "AB_smpl#_VWC_percent",
        "suction_col": "AB_smpl#_suction_kPa",
        "parameters": {
            'sand': 26.3,
            'clay': 10.1,
            'ksat': 1.5e-6,
            'rhosoil': 1590.0,  # kg/m3
        },
        "limits": OrderedDict(

            user_wfcdp=[0.5, 0.99],  # of Wsat
            ksat=[4.05E-07, 4.61E-06],
            rhosoil=[1540, 1646]
        )

    },

    "BC": {
        "vol_swc_col": "BC_smpl#_VWC_percent",
        "suction_col": "BC_smpl#_suction_kPa",
        "parameters": {
            'sand': 84.3,
            'clay': 3.6,
            'ksat': 1.8e-5,
            'rhosoil': 1635.0,  # kg/m3
        },
        "limits": OrderedDict(

            user_wfcdp=[0.5, 0.99],  # of Wsat
            ksat=[4.40E-06, 4.93E-05],
            rhosoil=[1555, 1716]
        )
    }
}

# name of the soil parameters to be set in the MEASH_parameters.txt file
soil_params = [
    "sand", "clay", "ksat", "wsat", "wfc", "wwilt",
    'bcoef', 'psisat', "rhosoil"
]


enclosures_info = {

    "E1": {
        "layering": "15:2.5:CM + 160:5:CM + 15:2.5:CM"
    },

    "E2": {
        "layering": "15:2.5:CM + 15:5:CM + 80:5:BC + 65:5:CM + 15:2.5:CM"
    },

    "E3": {
        "layering": "15:2.5:CM + 95:5:AB + 65:5:CM + 15:2.5:CM"
    }
}
# guide on how to create the layering string:
# `15:2.5:CM` means 15 cm of cover material divided into 2.5 cm layers
# `+` means add another block of soil
# `95:5:AB` means 95 cm of AB soil divided into 5 cm layers

# parameters describing the site, identical for all enclosures
site_params = {
    "deglat": 45.82,
    "deglng": -72.37,
    "slop": 0.02,
    "zusl": 10.0,
    "ztsl": 1.5,
    "observed_forcing": "height",
    "draindens": 0.0,
    "vf_type": 13  # vegetation type: 13 = grass
}

# SVS internal parameters
svs_params = {
    "SCHMSOL": "SVS",
    "lsoil_freezing_svs1": ".TRUE.",
    "soiltext": "NIL",
    "read_user_parameters": 1,
    "save_minutes_csv": 0,
    "water_ponding": 0,
    "KFICE": 0,
    "OPT_SNOW": 2,
    "OPT_FRAC": 1,
    "OPT_LIQWAT": 1,
    "WAT_REDIS": 0,

    "tperm": 280.15,  # deep soil temperature
    "user_wfcdp": 0.15  # breakthroug water content for lower boundary condition
}

# ensure that Wfc is always greater than Wwilt by this amount; for CM soil type
CM_FC_WILT_DELTA = 0.025

# suction at which the soil is considered at field capacity
FIELD_CAP_SUCTION = 33  # kPa

# treat air entry suction as an empirical parameter for curve fitting
EMP_AIR_ENTRY = False

# meteorological variables to compare between the field weather station data
# and the main meteo data
# name of the variables must be the same in both datasets: check whether the
# variable names are the same in the field weather station data and the main
# meteo data
# `additive` indicates whether the perturbation is additive or multiplicative
compare_vars = {

    "Air temperature (degC)": {
        "additive": True
    },
    "Dew point temperature (degC)": {
        "additive": True
    },
    "Relative humidity (%)": {
        "additive": False,
        "var_max": 100
    },
    "Atmospheric pressure (kPa)": {
        "additive": True
    },
    "Shortwave radiation (W.m-2)": {
        "additive": False
    },

    "Wind speed (m.s-1)": {
        "additive": False
    },
}

# variables that are not perturbed
not_perturbed = ["Longwave radiation (W.m-2)", 'Wind speed (m.s-1)']

# perturbation factor for precipitation - Uniform [0.95, 1.05]
PREC_FACT = 0.05

# mapping of the required meteo variables to the labels in the `main_meteo.csv`
# needed for the `save_csv_as_met` function
meteo_cols = {
    "utc_dtime": "datetime_utc",
    "air_temperature": "Air temperature (degC)",
    "precipitation": "Precipitation (mm) corrected",
    # "precipitation": "Precipitation (mm)",
    "wind_speed": "Wind speed (m.s-1)",
    "atmospheric_pressure": "Atmospheric pressure (Pa)",
    "shortwave_radiation": "Shortwave radiation (W.m-2)",
    "longwave_radiation": "Longwave radiation (W.m-2)",
    "specific_humidity": "Specific humidity (kg.kg-1)",
    "relative_humidity": "Relative humidity (%)"
}
# Note: if you use `Precipitation (mm) corrected` as the precipitation column,
# it means that you want the precipitation to be corrected for the wind

# label for the datetime column in the `field_weather.csv`
FIELD_WEATHER_DT_LABEL = "datetime_local"

# tz of the field weather station data
FIELD_WEATHER_DT_TZ = "UTC-4"

# tz of the main meteo data
MAIN_MET_TZ = "UTC+0"

# label for the dew point temp in the `main_meteo.csv`
MAIN_METEO_DEW_POINT = "Dew point temperature (degC)"

# simulation start_date
START_DATE = "2018-182-05-00"

# spinup end date
SPINUP_END_DATE = "2019-07-01 00:00:00"

# spinup date timezone
SPINUP_END_DATE_TZ = "America/Montreal"

# ___________________________________________________________________ <<< >>>

# <<< main >>> ----------------------------------------------------------------
# create the MOAI dataclass


@dataclass
class MOAI:
    """
    This is the dataclass that contains all the information required for the
    simulations.
    """

    # number of ensemble members
    n_ens: int = N_ENS

    # number of CPUs to use
    n_cpus: int = N_CORES

    # dictionary of paths
    paths_dict: dict = field(default_factory=lambda: deepcopy(paths_dict))

    # soil types information
    soil_types_dict: dict = field(
        default_factory=lambda: deepcopy(soil_types_dict))

    # suction at which the soil is considered at field capacity
    field_cap_suction: int = FIELD_CAP_SUCTION

    # treat air entry suction as an empirical parameter for curve fitting
    emp_air_entry: bool = EMP_AIR_ENTRY

    # meteorological variables to compare between the field weather station data
    # and the main meteo data
    compare_vars: dict = field(default_factory=lambda: deepcopy(compare_vars))

    # variables that are not perturbed
    not_perturbed: list = field(
        default_factory=lambda: deepcopy(not_perturbed))

    # perturbation factor for precipitation - Uniform [1 - prec_fact, 1 + prec_fact]
    prec_fact: float = PREC_FACT

    # mapping of the required meteo variables to the labels in the `main_meteo.csv`
    meteo_cols: dict = field(default_factory=lambda: deepcopy(meteo_cols))

    # label for the datetime column in the `field_weather.csv`
    field_weather_datetime: str = FIELD_WEATHER_DT_LABEL

    # tz of the field weather station data
    field_weather_tz: str = FIELD_WEATHER_DT_TZ

    # tz of the main meteo data
    main_meteo_tz: str = MAIN_MET_TZ

    # label for the dew point temp in the `main_meteo.csv`
    main_meteo_dew_point: str = MAIN_METEO_DEW_POINT

    # svs parameters
    svs_params: dict = field(default_factory=lambda: deepcopy(svs_params))

    # site parameters
    site_params: dict = field(default_factory=lambda: deepcopy(site_params))

    # enclosure info
    enclosures_info: dict = field(
        default_factory=lambda: deepcopy(enclosures_info)
    )

    # soil parameters
    soil_params: dict = field(default_factory=lambda: deepcopy(soil_params))

    # ensure that Wfc is always greater than Wwilt by this amount; for CM soil type
    cm_fc_wilt_delta: float = CM_FC_WILT_DELTA

    # random seed
    random_seed: int = RANDOM_SEED

    # simulation start_date
    start_date: str = START_DATE

    # spinup end date
    spinup_end_date: str = SPINUP_END_DATE

    # spinup date timezone
    spinup_end_date_tz: str = SPINUP_END_DATE_TZ

    # a dict to store the bias/std of the meteo variables
    bias_std: dict = field(default_factory=dict)

    def __post_init__(self):

        # check if the path to the HYPROP data exists
        if self.paths_dict["hyprop_data"] is not None:
            assert_file_exist(self.paths_dict["hyprop_data"])


# a class representing the article
class SVSPaper:
    """
    This class represents the SVS paper; it can be used to run the simulations.
    """

    def __init__(self):
        self.all_info = MOAI()

        # list of the perturbed meteo dataframes
        self.met_pert_list = []

        # the field weather data
        self.df_field_weather = pd.DataFrame()

        # the main meteo data
        self.df_main_meteo = pd.DataFrame()

        # SVS instances
        self.svs_instances = {}

        # the base .met file path
        self.base_dotmetfile_path = Path("NoWhereYet")

        # the ensemble outputs
        self.ens_outputs = pd.DataFrame()

    def get_swrc_params(self):
        """
        Read HYPROP data and fit SWRC equation to the data to get the
        parameters of the SWRC equation and other related soil hydraulic
        properties.
        """

        # if no path to HYPROP data is provided, raise an error and let the user
        # know that this means that this method cannot be used
        if self.all_info.paths_dict["hyprop_data"] is None:
            raise Exception(
                "No path to HYPROP data is provided. This means that the "
                "method `get_swrc_params` cannot be used."
            )

        # fill the soil types dictionary with the HYPROP data
        self.all_info.soil_types_dict = fill_hyprop_dict(
            self.all_info.soil_types_dict,
            self.all_info.paths_dict["hyprop_data"]
        )

        # fit SWRC equation to the HYPROP data and get the parameters of the
        # SWRC equation and other related soil hydraulic properties
        self.all_info.soil_types_dict = fit_swrc(
            self.all_info.soil_types_dict,
            self.all_info.field_cap_suction,
            self.all_info.emp_air_entry
        )

        # find the best sample for each soil type: sort the samples based on
        # their MAE and R2
        self.all_info.soil_types_dict = find_best_sample(
            self.all_info.soil_types_dict
        )

        # assign to each soil types, the parameters of the best sample for that
        # soil type
        self.assign_best_sample_params()

    def assign_best_sample_params(self):
        """
        Assign to each soil types, the parameters of the best sample for that
        soil type.

        Best sample is the sample with the highest `sample_score` value in the
        `samples_info` dataframe stored as a key in the `soil_types_dict`.
        """

        # loop over the soil types
        for stype_dict in self.all_info.soil_types_dict.values():

            # get the `samples_info` dataframe
            samples_info = stype_dict["samples_info"].copy()

            if len(samples_info) > 1:

                # get the index of the row with the highest `sample_score` value
                best_sample_idx = samples_info["sample_score"].idxmax()

                # get the row with the highest `sample_score` value
                best_sample = samples_info.loc[best_sample_idx]

            elif len(samples_info) == 1:
                best_sample = samples_info.iloc[0]

            else:
                raise Exception(
                    "The `samples_info` dataframe is empty. This means that "
                    "no samples were found for this soil type."
                )

            # assign the best sample parameters to the soil type dictionary
            stype_dict["parameters"]["wsat"] = best_sample["Wsat_percent"] / 100
            stype_dict["parameters"]["psisat"] = np.round(best_sample["psi_ae_kPa"], 4)
            stype_dict["parameters"]["wwilt"] = best_sample["Wwilt_percent"] / 100
            stype_dict["parameters"]["wfc"] = best_sample["Wfc_percent"] / 100
            stype_dict["parameters"]["bcoef"] = np.round(best_sample["bcoef"], 2)

        return soil_types_dict

    def create_perturbed_metfiles(self, use_previous: bool = False):
        """
        Create perturbed meteo files for the ensemble members.
        """

        # a dictionary to store the bias and standard deviation of the difference between
        # the two dataframes for each variable
        bias_std = {}

        # a dictionary to store the values of the variables for each source and each
        # variable and the associated biases
        var_info = {}

        # read field station data
        self.df_field_weather = read_field_meteo(
            self.all_info.paths_dict["field_weather"],
            self.all_info.field_weather_datetime,
            self.all_info.field_weather_tz
        )

        # read saint-g station + ERA5 data
        self.df_main_meteo = read_main_met_source(
            # paths_dict["stg_met"], meteo_cols["utc_dtime"], "UTC+0"
            self.all_info.paths_dict["main_meteo"],
            self.all_info.meteo_cols["utc_dtime"],
            self.all_info.main_meteo_tz
        )

        # apply precipitation catch ratio correction
        self.df_main_meteo = prec_cr(
            self.df_main_meteo,
            u_lab=self.all_info.meteo_cols["wind_speed"],
            t_lab=self.all_info.meteo_cols["air_temperature"]
        )

        # add specific humidity column
        self.df_main_meteo[self.all_info.meteo_cols["specific_humidity"]] = (
            calc_specific_humidity(
                self.df_main_meteo.loc[:, self.all_info.meteo_cols["atmospheric_pressure"]] / 1000, # to kPa
                self.df_main_meteo.loc[:, self.all_info.main_meteo_dew_point],
            )
        )

        # save it in the `svs_runs` folder
        self.base_dotmetfile_path = (
            self.all_info.paths_dict["svs_runs"] / "basin_forcing.met"
        )
        save_csv_as_met(
            self.base_dotmetfile_path, self.all_info.meteo_cols, self.df_main_meteo
        )

        # get bias and standard deviation of the difference between the two dataframes
        bias_std, var_info = calc_bias_std(
            list(self.all_info.compare_vars),
            self.df_field_weather, self.df_main_meteo
        )

        # fitting the AR1 model to the bias of each variable
        # add a `phi` key to each dict nested in `bias_std`
        for var_name in bias_std:
            phi = get_ar1_param(var_info[var_name]['var_2'])
            bias_std[var_name]['phi'] = phi

        # own the `bias_std` dict to the `all_info` object
        self.all_info.bias_std = deepcopy(bias_std)

        # create perturbed dataframes
        for mem_i in range(N_ENS):

            # no need to create the perturbed dataframes if they already exist
            if use_previous:
                break

            df_pert = create_perturbed_data(
                self.df_main_meteo,
                self.all_info.compare_vars,
                self.all_info.not_perturbed,
                self.all_info.prec_fact, bias_std, mem_i,
                dew_point_label=self.all_info.main_meteo_dew_point,
                precipitation_label=self.all_info.meteo_cols["precipitation"]
            )

            # a fix for my specific case: should be handled in a more general way!
            # replace "Atmospheric pressure (kPa)" with "Atmospheric pressure (Pa)"
            df_pert[self.all_info.meteo_cols["atmospheric_pressure"]] = (
                df_pert['Atmospheric pressure (kPa)'] * 1000
            )
            df_pert.drop('Atmospheric pressure (kPa)', axis=1, inplace=True)

            # create a datetime column
            df_pert[self.all_info.meteo_cols["utc_dtime"]] = df_pert.index

            # save_path for dfmet
            save_path = self.all_info.paths_dict["met_ens"] / \
                f"basin_forcing_{mem_i}.met"
            save_csv_as_met(save_path, self.all_info.meteo_cols, df_pert)

            # # save the perturbed dataframe
            # df_pert.to_csv(
            #     self.all_info.paths_dict["met_ens"] / f"met_ens_{mem_i}.csv")

            # add the perturbed dataframe to the list
            self.met_pert_list.append(df_pert)

    def create_base_params_for_enclosures(self, create_svs=False):
        """
        Create the base parameters for the enclosures.

        Parameters
        ----------
        create_svs : bool, optional
            Whether to create an instance of SVS model for each enclosure

        """

        for enc, enc_dict in self.all_info.enclosures_info.items():

            # create a `soil_types` dictionary from the `soil_types_dict` dictionary
            # where each key is the soil type name and each value is the values
            # in the `parameters` key of the corresponding soil type dictionary
            # `stype_dict["parameters"]` is a dictionary where each key is the
            # name of the parameter and each value is the value of the parameter
            soil_types = {
                stype_name: stype_dict["parameters"]
                for stype_name, stype_dict in self.all_info.soil_types_dict.items()
            }

            # create a `enc_layering` dict from the self.all_info.enclosures_info
            # where each key is the name of the enclosure and each value is the
            # layering of the enclosure
            enc_layering = {
                enc_name: enc_dict["layering"]
                for enc_name, enc_dict in self.all_info.enclosures_info.items()
            }

            # create the dataframe containing the soil parameters for each layer
            dfenc = create_enclosure_dataframe(
                enc_layering, enc, self.all_info.soil_params, soil_types
            )

            # add the site parameters to the dataframe
            for key, value in self.all_info.site_params.items():
                dfenc[key] = value

            # input data for the SVS model
            required_data = ModelInputData(
                work_dir_path=self.all_info.paths_dict["svs_runs"],
                lyinfo=dfenc,
                metfile_path=self.all_info.paths_dict["main_meteo"],
                exec_file_path=self.all_info.paths_dict["svs_executable"],
                host_dir_name=enc,
                meteo_col_names=self.all_info.meteo_cols,
                param_col_names={
                    **{k: k for k in self.all_info.soil_params},
                    **{k: k for k in self.all_info.site_params}
                },
                model_params=self.all_info.svs_params,
                copy_metfile=self.base_dotmetfile_path,
                start_date=self.all_info.start_date,
                spinup_end_date=self.all_info.spinup_end_date,
                time_zone=self.all_info.spinup_end_date_tz
            )

            # assign the `required_data` object to the `enclosures_info` dictionary
            self.all_info.enclosures_info[enc]["input_data"] = required_data

            # create an instance of SVS model for each enclosure
            if create_svs:
                self.svs_instances[enc] = SVSModel(required_data, True)

    def create_ensembles(self):
        """
        Create the ensemble scenarios for each enclosure.

        """

        self.all_info.enclosures_info, self.all_info.soil_types_dict = (
            create_ensemble_scenarios(
                self.all_info.enclosures_info, self.all_info.soil_types_dict,
                self.all_info.n_ens, self.all_info.random_seed,
                self.all_info.cm_fc_wilt_delta,
                class_method=True
            )
        )

    def run_ensembles(self, save_name="ensemble_outputs.csv"):
        """
        Run the ensembles in parallel.

        """

        for enc, enc_dict in self.all_info.enclosures_info.items():

            message = "\n" + "-" * 80 + "\n"
            message += F"Running the ensemble for enclosure {enc}...\n"
            print(message)


            pert_enc = PerturbAndRun(
                enc_dict["input_data"], enc_dict["scenarios"],
                self.all_info.paths_dict["met_ens"],
                njobs=self.all_info.n_cpus
            )

            pert_enc.run_all_parallel()

            enc_dict["ensemble_run"] = pert_enc

        # get the outputs
        self.get_outputs(save_name)

    def get_outputs(self, save_name):
        """
        Get the outputs from the ensemble runs and save them in a csv file.

        """

        for enc, enc_dict in self.all_info.enclosures_info.items():

            # get the ensemble run object
            ens_run = enc_dict["ensemble_run"]

            # get the output dataframe
            dfout = ens_run.dfoutput.copy()

            # add a column with the enclosure name
            dfout["enclosure"] = enc

            # add the output dataframe to the `ens_outputs` dataframe
            self.ens_outputs = pd.concat([self.ens_outputs, dfout], axis=0)

        # save it in `svs_runs` directory
        self.ens_outputs.to_csv(
            self.all_info.paths_dict["svs_runs"] / save_name
        )

        # report
        print("\n" + "-" * 80 + "\n")
        print("Done with the ensemble runs.\n")
        print("The outputs are saved in the `svs_runs` directory.\n")



# ________________________________________________________________ <<< main >>>
