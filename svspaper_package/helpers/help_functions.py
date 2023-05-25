'''
A collection of functions used in the notebooks.
'''

import numpy as np
import re
import pandas as pd
from pandas.core.frame import DataFrame
from pathlib import Path
from pytz import UTC
import datetime
from collections import OrderedDict
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from copy import deepcopy
from scipy.optimize import curve_fit
from scipy.stats import pearsonr


from statsmodels.tsa.arima.model import ARIMA
from smt.sampling_methods import LHS
from metpy.calc import specific_humidity_from_dewpoint
from metpy.units import units
from sklearn.metrics import mean_absolute_error
# from bias_correction import BiasCorrection

from ..svs.prep_svs import check_metvar_units, create_dtime_vars, process_pcpn
from ..svs.prep_svs import order_dfmet, check_meteo_colnames_dict, calc_psnow


def swrc_svs(theta_sat: float, psi_sat: float, estimate_psi: bool = False):
    '''
    This function returns a function that estimates swc using the SWRC model
    implemented in the SVS model. Suitable to use with `scipy.optimize.curve_fit`.

    Parameters
    ----------
    theta_sat : the volumetric soil water content at saturation (%)

    psi_sat : the air-entery suction (kPa)

    estimate_psi : if True, the function returns a function that also requires
                    `psi_sat` as an input. If False, the function returns a
                    function that only requires `bcoef` as an input parameter to
                    estimate the swc.

    Returns
    -------
    out : a function that accepts the independent variable as its first arg
          and the fitting parameter as the second.
    '''

    def suction_function(swc, bcoef):
        '''
        SVS SWRC equation to estimate suction based on SWC.

        Parameters
        ----------
        swc : the volumetric soil water content (%).

        bcoef : fitting parameter

        Returns
        -------
        out : suction (kPa)
        '''

        nonlocal theta_sat, psi_sat

        suction_kpa = psi_sat * pow((swc / theta_sat), -bcoef)

        return suction_kpa

    def suction_function_psi(swc, bcoef, psi_sat):
        '''
        SVS SWRC equation to estimate suction based on SWC.

        Parameters
        ----------
        swc : the volumetric soil water content (%).

        bcoef : fitting parameter

        psi_sat : the air-entery suction (kPa) - fitting parameter

        Returns
        -------
        out : suction (kPa)
        '''

        nonlocal theta_sat

        suction_kpa = psi_sat * pow((swc / theta_sat), -bcoef)

        return suction_kpa

    if estimate_psi:

        return suction_function_psi

    else:

        return suction_function


def swc_function(suction_kpa, bcoef: float, theta_sat: float, psi_sat: float):
    '''
    SVS SWRC equation to estimate swc based on suction.

    Parameters
    ----------
    suction_kpa : the suction (kPa).

    bcoef : the fitting parameter to be estimated.

    theta_sat : the volumetric soil water content at saturation (%)

    psi_sat : the suction at air-entery (kPa)

    Returns
    -------
    out : Vol SWC (%)
    '''

    swc = (suction_kpa / psi_sat)**(-1/bcoef) * theta_sat
    return swc


def fill_hyprop_dict(hyprop_dict: dict, hyprop_data_path: Path):
    """
    This function reads the HYPROP data and fills the `hyprop_dict` with the
    data for each soil type.

    Parameters
    ----------
    hyprop_dict : dict
        A dict that each key is a soil type and the value is a dict that contains
        required info to read the data from the csv file

    hyprop_data_path : Path
        Path to the HYPROP data

    Returns
    -------
    hyprop_dict : dict
    """

    # copy the dict
    hyprop_dict = deepcopy(hyprop_dict)

    # read the data
    hyprop_data = pd.read_csv(hyprop_data_path)

    # loop over the soil types
    for soil_type_dict in hyprop_dict.values():

        # find the samples for this soil type by searching for the column names
        # that contains `vol_swc_col`, replace the `#` with a regex that matches
        # one or two digits
        soil_type_labels = []

        for temp_label in ["vol_swc_col", "suction_col"]:
            soil_type_labels += [
                st1 for st1 in hyprop_data.columns.values
                if re.findall(
                    soil_type_dict[temp_label].replace("#", "[0-9]"),
                    st1
                ) or re.findall(
                    soil_type_dict[temp_label].replace("#", "[0-9][0-9]"),
                    st1
                )
            ]

        soil_type_dict["nsamples"] = len(soil_type_labels) // 2
        soil_type_dict["samples"] = hyprop_data.loc[:, soil_type_labels]

    return hyprop_dict

def fit_swrc(
        hyprop_dict: dict, suction_fc: float = 33,
        treat_aes_as_param: bool = False
    ) -> dict:
    """
    This function fits the soil water retention curve (SWRC) for each sample
    in the HYPROP data.

    Parameters
    ----------
    hyprop_dict : dict
        A dict that each key is a soil type and the value is a dict that contains

    suction_fc : float, optional
        Suction at field capacity (kPa), by default 33

    treat_aes_as_param : bool, optional
        Treat air entry suction as a parameter or not, by default False

    Returns
    -------
    hyprop_dict : dict
    """

    # copy the dict
    hyprop_dict = deepcopy(hyprop_dict)

    # loop over the soil types
    for soil_type_dict in hyprop_dict.values():

        # dataframe to store the fit results
        df_fit = pd.DataFrame(index=range(soil_type_dict["nsamples"]))

        # loop over the samples
        for smpl_i in range(soil_type_dict["nsamples"]):

            # get the vol swc and suction values for this sample
            vol_swc = soil_type_dict["samples"].loc[
                :, soil_type_dict["vol_swc_col"].replace("#", str(smpl_i + 1))
            ].values
            suction = soil_type_dict["samples"].loc[
                :, soil_type_dict["suction_col"].replace("#", str(smpl_i + 1))
            ].values

            # skip this sample if there is no data for it
            if (len(vol_swc) * len(suction)) == 0:
                continue
            else:
                dfsample = pd.DataFrame(
                data=np.vstack((vol_swc, suction)).T,
                columns=["vol_swc", "suction"]
                ).dropna()

            # Wsat is the largest value of the vol swc values
            Wsat = dfsample["vol_swc"].max()

            # Wwilt is the smallest value of the vol swc values
            Wwilt = dfsample["vol_swc"].min()

            # suction at air entry value
            # this is the suction value at which the vol swc is `Wsat`
            psi_ae = dfsample.loc[
                dfsample["vol_swc"] == Wsat, "suction"
            ].values[0]

            # suction at Wwilt
            psi_wilt = dfsample.loc[
                dfsample["vol_swc"] == Wwilt, "suction"
            ].values[0]

            # Wfc is the value of the vol swc values at which the suction is
            # closest to `suction_fc`
            Wfc = dfsample.loc[
                np.abs(dfsample["suction"] - suction_fc).idxmin(), "vol_swc"
            ]

            # suction at Wfc
            psi_fc = dfsample.loc[
                dfsample["vol_swc"] == Wfc, "suction"
            ].values[0]

            # curve fitting - Suction ~ VolSWC
            xtrue = dfsample["vol_swc"].values / 100
            ytrue = dfsample["suction"].values

            curve_fit_results = curve_fit(
                swrc_svs(Wsat/100, psi_ae, treat_aes_as_param),
                method="lm",
                xdata=xtrue,
                ydata=ytrue,
                maxfev = 100000
            )

            if treat_aes_as_param:
                bcoef, psi_ae = curve_fit_results[0]
                ypreds = swrc_svs(Wsat/100, psi_ae, treat_aes_as_param)(xtrue, bcoef, psi_ae)
            else:
                bcoef = curve_fit_results[0][0]
                ypreds = swrc_svs(Wsat/100, psi_ae, treat_aes_as_param)(xtrue, bcoef)



            # calc metrics
            r2 = pearsonr(ytrue, ypreds)[0]**2
            mae = mean_absolute_error(ytrue, ypreds)

            # store the results in `df_fit`
            df_fit.loc[smpl_i, "Wsat_percent"] = Wsat
            df_fit.loc[smpl_i, "psi_ae_kPa"] = psi_ae
            df_fit.loc[smpl_i, "Wwilt_percent"] = Wwilt
            df_fit.loc[smpl_i, "psi_wilt_kPa"] = psi_wilt
            df_fit.loc[smpl_i, "Wfc_percent"] = Wfc
            df_fit.loc[smpl_i, "psi_fc_kPa"] = psi_fc
            df_fit.loc[smpl_i, "bcoef"] = bcoef
            df_fit.loc[smpl_i, "r2"] = r2
            df_fit.loc[smpl_i, "mae_kPa"] = mae

        # store a copy of `df_fit` in the dict
        soil_type_dict["samples_info"] = df_fit.copy()

    return hyprop_dict

# find best sample for each soil type
def find_best_sample(hyprop_dict: dict) -> dict:
    """
    This function finds the best sample for each soil type.

    Parameters
    ----------
    hyprop_dict : dict
        A dict that each key is a soil type and the value is a dict that contains

    Returns
    -------
    hyprop_dict : dict
    """

    # copy the dict
    hyprop_dict = deepcopy(hyprop_dict)

    # loop over the soil types
    for soil_type_dict in hyprop_dict.values():

        # normalize the `mae` values to the range [0, 1]
        soil_type_dict["samples_info"]["mae_norm"] = (
            soil_type_dict["samples_info"]["mae_kPa"] -
            soil_type_dict["samples_info"]["mae_kPa"].min()
        ) / (
            soil_type_dict["samples_info"]["mae_kPa"].max() -
            soil_type_dict["samples_info"]["mae_kPa"].min()
        )

        # `sample_score`: `r2` - `mae_norm`
        soil_type_dict["samples_info"]["sample_score"] = (
            soil_type_dict["samples_info"]["r2"] -
            soil_type_dict["samples_info"]["mae_norm"]
        )

        # sort
        soil_type_dict["samples_info"] = soil_type_dict["samples_info"].sort_values(
            by=["sample_score"], ascending=False
        ).reset_index(drop=True)

    return hyprop_dict

def read_field_meteo(data_path: Path, dt_label: str, tz="UTC-4") -> pd.DataFrame:
    """
    Read the field station meteo data and return a dataframe with the
    datetime index.

    Parameters
    ----------
    data_path : Path
        Path to the field station meteo data.

    dt_label : str
        Label of the datetime column.

    tz : str, optional
        Timezone of the data, by default "UTC+4"

    Returns
    -------
    df_field : pd.DataFrame
        Dataframe with the datetime index.
    """

    # read field station data
    df_field = pd.read_csv(data_path)

    # remove duplicate rows
    df_field = df_field[~df_field.duplicated()]

    # convert the index to datetime
    try:
        df_field.index = pd.to_datetime(
                df_field[dt_label], format='%Y-%m-%d %H:%M:%S'
        )
    except ValueError:
        df_field.index = pd.to_datetime(
                df_field[dt_label], format='%Y-%m-%d %H:%M'
        )


    # if tz is given as a string `UTC+#`
    if tz.startswith("UTC"):

        # remove duplicate indices
        df_field = df_field[~df_field.index.duplicated(keep='first')]

        # add/subtract x hours to the time index
        # (In my case, This is done because the datetime is stored in UTC-4,
        #  which is local time for the field station without daylight saving time)
        if tz[3] == "-":
            df_field.index = df_field.index + (pd.Timedelta(hours=int(tz[4:])))
        elif tz[3] == "+":
            df_field.index = df_field.index - (pd.Timedelta(hours=int(tz[4:])))

        # localize the index to UTC
        df_field.index = df_field.index.tz_localize('UTC')

    else:

        # localize the index to the given timezone
        df_field.index = df_field.index.tz_localize(tz)

        # convert the index to UTC
        df_field.index = df_field.index.tz_convert('UTC')

    # convert wind speed from km/h to m/s
    if (
        "Wind speed (m.s-1)" not in df_field.columns and
        "Wind speed (km.h-1)" in df_field.columns
    ):
        df_field["Wind speed (m.s-1)"] = df_field["Wind speed (km.h-1)"] / 3.6

        # replace zeros with NaNs
        df_field["Wind speed (m.s-1)"] = df_field["Wind speed (m.s-1)"].replace(
            0, np.nan
        )

    return df_field

def read_main_met_source(
        data_path: Path, dt_label: str, tz="UTC"
) -> pd.DataFrame:
    """
    Read the main meteo data source and return a dataframe with the
    datetime index.

    Parameters
    ----------
    data_path : Path
        Path to the main meteo data source.

    dt_label : str
        Label of the datetime column.

    tz : str, optional
        Timezone of the data, by default "UTC"

    Returns
    -------
    df_met : pd.DataFrame
        Dataframe with the datetime index.
    """

    # read main meteo data source
    df_met = pd.read_csv(data_path)

    # remove duplicate rows
    df_met = df_met[~df_met.duplicated()]

    # convert the index to datetime
    try:
        df_met.index = pd.to_datetime(
                df_met[dt_label], format='%Y-%m-%d %H:%M:%S'
        )
    except ValueError:
        df_met.index = pd.to_datetime(
                df_met[dt_label], format='%Y-%m-%d %H:%M'
        )

    # if tz is given as a string `UTC+#`
    if tz.startswith("UTC"):

        # remove duplicate indices
        df_met = df_met[~df_met.index.duplicated(keep='first')]

        # add/subtract x hours to the time index
        # (In my case, This is done because the datetime is stored in UTC-4,
        #  which is local time for the field station without daylight saving time)
        if tz[3] == "-":
            df_met.index = df_met.index + (pd.Timedelta(hours=int(tz[4:])))
        elif tz[3] == "+":
            df_met.index = df_met.index - (pd.Timedelta(hours=int(tz[4:])))

        # localize the index to UTC
        df_met.index = df_met.index.tz_localize('UTC')

    else:

        # localize the index to the given timezone
        df_met.index = df_met.index.tz_localize(tz)

        # convert the index to UTC
        df_met.index = df_met.index.tz_convert('UTC')



    return df_met

def prec_cr(
    data, u_lab='Wind speed (m.s-1)', t_lab='Air temperature (°C)',
    coeff_a=0.021, coeff_b=0.74, coeff_c=0.66, u_max=8, t_max=2
):
    """
    Precipitation catch ratio correction based on wind speed and air temperature.

    Parameters
    ----------
    data : pandas.DataFrame
        Dataframe containing the precipitation and meteorological variables.

    u_lab : str, optional
        Label for wind speed, by default 'Wind speed (m.s-1)'

    t_lab : str, optional
        Label for air temperature, by default 'Air temperature (°C)'

    coeff_a : float, optional
        Parameter a, by default 0.021 for a double alter rain gauge

    coeff_b : float, optional
        Parameter b, by default 0.74 for a double alter rain gauge

    coeff_c : float, optional
        Parameter c, by default 0.66 for a double alter rain gauge

    u_max : float, optional
        Maximum wind speed, by default 6 m.s-1

    t_max : float, optional
        Maximum air temperature, by default 2 °C

    Returns
    -------
    pandas.DataFrame
        Dataframe with a new column for corrected precipitation.

    Reference
    ----------
    `The quantification and correction of wind-induced precipitation measurement
      errors`, Kochendorfer et al. (2017): check the corrigendum
    """

    data = data.copy()

    # clip wind speed to 6 m.s-1
    data[u_lab] = data[u_lab].clip(upper=u_max)

    # calculate correction factor
    data['prec_cr'] = np.exp(
        -coeff_a * data[u_lab] *
        (1 - np.arctan(coeff_b * data[t_lab]) + coeff_c)
    )

    # set CR to 1 if temperature is above 2 °C
    data.loc[data[t_lab] > t_max, 'prec_cr'] = 1

    # a new column for corrected precipitation
    data['Precipitation (mm) corrected'] = data['Precipitation (mm)'] / \
        data['prec_cr']

    return data

def qm_bias_correction(
        variable_label: str, df_source_1: pd.DataFrame,
        df_source_2: pd.DataFrame
) -> pd.DataFrame:
    """
    Quantile mapping bias correction.

    Parameters
    ----------
    variable_label : str
        Label of the variable to correct.

    df_source_1 : pd.DataFrame
        Dataframe containing the variable to correct.

    df_source_2 : pd.DataFrame
        Dataframe containing the variable to use for the bias correction.

    Returns
    -------


    Conditions
    ----------
    The two dataframes must have a datetime index with a UTC timezone and
    the same frequency.

    """

    assert df_source_1.index.tz in ('UTC', datetime.timezone.utc), 'df_source_1 must have a UTC timezone'
    assert df_source_2.index.tz in ('UTC', datetime.timezone.utc), 'df_source_2 must have a UTC timezone'
    assert (df_source_1.index.freq == df_source_2.index.freq), (
        'df_source_1 and df_source_2 must have the same frequency'
    )
    assert (variable_label in df_source_1.columns), (
        f'{variable_label} is not a column of df_source_1'
    )
    assert (variable_label in df_source_2.columns), (
        f'{variable_label} is not a column of df_source_2'
    )

    # initiate the bias correction class

    # not developed yet
    raise NotImplementedError

def calc_bias_std(
    compare_variables: list, df_source_1: DataFrame, df_source_2: DataFrame,
) -> dict:
    '''
    Calculate the bias and standard deviation of the difference between two
    dataframes.

    bias = (df_source_1 - df_source_2).mean() for the variables
    std = (df_source_1 - df_source_2).std() for the variables

    Parameters
    ----------
    compare_variables : list
        List of variables to compare.

    df_source_1 : pandas.DataFrame
        Dataframe containing the variables to compare.

    df_source_2 : pandas.DataFrame
        Dataframe containing the variables to compare.

    Returns
    -------
    dict
        Dictionary containing the bias and standard deviation of the difference
        between the two dataframes for each variable.

    Conditions
    ----------
    The two dataframes must have a datetime index with a UTC timezone and
    the same frequency.
    '''

    # assertions --------------------------------------------------------------
    assert df_source_1.index.tz in ('UTC', datetime.timezone.utc), 'df_source_1 must have a UTC timezone'
    assert df_source_2.index.tz in ('UTC', datetime.timezone.utc), 'df_source_2 must have a UTC timezone'
    assert (df_source_1.index.freq == df_source_2.index.freq), (
        'df_source_1 and df_source_2 must have the same frequency'
    )
    assert all([var in df_source_1.columns for var in compare_variables]), (
        'Some variables to compare are not in df_source_1'
    )
    assert all([var in df_source_2.columns for var in compare_variables]), (
        'Some variables to compare are not in df_source_2'
    )
    # -------------------------------------------------------------------------

    # collect bias and std in a dictionary for each variable
    bias_std = {key_var: {'bias': None, 'std': None}
                for key_var in compare_variables}

    # collect the variable values and the associated bias
    var_info = {}

    for kvar in compare_variables:

        # extract the variables from the respective dataframes
        var_1 = df_source_1[kvar]
        var_2 = df_source_2[kvar]

        # merge the two variables into a dataframe - mutual index
        merged_df = pd.merge(var_1, var_2, how="inner",
                             left_index=True, right_index=True)
        # remove NaNs
        merged_df = merged_df.dropna()

        # calculate the bias and std
        bias_std[kvar]['bias'] = (
            merged_df[kvar + '_x'] - merged_df[kvar + '_y']).mean()
        bias_std[kvar]['std'] = (
            merged_df[kvar + '_x'] - merged_df[kvar + '_y']).std()

        # collect the variable values and the associated bias
        var_info[kvar] = {
            'var_1': merged_df[kvar + '_x'],
            'var_2': merged_df[kvar + '_y'],
            'bias': merged_df[kvar + '_x'] - merged_df[kvar + '_y'],
        }

    return bias_std, var_info


def plot_var_bias(var_info, compare_vars):
    '''
    Plot the bias vs the variable value for each variable in `compare_vars`.
    Print the correlation between the variable and the bias.

    Parameters
    ----------
    var_info : dict
        A dictionary containing the values of the variables for each source and
        each variable and the associated biases.

    compare_vars : list
        A list containing the names of the variables to compare between the two
    '''

    for varname in compare_vars:

        # create a dataframe with the variable value and the bias
        df1 = pd.DataFrame()
        df1['var_1'] = var_info[varname]['var_1']
        df1['bias'] = var_info[varname]['bias']

        # print the correlation between the variable and the bias
        print(
            f"Correlation between {varname} and bias: {df1.corr().iloc[0,1]}")

        # plot var_1 vs bias and the regression line (linear fit)
        sns.regplot(
            x='var_1', y='bias', data=df1, scatter_kws={'alpha': 0.2},
              line_kws={'color': 'red'}
        )
        plt.title(f"{varname} vs bias")
        plt.xlabel(f"{varname}")
        plt.ylabel("Bias")
        plt.show()


def get_ar1_param(bias):
    """
    Fit an AR(1) model to the bias and return the value of the AR(1) parameter.

    Parameters
    ----------
    bias : pandas.Series
        A pandas Series containing the bias.

    Returns
    -------
    float
        The value of the AR(1) parameter.

    Assumptions
    -----------
    The bias is a pandas Series with a datetime index and a UTC timezone and
    a frequency of 1 hour.

    If the time period is not fully covered by the
    bias, the bias is reindexed with a full range of time.

    Constraints
    -----------
    The intercept of the AR(1) model is fixed to 0.
    """

    # assertions --------------------------------------------------------------
    assert isinstance(bias, pd.Series), 'bias must be a pandas Series'
    assert bias.index.tz in ('UTC', datetime.timezone.utc), 'bias must have a UTC timezone'
    # -------------------------------------------------------------------------

    # make a copy of the bias
    bias = bias.copy()

    # get the start and end time of the bias
    t1 = bias.index[0]
    t2 = bias.index[-1]

    # create a full range of time
    full_range = pd.date_range(t1, t2, freq='H')

    # reindex the bias with the full range of time
    bias = bias.reindex(full_range, fill_value=np.nan)

    # fit AR(1) model- fix the intercept to 0
    autor1 = ARIMA(bias, order=(1, 0, 0), freq='H')

    with autor1.fix_params({'const': 0}):
        autor1_fit = autor1.fit()

    return autor1_fit.params['ar.L1']


def perturb_var(
    data, var_label, phi, var_std, additive,
    plimit_u=0.2, plimit_l=-0.2,
    seed=None, var_max=None, var_min=None,
    use_bias=False, var_bias=None
):
    """
    Create a perturbed variable from the original variable using the AR1 model


    Parameters
    ----------
    data : pandas.DataFrame
        The dataframe containing the original variable.

    var_label : str
        The column label of the variable.

    phi : float
        The auto-regressive (AR1) model parameter.

    var_std : float
        The standard deviation of the differences between the original variable
        and the variable from the other source.

    additive : bool
        True if the perturbation is additive, False if it is multiplicative.

    plimit_u : float, optional
        The upper limit for the multiplicative perturbation. The default is 0.25.

    plimit_l : float, optional
        The lower limit for the multiplicative perturbation. The default is -0.25.

    seed : int, optional
        The seed for the random number generator. The default is None.

    var_max : float, optional
        The maximum value of the variable. The default is None.

    var_min : float, optional
        The minimum value of the variable. The default is None.

    use_bias : bool, optional
        True if the bias should be used in the perturbation. The default is False.


    var_bias : float, optional
        The bias of the variable. The default is 0.

    Returns
    -------
    addvar_pert : pandas.Series
        The perturbed variable.

    """

    # calculate sigma-2
    sigma2 = (1 - phi**2) * var_std

    # create a random number generator
    drng = np.random.default_rng(seed)

    # generate the white noise
    # if not use_bias:
    #     white_noise = drng.normal(0, sigma2**0.5, len(data))
    # else:
    #     white_noise = drng.normal(var_bias, sigma2**0.5, len(data))

    white_noise = drng.normal(0, sigma2**0.5, len(data))


    # create the perturbation values for the variable
    perturb = np.zeros(len(data))

    for i in range(1, len(data)):
        perturb[i] = (phi * perturb[i-1]) + white_noise[i]

    # create the perturbed variable
    if additive:
        addvar_pert = data[var_label] + perturb
    else:
        # limit the perturbation
        perturb = np.clip(perturb, plimit_l, plimit_u)
        addvar_pert = data[var_label] * (1 + perturb)

    if var_max is not None:
        addvar_pert = np.clip(addvar_pert, None, var_max)

    if var_min is not None:
        addvar_pert = np.clip(addvar_pert, var_min, None)

    return addvar_pert


def calc_specific_humidity(pressure, dewpoint):
    """
    Calculate the specific humidity from the pressure and dewpoint temperature.

    Parameters
    ----------
    pressure : pandas.Series
        The pressure.

    dewpoint : pandas.Series
        The dewpoint temperature.

    Returns
    -------
    pandas.Series
        The specific humidity.

    Assumptions
    -----------
    The pressure and dewpoint temperature are pandas Series in units of kPa and
    degC, respectively.

    """

    # make them unit aware
    pressure = pressure.to_numpy() * units.kPa
    dewpoint = dewpoint.to_numpy() * units.degC

    # calculate the specific humidity
    spech = specific_humidity_from_dewpoint(pressure, dewpoint)

    return spech.magnitude


def create_perturbed_data(
        original_data, compare_vars, not_perturbed, PREC_FACT, bias_std,
        seed=None,
        precipitation_label="Precipitation (mm) corrected",
        pressure_label="Atmospheric pressure (kPa)",
        dew_point_label="Dew point temperature (°C)"
):
    """
    Creates a perturbed dataframe based on the original data.

    Parameters
    ----------
    original_data : pandas.DataFrame
        The original dataframe.

    compare_vars : dict
        A dictionary of the variables to compare between the original dataframe
        and the perturbed dataframe. The keys are the variable names and the
        values are dictionaries containing the perturbation parameters.

    not_perturbed : list
        A list of the variables that are not perturbed.

    PREC_FACT : float
        The perturbation factor for precipitation. The perturbation is
        calculated as a uniform distribution between 1-PREC_FACT and
        1+PREC_FACT.

    bias_std : dict
        A dictionary containing the bias and standard deviation of the
        difference between the original dataframe and the perturbed dataframe
        for each variable.
        It also contains the AR1 parameter for each variable.

    seed : int, optional
        The seed to use for the random number generator. If None, the seed is
        not set.

    precipitation_label : str, optional
        The label of the precipitation variable in the original dataframe.

    pressure_label : str, optional
        The label of the atmospheric pressure variable in the original
        dataframe.

    dew_point_label : str, optional
        The label of the dew point temperature variable in the original
        dataframe.

    Returns
    -------
    df_pert : pandas.DataFrame
        The perturbed dataframe.
    """

    # make a copy of the original dataframe with the `not_perturbed` variables
    df_pert = original_data.loc[:, not_perturbed].copy()

    # perturb the precipitation data with a uniform distribution
    dfrng = np.random.default_rng(seed=seed)
    pert_factors = dfrng.uniform(1-PREC_FACT, 1+PREC_FACT, size=len(df_pert))
    df_pert[precipitation_label] = (
        original_data[precipitation_label] * pert_factors
    )

    # perturb the other variables
    for var_name, var_dict in compare_vars.items():
        df_pert[var_name] = perturb_var(
            data=original_data, var_label=var_name, phi=bias_std[var_name]['phi'],
            var_std=bias_std[var_name]['std'], seed=seed,
            use_bias=True, var_bias=bias_std[var_name]["bias"],
            **var_dict
        )

    # calculate specific humidity
    df_pert['Specific humidity (kg.kg-1)'] = calc_specific_humidity(
        df_pert.loc[:, pressure_label],
        df_pert.loc[:, dew_point_label]
    )

    return df_pert


def decode_enclosure_layering(enclosure_layering: dict) -> dict:
    """
    Decode the `enclosure_layering` dictionary.

    Parameters
    ----------
    enclosure_layering : dict
        The dictionary with the layering of the enclosures.

    Returns
    -------
    dict
        The decoded dictionary with the layering of the enclosures.
    """

    # create the decoded dictionary
    enclosure_layering_decoded = {}

    # loop over the enclosures
    for enclosure, layering in enclosure_layering.items():

        # split the string
        layering_split = layering.split("+")

        # create a list of dictionaries for each layer
        layering_decoded = []

        # loop over the layers
        for layer in layering_split:

            # split the layer string
            layer_split = layer.split(":")

            # strip the strings
            layer_split = [x.strip() for x in layer_split]

            # create a dictionary for the layer
            try:
                layer_decoded = {
                    "total_thickness": float(layer_split[0]),
                    "per_layer_thickness": float(layer_split[1]),
                    "soil_type": layer_split[2]
                }
            except ValueError:
                print("Make sure the layering string is correct.")
                print("The layering string should be in the format:")
                print("thickness of the stratum:layer_thickness:soil_type")
                print("e.g. `30:5:CM`")
                raise

            # append the layer dictionary to the list
            layering_decoded.append(layer_decoded)

        # add the list of dictionaries to the enclosure dictionary
        enclosure_layering_decoded[enclosure] = layering_decoded

    return enclosure_layering_decoded


def create_enclosure_dataframe(
    enclosure_layering: dict, enclosure_name: str,
    soil_params: list, soil_types: dict
):
    """
    Creates a dataframe with the soil profile for each enclosure.

    Parameters
    ----------
    enclosure_layering : dict
        A dictionary with the soil profile for each enclosure.

    enclosure_name : str
        The name of the enclosure; a key in the `enclosure_layering_decoded` dictionary.

    soil_params : list
        A list of the soil parameters.

    soil_types : dict
        A dictionary with the soil parameters for each soil type.

    Returns
    -------
    en_df : pd.DataFrame
        A dataframe with the soil profile for the given enclosure.
    """

    # decode the `enclosure_layering` dictionary
    enclosure_layering_decoded = decode_enclosure_layering(enclosure_layering)

    en_df = pd.DataFrame(
        columns=["depth", "thickness", "soil_type"]
    )

    start_depth = 0
    for layer in enclosure_layering_decoded[enclosure_name]:

        column_1 = np.arange(
            start_depth + layer["per_layer_thickness"],
            start_depth + layer["total_thickness"] +
            layer["per_layer_thickness"],
            layer["per_layer_thickness"]
        )
        column_2 = np.repeat(layer["per_layer_thickness"], len(column_1))
        column_3 = np.repeat(layer["soil_type"], len(column_1))

        en_df = pd.concat(
            [
                en_df,
                pd.DataFrame(
                    data=np.array([column_1, column_2, column_3]).T,
                    columns=["depth", "thickness", "soil_type"]
                )
            ],
            ignore_index=True
        )

        start_depth = column_1[-1]

    # make `depth` and `thickness` columns numerical
    en_df["depth"] = pd.to_numeric(en_df["depth"])
    en_df["thickness"] = pd.to_numeric(en_df["thickness"])

    # populate the dataframe with the soil properties
    for prop in soil_params:
        en_df[prop] = en_df.loc[:, "soil_type"].map(
            soil_types.get).map(lambda x: x.get(prop))

    return en_df


# functions used in `create_ensemble_members.ipynb` ___________________________
def round_params(p_dict):
    '''
    Round the values of the parameters in the `p_dict`.

    Parameters
    ----------
    p_dict : dict
        A dictionary with the parameters.

    Returns
    -------
    p_dict : dict
        A dictionary with the rounded parameters.
    '''

    for key, val in p_dict.items():
        match key:
            case "wsat" | "wfc" | "wwilt" | "wunfrz" | "user_wfcdp" | "psisat":
                p_dict[key] = np.round(val, 4)

            case "tperm" | "z0v" | "d95" | "d50" | "conddry" | "condsld" | "bcoef":
                p_dict[key] = np.round(val, 2)

            case "sand" | "clay":
                p_dict[key] = np.round(val, 1)

            case "ksat":
                p_dict[key] = [float(F"{kval:.2e}") for kval in val]

    return p_dict


def load_input_data(enclosures: dict) -> None:
    """
    Load the input data pickle files for each enclosure and store them in the
    `enclosures` dictionary.

    Parameters
    ----------
    enclosures : dict
        A dictionary containing the information for each enclosure.

    """

    for enc_dict in enclosures.values():
        with open(enc_dict["path"], "rb") as f:
            enc_dict["input_data"] = pickle.load(f)


def pop_all_limits(soil_types: dict) -> OrderedDict:
    """
    Populate the all_limits dictionary; prefix the soil type to the parameter
    name.

    Used inside the `create_ensemble_members` notebook.

    Parameters
    ----------
    soil_types : dict
        A dictionary containing the limits for each soil type

    """

    all_limits = OrderedDict()

    for soil_type, soil_type_dict in soil_types.items():
        for param, limits in soil_type_dict["limits"].items():
            all_limits[f"{soil_type}_{param}"] = limits

    return all_limits


def pop_all_samples(
    all_limits: OrderedDict, n_s: int, r_seed: int
) -> OrderedDict:
    """
    Populate the all_samples dictionary.
    (used in `create_ensemble_members.ipynb`)

    Parameters
    ----------
    all_limits : OrderedDict
        A dictionary containing the limits for each soil type.

    n_s : int
        Number of samples.

    r_seed : int
        Random seed.


    Returns
    -------
    all_samples : OrderedDict
        A dictionary containing the samples for the parameters to be perturbed
        for each soil type.
    """

    all_samples = OrderedDict()

    # create the sampling object
    sampling = LHS(
        xlimits=np.array(list(all_limits.values())),
        random_state=r_seed
    )

    # create the samples
    samples = sampling(n_s)

    # populate the all_samples dictionary
    for i, param in enumerate(all_limits.keys()):
        all_samples[param] = samples[:, i]

    return all_samples


def store_samples(soil_types, all_samples):
    """
    Store the samples in the `soil_types` dictionary. Check if the samples are
    within the limits.

    (used in `create_ensemble_members.ipynb`)

    Parameters
    ----------
    soil_types : dict
        A dictionary containing the information for each soil type.

    all_samples : OrderedDict
        A dictionary containing the samples for the parameters to be perturbed
        for each soil type.
    """

    # store the samples in the `soil_types` dictionary
    for soil_type, soil_type_dict in soil_types.items():
        soil_type_dict["samples"] = OrderedDict()
        for param in soil_type_dict["limits"].keys():
            soil_type_dict["samples"][param] = all_samples[f"{soil_type}_{param}"]

    # for each soil type, check if the samples are within the limits
    for soil_type, soil_type_dict in soil_types.items():

        for param, param_samples in soil_type_dict["samples"].items():

            min_val, max_val = soil_type_dict["limits"][param]

            if not np.all((param_samples >= min_val) & (param_samples <= max_val)):
                raise ValueError(
                    f"Samples for {param} are not within the limits: "
                    f"{min_val} - {max_val}"
                )


def create_scnearios(
    enclosures: dict, soil_types: dict, n_ens: int, fc_delta_wilt_cm: float
) -> None:
    """
    Create scenarios for each enclosure and store them in the `enclosures` dict.
    (used in `create_ensemble_members.ipynb)

    Parameters
    ----------
    enclosures : dict
        A dictionary containing the information for each enclosures.

    soil_types : dict
        A dictionary containing the limits for each soil type.

    n_ens : int
        Number of ensemble members to create.

    fc_d_wilt : float
        The minimum difference between the field capacity and wilting point
        for the CM soil type.

    """

    # now we ought to make scenarios for each enclosure
    for enc_dict in enclosures.values():

        # create a key for the scenarios in the enclosure dictionary associated with
        # the current enclosure
        enc_dict["scenarios"] = OrderedDict()

        # create `N_ENS` scenarios for the current enclosure
        for scn_i in range(n_ens):

            # copy the base input data for the current enclosure
            df_base = enc_dict["input_data"].lyinfo.copy()

            # create a key for the current scenario
            enc_dict["scenarios"][f"scn_{scn_i}"] = dict()

            # change the soil parameters for each soil type according to the
            # samples
            for stype in df_base.soil_type.unique():

                for param in soil_types[stype]["samples"]:

                    if param != "user_wfcdp":
                        df_base.loc[
                            df_base.soil_type == stype, param
                        ] = soil_types[stype]["samples"][param][scn_i]

                        enc_dict["scenarios"][f"scn_{scn_i}"][param] = df_base.loc[
                            :, param
                        ].values

                    else:
                        # for the user_wfcdp parameter, we need to multiply the
                        # sample value by the Wsat value of the last layer in the
                        # soil profile
                        user_wfcdp = (
                            soil_types[stype]["samples"][param][scn_i] *
                            df_base.loc[:, "wsat"].values[-1]
                        )
                        enc_dict["scenarios"][f"scn_{scn_i}"][param] = user_wfcdp

                    # for the CM soil type, we need to ensure that Wfc is always
                    # greater than Wwilt by a certain amount
                    if stype == "CM":
                        if df_base.loc[
                            df_base.soil_type == stype, "wfc"
                        ].values[0] < (
                            df_base.loc[
                                df_base.soil_type == stype, "wwilt"
                            ].values[0] + fc_delta_wilt_cm
                        ):
                            df_base.loc[
                                df_base.soil_type == stype, "wfc"
                            ] = (
                                df_base.loc[
                                    df_base.soil_type == stype, "wwilt"
                                ].values[0] + fc_delta_wilt_cm
                            )

                            enc_dict["scenarios"][f"scn_{scn_i}"]["wfc"] = df_base.loc[
                                :, "wfc"
                            ].values

            # round the scenario
            enc_dict["scenarios"][f"scn_{scn_i}"] = round_params(
                enc_dict["scenarios"][f"scn_{scn_i}"]
            )


def create_ensemble_scenarios(
        enclosures: dict, soil_types: dict, n_ens: int, r_seed: int,
        fc_delta_wilt_cm: float, class_method: bool = False
) -> dict:
    """
    Create the ensemble scenarios.

    Parameters
    ----------
    enclosures : dict
        A dictionary containing the information for each enclosures.

    soil_types : dict
        A dictionary containing the limits for each soil type.

    n_ens : int
        Number of ensemble members to create.

    r_seed : int
        Random seed.

    fc_delta_wilt_cm : float
        The minimum difference between the field capacity and wilting point
        for the CM soil type.

    Returns
    -------
    enclosures : dict
        A dictionary containing the information for each enclosures.

    soil_types : dict
        A dictionary containing the limits for each soil type.
    """

    # create a dictionary to store the limits for all soil types:
    # this is to ensure that, for instance, the CM-Sand for enclosure 1 has the
    # same value as the CM-Sand for enclosure 2 for scenario #1
    all_limits = OrderedDict()

    # create a dictionary to store the samples for all soil types
    all_samples = OrderedDict()

    if not class_method:
        # make a copy of the enclosures/soil_types dictionaries
        enclosures = deepcopy(enclosures)
        soil_types = deepcopy(soil_types)

        # load the input data for each enclosure
        load_input_data(enclosures)


    # populate the all_limits dictionary; prefix the soil type to the parameter name
    all_limits = pop_all_limits(soil_types)

    # populate the all_samples dictionary;
    all_samples = pop_all_samples(all_limits, n_ens, r_seed)

    # store the samples in the `soil_types` dictionary
    store_samples(soil_types, all_samples)

    # now we ought to make scenarios for each enclosure
    create_scnearios(enclosures, soil_types, n_ens, fc_delta_wilt_cm)

    return enclosures, soil_types


# a function to save the CSV hourly meteo data as .met file for SVS
def save_csv_as_met(save_path: str, col_names, dfmet):
    '''
    Create `basin_forcing.met` file for SVS by having the hourly meteo data
    in a CSV file.

    Parameters
    ----------
    save_path : str
        The full path to save the `basin_forcing.met` file, including the
        file name.

    col_names : dict
        A dictionary containing the names of the columns in the CSV file.

    dfmet : pandas.DataFrame
        A dataframe containing the hourly meteo data.

    Returns
    -------
    None

    '''

    # copy the dataframe
    dfmet = dfmet.copy()

    # check whether the required columns are present in the dictionary
    check_meteo_colnames_dict(col_names)

    # check units of the variables
    check_metvar_units(dfmet, col_names)

    # add the required time columns in the dataframes
    dfmet, not_used, not_used_2 = create_dtime_vars(
        dfmet, col_names["utc_dtime"]
    )

    # process the precipitation for the dataframe
    dfmet = process_pcpn(dfmet, col_names["precipitation"], calc_psnow,
                            **dict(air_temp=dfmet.loc[:, col_names["air_temperature"]],
                                rel_humidity=dfmet.loc[:,
                                                        col_names["relative_humidity"]]
                                )
                            )

    # order the dataframe columns
    dfmet = order_dfmet(dfmet, col_names) 

    # let the instance own the text file
    met_file = dfmet.to_string(
        header=False, col_space=4, index=False)

    # save the file
    with open(save_path, "w", encoding="utf-8") as mfile:
        mfile.write(met_file)
