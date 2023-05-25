# <<< doc >>> -----------------------------------------------------------------
'''
This is the module in which I keep the definition of the helper functions that
are used in other modules.

'''
# _________________________________________________________________ <<< doc >>>

# <<< imports >>> -------------------------------------------------------------
import os
from pathlib import Path
import re
from typing import Union, Callable

import numpy as np
from numpy.typing import ArrayLike
import pandas as pd
from pandas.core.frame import DataFrame

from sklearn.metrics import make_scorer, mean_absolute_error
import cdsapi
# _____________________________________________________________ <<< imports >>>

# <<< main >>> ----------------------------------------------------------------

def send_request(
    name: str, originating_center: str, system: str, variable: list,
    year: list, month: list, day: list, leadtime_hour: list, area: list,
    file_format: str, save_path: Union[Path, str]
):
    '''
    Sends the data download request to CDS API.

    Parameters
    ----------
    name : str
        The name of the dataset to download.

    originating_center : str
        The name of the originating center.

    system : str
        The name of the system.

    variable : list
        A list of the variables to download.

    year : list
        A list of the years to download.

    month : list
        A list of the months to download.

    day : list
        A list of the days to download.

    leadtime_hour : list
        A list of the leadtimes to download.

    area : list
        A list of the bounding box coordinates to download.

    file_format : str, default='netcdf'
        The format of the data to download.

    save_path : str
        The path to save the downloaded data. The path must include the name
        of the file to save the data to.

    '''

    # the request
    request = {
        'originating_center': originating_center,
        'system': system,
        'variable': variable,
        'year': year,
        'month': month,
        'day': day,
        'leadtime_hour': leadtime_hour,
        'area': area,
        'format': file_format
    }

    # send the request
    cdsapi.Client().retrieve(name, request, save_path)






def calc_nmbe(ytrue: ArrayLike, ypreds: ArrayLike) -> float:
    '''
    Calculate the normalized mean-absolute-error (NMAE).

    Parameters
    ----------
    ytrue : ArrayLike
        The vector of true or observation values.

    ypdres : ArrayLike
        The vector of estimated values.

    Returns
    -------
    nmbe : float
        A scalar value representing the NMBE between the observations and estimates.
    '''

    # calc the normalized mean-absolute-error
    nmbe = 100 * np.mean(ypreds - ytrue) / np.mean(ytrue)

    return nmbe


def calc_nmae(ytrue: ArrayLike, ypreds: ArrayLike) -> float:
    '''
    Calculate the normalized mean-absolute-error (NMAE).

    Parameters
    ----------
    ytrue : ArrayLike
        The vector of true or observation values.

    ypdres : ArrayLike
        The vector of estimated values.

    Returns
    -------
    nmae : float
        A scalar value representing the NMAE between the observations and estimates.
    '''

    # calc the normalized mean-absolute-error
    nmae_ = 100 * mean_absolute_error(ytrue, ypreds) / np.mean(ytrue)

    return nmae_


def nmae(calc_nmae_func: Callable = calc_nmae) -> Callable:
    '''
    Turn the function that calculates the normalized mean-absolute-error (NMAE)
    intoa sklearn scorer.

    Parameter
    --------
    calc_nmae : Callable
        The function that accepts 'ytrue' and 'ypreds' to calculate the NMAE.

    Returns
    -------
    nmae : Callable
    '''

    nmae_ = make_scorer(calc_nmae_func, greater_is_better=False)
    return nmae_


def create_spatial_cv_folds(
    file_path: Path, nresamples: int, loc_column: str, n_analysis: int,
    rseed: int
) -> list:
    '''
    This function creates indices for resamples that are going to be created
    by spatially separating analysis and assessment data.

    Parameters:
    ----------
    file_path : Path
        The path to the training data used for the tuning/selection scheme.

    nresamples : int
        The number of resample to create for the cross-validation.

    loc_column : str
        The column's label that should be used as the location identifier for
        the rows. In my case, its ID of the FLUXNET site. It could be lat/lon
        or etc in other cases.

    n_analysis : int
        The number of distinct location that the analysis data will come from.
        In my case, its the number of FLXUNET sites whose data will make up the
        anayslis data.

    rseed : int
        The random number seed to be used in the random number generator.

    Returns:
    --------
    list_resample_indices : list
        A list of size-2 tuples, where first member of each tuple is indices needed
        to create the analysis data, and the second member contains the indices
        corresponding to the assessment data for the split.
    '''

    # a list that contains the indices that will be used to select the data
    # from `data` for each resample. Every member of this list is a tuple.
    # The first element of this tuple is a list of indices for the analysis data.
    # the second element of the tuple is a list of indices for the assessment data.
    list_resample_indices = []

    # hold the training data
    training_data = None

    # the random number generator
    drng = np.random.default_rng(rseed)

    # counter to iterate as many as `nresamples`
    res_cntr = 0

    # read the data
    assert_file_exist(file_path)
    training_data = pd.read_csv(file_path)

    assert (loc_column in list(training_data.columns)), (
        F"the column: {loc_column} does not exist in the `data` dataframe."
    )

    # determine the total number of locations
    ntot = training_data[loc_column].unique().size

    # create the resamples one at a time
    while res_cntr < nresamples:

        # choose the analysis location from all of the location randomly
        # mind that these are just labels for now.
        an_locs = training_data[loc_column].unique()[
            drng.choice(range(ntot), n_analysis, False)
        ]

        # now we can have the assessment locations too
        as_locs = set(training_data[loc_column].unique()) - set(an_locs)

        # the indices (inside `training_data`) for the analyses and assessment sites
        an_inds = training_data.loc[training_data[loc_column].isin(
            an_locs), ].index.values
        as_inds = training_data.loc[training_data[loc_column].isin(
            as_locs), ].index.values

        # cache the indices in the list
        list_resample_indices.append((an_inds, as_inds))

        # increment the resample counter
        res_cntr += 1

    print(
        F"Done! Created {nresamples} resamples via spatial separation. "
        F"The column: `{loc_column}` was used as the location identifier.\n"
        F"For each resample, the analysis data come from {n_analysis} locations, and the "
        F"assessment data come from {ntot - n_analysis} different locations.\n"
    )

    # return the list
    return list_resample_indices


def cvresults_to_df(
    dict_cv_results: dict, metric_name: str, metric_pos: bool = True
) -> DataFrame:
    '''
    Converts the `cv_results_` attributes of the `GridSearchCV` object to
     a pandas dataframe sorted based on the test scores.

    Parameters
    ----------
    dict_cv_results : dict
        An attribute of a fitted `GridSearchCV` object.

    metric_name : str
        The name of the metric passed to the `GridSearchCV`'s `scoring`
        parameter.

    metric_pos : bool
        A boolean indicating whether the metric values are supposed to be
        positive. Set to True, this will call `abs()` on the test and train
        scores stored in `dict_cv_results`.

    Returns
    -------
    dfresults : Dataframe
        The dataframe containing train and test scores for each set of
        hyperparameters tried within a `GridSearchCV`. It also contains the
        corresponding hyperparameters.
    '''

    # the cv results
    dfresults = pd.DataFrame(dict_cv_results)

    # only select these columns from `dfresults`
    sel_cols = [F"mean_train_{metric_name}", F"mean_test_{metric_name}"]

    # make the metric values positive if:
    if metric_pos:
        dfresults.loc[:, sel_cols] *= -1

    # find columns starting with `param_`, these correspond to the
    #  hyperparameter of the ML model
    the_params_cols = list(
        filter(re.compile("^param_(.*$)").findall, dfresults.columns)
    )

    # add the above cols to the `sel_cols`
    sel_cols += the_params_cols

    # now limit the dataframe and sort it
    dfresults = dfresults.loc[:, sel_cols]
    dfresults = dfresults.sort_values(by=F"mean_test_{metric_name}")

    return dfresults


def assert_file_exist(the_path: str) -> None:
    '''
    A function that performs an assertion to check whether the provided
    path points to an existing file.

    Parameters
    ----------
    the_path : str, path object
        the path to the file.
    '''

    assert_message = F"`{the_path}` does not point to an existing file!"
    assert os.path.isfile(the_path), assert_message


def assert_dir_exist(abs_path: str) -> None:
    '''
    A function that performs an assertion to check whether the provided
    absolute path points to an existing directory.

    Parameters
    ----------
    abs_path : str
        the absolute path to the directory.
    '''

    assert_message = F"`{abs_path}` does not point to an existing directory!"
    assert os.path.isdir(abs_path), assert_message


def print_dict_nice(
    the_dict: dict, report_header: str = "", kv_gap: int = 20
) -> None:
    '''
    Prints the key, value pair of a dict in a neat way.

    Parameters
    ----------
    the_dict : dict
        The dictionary that we want to print its content.

    report_header : str
        The header that needs to printed before `the_dict` content.

    kv_gap : int
        The gap between the key-value pairs printed out.

    '''

    if report_header:
        print("\n" + report_header + "\n")

    for key, value in the_dict.items():
        message = F"\t{key}: "
        message += (kv_gap - len(key)) * "-" + "> " + F"{value}\n"
        print(message)


def get_optimal_set(gsearch_results: DataFrame, print_set: bool = True) -> dict:
    '''
    Returns the optimal values for the tuning hyperparameters.

    Parameters
    ----------
    gsearch_results : DataFrame
        The sorted dataframe that contains the cv_results_ of a GridSearchCV.

    print_set : bool
        If True, it will print the hyperparameters and their optimal values.
    '''

    # a dict to hold the name and optimal value for the h parameters
    dict_best_params = dict()

    # list of the column-labels corresponding to the tuning hparameters
    list_param_cols = list(
        filter(re.compile("^param_(.*$)").findall, gsearch_results.columns)
    )

    for ith_param in list_param_cols:
        dict_best_params[ith_param.removeprefix(
            "param_")] = gsearch_results.loc[:, ith_param].iloc[0]

    if print_set:
        hdr_message = "The optimal values for the tuning hyperparameters"
        print_dict_nice(dict_best_params, hdr_message)

    return dict_best_params


def remove_file_dir(path: Union[str, Path]) -> None:
    '''
    Removes a file, or a directory and its content.

    Parameters
    ----------
    path : Union[str, Path]
        The absolute path to the directory that should be removed.

    Returns
    -------
    out : None
    '''

    # change it from str to Path if necessary
    if isinstance(path, str):
        path = Path(path)

    # if `path` points to a file, then remove the file.
    if path.is_file():
        try:
            path.unlink()
        except FileNotFoundError:
            pass

        return

    # if `path` points to an existing dir
    if os.path.isdir(path):
        # next line is necessary for Windows users
        os.chmod(path, 0o777)

        # remove the contents (files) of this dir by recursion
        for file_i in path.iterdir():
            remove_file_dir(file_i)

        # now remove the empty dir
        try:
            path.rmdir()
        except FileNotFoundError:
            pass


def check_start_end_date_format(date: str) -> None:
    '''
    Check if the date string is in the desired format.
    The desired format is suitable for `MESH_input_run_options`.
    It has four digits for the year, up to three digits for the day of the year,
    two digits for the hour and miniute of the time step; all separated by '-'.
    Example:
    >>> start_date = "2018-183-04-00"

    Parameters
    ----------
    date : str
        The date value as a string.
    '''

    try:
        pd.to_datetime(date, format="%Y-%j-%H-%M", utc=True)
    except Exception as exc:
        raise ValueError(
            F"Make sure your date (i.e. `{date}`) is a string with the following format: "
            F"%Y-%m-%d %H-%M-%S\n"
        ) from exc


def check_spinup_end_date_format(date: str) -> None:
    '''
    Check if the simulation date string is in the desired format.
    This is the date that marks the beggining of the simulation and
    end of the spin-up period.

    Example:
    >>> sim_date = "2018-07-01 13:00:00"

    Parameters
    ----------
    date : str
        The date value as a string.
    '''

    try:
        pd.to_datetime(date, format="%Y-%m-%d %H:%M:%S", utc=True)

    except Exception as exc:
        raise ValueError(
            F"Make sure your date (i.e. `{date}`) is a string with the following format: "
            F"%Y-%m-%d %H-%M-%S\n"
        ) from exc

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

            case "ksat":
                p_dict[key] = [float(F"{kval:.2e}") for kval in val]

    return p_dict
# ________________________________________________________________ <<< main >>>
