__author__ = ['Gregory A. Greene, map.n.trowel@gmail.com']

import os
import pandas as pd
import numpy as np
from typing import Union, Optional
import unittest

# EQUATION PARAMETERS FOR BIOMASS MODELING
biomass_df = pd.read_csv(os.path.join(os.path.dirname(__file__),
                                      r'Supplementary_Data\Biomass_EquationParameters.csv'),
                         header=0,
                         index_col=False)
biomass_df.set_index('Species', inplace=True)

# GENERAL LIST OF CONIFER SPECIES IN BC
softwood_list = biomass_df[biomass_df['TreeType'] == 'softwood'].index.to_list()

# VARIABLES FOR BIOMASS MODELING
# Tree components list
tree_components_list = ['wood', 'bark', 'branches', 'foliage']

# Tree components dictionary
tree_components_dict = {
    'wood': 0,
    'bark': 1,
    'branches': 2,
    'foliage': 3
}

# Proportions of softwood tree components remaining based on decay classes
softwood_decay_dict = {
    1: [1, 1, 1, 1],
    2: [1, 1, 1, 1],
    3: [1, 1, 1, 0],
    4: [0.95, 0.75, 0.5, 0],
    5: [0.80, 0.5, 0, 0],
    6: [0.66, 0.25, 0, 0],
    7: [0.5, 0.1, 0, 0],
    8: [0.33, 0.05, 0, 0],
    9: [0.05, 0, 0, 0]
}

# Proportions of hardwood tree components remaining based on decay classes
hardwood_decay_dict = {
    1: [1, 1, 1, 1],
    2: [1, 1, 1, 1],
    3: [1, 1, 1, 0],
    4: [0.95, 0.75, 0.5, 0],
    5: [0.66, 0.25, 0, 0],
    6: [0.05, 0, 0, 0]
}


### CALCULATE BIOMASS COMPONENTS OF ENTIRE TREE FOR ALL SPECIES
def getTreeBiomass(spp: str,
                   decayclass: int,
                   components: Union[str, list],
                   dbh: float,
                   height: float = None) -> Union[float, tuple]:
    """
    Function returns species-specific biomass values for indiviudal tree components using
    the Canadian National Biomass equations (Lambert et al. 2005, Ung et al. 2008).
    If only DBH is entered, the "DBH-only" model will be used.
    If both DBH and Height are entered, the DBH-HEIGHT model will be used.
    :param spp: tree species code (uses BC two-letter codes)
    :param decayclass: softwood = 1-9, hardwood = 1-6
    :param components: string (individual) or list (multiple) of tree components ('wood', 'bark', 'branches', 'foliage')
    :param dbh: tree diameter at breast height (cm; value > 0 cm) - REQUIRED INPUT
    :param height: tree height (m; value > 1.3 m)
    :return: biomass (kg) of tree component(s) as either (1) a single value or (2) a tuple, as follows...
        (1) A single value is returned if the input "components" parameter is a string type.
        (2) A list of values is returned if the input "components" parameter is a list.
            Values are returned in the order they are received
    """
    # ### VERIFY INPUTS
    # Check if components object is string or list
    if isinstance(components, str):
        # If it is a string list, make it a list
        component_list = [components]
    elif not isinstance(components, list):
        raise TypeError('"Components" input parameter is not a string or list')
    else:
        component_list = components
    # Check if any input components are invalid
    invalid_components = [comp for comp in component_list if comp not in tree_components_list]
    if invalid_components:
        raise ValueError(f'The following tree components are invalid: {invalid_components}')

    # Bring in global variables
    global biomass_df, softwood_list
    global tree_components_list, tree_components_dict
    global softwood_decay_dict, hardwood_decay_dict

    # Create list to store biomass values
    biomass_list = []

    # Assign decay dictionary based on tree type (softwood or hardwood)
    if spp in softwood_list:
        decay_dict = softwood_decay_dict
    else:
        decay_dict = hardwood_decay_dict

    # Calculate biomass
    for component in component_list:
        if height is None:
            # Biomass for DBH only equations
            biomass = (decay_dict.get(decayclass)[tree_components_dict.get(component)] *
                       biomass_df.loc[spp, f'BIOMASS_DBH-HT_B{component}1'] *
                       pow(dbh, biomass_df.loc[spp, f'BIOMASS_DBH_B{component}2']))
        else:
            # Biomass for DBH + HEIGHT equations
            biomass = (decay_dict.get(decayclass)[tree_components_dict.get(component)] *
                       biomass_df.loc[spp, f'BIOMASS_DBH-HT_B{component}1'] *
                       pow(dbh, biomass_df.loc[spp, f'BIOMASS_DBH-HT_B{component}2']) *
                       pow(height, biomass_df.loc[spp, f'BIOMASS_DBH-HT_B{component}3']))

        if isinstance(components, str):
            return biomass
        else:
            biomass_list += [biomass]

    return tuple(biomass_list)


# GET BIOMASS VALUES FOR PHOTOLOAD PLANT PHOTO SERIES
def getPhotoloadBiomass(pl_code: str,
                        pct_cvr: float,
                        height: float = None) -> float:
    """
    Function to estimate plant biomass using equations derived for the Photoload plant species
    :param pl_code: photoload species code
    :param pct_cvr: percent cover of species (0-100)
    :param height: height of species (cm) - if no height entered, default Photoload heights will be used
    :return: biomass (kg/m2) - returns value of 0 if pl_code is invalid
    """
    # Return value of 0 if either percent cover or height are 0
    if (pct_cvr == 0) or (height == 0):
        return 0

    # Dictionary of default Photoload plant species heights
    height_dict = {
        'AMAL': 35.56,
        'BERE': 10.16,
        'SYAL': 35.56,
        'VAGL': 35.56,
        'VASC': 17.78,
        'ARLA': 30.48,
        'CARU': 20.3,
        'FESC': 30.5,
        'XETE': 25.4
    }

    # Use default photoload plant height if height is missing or not a finite number
    # Assign height = 0 if species code is invalid
    if (not np.isfinite(height)) or (height is None):
        height = height_dict.get(pl_code, None)  # Returns NONE if species code is invalid

    # If height is NONE, assign height as 0 and throw invalid species code warning
    if height is None:
        import warnings
        warnings.warn(f'Photoload species code "{pl_code}" is invalid. Biomass returned as 0.')
        return 0

    # Dictionary of photoload equations for each species code
    photoload_dict = {
        'AMAL': (height / 35.56) * 0.0148 * np.exp(0.0454 * pct_cvr),
        'BERE': (height / 10.16) * 0.0013 * pct_cvr,  # BERE is an old species code
        'MARE': (height / 10.16) * 0.0013 * pct_cvr,  # This is the replacement species code for BERE
        'SYAL': (height / 35.56) * 0.0107 * np.exp(0.0442 * pct_cvr),
        'VAGL': (height / 35.56) * 0.0052 * pct_cvr,
        'VASC': (height / 17.78) * 0.0135 * np.exp(0.035 * pct_cvr),
        'ARLA': (height / 30.48) * (0.000008 * np.power(pct_cvr, 2) + (0.0006 * pct_cvr) + 0.0022),
        'CARU': (height / 20.30) * 0.0194 * np.exp(0.0282 * pct_cvr),
        'FESC': (height / 30.5) * 0.0188 * np.exp(0.0298 * pct_cvr),
        'XETE': (height / 25.4) * 0.0154 * np.exp(0.0444 * pct_cvr)
    }
    return photoload_dict.get(pl_code, np.nan)


def getDuffLitterBiomass(spp: Union[str, list],
                         pct_list: Optional[list] = None,
                         return_type: str = 'bulk_density',
                         duff_depth: Optional[Union[int, float]] = None,
                         litter_depth: Optional[Union[int, float]] = None) -> Union[float, tuple]:
    """
    Function to return bulk density values or calculate biomass for duff and litter.

    Bulk density values are derived from Woodall and Monleon (2008). Species codes provided to this function
    are associated with the various forest types provided in Woodall and Monleon. The general approach should
    be to obtain a weighted average value of the duff/litter bulk density values based on the proportion
    of each species in a stand. Then calculate biomass by multiplying the weighted average bulk density
    value by the total duff/litter depth.
    Biomass is calculated as: bulk_density * depth
    :param spp: a single tree species code (uses BC two-letter codes), or a list of codes
    :param pct_list: if "spp" is a list of species codes, "pc_list" must be a list of
        percentages of each species in the list of species codes (range: 0-100). The sum of values must = 100.
    :param return_type: type of value returned (string: "bulk_density" or "biomass")
    :param duff_depth: depth of duff (cm)
    :param litter_depth: depth of litter (cm)
    :return: If return_type = "bulk_density", return duff and litter bulk density (kg/m3) as tuple (duff, litter).
        If return_type = "biomass", return duff and/or litter biomass (kg/m2).
            If both duff and litter depths are provided, biomass will be returned as a tuple (e.g., (duff, litter)).
            If a single duff or litter depth is provided, a single biomass value will be returned.
    """
    # ### VERIFY INPUTS
    # spp
    if isinstance(spp, list):
        for s in spp:
            if not isinstance(s, str):
                raise TypeError('Each species code in the "spp" list must be a string data type')
        if pct_list is None:
            raise ValueError('The "pct_list" parameter must contain a list of of percentages for each '
                             'species provided in the "spp" list')
    # pct_list
    if isinstance(pct_list, list):
        for pct in pct_list:
            if not isinstance(pct, (float, int)):
                raise TypeError('Each value in "pct_list" must be int or float data type')
        if round(sum(pct_list), 0) != 100:
            raise ValueError('The values in "pct_list" must sum to 100')
    # return_type
    if return_type not in ['bulk_density', 'biomass']:
        raise ValueError('The "return_type" parameter must be either "bulk_density" or "biomass"')
    elif return_type == 'biomass':
        if all(depth is None for depth in (duff_depth, litter_depth)):
            raise ValueError('Either "duff_depth", "litter_depth", or both parameters must be provided '
                             'when return_type = "biomass"')
    # duff_depth
    if duff_depth is not None:
        if not isinstance(duff_depth, (int, float)):
            raise TypeError('The "duff_depth" parameter must be either int or float type, and in units of cm')
        duff_depth /= 100  # Convert from cm to meters
    # litter_depth
    if litter_depth is not None:
        if not isinstance(litter_depth, (int, float)):
            raise TypeError('The "litter_depth" parameter must be either int or float type, and in units of cm')
        litter_depth /= 100  # Convert from cm to meters

    global biomass_df

    # Get bulk density values from dataframe
    if isinstance(spp, list):
        # Precompute biomass and prop lists
        bd_list = biomass_df.loc[spp, ['DUFF_BD', 'LITTER_BD']].to_numpy()
        d_bd_total, l_bd_total = 0, 0

        # Iterate over precomputed lists with zip
        for (d_bd, l_bd), prop in zip(bd_list, pct_list):
            d_bd_total += d_bd * prop / 100
            l_bd_total += l_bd * prop / 100
    else:
        d_bd_total, l_bd_total = biomass_df.loc[spp, ['DUFF_BD', 'LITTER_BD']].values

    # Return values
    if return_type == 'bulk_density':
        return float(d_bd_total), float(l_bd_total)
    elif return_type == 'biomass':
        # Check if components object is string or list
        if duff_depth and litter_depth:
            biomass = tuple(np.multiply(
                [d_bd_total, l_bd_total],
                [duff_depth, litter_depth]).tolist())
        elif duff_depth:
            biomass = d_bd_total * duff_depth
        elif litter_depth:
            biomass = l_bd_total * litter_depth
        else:
            raise Exception('Neither duff or litter depths were provided.')

    return biomass


if __name__ == '__main__':
    # Test getDuffLitterBiomass
    print('Testing the getDuffLitterBiomass() function...')
    print('Bulk Density:', getDuffLitterBiomass(spp='PY', return_type='bulk_density'))
    print('Bulk Density:', getDuffLitterBiomass(spp=['PY', 'FDI', 'PLI'],
                                                 pct_list=[60, 30, 10],
                                                 return_type='bulk_density'))
    print('Biomass:', getDuffLitterBiomass(spp='PY', return_type='biomass', litter_depth=0.7))
    print('Biomass:', getDuffLitterBiomass(spp='PY', return_type='biomass', duff_depth=3.4))
    print('Biomass:', getDuffLitterBiomass(spp='PY', return_type='biomass', litter_depth=0.7, duff_depth=3.4))
    print('Biomass:', getDuffLitterBiomass(spp=['PY', 'FDI', 'PLI'],
                                            pct_list=[60, 30, 10],
                                            return_type='biomass',
                                            litter_depth=0.7,
                                            duff_depth=3.4))
    try:
        print('Testing scenario where pct_list does not sum to 100')
        getDuffLitterBiomass(spp=['PY', 'FDI', 'PLI'],
                             pct_list=[60, 30, 1],
                             return_type='biomass',
                             litter_depth=0.7,
                             duff_depth=3.4)
    except ValueError as err:
        print(f'\tValueError: {err}')
        pass
