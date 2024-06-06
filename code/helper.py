
from cartopy.util import add_cyclic_point
import xarray as xr
import numpy as np

def find_varname_from_attribute(ds, attribute, pattern):
    """
    Find the variable name in a xarray dataset that has a specific unit.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset to search for variable name.
    unit : str
        Unit to search for.

    Returns
    -------
    str
        Variable name with the specific unit.
    """
    var_name = None
    for var_name, var_data in ds.variables.items():
        if var_data.attrs.get(attribute) == pattern:
            return var_name
    return var_name

def find_varname_from_keywords(ds, keywords):
    """
    Find the variable name in a xarray dataset that contains specific keywords in its long name.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset to search for variable name.
    keywords : list of str
        Keywords to search for in the variable's long name.

    Returns
    -------
    str
        Variable name where the long name contains all the specified keywords, or None if no such variable exists.
    """
    for var_name, var_data in ds.variables.items():
        # Retrieve the long name attribute of the variable
        long_name = var_data.attrs.get('long_name', '').lower()
        standard_name = var_data.attrs.get('standard_name', '').lower()

        # Check if any keyword is present in the standard name
        if any(keyword.lower() in standard_name for keyword in keywords):
            return var_name
        # Check if any keyword is present in the long name
        if any(keyword.lower() in long_name for keyword in keywords):
            return var_name

    return None  # Return None if no variable matches all keywords


def find_geo_coords(ds):
    """
    Find the latitude and longitude variable names in a xarray dataset.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset to search for latitude and longitude variable names.

    Returns
    -------
    tuple
        A tuple containing the longitude and latitude variable names.
    
    Raises
    ------
    ValueError
        If the latitude or longitude variable names cannot be determined automatically.
    """
    # Common names for latitude and longitude variables
    possible_lat_names = ['latitude', 'lat', 'LAT', 'Latitude']
    possible_lon_names = ['longitude', 'lon', 'LON', 'Longitude']

    lat_name = None
    lon_name = None

    # Find latitude and longitude in dataset
    for var in ds.coords:
        if var in possible_lat_names:
            lat_name = var
        if var in possible_lon_names:
            lon_name = var

    if lat_name is None or lon_name is None:
        raise ValueError("Could not automatically determine the latitude or longitude variable names.")
    
    return lon_name, lat_name

def area_weighted_global_mean(da):
    """
    Calculate the area-weighted global mean value of a 2D (lat x lon) xarray DataArray.

    Parameters
    ----------
    da : xarray.DataArray
        2D DataArray with latitude and longitude coordinates.

    Returns
    -------
    float
        Area-weighted global mean value.
    """
    
    lon_name, lat_name = find_geo_coords(da)
    
    # Convert latitude to radians
    lat_radians = np.deg2rad(da[lat_name])
    
    # Calculate the weights (area of each grid cell)
    weights = np.cos(lat_radians)
    
    # Ensure the weights sum to 1
    weights /= weights.sum(dim=lat_name)
    
    # Use the xr.DataArray.weighted method to apply weights and compute mean
    weighted_mean = da.weighted(weights).mean(dim=[lat_name, lon_name]).values
    
    return weighted_mean

def calc_tpw(q):

    # get dimension names
    vert_name = find_varname_from_attribute(q, 'axis', 'Z')

    g = 9.81  # Gravity in m/s^2

    # Calculate the pressure differences between successive levels using xarray's diff method
    dp = q[vert_name].diff(vert_name)

    if (q[vert_name].attrs.get('units') in ['hPa','mbar']):
        dp = dp * 100

    # Calculate the average specific humidity at the interfaces between pressure levels using xarray's rolling method
    q_avg = q['hus'].rolling({vert_name:2}, center=True).mean()

    # Truncate q_avg to match the size of dp (last value of q_avg will be NaN because there's no subsequent level)
    q_avg = q_avg.isel({vert_name:slice(None, -1)})


    # Integrate vertically
    # Correct the integration direction for increasing altitude with -1
    TPW = ((q_avg * dp / g).sum(vert_name)) * -1
    return TPW

def add_cyclic_data(data):
    lon_name, lat_name = find_geo_coords(data)
    # Add a cyclic point along the longitude axis and update the longitude coordinate
    data_cyclic, lon_cyclic = add_cyclic_point(data.values, coord=data[lon_name])
    # Update the DataArray with the new cyclic data and coordinate
    return xr.DataArray(data_cyclic, dims=data.dims, coords={**data.coords, lon_name: lon_cyclic})