
import pandas as pd

from pathlib import Path
from typing import Dict, Any

# ======| File I/O Utilities |==================================================

def load_sample(file: Path) -> pd.DataFrame:
    """
    Loads the raw DNA methylation data of a given sample, provided as a 
    `.parquet` with CpG probe ID as the index and `beta_value` as the column 
    name.

    Parameters
    ----------
    file : Path
        Path to a .parquet file containing the beta values of a sample.

    Returns
    -------
    pd.DataFrame
        Beta values of a sample loaded as a DataFrame.

    Raises
    ------
    FileNotFoundError
        If the file path does not exist or is not `.parquet`.
    """

    # Verify the file exists and is a `.parquet` file
    if not file.exists():
        raise FileNotFoundError(f"File was not found: {file}")
    if file.suffix != ".parquet":
        raise FileNotFoundError(f"File must be a `.parquet`: {file}")
    
    return pd.read_parquet(file)


# ======| Dictionary Utilities |================================================

def merge_dicts(base: Dict[str, Any], override: Dict[str, Any]):
    """
    Recursively merge two dictionaries. Values from `override` take precedence.

    Parameters
    ----------
    base : dict
        Base dictionary (e.g., defaults).
    override : dict
        Dictionary with overriding values (e.g., user config).

    Returns
    -------
    merged : dict
        Merged dictionary with user values overriding defaults.
    """

    result = base.copy()
    for key, value in override.items():
        if (key in result and isinstance(result[key], dict) 
                          and isinstance(value, dict)):
            result[key] = merge_dicts(result[key], value)
        else:
            result[key] = value
    return result


def check_dict(default: Dict[str, Any], 
               user: Dict[str, Any], 
               path: str = "") -> None:
    """
    Recursively verify that a user-provided configuration dictionary matches
    the structure, types, and constraints of a default configuration dictionary.

    Parameters
    ----------
    default : dict
        The default configuration dictionary serving as the schema. Keys and
        value types define what is valid.
    user : dict
        The user-provided configuration dictionary to validate.
    path : str, optional
        Dot-separated path of keys used internally to indicate the location
        in the nested dictionary for informative error messages. Default is "".

    Raises
    ------
    KeyError
        If a required key from `default` is missing in `user`.
    TypeError
        If the type of a value in `user` does not match the type in `default`.
    ValueError
        If a value in `user` violates a constraint (e.g., allowed enum values).
    """

    for key, val in default.items():
        full_key = f"{path}.{key}" if path else key

        # Verify the key exist in the user configurations
        if key not in user:
            raise KeyError(f"Missing key in configuration: {full_key}")
        
        user_val = user[key]

        # If the value itself is a dictionary in the default, recurse
        if isinstance(val, dict):
            if not isinstance(user_val, dict):
                raise TypeError(f"Key '{full_key}' must be a dictionary")
            
            # Do not check metadata as fields will vary
            if key != "metadata":
                check_dict(val, user_val, path = full_key)
                
        else:

            # Verify the value type
            if not isinstance(user_val, type(val)):
                raise TypeError(f"Key '{full_key}' must be of type "
                                f"{type(val).__name__}")
            
            # Verify value ranges individually
            if key == "split":
                if sum(user_val) != 1:
                    raise ValueError(f"Key '{full_key}' must sum to 1.0")
                
            elif key == "missing_threshold" or key == "outlier_threshold":
                if user_val < 0:
                    raise ValueError(f"Key '{full_key}' must be "
                                     "greater or equal to zero")
            
            elif key == "clip_values":
                if user_val[0] < 0 or user_val[1] > 1:
                    raise ValueError(f"Key '{full_key}' must range between "
                                     "zero and one")
