import logging
from pathlib import Path
import itertools
import random

import numpy as np
import pandas as pd
from osgeo import gdal
import matplotlib.colors

import pygeoprocessing

from .config import (
    DEFAULT_GTIFF_CREATION_TUPLE_OPTIONS,
    DEFAULT_LARGEST_BLOCK_SIZE,
)

logger = logging.getLogger(__name__)


class CombineRasterCategories(object):
    """Combine all rasters based on their categorical data, using the minimum number of integers for
    each raster (i.e., XXXyyZ for rasters X with 3 digits of classification, y with 2, and Z with 1).
    Assumes the same nodata value for all raster inputs.
    """

    def __init__(self, nodata_val, *array_values_list, target_nodata_val=None):
        self.nodata = nodata_val
        self.array_values = array_values_list
        self.target_nodata = (
            target_nodata_val if target_nodata_val is not None else nodata_val
        )
        self.reclassification_dict = {
            sum([k * (10**j) for j, k in enumerate(combo)]): i + 1
            for i, combo in enumerate(itertools.product(*array_values_list))
        }

    def __call__(self, *array_list):
        """Create multiple categorical raster overlay with minimum raster int type."""
        valid_mask = np.zeros(array_list[0].shape, dtype=bool)
        temp = np.empty_like(array_list[0], dtype=np.dtype("int16"))
        temp[:] = 0
        num_digits = 0
        for i, array in enumerate(array_list):
            local_valid_mask = ~pygeoprocessing.array_equals_nodata(array, self.nodata)
            temp[local_valid_mask] += array.astype("int")[local_valid_mask] * (
                10**num_digits
            )
            valid_mask &= local_valid_mask  # ANY nodata makes the result nodata
            num_digits += len(str(int(np.max(array[local_valid_mask]))))
        temp[~valid_mask] = self.nodata

        result = np.empty_like(temp)
        result[:] = self.target_nodata
        for key, value in self.reclassification_dict.items():
            index = np.where(temp == key)
            result[index] = value

        return result


def nodata_validation(dtype: np.dtype, nodata_val: int | float) -> int | float:
    """Validate nodata value for numpy dtype, assigning new value if necessary.

    Parameters:
        dtype (np.dtype): The numpy dtype to validate against.
        nodata_val (int | float): The nodata value to validate.

    Returns:
        int | float: The validated nodata value.
    """

    if not np.can_cast(np.min_scalar_type(nodata_val), dtype, casting="same_kind"):
        new_nodata_val = pygeoprocessing.choose_nodata(dtype)
        logger.warning(
            f"Nodata value {nodata_val} cannot be cast to dtype {dtype}, assigning new nodata value ({new_nodata_val})."
        )
        return new_nodata_val
    else:
        return nodata_val


def bin_raster(
    raster_path: str | Path,
    output_raster_path: str | Path,
    category_bins: list[float],
    category_integer_labels: list[int],
    raster_driver_creation_tuple: tuple = DEFAULT_GTIFF_CREATION_TUPLE_OPTIONS,
    largest_block: int = DEFAULT_LARGEST_BLOCK_SIZE,
) -> None:
    """Bin raster into a smaller number of integer classes.


    Parameters:
        raster_path (str): path to input raster
        output_raster_path (str): path to output raster
        category_bins (list[float]): list of bin edges
        category_integer_labels (list[int]): integer bin labels for each bin in category_bins
        raster_driver_creation_tuple (tuple, optional): A (driver_name, creation_options) tuple
            defining the output raster format. Defaults to DEFAULT_GTIFF_CREATION_TUPLE_OPTIONS
        largest_block (int, optional): Largest block size to use for raster calculations. Defaults to DEFAULT_LARGEST_BLOCK_SIZE.

    Returns:
        None
    """

    input_nodata_val = pygeoprocessing.get_raster_info(str(raster_path))["nodata"][0]
    dtype = np.result_type(
        *[np.min_scalar_type(v) for v in category_integer_labels]
    ).type
    output_nodata_val = nodata_validation(dtype, input_nodata_val)
    if output_nodata_val in category_integer_labels:
        logger.warning("Nodata value does not fit within raster value map.")
    gdal_dtype, type_creation_options = (
        pygeoprocessing.geoprocessing._numpy_to_gdal_type(dtype)
    )
    # Check if input raster nodata value is within the range of the output raster type and adapt if necessary

    # Create constraining raster calculator function
    def _bin_op(array, bins):
        """Bin and label rasters."""
        valid_mask = ~pygeoprocessing.array_equals_nodata(array, input_nodata_val)
        result = np.empty_like(array)
        result[:] = output_nodata_val
        result[valid_mask] = np.array(
            pd.cut(array[valid_mask], bins, labels=category_integer_labels)
        )
        return result

    # Replace values in intermediate raster with original raster values where constrained
    logger.debug(f"Binning raster {Path(raster_path).stem}")
    pygeoprocessing.raster_calculator(
        [
            (str(raster_path), 1),
            (category_bins, "raw"),
        ],
        _bin_op,
        str(output_raster_path),
        gdal_dtype,
        output_nodata_val,
        raster_driver_creation_tuple=raster_driver_creation_tuple,
        largest_block=largest_block,
    )


def assign_raster_labels_and_hex_colors(
    raster_path: str | Path,
    category_integer_labels: list[int],
    category_string_labels: list[str],
    category_hex_colors: list[str],
    raster_band_label: str,
) -> None:
    """Assign labels and colors to raster categories.

    Parameters:
        raster_path (str): path to input raster
        category_integer_labels (list): list of category integer labels (must be in ascending order!)
        category_string_labels (list): list of labels for each category value
        category_hex_colors (list): list of colors for each category value
        raster_band_label (int): Band label

    Returns:
        None
    """

    assert len(category_integer_labels) == len(
        category_string_labels
    ), "category_integer_labels and category_string_labels must be the same length"
    assert len(category_integer_labels) == len(
        category_hex_colors
    ), "category_integer_labels and category_hex_colors must be the same length"

    # Get raster info
    raster_info = pygeoprocessing.get_raster_info(str(raster_path))

    # Read raster band
    raster = gdal.OpenEx(str(raster_path), gdal.GA_Update | gdal.OF_RASTER)
    band = raster.GetRasterBand(1)

    # # Raster Attribute Table
    rat = gdal.RasterAttributeTable()

    rat.CreateColumn("Value", gdal.GFT_Integer, gdal.GFU_MinMax)
    rat.CreateColumn("Label", gdal.GFT_String, gdal.GFU_Name)

    for i, (int_label, str_label) in enumerate(
        zip(category_integer_labels, category_string_labels)
    ):
        rat.SetValueAsInt(i, 0, int(int_label))
        rat.SetValueAsString(i, 1, str(str_label))

    band.SetDefaultRAT(rat)

    # GDAL Color Table
    color_table = gdal.ColorTable()

    for int_label, color in zip(category_integer_labels, category_hex_colors):
        color_table.SetColorEntry(
            int(int_label),
            tuple(int(c * 255) for c in matplotlib.colors.to_rgba(color)),
        )

    # set color table and color interpretation
    band.SetRasterColorTable(color_table)
    band.SetRasterColorInterpretation(gdal.GCI_PaletteIndex)

    # set category names
    dtype = raster_info["numpy_type"]
    if dtype != np.uint8:
        logger.warning(
            "Raster data type is not uint8. Assigning category names may increase raster size."
        )

    dtype_min = (
        np.iinfo(dtype).min if np.issubdtype(dtype, np.integer) else np.finfo(dtype).min
    )
    dtype_max = (
        np.iinfo(dtype).max if np.issubdtype(dtype, np.integer) else np.finfo(dtype).max
    )

    label_dict = {
        int_label: str_label
        for int_label, str_label in zip(category_integer_labels, category_string_labels)
    }
    full_labels = [
        label_dict[value] if value in label_dict else str(value)
        for value in range(dtype_min, dtype_max + 1)
    ]
    band.SetRasterCategoryNames(full_labels)

    # set band name
    band.SetDescription(raster_band_label)

    del band, raster


def assign_raster_labels_and_colors_from_csv(
    raster_path: str | Path,
    clr_path: str | Path | pd.DataFrame,
    raster_band_label: str,
    index_column: str = "INDEX",
    r_column: str = "R",
    g_column: str = "G",
    b_column: str = "B",
    alpha_column: str = "ALPHA",
    label_column: str = "LABEL",
) -> None:
    """Assign labels and colors to raster categories from a CSV file.

    Parameters:
        raster_path (str | Path): path to input raster
        clr_path (str | Path | pd.DataFrame): path to CSV file with INDEX, R, G, B, ALPHA, and LABEL,
            or a DataFrame containing these columns
        raster_band_label (str): Band label to assign to raster band
        index_column (str, optional): Name of the index column in the CSV. Defaults to 'INDEX'.
        r_column (str, optional): Name of the red channel column in the CSV. Defaults to 'R'.
        g_column (str, optional): Name of the green channel column in the CSV. Defaults to 'G'.
        b_column (str, optional): Name of the blue channel column in the CSV. Defaults to 'B'.
        alpha_column (str, optional): Name of the alpha channel column in the CSV. Defaults to 'ALPHA'.
        label_column (str, optional): Name of the label column in the CSV. Defaults to 'LABEL'.

    Returns:
        None
    """
    expected_columns = [
        index_column,
        r_column,
        g_column,
        b_column,
        alpha_column,
        label_column,
    ]

    if type(clr_path) is pd.DataFrame:
        clr_df = clr_path
    else:
        clr_df = pd.read_csv(clr_path)
    for col in expected_columns:
        assert (
            col in clr_df.columns
        ), f"CSV File must have column: {col} (provided: {', '.join(clr_df.columns)})"

    # Get raster info
    raster_info = pygeoprocessing.get_raster_info(str(raster_path))

    # Read raster band
    raster = gdal.OpenEx(str(raster_path), gdal.GA_Update | gdal.OF_RASTER)
    band = raster.GetRasterBand(1)

    # Raster Attribute Table
    rat = gdal.RasterAttributeTable()

    rat.CreateColumn("Value", gdal.GFT_Integer, gdal.GFU_MinMax)
    rat.CreateColumn("Label", gdal.GFT_String, gdal.GFU_Name)
    rat.CreateColumn("Red", gdal.GFT_String, gdal.GFU_Red)
    rat.CreateColumn("Green", gdal.GFT_String, gdal.GFU_Green)
    rat.CreateColumn("Blue", gdal.GFT_String, gdal.GFU_Blue)
    rat.CreateColumn("Alpha", gdal.GFT_String, gdal.GFU_Alpha)

    for row in clr_df.itertuples():
        rat.SetValueAsInt(row.Index, 0, int(getattr(row, index_column)))
        rat.SetValueAsString(row.Index, 1, str(getattr(row, label_column)))
        rat.SetValueAsInt(row.Index, 2, int(getattr(row, r_column)))
        rat.SetValueAsInt(row.Index, 3, int(getattr(row, g_column)))
        rat.SetValueAsInt(row.Index, 4, int(getattr(row, b_column)))
        rat.SetValueAsInt(row.Index, 5, int(getattr(row, alpha_column)))

    band.SetDefaultRAT(rat)

    # GDAL Color Table
    color_table = gdal.ColorTable()

    for row in clr_df.itertuples():
        color_table.SetColorEntry(
            int(getattr(row, "INDEX")),
            # getattr(row, "LABEL"),
            (
                int(getattr(row, "R")),
                int(getattr(row, "G")),
                int(getattr(row, "B")),
                int(getattr(row, "ALPHA")),
            ),
        )

    # set color table and color interpretation
    band.SetRasterColorTable(color_table)
    band.SetRasterColorInterpretation(gdal.GCI_PaletteIndex)

    # set category names
    dtype = raster_info["numpy_type"]
    if dtype != np.uint8:
        logger.warning(
            "Raster data type is not uint8. Assigning category names may increase raster size."
        )

    dtype_min = (
        np.iinfo(dtype).min if np.issubdtype(dtype, np.integer) else np.finfo(dtype).min
    )
    dtype_max = (
        np.iinfo(dtype).max if np.issubdtype(dtype, np.integer) else np.finfo(dtype).max
    )

    label_dict = {
        int_label: str_label
        for int_label, str_label in zip(clr_df["INDEX"], clr_df["LABEL"])
    }
    full_labels = [
        label_dict[value] if value in label_dict else str(value)
        for value in range(dtype_min, dtype_max + 1)
    ]
    band.SetRasterCategoryNames(full_labels)

    # set band name
    band.SetDescription(raster_band_label)

    del band, raster



def categorize_and_bin_raster_hex(
    raster_path: str | Path,
    output_raster_path: str | Path,
    category_bins: list[float],
    category_string_labels: list[str],
    category_hex_colors: list[str],
    raster_band_label: str,
    category_integer_labels: list[int] = None,
) -> None:
    """Bin raster into a smaller number of integer classes and assign labels and colors to each class.


    Parameters:
        raster_path (str): path to input raster
        output_raster_path (str): path to output raster
        category_bins (list): list of bin edges
        category_string_labels (list): list of string labels for each bin
        category_hex_colors (list): list of colors for each bin
        raster_band_label (str): label to assign to the raster band
        category_integer_labels (list): (optional) list of integer labels for each bin

    Returns:
        None
    """

    if category_integer_labels is None:
        # Default to categories starting at 1
        category_integer_labels = [
            i + 1 for i in list(range(len(category_string_labels)))
        ]

    bin_raster(raster_path, output_raster_path, category_bins, category_integer_labels)

    assign_raster_labels_and_hex_colors(
        output_raster_path,
        category_integer_labels,
        category_string_labels,
        category_hex_colors,
        raster_band_label,
    )


def random_hex_color(num_colors: int = 1) -> list:
    """Generate a number of random hex colors.

    Parameters:
        num_colors (int): Number of colors to generate, defaults to 1.

    Returns:
        list: A list of random hex color strings.
    """

    return [
        "#" + "".join([random.choice("0123456789ABCDEF") for x in range(6)])
        for i in range(num_colors)
    ]
