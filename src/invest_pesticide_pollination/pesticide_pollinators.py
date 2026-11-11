import collections
import itertools
import logging
import os
from pathlib import Path
import importlib

import matplotlib
import yaml
import numpy as np
import pandas as pd
from osgeo import gdal
import pygeoprocessing
import natcap.invest.pollination
import natcap.invest.utils
from natcap.invest import validation

import nci_pollination as nci_pollination

from .model_spec import MODEL_SPEC
from .utils import (
    assign_raster_labels_and_colors_from_csv,
    assign_raster_labels_and_hex_colors,
    categorize_and_bin_raster_hex,
    random_hex_color,
)

logger = logging.getLogger(__name__)

pathlike = str | Path

_INDEX_NODATA = natcap.invest.pollination._INDEX_NODATA

SIGNED_NODATA_VAL = -1
UNSIGNED_NODATA_VAL = 0

# These patterns are expected in the biophysical table
_PESTICIDE_APPLICATION_RATE_FIELD = "pesticide_application_rate"

# These patterns are expected in the guild table
_TOXICITY_THRESHOLD_FIELD = "toxicity_threshold"
_POLLEN_CONSUMPTION_RATE_FIELD = "pollen_consumption_rate"

# These are patterns expected in the decision table
_DECISION_CODE_HEADER = "decision_code"
_DECISION_LABEL_HEADER = "decision_label"
_DECISION_COLOR_HEADER = "decision_color"
_EXPECTED_DECISION_HEADERS = [
    _DECISION_CODE_HEADER,
    _DECISION_LABEL_HEADER,
    _DECISION_COLOR_HEADER,
]


# Raster Labels
_POLLINATOR_MV_DICT = {
    "thresholds": (0.02, 0.04),  # TODO Check these values
    "values": [1, 2, 3],
    "labels": [
        "Low Restoration Potential",
        "Medium Restoration Potential",
        "High Restoration Potential",
    ],
    "colors": ["#e5f5e0", "#a1d99b", "#31a354"],
}
_POLLINATOR_LANDSCAPE_SCORE_DICT = {
    "thresholds": (0.1, 0.3),
    "values": [1, 2, 3],
    "labels": ["Low Habitat Quality", "Medium Habitat Quality", "High Habitat Quality"],
    "colors": ["#e5f5e0", "#a1d99b", "#31a354"],
}
_TOX_DICT = {
    "values": [1, 2, 3],
    "labels": ["Low Toxicity", "Medium Toxicity", "High Toxicity"],
    "colors": ["#fee0d2", "#fc9272", "#de2d26"],
}

COMBINATION_DICT = {
    "categories": [],
    "complex_values": [],
    "values": [],
}
for i, j in enumerate(
    itertools.product(
        _POLLINATOR_MV_DICT["values"],
        _POLLINATOR_LANDSCAPE_SCORE_DICT["values"],
        _TOX_DICT["values"],
    )
):
    COMBINATION_DICT["categories"].append(j)
    COMBINATION_DICT["values"].append(i + 1)
    COMBINATION_DICT["complex_values"].append(
        sum([l * (10 ** (len(j) - k - 1)) for k, l in enumerate(j)])
    )

COMBINATION_DICT["labels"] = [
    " | ".join(str(j) for j in i)
    for i in itertools.product(
        _POLLINATOR_MV_DICT["labels"],
        _POLLINATOR_LANDSCAPE_SCORE_DICT["labels"],
        _TOX_DICT["labels"],
    )
]

# Determine appropriate raster dtype for combination raster
DTYPE = np.result_type(
    *[np.min_scalar_type(v) for v in COMBINATION_DICT["values"]]
).type
if DTYPE in [np.float16]:
    DTYPE = np.float32  # NOTE: GDAL<3.11 doesn't support float16 rasters
GDAL_TYPE, _ = pygeoprocessing.geoprocessing._numpy_to_gdal_type(DTYPE)


def _parse_scenario_variables(
    args,
):
    """Parse out scenario variables from input parameters.

    This function parses through the chemical compound table to generate

    Parameter:
        args (dict): this is the args dictionary passed in to the `execute`
            function, requires a 'guild_table_path' key.

    Returns:
        A dictionary with the keys:
            * season_list (list of string)
            * substrate_list (list of string)
            * species_list (list of string)
            * alpha_value[species] (float)
            * landcover_substrate_index[substrate][landcover] (float)
            * landcover_floral_resources[season][landcover] (float)
            * species_abundance[species] (string->float)
            * species_foraging_activity[(species, season)] (string->float)
            * species_substrate_index[species][substrate] (string->float)
            * foraging_activity_index[(species, season)] (tuple->float)
    """

    scenario_variables = nci_pollination.nci_pollination._parse_scenario_variables(args)

    # Create reclass dicts for each of the yield-related fields, mapping LULC codes to the appropriate values from the biophysical table
    landcover_biophysical_df = MODEL_SPEC.get_input(
        "landcover_biophysical_table_path"
    ).get_validated_dataframe(args["landcover_biophysical_table_path"], args=args)

    # Add pesticide application rate to scenario variables dict
    scenario_variables[_PESTICIDE_APPLICATION_RATE_FIELD] = collections.defaultdict(
        dict
    )
    for landcover_id, row in landcover_biophysical_df.iterrows():
        scenario_variables[_PESTICIDE_APPLICATION_RATE_FIELD][landcover_id] = row[
            _PESTICIDE_APPLICATION_RATE_FIELD
        ]

    guild_df = MODEL_SPEC.get_input("guild_table_path").get_validated_dataframe(
        args["guild_table_path"],
        args=args,
    )

    # Ingest fields from guild table for toxicity thresholds and exposure rates, mapping species to the appropriate values from the guild table
    for field in [
        _TOXICITY_THRESHOLD_FIELD,
        _POLLEN_CONSUMPTION_RATE_FIELD,
    ]:
        scenario_variables[field] = collections.defaultdict(dict)
        for species in scenario_variables["species_list"]:
            scenario_variables[field][species] = guild_df[field][species]

    # Create reclass dicts for decision table
    decision_csv = str(
        importlib.resources.files("invest_pesticide_pollination") / "decision_table.csv"
    )
    decision_df = MODEL_SPEC.get_input("decision_table_path").get_validated_dataframe(
        decision_csv,
        args=args,
    )

    for field in _EXPECTED_DECISION_HEADERS:
        scenario_variables[field] = collections.defaultdict(dict)
        for landcover_id, row in decision_df.iterrows():
            scenario_variables[field][landcover_id] = row[field]

    # * Decision reclass dict (dict of int->string)
    scenario_variables["decision_reclass_dict"] = {
        int(i): None for i in sorted(decision_df.index)
    }
    for decision_code, row in decision_df.iterrows():
        scenario_variables["decision_reclass_dict"][int(decision_code)] = int(
            row[_DECISION_CODE_HEADER]
        )
    scenario_variables["decision_df"] = (
        decision_df.copy()[_EXPECTED_DECISION_HEADERS]
        .drop_duplicates()
        .sort_values(by=[_DECISION_CODE_HEADER])
        .reset_index(drop=True)
    )
    scenario_variables["decision_df"][_DECISION_CODE_HEADER] = pd.to_numeric(
        scenario_variables["decision_df"][_DECISION_CODE_HEADER], downcast="integer"
    )

    return scenario_variables


def execute(args: dict):
    # TODO Fix the docstrings
    """Plug-in executable expanding the InVEST Pollination model to assess pesticide exposure
    risk to pollinators alongside their current and potential habitat quality. Calculates
    current pollinator habitat quality, potential restoration impacts on habitat quality, and
    pesticide exposure risk to pollinators, creating a synthesis layer of recommended conservation
    actions based on the combination of these factors.

    Parameters:
         args (dict): a dictionary containing the following keys and values:

        **InVEST args:**

        args['workspace_dir'] (string): a path to the output workspace folder.
            Will overwrite any files that exist if the path already exists.
        args['results_suffix'] (string): string appended to each output
            file path.
        args['landcover_raster_path'] (string): file path to a landcover
            raster.
        args['guild_table_path'] (string): file path to a table indicating
            the bee species to analyze in this model run.  Table headers
            must include:

                * 'species': a unique name or identifier for each pollinator
                    species or guild of interest. Column string names will
                    be referred to in other tables and the model will output
                    analyses per species.
                * one or more columns matching
                    'nesting_suitability_[SUBSTRATE]_index' with values in
                    the range [0, 1] indicating the utilization of the given
                    substrate by this species, where 1 indicates the nesting
                    substrate is fully utilized and 0 indicates it is not
                    utilized at all. Replace [SUBSTRATE] with substrate
                    names matching those in the biophysical table.
                * one or more columns matching
                    'foraging_activity_[SEASON]_index' with values in the
                    range [0, 1] indicating the relative level of pollinator
                    activity for this species/guild during a particular
                    season. Replace [SEASON] with season names matching
                    those in the biophysical table.
                * 'alpha': the average distance, in meters, that this
                    species or guild travels to forage on flowers.
                * 'relative_abundance': the proportion of total pollinator
                    abundance that consists of this species/guild.
                * 'toxicity_threshold': the exposure to the pesticide that
                    results in toxicity (LD50, LOAEL, etc) for this species
                    or guild, in mg/kg. Required if
                    ``args['pesticide_analysis']`` is True.
                * 'pollen_consumption_rate': how much pollen this
                    species/guild eats per day, in mg/day. Required if
                    ``args['pesticide_analysis']`` is True.

        args['landcover_biophysical_table_path'] (string): path to a table
            mapping landcover codes in `args['landcover_raster_path']` to
            indexes of nesting availability for each nesting substrate
            referenced in the guild table, indexes of abundance of floral
            resources on that landcover type per season referenced in the
            guild table, crop pollination dependence and value where
            appropriate, and the pesticide application rate for that
            landcover type.

            Columns in the table must include:
                * 'lucode': representing all the unique landcover codes in
                    the raster at `args['landcover_raster_path']`.
                * one or more columns matching
                    'nesting_[SUBSTRATE]_availability_index' with values in
                    the range [0, 1], one for each substrate referenced in
                    the guild table.
                * one or more columns matching
                    'floral_resources_[SEASON]_index' with values in the
                    range [0, 1], one for each season referenced in the
                    guild table.
                * 'crop_pollinator_dependence_index': proportion of annual
                    crop yield dependent on pollinators, in the range
                    [0, 1]. Required if ``args['crop_yield']`` is True.
                * 'half_saturation_coefficient': half saturation coefficient
                    for crop yield modeling, in the range [0, 1]. Required
                    if ``args['crop_yield']`` is True.
                * 'crop_value_per_hectare': annual crop yield value per
                    hectare. Required if ``args['crop_yield']`` is True.
                * 'pesticide_application_rate': annual pesticide
                    application rate, in kg/ha/year. Required if
                    ``args['pesticide_analysis']`` is True.

        **Marginal Value args:**

        args['scenario_analysis'] (boolean): whether to perform scenario analysis.
            If False (or absent), the model returns after running the baseline InVEST
            model. If True, the model will run the baseline InVEST model and then the
            marginal value model for the scenario raster provided.
        args['scenario_landcover_raster_path'] (string): file path to a landcover
            raster.
        args['scenario_suffix'] (string): suffix to add to output files for the
            scenario being evaluated. Appends after the 'results_suffix' from
            the original model.

        **Pesticide model args:**

        args['pesticide_analysis'] (boolean): whether to calculate
            pesticide risk for pollinators from the scenario raster. If
            False (or absent), the model returns after running the NCI
            Pollination marginal value model.
        args['pesticide_name'] (string): name of the pesticide to use in
            the analysis. Required if ``args['pesticide_analysis']`` is
            True.
        args['drift_distance_m'] (float): average distance, in meters,
            that the pesticide drifts via wind. Required if
            ``args['pesticide_analysis']`` is True.

        **Landscape Decision Analysis args:**

        args['landscape_decision_analysis'] (boolean): whether to perform
            landscape decision analysis for pollinators from the scenario
            raster, producing a combined decision raster from habitat
            quality, restoration potential, and pesticide toxicity.
            Requires ``args['pesticide_analysis']`` to be True.

    """

    # Run NCI Pollination Marginal Value model
    args, file_registry, _ = MODEL_SPEC.setup(args)

    # Get the file registry for the baseline pollination model and the marginal value
    _, pollination_file_registry, _ = nci_pollination.MODEL_SPEC.setup(args)

    # Get the file registry for the restoration scenario
    mv_args = nci_pollination.combine_results_and_scenario_suffixes(args)
    pollination_mv_args, pollination_mv_file_registry, _ = (
        nci_pollination.MODEL_SPEC.setup(mv_args)
    )
    del mv_args

    # Test if pollination model has already been run for this species
    pol_scenario_variables = nci_pollination.nci_pollination._parse_scenario_variables(
        pollination_mv_args
    )

    required_files = [
        pollination_mv_file_registry[key, species]
        for key in [
            "pollinator_supply_[SPECIES]",
            "d_pollinator_supply_[SPECIES]",
        ]
        for species in pol_scenario_variables["species_list"]
    ] + [
        pollination_file_registry[
            "kernel_[ALPHA]",
            f"{(pol_scenario_variables['alpha_value'][species] / pol_scenario_variables['landcover_mean_pixel_size']):.6f}",
        ]
        for species in pol_scenario_variables["species_list"]
    ]

    if not all([os.path.exists(f) for f in required_files]):
        logger.info("Pollination model outputs not found, running model")
        for f in required_files:
            if not os.path.exists(f):
                logger.debug(f"Missing file: {f}")
        nci_pollination.execute(args)

    # If pesticide analysis is not requested, we can exit early after running the
    # NCI Pollination model. Otherwise, we continue to process the pesticide analysis
    if "pesticide_analysis" not in args or not args["pesticide_analysis"]:
        return

    pesticide = args["pesticide_name"]

    # If decision analysis is requested, add relevant arguments to args
    if "landscape_decision_analysis" in args and args["landscape_decision_analysis"]:
        # Hard code pollination thresholds for decision analysis
        args["pollination_low_threshold"] = (
            _POLLINATOR_LANDSCAPE_SCORE_DICT["thresholds"][0]
            if (
                "pollination_low_threshold" not in args
                or args["pollination_low_threshold"] is None
            )
            else args["pollination_low_threshold"]
        )
        args["pollination_high_threshold"] = (
            _POLLINATOR_LANDSCAPE_SCORE_DICT["thresholds"][1]
            if (
                "pollination_high_threshold" not in args
                or args["pollination_high_threshold"] is None
            )
            else args["pollination_high_threshold"]
        )

        # Hard code restoration thresholds for decision analysis
        args["restoration_low_threshold"] = (
            _POLLINATOR_MV_DICT["thresholds"][0]
            if (
                "restoration_low_threshold" not in args
                or args["restoration_low_threshold"] is None
            )
            else args["restoration_low_threshold"]
        )
        args["restoration_high_threshold"] = (
            _POLLINATOR_MV_DICT["thresholds"][1]
            if (
                "restoration_high_threshold" not in args
                or args["restoration_high_threshold"] is None
            )
            else args["restoration_high_threshold"]
        )

    # Get scenario variables such as drift distance, toxicity thresholds, and exposure rates
    scenario_variables = _parse_scenario_variables(args)

    # Create pesticide application raster
    reclass_error_details = {
        "raster_name": "LULC",
        "column_name": "lucode",
        "table_name": "Biophysical",
    }
    logger.debug(f"Reclassifying landcover to pesticide application rate")
    natcap.invest.utils.reclassify_raster(
        (args["landcover_raster_path"], 1),
        scenario_variables[_PESTICIDE_APPLICATION_RATE_FIELD],
        file_registry["pesticide_application_rate_[PESTICIDE]", pesticide],
        gdal.GDT_Float32,
        _INDEX_NODATA,
        reclass_error_details,
    )

    # Calculate pesticide application per pixel
    # NOTE this is not in Eric's original code, we should decide if we want per-ha or per-pixel
    logger.debug(f"Calculating pesticide application per pixel")
    lulc_raster_info = pygeoprocessing.get_raster_info(args["landcover_raster_path"])
    cell_area_ha = abs(np.prod(lulc_raster_info["pixel_size"])) * 0.0001
    pygeoprocessing.raster_map(
        op=lambda application_rate: application_rate * cell_area_ha,
        rasters=[
            file_registry["pesticide_application_rate_[PESTICIDE]", pesticide],
        ],
        target_path=file_registry["pesticide_application_[PESTICIDE]", pesticide],
        target_dtype=np.float32,
        target_nodata=_INDEX_NODATA,
    )

    # Calculate pesticide drift
    logger.debug(
        f"Convolving pesticide application with drift kernel to calculate pesticide drift"
    )
    drift_distance_pixels = (
        args["drift_distance_m"] / scenario_variables["landcover_mean_pixel_size"]
    )
    pygeoprocessing.kernels.exponential_decay_kernel(
        target_kernel_path=file_registry["drift_kernel_[PESTICIDE]", pesticide],
        expected_distance=drift_distance_pixels,
        max_distance=drift_distance_pixels * 5,
    )
    pygeoprocessing.convolve_2d(
        signal_path_band=(
            file_registry["pesticide_application_[PESTICIDE]", pesticide],
            1,
        ),
        kernel_path_band=(file_registry["drift_kernel_[PESTICIDE]", pesticide], 1),
        target_path=file_registry["pesticide_drift_[PESTICIDE]", pesticide],
        target_datatype=gdal.GDT_Float32,
        target_nodata=_INDEX_NODATA,
    )

    # Calculate pesticide exposure per species based on drift and species-specific exposure rates
    for species in scenario_variables["species_list"]:
        logger.debug(f"Calculating initial pesticide exposure via foraging ({species})")
        alpha = (
            pol_scenario_variables["alpha_value"][species]
            / pol_scenario_variables["landcover_mean_pixel_size"]
        )
        # NOTE The original equation was weighted by floral resources, but here we simplify to drift by species foraging radius
        # Eric's original equation:
        # exposure = sumproduct(drift, floral_resources, kernel) / (supply / nesting) * realized_exposure

        pygeoprocessing.convolve_2d(
            signal_path_band=(
                file_registry["pesticide_drift_[PESTICIDE]", pesticide],
                1,
            ),
            kernel_path_band=(
                pollination_file_registry["kernel_[ALPHA]", f"{alpha:6f}"],
                1,
            ),
            target_path=file_registry[
                "landscape_pesticide_load_[PESTICIDE]_[SPECIES]", pesticide, species
            ],
            target_datatype=gdal.GDT_Float32,
            target_nodata=_INDEX_NODATA,
        )

        logger.debug(
            f"Calculating expected pesticide exposure using pollination consumption ({species})"
        )

        def exposure_op(
            initial_exposure,
            pollen_consumption_rate=scenario_variables["pollen_consumption_rate"][
                species
            ],
            exposure_coef=0.31,  # NOTE Exposure coefficient may eventually be by chemcial, but for now is a single hard-coded value
            load_intercept=1,
        ):
            pollen_concentration = (
                initial_exposure + load_intercept
            ) ** exposure_coef - load_intercept
            pesticide_consumption = pollen_concentration * pollen_consumption_rate
            return pesticide_consumption

        pygeoprocessing.raster_map(
            op=exposure_op,
            rasters=[
                file_registry[
                    "landscape_pesticide_load_[PESTICIDE]_[SPECIES]", pesticide, species
                ]
            ],
            target_path=file_registry[
                "pesticide_exposure_[PESTICIDE]_[SPECIES]", pesticide, species
            ],
            target_dtype=np.float32,
            target_nodata=_INDEX_NODATA,
        )

        logger.debug(f"Calculating expected pesticide hazard for {species}")
        pygeoprocessing.raster_map(
            op=lambda exposure: exposure
            / scenario_variables[_TOXICITY_THRESHOLD_FIELD][species],
            rasters=[
                file_registry[
                    "pesticide_exposure_[PESTICIDE]_[SPECIES]", pesticide, species
                ],
            ],
            target_path=file_registry[
                "pesticide_hazard_[PESTICIDE]_[SPECIES]", pesticide, species
            ],
            target_dtype=np.float32,
            target_nodata=_INDEX_NODATA,
        )

    # If the landscape decision matrix step is not requested, skip the rest of this section
    if (
        "landscape_decision_analysis" not in args
        or not args["landscape_decision_analysis"]
    ):
        return

    # Create cmap table for eventual decision raster
    cmap_dict = {
        "INDEX": [],
        "LABEL": [],
        "Description": [],
        "R": [],
        "G": [],
        "B": [],
        "ALPHA": [],
    }
    for index, label, color in zip(
        scenario_variables["decision_df"][_DECISION_CODE_HEADER],
        scenario_variables["decision_df"][_DECISION_LABEL_HEADER],
        scenario_variables["decision_df"][_DECISION_COLOR_HEADER],
    ):
        r, g, b, a = tuple(int(c * 255) for c in matplotlib.colors.to_rgba(color))
        cmap_dict["INDEX"].append(index)
        cmap_dict["LABEL"].append(label)
        cmap_dict["Description"].append("")
        cmap_dict["R"].append(r)
        cmap_dict["G"].append(g)
        cmap_dict["B"].append(b)
        cmap_dict["ALPHA"].append(a)
    cmap_df = pd.DataFrame(cmap_dict)

    cmap_df.to_csv(file_registry["decision_matrix_colormap"], index=False)

    # Iterate through species and categorize rasters for each species
    for species in scenario_variables["species_list"]:

        # Categorize baseline InVEST outputs into High, Medium, and Low Quality pollinator habitat
        logger.debug(f"{species} | Categorizing habitat quality")

        categorize_and_bin_raster_hex(
            pollination_file_registry["pollinator_supply_[SPECIES]", species],
            file_registry["habitat_quality_[SPECIES]", species],
            [
                -float("inf"),
                args["pollination_low_threshold"],
                args["pollination_high_threshold"],
                float("inf"),
            ],
            _POLLINATOR_LANDSCAPE_SCORE_DICT["labels"],
            _POLLINATOR_LANDSCAPE_SCORE_DICT["colors"],
            "Habitat Quality",
            category_integer_labels=_POLLINATOR_LANDSCAPE_SCORE_DICT["values"],
        )

        # Categorize Marginal Value outputs into High, Medium, and Low Potential pollinator habitat
        logger.debug(f"{species} | Categorizing pollination marginal value")

        categorize_and_bin_raster_hex(
            pollination_mv_file_registry["d_pollinator_supply_[SPECIES]", species],
            file_registry["marginal_pollinator_value_[SPECIES]", species],
            [
                -float("inf"),
                args["restoration_low_threshold"],
                args["restoration_high_threshold"],
                float("inf"),
            ],
            _POLLINATOR_MV_DICT["labels"],
            _POLLINATOR_MV_DICT["colors"],
            "Habitat Quality",
            category_integer_labels=_POLLINATOR_MV_DICT["values"],
        )

        # Calculate pesticide toxicity
        logger.debug(f"{species} | Categorizing pesticide toxicity")

        categorize_and_bin_raster_hex(
            file_registry[
                "pesticide_exposure_[PESTICIDE]_[SPECIES]", pesticide, species
            ],
            file_registry[
                "pesticide_toxicity_[PESTICIDE]_[SPECIES]", pesticide, species
            ],
            [
                -float("inf"),
                scenario_variables[_TOXICITY_THRESHOLD_FIELD][species]
                / 2,  # NOTE Using half of the toxicity threshold as the low threshold
                scenario_variables[_TOXICITY_THRESHOLD_FIELD][species],
                float("inf"),
            ],
            _TOX_DICT["labels"],
            _TOX_DICT["colors"],
            "Pesticide Exposure",
            category_integer_labels=_TOX_DICT["values"],
        )

        logger.debug(
            f"{species} | Creating combined pesticide-pollination-habitat raster"
        )

        value_map = dict(
            zip(COMBINATION_DICT["complex_values"], COMBINATION_DICT["values"])
        )
        keys = sorted(np.array(list(value_map.keys())))
        values = np.array([value_map[x] for x in keys])

        def _combine_raster_op(*array_list):
            result = np.empty_like(array_list[0], dtype=DTYPE)
            result[:] = 0
            num_digits = 0
            for i, array in enumerate(array_list):
                result += array.astype(DTYPE) * (10**num_digits)
                num_digits += len(str(int(np.max(array))))

            # Reclassify to simpler format
            out_array = np.full(result.shape, UNSIGNED_NODATA_VAL, dtype=DTYPE)
            index = np.digitize(result, keys, right=True)
            out_array = values[index]
            return out_array

        # TODO This is outputting the multi-digit format 121 etc, fix and simplify
        pygeoprocessing.raster_map(
            _combine_raster_op,
            [
                file_registry["marginal_pollinator_value_[SPECIES]", species],
                file_registry["habitat_quality_[SPECIES]", species],
                file_registry[
                    "pesticide_toxicity_[PESTICIDE]_[SPECIES]", pesticide, species
                ],
            ],
            file_registry[
                "combination_matrix_[PESTICIDE]_[SPECIES]", pesticide, species
            ],
            target_nodata=UNSIGNED_NODATA_VAL,
            target_dtype=DTYPE,
        )

        # assign colors
        assign_raster_labels_and_hex_colors(
            file_registry[
                "combination_matrix_[PESTICIDE]_[SPECIES]", pesticide, species
            ],
            COMBINATION_DICT["values"],
            COMBINATION_DICT["labels"],
            random_hex_color(len(COMBINATION_DICT["values"])),
            "Combination",
        )

        # Reclassify combination matrix raster to decision matrix
        pygeoprocessing.reclassify_raster(
            (
                file_registry[
                    "combination_matrix_[PESTICIDE]_[SPECIES]", pesticide, species
                ],
                1,
            ),
            scenario_variables["decision_reclass_dict"],
            file_registry["decision_matrix_[PESTICIDE]_[SPECIES]", pesticide, species],
            gdal.GDT_Byte,
            UNSIGNED_NODATA_VAL,
        )

        assign_raster_labels_and_colors_from_csv(
            raster_path=file_registry[
                "decision_matrix_[PESTICIDE]_[SPECIES]", pesticide, species
            ],
            clr_path=file_registry["decision_matrix_colormap"],
            raster_band_label="Decision",
        )


@validation.invest_validator
def validate(args, limit_to=None):
    return validation.validate(args, MODEL_SPEC)


if __name__ == "__main__":
    config = Path("configs/configs.yaml")
    with open(config) as yaml_data_file:
        args = yaml.load(yaml_data_file, Loader=yaml.FullLoader)
    execute(args)
