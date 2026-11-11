from natcap.invest import gettext
from natcap.invest import spec
from natcap.invest.spec import u

from .nci_pollination.model_spec import MODEL_SPEC as input_model_spec

# Define component parts of the model spec prior to creation, so as to use the update= argument correctly

# 1. Define basic model input parameters

model_id = "pesticide_pollination"
model_title = "Pesticide Pollinator Risk"


# 2. Create inputs and input field order

inputs = input_model_spec.inputs

# Get index of the landcover table
for i, input_spec in enumerate(inputs):
    if input_spec.id == "landcover_biophysical_table_path":
        _landcover_table_index = i
        break

# Get index of the guild table
for i, input_spec in enumerate(inputs):
    if input_spec.id == "guild_table_path":
        _guild_table_index = i
        break


# Replace the landcover table spec with the new one that has the additional columns
# for pixel-level crop yield modeling and pesticide application rate
inputs[_landcover_table_index] = spec.CSVInput(
    id="landcover_biophysical_table_path",
    name=gettext("biophysical table"),
    about=gettext(
        "A table mapping each LULC class to nesting availability and floral"
        " abundance data for each substrate and season in that LULC class,"
        " as well as crop pollination dependency and value where appropriate."
        " All values in the LULC raster must have corresponding entries in this"
        " table."
    ),
    columns=[
        spec.LULC_TABLE_COLUMN,
        spec.RatioInput(
            id="nesting_[SUBSTRATE]_availability_index",
            about=gettext(
                "Index of availability of the given substrate in this LULC class."
                " Replace [SUBSTRATE] with substrate names matching those in the"
                " Guild Table, so that there is a column for each substrate."
            ),
            units=None,
        ),
        spec.RatioInput(
            id="floral_resources_[SEASON]_index",
            about=gettext(
                "Abundance of flowers during the given season in this LULC class."
                " This is the proportion of land area covered by flowers,"
                " multiplied by the proportion of the season for which there is"
                " that coverage. Replace [SEASON] with season names matching"
                " those in the Guild Table, so that there is a column for each"
                " season."
            ),
            units=None,
        ),
        spec.RatioInput(
            id="crop_pollinator_dependence_index",
            about=gettext("Proportion of annual crop yield dependent on pollinators."),
            required="crop_yield",
            units=None,
        ),
        spec.RatioInput(
            id="half_saturation_coefficient",
            about=gettext(
                "Half saturation coefficient for crop yield modeling."  # TODO add more explanation of what this means
            ),
            required="crop_yield",
            units=None,
        ),
        spec.NumberInput(
            id="crop_value_per_hectare",
            about=gettext("Annual crop yield value per hectare."),
            required="crop_yield",
            units=u.currency / u.hectare / u.year,
        ),
        spec.NumberInput(
            id="pesticide_application_rate",
            about=gettext("Annual pesticide application rate per hectare."),
            required="pesticide_analysis",
            units=u.kilograms / u.hectare / u.year,
        ),
    ],
    index_col="lucode",
)


# Replace the guild table spec with the new one that has the additional columns
# for pesticide exposure modeling
inputs[_guild_table_index] = spec.CSVInput(
    id="guild_table_path",
    name=gettext("Guild Table"),
    about=gettext(
        "A table mapping each pollinator species or guild of interest to its"
        " pollination- and pesticide-related parameters."
    ),
    columns=[
        spec.StringInput(
            id="species",
            about=gettext(
                "Unique name or identifier for each pollinator species or guild"
                " of interest."
            ),
            regexp=None,
        ),
        spec.RatioInput(
            id="nesting_suitability_[SUBSTRATE]_index",
            about=gettext(
                "Utilization of the substrate by this species, where 1 indicates"
                " the nesting substrate is fully utilized and 0 indicates it is"
                " not utilized at all. Replace [SUBSTRATE] with substrate names"
                " matching those in the Biophysical Table, so that there is a"
                " column for each substrate."
            ),
            units=None,
        ),
        spec.RatioInput(
            id="foraging_activity_[SEASON]_index",
            about=gettext(
                "Pollinator activity for this species/guild in each season. 1"
                " indicates maximum activity for the species/guild, and 0"
                " indicates no activity. Replace [SEASON] with season names"
                " matching those in the biophysical table, so that there is a"
                " column for each season."
            ),
            units=None,
        ),
        spec.NumberInput(
            id="alpha",
            about=gettext(
                "Average distance that this species or guild travels to forage on"
                " flowers."
            ),
            units=u.meter,
        ),
        spec.RatioInput(
            id="relative_abundance",
            about=gettext(
                "The proportion of total pollinator abundance that consists of"
                " this species/guild."
            ),
            units=None,
        ),
        spec.NumberInput(
            id="toxicity_threshold",
            about=gettext(
                "The exposure to the pesticide that results in toxicity (LD50, LOAEL, etc) for this species or guild."
            ),
            required="pesticide_analysis",
            units=u.mg / u.kg,  # TODO confirm units for toxicity threshold
        ),
        spec.NumberInput(  # TODO Check if this should be a ratio (0-1) or number (e.g. mg of pesticide that adheres to the pollinator)
            id="pollen_consumption_rate",
            about=gettext("How much pollen bees eat per day in mg/day"),
            required="pesticide_analysis",
            units=u.mg / u.day,
        ),
    ],
    index_col="species",
)

# Add scenario LULC raster and suffix
input_field_order = input_model_spec.input_field_order
input_field_order.extend(
    [
        ["pesticide_analysis", "pesticide_name", "drift_distance_m"],
        [
            "landscape_decision_analysis",
        ],
    ]
)
inputs.extend(
    [
        spec.BooleanInput(
            id="pesticide_analysis",
            name=gettext("calculate pesticide risk"),
            about=gettext(
                "Calculate pesticide risk for pollinators from a scenario raster."
            ),
        ),
        spec.StringInput(
            id="pesticide_name",
            name=gettext("pesticide name"),
            about=gettext("Name of the pesticide to use in the analysis."),
            required="pesticide_analysis",
        ),
        spec.NumberInput(
            id="drift_distance_m",
            name=gettext("drift distance (m)"),
            about=gettext("Average distance that the pesticide drifts via wind."),
            required="pesticide_analysis",
            allowed="pesticide_analysis",
            units=u.meters,
        ),
        spec.BooleanInput(
            id="landscape_decision_analysis",
            name=gettext("calculate landscape decision matrix"),
            about=gettext(
                "Perform landscape decision analysis for pollinators from a scenario raster."
            ),
            allowed="scenario_analysis and pesticide_analysis",
        ),
        # Hidden inputs, these are currently hard-coded into the package.
        spec.CSVInput(
            id="decision_table_path",
            name=gettext("decision table"),
            about=gettext("Path to the decision table CSV file."),
            required=False,
            hidden=True,
            columns=[
                spec.IntegerInput(
                    id="code",
                    about=gettext(
                        "Unique code from the combination of baseline pollination, restoration potential, and pesticide data."
                    ),
                    units=None,
                ),
                spec.IntegerInput(
                    id="decision_code",
                    about=gettext(
                        "Unique code for each decision in the decision table."
                    ),
                    units=None,
                ),
                spec.StringInput(
                    id="decision_label",
                    about=gettext("Label for each decision in the decision table."),
                    units=None,
                ),
                spec.StringInput(
                    id="decision_color",
                    about=gettext("Hex color for each decision in the decision table."),
                    units=None,
                ),
            ],
            index_col="code",
        ),
        spec.RatioInput(
            id="pollination_low_threshold",
            name=gettext("pollination threshold (low)"),
            about=gettext(
                "The threshold of pollination supply below which pollination is categorized as 'Low'."
            ),
            allowed="pesticide_analysis",
            required=False,
            hidden=True,
        ),
        spec.RatioInput(
            id="pollination_high_threshold",
            name=gettext("pollination threshold (high)"),
            about=gettext(
                "The threshold of pollination supply above which pollination is categorized as 'High'."
            ),
            allowed="pesticide_analysis",
            required=False,
            hidden=True,
        ),
        spec.RatioInput(
            id="restoration_low_threshold",
            name=gettext("restoration threshold (low)"),
            about=gettext(
                "The threshold of pollination restoration potential below which restoration is categorized as 'Low'."
            ),
            allowed="pesticide_analysis",
            required=False,
            hidden=True,
        ),
        spec.RatioInput(
            id="restoration_high_threshold",
            name=gettext("restoration threshold (high)"),
            about=gettext(
                "The threshold of pollination restoration potential above which restoration is categorized as 'High'."
            ),
            allowed="pesticide_analysis",
            required=False,
            hidden=True,
        ),
    ]
)

# 3. Create outputs
outputs = input_model_spec.outputs

# Add marginal value raster outputs
outputs.extend(
    [
        spec.SingleBandRasterOutput(
            id="pesticide_application_rate_[PESTICIDE]",
            path="pesticide_application_rate_[PESTICIDE].tif",
            about=gettext(
                "Raster representing the pesticide application rate per hectare."
            ),
            data_type=float,
            units=u.kilograms / u.hectare / u.year,
            created_if="pesticide_analysis",
        ),
        spec.SingleBandRasterOutput(
            id="pesticide_application_[PESTICIDE]",
            path="intermediate_outputs/pesticide_application_[PESTICIDE].tif",
            about=gettext(
                "Raster representing the total pesticide application per pixel"
            ),
            data_type=float,
            units=u.kilograms / u.year,
            created_if="pesticide_analysis",
        ),
        spec.SingleBandRasterOutput(
            id="drift_kernel_[PESTICIDE]",
            path="intermediate_outputs/drift_kernel_[PESTICIDE].tif",
            about=gettext("Exponential decay kernel for pesticide drift"),
            data_type=float,
            units=None,
            created_if="pesticide_analysis",
        ),
        spec.SingleBandRasterOutput(
            id="pesticide_drift_[PESTICIDE]",
            path="intermediate_outputs/pesticide_drift_[PESTICIDE].tif",
            about=gettext("Raster representing the pesticide drift per pixel"),
            data_type=float,
            units=u.kilograms / u.year,  # TODO Check units
            created_if="pesticide_analysis",
        ),
        spec.SingleBandRasterOutput(
            id="landscape_pesticide_load_[PESTICIDE]_[SPECIES]",
            path="intermediate_outputs/landscape_pesticide_load_[PESTICIDE]_[SPECIES].tif",
            about=gettext(
                "Raster representing the initial pesticide exposure per pixel for the given species."
            ),
            data_type=float,
            units=u.kilograms / u.year,  # TODO Check units
            created_if="pesticide_analysis",
        ),
        spec.SingleBandRasterOutput(
            id="pesticide_exposure_[PESTICIDE]_[SPECIES]",
            path="intermediate_outputs/pesticide_exposure_[PESTICIDE]_[SPECIES].tif",
            about=gettext(
                "Raster representing the pesticide exposure per pixel for the given species."
            ),
            data_type=float,
            units=u.kilograms / u.year,  # TODO Check units
        ),
        spec.SingleBandRasterOutput(
            id="pesticide_hazard_[PESTICIDE]_[SPECIES]",
            path="pesticide_hazard_[PESTICIDE]_[SPECIES].tif",
            about=gettext(
                "Raster representing the pesticide hazard per pixel for the given species."
            ),
            data_type=float,
            units=None,  # TODO Check units
            created_if="pesticide_analysis",
        ),
        spec.SingleBandRasterOutput(
            id="habitat_quality_[SPECIES]",
            path="habitat_quality_[SPECIES].tif",
            about=gettext(
                "Categorization of habitat quality for the given species, based on InVEST"
                " Pollinator Supply output and binned into Low, Medium, and High categories."
            ),
            data_type=int,
            units=None,
            created_if="landscape_decision_analysis",
        ),
        spec.SingleBandRasterOutput(
            id="marginal_pollinator_value_[SPECIES]",
            path="marginal_pollinator_value_[SPECIES].tif",
            about=gettext(
                "Categorization of marginal pollinator value for the given species under the scenario, based on InVEST"
                " Pollinator Supply output and binned into Low, Medium, and High categories."
            ),
            data_type=int,
            units=None,
            created_if="landscape_decision_analysis",
        ),
        spec.SingleBandRasterOutput(
            id="pesticide_toxicity_[PESTICIDE]_[SPECIES]",
            path="pesticide_toxicity_[PESTICIDE]_[SPECIES].tif",
            about=gettext(
                "Categorization of pesticide toxicity for the given species under the scenario,"
                " based on pesticide exposure and binned into Not Toxic and Toxic categories."
            ),
            data_type=int,
            units=None,
            created_if="landscape_decision_analysis",
        ),
        spec.SingleBandRasterOutput(
            id="combination_matrix_[PESTICIDE]_[SPECIES]",
            path="intermediate_outputs/combination_matrix_[PESTICIDE]_[SPECIES].tif",
            about=gettext(
                "Combined matrix of pollinator value, habitat quality, and pesticide toxicity for the given species."
            ),
            data_type=int,
            units=None,
            created_if="landscape_decision_analysis",
        ),
        spec.CSVOutput(
            id="decision_matrix_colormap",
            path="intermediate_outputs/decision_matrix_colormap.csv",
            about=gettext("Colormap used for the decision_matrix raster."),
            created_if="landscape_decision_analysis",
        ),
        spec.SingleBandRasterOutput(
            id="decision_matrix_[PESTICIDE]_[SPECIES]",
            path="decision_matrix_[PESTICIDE]_[SPECIES].tif",
            about=gettext(
                "Decision matrix for the given species, indicating the recommended action based on pollinator value, habitat quality, and pesticide toxicity."
            ),
            data_type=int,
            units=None,
            created_if="landscape_decision_analysis",
        ),
    ]
)


# 4. Create model spec from the base
MODEL_SPEC = input_model_spec.model_copy(
    deep=True,
    update={
        "model_id": model_id,
        "model_title": model_title,
        "input_field_order": input_field_order,
        "inputs": inputs,
        "outputs": outputs,
    },
)
