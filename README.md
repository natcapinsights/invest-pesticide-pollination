# NCI Pesticide Pollination

An [InVEST](https://naturalcapitalproject.stanford.edu/software/invest) plugin that extends the InVEST Pollination model to assess pesticide exposure risk to pollinators alongside their current and potential habitat quality.

## What it does

Starting from a land-use/land-cover (LULC) raster, this model chains together three analyses:

1. **Pollinator habitat quality and marginal value** — runs the [`nci_pollination`](https://github.com/natcapinsights/nci-pollination) marginal-value model (itself built on the InVEST Pollination model) to estimate current pollinator supply and the marginal value of restoring each pixel, per species/guild.
    - Requires the typical InVEST Pollinator model input data, plus an additional LULC raster representing the potential restoration scenario.
2. **Pesticide exposure and hazard** — reclassifies the landscape by pesticide application rate, spreads it via a wind-drift kernel, convolves it with each species' foraging range, and converts the resulting pollen concentration into a per-species pesticide hazard map (exposure relative to a toxicity threshold).
    - Requires an additional `pesticide_application_rate` field in the biophysical table and additional `toxicity_threshold` and `pollen_consumption_rate` fields in the guild table, along with new model arguments of `pesticide_name` and `drift_distance_m`.
3. **Landscape decision matrix** — bins habitat quality, restoration potential, and pesticide toxicity into Low/Medium/High categories, combines them, and reclassifies the combination into a recommended conservation action (e.g. *Protect*, *Restore*, *Reduce*, and combinations therein) using a configurable decision table.

The output is a set of rasters that can be used to prioritize where to protect existing pollinator habitat, where restoration would have the most impact, and where pesticide use poses the greatest risk for each pollinator species/guild of interest.

## Installation

```bash
pip install .
```

Requires Python ≥3.10 and GDAL. Depends on `natcap.invest` and `pygeoprocessing`.

## Usage

This package is an InVEST plugin (`package_name = invest_invest_pesticide_pollination` in [pyproject.toml](pyproject.toml)) and can be run through the `natcap.invest` plugin interface.

## New data requirements

There are two additional sub-modules with optionally-required input data, each toggled by a boolean argument (`pesticide_analysis`, `landscape_decision_analysis`).

### Pesticide Analysis

The following input data are required if the `pesticide_analysis` modeling option is selected:

- `pesticide_name`: name of the pesticide to use in the analysis
- `drift_distance_m`: average distance that the pesticide drifts via wind

The following columns are newly required in the biophysical and guild tables if the `pesticide_analysis` modeling option is selected:

#### Biophysical Table

- `pesticide_application_rate`: the annual pesticide application rate for the chemical in question in `kg/ha`

#### Guild Table

- `toxicity_threshold`: the exposure rate to the pesticide that results in toxicity (LD50, LOAEL, etc) for this species or guild, in `mg pollen / kg chemical`
- `pollen_consumption_rate`: how much pollen this species or guild eats per day in `mg/day`

### Landscape Decision Analysis

There are no additional arguments required beyond the boolean toggle `landscape_decision_analysis`.

## New output rasters

### Core pollination outputs

| Output | Description | Filename |
| --- | --- | --- |
| Pollinator supply | Total pollinator supply across all species. This is considered the primary biophysical output and can be considered an index of pollinator habitat quality | `pollinator_supply_index.tif` |
| Pollinator visitation | Total pollinator visitation across all species, replacing InVEST's floral-weighted visitation measure | `pollinator_visitation_index.tif` |
| Habitat nesting density | Nesting density by species, used in marginal value calculations | `intermediate_outputs/habitat_nesting_density_[SPECIES].tif` |

### Scenario analysis outputs (`scenario_analysis` option)

| Output | Description | Filename |
| --- | --- | --- |
| Attributed change in pollinator supply | Changes in pollinator supply attributed back to the scenario locations that caused it, across all species | `d_pollinator_supply_index.tif` |
| Attributed change in pollinator supply by species | Attributed changes in pollinator supply for the given species | `intermediate_outputs/d_pollinator_supply_[SPECIES].tif` |

### Pesticide Analysis (`pesticide_analysis` option)

| Output | Description | Filename |
| --- | --- | --- |
| Pesticide application rate per hectare | Amount of pesticide applied per hectare per crop for the study area | `pesticide_application_rate_[PESTICIDE].tif` |
| Pesticide hazard to SPECIES | An index of pesticide risk to each SPECIES from the analyzed chemical | `pesticide_hazard_[PESTICIDE]_[SPECIES].tif` |
| Species-specific pesticide exposure | Pesticide exposure via pollen consumption by SPECIES, in `kg/yr` | `intermediate_outputs/pesticide_exposure_[PESTICIDE]_[SPECIES].tif` |

### Landscape Decision Analysis (`landscape_decision_analysis` option)

| Output | Description | Filename |
| --- | --- | --- |
| Habitat Quality | Categorized habitat quality measure into `Low`, `Medium`, and `High` | `habitat_quality_[SPECIES].tif` |
| Restoration Potential | Categorized restoration potential into `Low`, `Medium`, and `High` | `marginal_pollinator_value_[SPECIES].tif` |
| Pesticide Toxicity | Categorized pesticide hazard into `Low`, `Medium`, and `High` | `pesticide_toxicity_[PESTICIDE]_[SPECIES].tif` |
| Decision Matrix | Decision matrix for the given species, indicating the recommended action based on restoration value, habitat quality, and pesticide toxicity with a colormap and labels drawn from [`decision_table.csv`](invest_pesticide_pollination/decision_table.csv) | `decision_matrix_[PESTICIDE]_[SPECIES].tif` |

Most other intermediate calculation rasters are written to `intermediate_outputs/` and are not part of the primary output set.

See [`model_spec.py`](invest_pesticide_pollination/model_spec.py) for the full input/output specification.

## Project context

Developed as part of the USFWS Great Lakes Pollinators project.

## License

Apache-2.0
