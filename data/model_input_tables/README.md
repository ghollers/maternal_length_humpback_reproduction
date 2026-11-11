# Anonymized tag-derived model-input tables

The seven aggregated tables the empirical sub-models are fit from, with whale IDs
replaced by anonymous codes. These are per-animal (or per-animal-per-behaviour)
summaries

IDs are role-specific: mothers `mother_01..` (20), calves `calf_01..` (31), each
consistent across its two files. The frequency/chord/area models are fit on all
deployments (mothers, calves, escort males, subadults), so animals that appear
only in `tailbeat_frequency.csv` are labelled `other_01..`. The two rows per ID
in the speed and fluking tables are the "rest" and "other" behaviour categories.

| File | Fits | Columns |
|---|---|---|
| `fluke_chord.csv` | `ChordLength ~ BodyLength` | BodyLength, ChordLength |
| `fluke_area.csv` | `FlukePlanarArea ~ BodyLength` | BodyLength, FlukePlanarArea |
| `tailbeat_frequency.csv` | `log10(mean_freq) ~ log10(BodyLength)` | ID, BodyLength, mean_freq |
| `swim_speed_calf.csv` | `mean_speed ~ Behavior + (1\|ID)` | ID, Behavior, BodyLength, mean_speed |
| `swim_speed_mother.csv` | `mean_speed ~ Behavior + (1\|ID)` | ID, Behavior, BodyLength, mean_speed |
| `prop_fluking_mother.csv` | `prop_sec ~ Behavior + (1\|ID)` | ID, Behavior, prop_sec |
| `prop_fluking_calf.csv` | `prop_sec ~ Behavior + BodyLength + (1\|ID)` | ID, Behavior, BodyLength, prop_sec |
