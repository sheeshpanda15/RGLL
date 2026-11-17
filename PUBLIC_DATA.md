# Additional public-data RI experiments

This protocol was fixed before fitting the new data. It does not replace the
Student/NLSY79 experiments, the simulation, or their outputs.

## Inputs and provenance

- UCI Apartment for Rent Classified, 10K version, DOI 10.24432/C5X623:
  https://archive.ics.uci.edu/dataset/555/apartment+for+rent+classified
- UCI Beijing Multi-Site Air Quality, all 12 original station files,
  DOI 10.24432/C5RK5G:
  https://archive.ics.uci.edu/dataset/501/beijing+multi+site+air+quality+data
- Both official sources use CC BY 4.0. Cite the datasets in the paper.
- The archived inputs under `real_data/new_public/` retain the source URLs
  and SHA256 hashes.

## Fixed analysis

RI only; 50 paired repetitions; LM, KM, OBS, CPF, BLM and PGR. Reuse the existing
estimators without changing their tuning: lambda=0, tau=5(p+1), SA25, ML-EM300,
tolerance 1e-5, eight-point KM K grid, one MiniBatchKmeans initialization and
25 clustering iterations. ClusterR's internal default seed remains 1, matching
the earlier experiments. This does not test residual augmentation or alternative
initializations. CPF/BLM are the existing pseudo-group adaptations, not CPFO or
BLMX. OBS/CPF/BLM share 20 artificial groups. Original-label RI models are
additional-information benchmarks, not algorithms feasible without labels.

Training/reference/test are approximately 50/20/30 percent by splitting unit.
Only training observations estimate X scaling and y mean/SD. Reference responses
are not supplied to estimators; test responses only score predictions. All
methods share each split and feature matrix. No outcome-dependent trimming,
feature selection, tuning, split replacement, or changing failure thresholds.

Apartment response: natural log monthly rent in USD. Retain positive-price,
positive-area apartment advertisements with numeric bedrooms, bathrooms and
valid coordinates; remove exact duplicate records ignoring ID and timestamp.
Keep extreme positive rents rather than choosing outcome cutoffs. The raw
`null` string denotes missingness. Do not use price_display, text, listing ID,
address, source, city or state in X. Full X: log area, bedrooms, bathrooms,
latitude, longitude, amenities-recorded flag, and fixed indicators for AC,
Dishwasher, Parking, Patio/Deck, Pool, Washer Dryer, Gym and Refrigerator.
Excl removes latitude/longitude, so here Excl specifically tests loss of
location proxies, not the historical datasets' exclusion lists. Missing amenity
tokens mean unrecorded, not verified absence. Split whole address+city signatures;
where no address is recorded use a coordinate/floor-plan signature. These are
conservative linkage rules, not verified property IDs. RI grouping benchmarks:
state, city (nested in state), source, city by source. Unseen label groups use
marginal predictions, without dropping test rows.

Beijing response: log1p daily mean PM2.5. Require at least 18 observed hours for
the target and each selected weather covariate on each station-day; do not
interpolate missing measurements. Full X: daily TEMP, PRES, DEWP, RAIN, WSPM,
mean sine/cosine of wind direction, annual sine/cosine and year index. Excl
removes calendar terms. No other pollutant, station label, block ID or raw date
enters X. Same-day observed weather is assumed available: this is prediction
for withheld periods, not a weather-free future forecasting task. Partition
calendar days into fixed 28-day blocks from 2013-03-01, discard the first and
last two days within each block for all stations, then randomly split blocks.
All station records from a calendar block stay together. This leaves at least
four omitted days between adjacent retained blocks; it reduces but does not
eliminate temporal dependence. RI grouping benchmarks: station, season and
station by season. Random-intercept models may remain misspecified because
temporal/spatial dependence is not explicitly fitted.

## Archived results

The current table workflow for Beijing is described in
[the standalone experiment README](beijing_experiment/README.md).
Raw failures and finite nonconverged fits are retained and reported. Primary
comparisons use common successful repetitions across all six main methods.
Original-group tables use the subset where those group baselines also succeed.
Reported mean +/- 1.96 MCSE describes repeated-split variability conditional
on the fixed dataset, not a population confidence interval or proof of
robustness. Output is separate under `result/public_data_RI/` and
`output/public_data_RI/`.
