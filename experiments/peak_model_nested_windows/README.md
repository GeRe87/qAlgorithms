# Nested hyperbolic windows on TripleTOF profiles

Pair the 9- and 17-point windows from the same injection, scan and original
anchor index. Both fits use dimensionless index coordinates with step 0.35.
For paired fits, verify exact nesting of the raw and m/z arrays and exact
agreement of stored whole-line areas between the two input artifacts. Keep
cohort membership in the CSV and deduplicate across cohorts for pooled
and injection-level results.

The whole-line area change is A17 / A9 - 1. To isolate how the refitted curves
differ inside the same observed support, integrate both curves over the
9-point span, retaining each fit's original raw least-squares amplitude.
The 10% whole-line change threshold is descriptive. Real reference areas
are unavailable, so the results quantify sensitivity to window choice.

## Results

Of 19,803 distinct anchors with both window sizes recorded, 318 have both
matched fits. Within those 318, 110 have two positive finite whole-line
areas, 30 have a finite area only at 9 points, 68 only at 17 points, and
110 at neither size.

Among the 110 anchors with two areas, 55 change by more than 10%.
The signed changes at the 10th, 50th and 90th percentiles are -9.3%, +8.2%
and +22.9%. Within the common 9-point span, the median absolute change is
13.0% across the 318 matched fits.

## Exploratory diagnostic transfer

Use three features available from the 9-point fit: estimated whole-line area
fraction outside the window, RMS of log-fit residuals, and mean endpoint
intensity divided by peak window intensity. Endpoint intensity is an
observed proxy and can contain signal from the peak itself. Fit a linear
model on one injection using standardized features and absolute log area
ratio as the target. On the other injection, compare the Spearman
correlation of its predictions with that of each feature individually.

| Train / evaluate | N train / evaluate | Tail fraction | Log RMS | Endpoint intensity | Combined |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 / 19 | 38 / 72 | 0.350 | 0.142 | 0.357 | 0.291 |
| 19 / 1 | 72 / 38 | 0.292 | 0.356 | 0.272 | 0.368 |

The combined score improves on tail fraction in one direction and degrades
in the other. The comparison uses two injections and conditions on both
fits having a finite area, so the score's behavior at anchors with missing
areas is unknown. The current evidence supports recording window
sensitivity and area availability as diagnostics. A threshold for
production would require validation on additional independent injections.

## Reproduce

Use the segmentation evidence from Actions run 42999805753 and the area
evidence from run 43039277241. The script checks the commit in each input
artifact. With the dependencies pinned in the area audit installed, run:

```sh
python -B experiments/peak_model_nested_windows/nested_windows.py segments areas results
```
