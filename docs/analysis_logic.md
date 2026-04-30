# Analysis Logic

For each PAM `i`, the control abundance is `X_i` and experiment abundance is `Y_i`, both stored as percentages.

## Raw Fold Change

```text
FC_raw_i = Y_i / (X_i + epsilon)
```

## Corrected Fold Change

```text
correction_factor = median(FC_raw_i for all PAMs)
Y_corrected_i = Y_i / correction_factor
FC_corrected_i = Y_corrected_i / (X_i + epsilon)
```

The fixed workflow intentionally does not apply a low-abundance filter. This keeps strong depletion candidates from being removed by an abundance cutoff.

## Thresholds

For sigma thresholds:

```text
logFC_i = ln(FC_i)
lower = exp(mean(logFC) - sigma * std(logFC))
upper = exp(mean(logFC) + sigma * std(logFC))
```

For hard fold-change threshold:

```text
lower = 0.15
upper = 1 / 0.15
```

PAMs with `FC <= lower` are depleted. PAMs with `FC >= upper` are enriched.

## Interpretation

- A robust PAM signal should appear in `CasX` positive control.
- A raw-only global shift that disappears after correction suggests broad abundance shift, toxicity, bottlenecking, or compositional effects.
- A signal consistent across raw and corrected modes is more likely sequence-specific.
