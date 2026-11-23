# Addendum to FOLLOWUP.md — two reporting corrections

Scope: wording only. No numerical experiment was rerun for this addendum, and no
result file was modified. `out/validation_report.json`,
`out/validation_report_direct.json`, `out/validation_report_dor.json` and
`out/edit_implementation_comparison.json` are unchanged.

---

## A1. "Failed 62/64" was a garbled sentence

`FOLLOWUP.md` §1.2 says of the original implementation-choice criterion:
"**The original failed at 62/64.**" That sentence is wrong as written, and reads as
though the direct decoder-delta lost in 62 of 64 conditions. It did not. The facts,
all from `out/edit_implementation_comparison.json`:

| Quantity | Value |
|---|---|
| Non-zero edit conditions compared | 64 |
| Conditions where DD had the **lower** relative L2 error than DoR | **62** |
| Conditions where DD had the **higher** relative L2 error than DoR | **2** |
| The 2 exceptions | `O61866`, case B, `strongest_active` — DD worse by 1.2% and 6.2% |

The original check #37 asked a single all-or-nothing question: *is
`rel_err_DD <= rel_err_DoR` in **every** condition?* Because 2 conditions violated it,
that **one check** returned FAIL. "Failed 62/64" conflated the check's binary outcome
with the per-condition tally and inverted the latter.

**Correct wording, to be used in place of the original sentence:**

> The direct decoder-delta had the lower relative L2 error in 62 of 64 conditions and
> the higher error in 2. The original criterion required the ordering to hold in every
> condition, so those 2 conditions made that single check fail.

The claim in §3 of FOLLOWUP.md — "62/64 conditions" as a point in favour of the direct
form — was always the correct direction and is unaffected.

---

## A2. The 1.1× per-condition tolerance is retrospective

Check #38 ("no condition worse than 1.1×", observed worst ratio 1.066) was written
**after** the 96-condition comparison had been run and its worst ratio observed. The
threshold was chosen in the knowledge that 1.066 was the number it had to admit.

It is therefore a **retrospective, post-observation tolerance**, and this addendum
labels it as one. The consequences are specific:

- **It is not an independent test.** A bound fitted to the data it judges cannot fail
  on that data. Check #38 passing carries no evidence beyond the descriptive statement
  "the worst observed per-condition regression was 1.066×".
- **It should not be quoted as a passed criterion** in the same register as the
  prospective checks (no-op exactness, delta localisation, top-k support integrity,
  reference agreement), which were specified before their data existed.
- **It becomes a real test only prospectively** — if 1.1× is now frozen and applied
  unchanged to conditions not in this comparison (different proteins, features,
  strengths, or positions). Until that happens it is a description, not a bar.

The corresponding note in `edit_implementation_comparison.json` under
`impl_comparison.criterion_note` records the criterion change but does not carry this
label; the file was left untouched, so read that note together with this section.

Nothing about the recommendation changes. The case for the direct decoder-delta rests
on the prospectively-defined quantities — median relative L2 error 3.367e-07 against
1.554e-06, maximum 9.829e-07 against 7.000e-06, minimum significant digits 6.01
against 5.15, and the absence of catastrophic cancellation by construction — none of
which depend on check #38.
