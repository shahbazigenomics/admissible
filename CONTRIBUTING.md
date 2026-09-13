# Contributing

## Getting set up

```bash
git clone https://github.com/shahbazigenomics/admissible
cd admissible
pip install -e ".[dev]"
pytest          # 110 tests, a few seconds
ruff check src tests
```

Zero runtime dependencies is a deliberate constraint, not an accident — see
"Nothing raises" in the README. A pull request that adds a required runtime
dependency needs to argue for it.

## The rules this codebase is built on

If you change code, these are the invariants to preserve. Each exists because the
obvious alternative is actively misleading in a clinical context.

1. **No check may raise.** A tool whose purpose is to describe broken files must not
   die on one. Every check returns a typed result with a stated reason. `tests/
   test_robustness.py` fuzzes for this.

2. **`UNKNOWN` means the inputs cannot answer the question.** It never means "not
   implemented". If you add a check, it reports `UNKNOWN` with a reason, never a
   reassuring default.

3. **Never report zero where you mean "could not evaluate".** Zero candidates looks
   like evidence of absence and gets written into papers as such.

4. **Never report a count without its denominator.** Zero over 20% of the target and
   zero over 95% are different findings.

5. **State evidence, not cause.** The tool reports what is observable and lists the
   mechanisms compatible with it. It does not assert that a duplicate came from the
   wet lab rather than from file handling, because the same observation arises from
   both.

6. **Bundle no annotation database.** Everything is read from the user's own files.
   Where a check needs a population frequency or a gene assignment, it takes one the
   user names and degrades honestly without it.

7. **Thresholds that are conventions must say so wherever they are printed.**

## Testing

New behaviour needs a test that would fail without it. Two kinds are especially
welcome:

- **Public data with published truth.** `tests/test_public_ceph.py` is the model:
  it skips cleanly when the data is absent, and it caught a real false-positive bug
  that every synthetic fixture had missed.
- **Tests that assert a limitation.** `test_second_degree_is_correct_on_average_
  but_not_per_pair` documents a failure mode rather than hiding it. Those are worth
  as much as the passing ones.

Fixtures are synthetic and deterministic. **No patient data, ever** — not in tests,
not in issues, not in bug reports. If you need to report a bug against real data,
reduce it to a synthetic file that reproduces the behaviour.

## Reporting a bug

Include the `admissible --version`, the command you ran, the JSON report if you have
one, and a synthetic input that reproduces it. "It said UNKNOWN and I expected a
number" is a useful report — that usually means either your inputs genuinely cannot
answer the question, or the reason given was not clear enough, and both are worth
fixing.
