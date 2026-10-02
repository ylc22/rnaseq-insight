# Contributing

Thanks for helping improve RNAseq Insight.

## Workflow

1. Create a focused branch from `main`.
2. Make one logical change per pull request.
3. Run the test suite locally with `pytest`.
4. Update documentation when behavior changes.
5. Open a pull request with a clear summary and validation notes.

## Development setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pytest
```

## Pull request checklist

- [ ] Change is scoped and documented
- [ ] Tests pass locally
- [ ] New behavior is covered by tests where appropriate
- [ ] No private or proprietary data is included
