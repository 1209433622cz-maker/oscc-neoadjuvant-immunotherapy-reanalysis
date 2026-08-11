# Analysis workspace

This directory contains the public computational workflow supporting the study.

- `analysis/` contains discovery, external-validation, sensitivity, figure-source and numerical-audit scripts.
- `config/` contains response harmonisation rules, model specifications and public-data manifests.
- `env/` contains environment setup, data acquisition and workflow entry points.

Run commands from the repository root and pass `-Workspace $PWD` to PowerShell entry points. Generated data, cached objects, figures, tables and logs remain outside version control according to the root `.gitignore`.

The numbered prefixes group related computational stages; they are execution identifiers, not manuscript or software release numbers.
