# spectral3d-cli

A command-line front end for [`spectral3d`](https://github.com/motoZ-crypto/spectral3d): spectral **identity** and mesh **repair** for 3D models.

Give it an OBJ and it derives a stable fingerprint that survives moving, rescaling, or roughing up the model. Hand it back the same fingerprint later and it tells you whether a second OBJ is the same object. A separate `repair` command rebuilds a messy mesh into a watertight solid.

## What it does

- **Register** an OBJ into a short identity hash, plus a small helper record that pins the hash down.
- **Verify** any OBJ against that record. Pass means same object up to pose, scale, and noise. Fail means it isn't.
- **Repair** a broken mesh into a closed manifold, with a full report of every edit.

## Build

Needs a recent stable Rust toolchain.

```sh
git clone https://github.com/motoZ-crypto/spectral3d-cli
cd spectral3d-cli
cargo build --release
```

## Usage

```
spectral3d-cli register <OBJ>              # register, explicitly
spectral3d-cli verify <OBJ> <RECORD.json>  # check an OBJ against a record
spectral3d-cli repair <OBJ>                # rebuild into a closed solid
```

### Register

```sh
spectral3d-cli register model.obj
# identity 9f2c…
# output   model-register.json
```

Reads the OBJ, computes its identity, and writes `model-register.json` next to it.

### Verify

```sh
spectral3d-cli verify other.obj model-register.json
# PASS ✓
```

Recomputes the identity of `other.obj` and compares it to the stored hash. `PASS ✓` if they match, `FAIL ✗` otherwise.
Use it to confirm a re-exported, re-scaled, or lightly remeshed model is still the original.

### Repair

```sh
spectral3d-cli repair model.obj
# dropped 3 degenerate, welded 12, dropped 1 dup faces; …
# result is now a closed manifold
# model-repair.obj
```

Runs a deterministic clean-up pipeline and writes `model-repair.obj` alongside the input.
