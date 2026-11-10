# QCanopy: Batched, Reuse-Aware GPU Simulation for Large-Scale Quantum Circuit Equivalence Checking

The deterministic circuit-pair generator is available as `qcanopy_generate` and
the reusable `qcanopy::generator` C++17 library. It writes labeled equivalent and
non-equivalent pairs in OpenQASM 3.0 or portable JSON IR, with seeds, construction
proofs, and raw/normalized prefix-sharing measurements.

See [generator build, usage, formats, and validation](docs/test_case_generator.md).
