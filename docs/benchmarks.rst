Timing Benchmarks
=================

pyCSRML fingerprinting speed is measured on five molecule-size-stratified
benchmark sets extracted from the CLinventory. Each set contains 500 molecules;
timing is the median of 5 repetitions of :meth:`~pyCSRML.Fingerprinter.fingerprint_batch`.

Benchmark sets
--------------

Sets are generated from the CLinventory by ``scripts/create_size_benchmarks.py``
and stored in ``tests/test_data/size_benchmarks/``.

.. list-table::
   :header-rows: 1
   :widths: 20 25 15

   * - Set
     - Heavy-atom range
     - Molecules
   * - ``bench_tiny``
     - 1 – 10
     - 500
   * - ``bench_small``
     - 11 – 20
     - 500
   * - ``bench_medium``
     - 21 – 35
     - 500
   * - ``bench_large``
     - 36 – 60
     - 500
   * - ``bench_xlarge``
     - 61 +
     - 500


pyCSRML timing results
----------------------

Measured on **Snapdragon X Elite X1E78100** (ARM64, 12 cores, ~32 GB RAM),
Python 3.14.2, RDKit 2025.09.3, NumPy 2.3.5.

.. list-table::
   :header-rows: 1
   :widths: 20 20 25 20

   * - Set
     - Heavy atoms
     - ToxPrint v2 (ms/mol)
     - TxP_PFAS v1 (ms/mol)
   * - ``bench_tiny``
     - 1 – 10
     - 3.76
     - 0.73
   * - ``bench_small``
     - 11 – 20
     - 5.47
     - 1.01
   * - ``bench_medium``
     - 21 – 35
     - 8.23
     - 1.53
   * - ``bench_large``
     - 36 – 60
     - 12.32
     - 2.19
   * - ``bench_xlarge``
     - 61 +
     - 23.20
     - 4.46

The TxP_PFAS v1 fingerprinter (129 bits) is roughly 5× faster than ToxPrint v2
(729 bits) across all size bins.  Both fingerprinters scale approximately
linearly with heavy-atom count: ToxPrint v2 ranges from 3.76 ms/mol (tiny) to
23.2 ms/mol (xlarge), a ~6× increase over the full size range.

Baseline file: ``tests/test_data/size_benchmarks/pycsrml_timing_baseline.json``.


Reproducing the benchmarks
--------------------------

.. code-block:: bash

   # Create benchmark sets (one-time; requires CLinventory CSV)
   python scripts/create_size_benchmarks.py

   # Time pyCSRML (saves pycsrml_timing_baseline.json)
   python scripts/benchmark_pycsrml_timing.py          # 5 reps (default)
   python scripts/benchmark_pycsrml_timing.py --reps 3  # fewer reps, faster

   # Run regression tests (require zips from ChemoTyper)
   pytest tests/test_benchmark_regression.py -v -m slow

Timing regression tests (``tests/test_benchmark_regression.py``) fail if any
set runs more than 30 % slower than the saved baseline.  They skip gracefully
until ``pycsrml_timing_baseline.json`` exists.


System information
------------------

Full details are recorded in
``tests/test_data/size_benchmarks/SYSTEM_INFO.md``.

.. list-table::
   :header-rows: 1
   :widths: 30 45

   * - Property
     - Value
   * - Host
     - ZenbookA14
   * - OS
     - Windows 11
   * - CPU
     - Snapdragon X Elite X1E78100 — Qualcomm Oryon
   * - Architecture
     - ARM64
   * - Physical cores
     - 12
   * - RAM
     - ~32 GB
   * - Python
     - 3.14.2
   * - RDKit
     - 2025.09.3
   * - NumPy
     - 2.3.5
   * - pyCSRML
     - 0.1.0 (editable install)
