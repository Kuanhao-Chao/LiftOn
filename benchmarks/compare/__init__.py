"""LiftOn three-way benchmark comparison pipeline.

Runs Liftoff, miniprot, and LiftOn independently on chromosome-subset data for
six reference->target benchmarks and compares annotation completeness, protein
identity, and DNA identity. See ``run_compare.py`` for the driver and the plan at
``plans/`` (benchmark_comparison_report.md) for the design.
"""
