def pytest_configure(config):
    config.addinivalue_line(
        "markers", "perf: opt-in performance / memory profile harness"
    )
