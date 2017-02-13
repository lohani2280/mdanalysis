##Introduction
Software testing is extremely important to the MDAnalysis project, so that users of the package can be confident in the results of their analysis. With the current test package [nose, ceasing to be developed](http://nose.readthedocs.io/en/latest/#note-to-users) we have decided to move to using [py.test](http://doc.pytest.org/en/latest/).

## Current state
There are a total of 4731 test cases(that nose executes).
Running the testsuite, with pytest gives the following output:
`1789 failed, 2630 passed, 52 skipped, 1610 pytest-warnings, 26 error`

## How?

### Existing Plugins and how to import them
### Tests passing with pytests
### Tests failing currently
* Classes with `__init__` methods are not collected for test cases.

### Examples of failed tests
### Where can we replace nose with pytest



